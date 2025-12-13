#!/usr/bin/bash

# --- run_haslam_comparative_analysis.sh ---
# Comparative genomics analysis for Haslam samples
# Uses N2SKTQ_1_1 as the reference (ancestor) strain
# Runs Parsnp for core genome alignment and phylogeny

set -e

# --- CONDA SETUP ---
source ~/miniforge3/etc/profile.d/conda.sh

# --- CONFIGURATION ---

# Reference assembly (N2SKTQ_1_1 - your ancestor strain)
# Primary location (most likely)
REFERENCE_ASSEMBLY="N2SKTQ_results/N2SKTQ_1_1_hybrid_assembly/unicycler_output/assembly.fasta"

# Alternative location
if [ ! -f "$REFERENCE_ASSEMBLY" ]; then
    REFERENCE_ASSEMBLY="N2SKTQ_results/N2SKTQ_1_1_hybrid_assembly/N2SKTQ_1_1_hybrid_assembly.fasta"
fi

# Check if reference exists
if [ ! -f "$REFERENCE_ASSEMBLY" ]; then
    echo "ERROR: Reference assembly not found at either:"
    echo "  - N2SKTQ_results/N2SKTQ_1_1_hybrid_assembly/unicycler_output/assembly.fasta"
    echo "  - N2SKTQ_results/N2SKTQ_1_1_hybrid_assembly/N2SKTQ_1_1_hybrid_assembly.fasta"
    exit 1
fi

# Number of threads
THREADS=${1:-32}

# Minimum alignment length (bp)
MIN_LEN=${2:-1000}

# Output directory
OUTPUT_DIR="haslam_comparative_results"

# Haslam samples
HASLAM_SAMPLES=(
    "Control-1"
    "Control-2"
    "Control-3"
    "Control-4"
    "NIL-1"
    "NIL-2"
)

echo "============================================"
echo "   Haslam Comparative Genomics Analysis"
echo "============================================"
echo "Reference: N2SKTQ_1_1 (ancestor strain)"
echo "Samples: ${#HASLAM_SAMPLES[@]} new strains"
echo "Threads: $THREADS"
echo "Min alignment length: ${MIN_LEN} bp"
echo "Output: $OUTPUT_DIR"
echo "============================================"
echo ""

echo "✓ Reference assembly found: $REFERENCE_ASSEMBLY"
echo ""

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Create assemblies subdirectory
ASSEMBLIES_DIR="${OUTPUT_DIR}/assemblies"
mkdir -p "$ASSEMBLIES_DIR"

echo "Step 1: Collecting assemblies..."
echo ""

# Copy reference assembly
echo "  Copying reference (N2SKTQ_1_1)..."
cp "$REFERENCE_ASSEMBLY" "${ASSEMBLIES_DIR}/N2SKTQ_1_1.fasta"

# Copy Haslam assemblies
FOUND_COUNT=0
for sample in "${HASLAM_SAMPLES[@]}"; do
    assembly_file="${sample}_hybrid_assembly/${sample}_hybrid_assembly.fasta"

    if [ -f "$assembly_file" ]; then
        echo "  ✓ Found assembly: $sample"
        cp "$assembly_file" "${ASSEMBLIES_DIR}/${sample}.fasta"
        ((FOUND_COUNT++))
    else
        echo "  ✗ Missing assembly: $sample"
        echo "    Expected: $assembly_file"
    fi
done

echo ""
echo "Collected assemblies: $((FOUND_COUNT + 1)) total (1 reference + $FOUND_COUNT samples)"
echo ""

if [ $FOUND_COUNT -eq 0 ]; then
    echo "ERROR: No Haslam assemblies found"
    echo "Please run hybrid assemblies first: ./run_haslam_hybrid_assemblies.sh"
    exit 1
fi

# Activate assembly-qc environment for Parsnp
echo "Step 2: Running Parsnp core genome analysis..."
echo ""

conda activate assembly-qc

# Run Parsnp with N2SKTQ_1_1 as reference
echo "Running Parsnp..."
echo "  Reference: N2SKTQ_1_1"
echo "  Assemblies: ${ASSEMBLIES_DIR}/*.fasta"
echo ""

parsnp -r "${ASSEMBLIES_DIR}/N2SKTQ_1_1.fasta" \
       -d "$ASSEMBLIES_DIR" \
       -p $THREADS \
       -c \
       -o "$OUTPUT_DIR" \
       -v \
       2>&1 | tee "${OUTPUT_DIR}/parsnp.log"

echo ""
echo "✓ Parsnp analysis complete"
echo ""

# Generate SNP matrix if not already created by Parsnp
if [ -f "${OUTPUT_DIR}/parsnp.vcf" ]; then
    echo "Step 3: Generating SNP distance matrix..."

    # Use bcftools to extract SNP info if available
    if command -v bcftools &> /dev/null; then
        bcftools stats "${OUTPUT_DIR}/parsnp.vcf" > "${OUTPUT_DIR}/snp_stats.txt"
        echo "  ✓ SNP statistics saved to snp_stats.txt"
    fi
fi

conda deactivate

echo ""
echo "============================================"
echo "   Comparative Analysis Complete!"
echo "============================================"
echo ""
echo "Output directory: $OUTPUT_DIR"
echo ""
echo "Key output files:"
echo "  - parsnp.tree          : Phylogenetic tree (Newick format)"
echo "  - parsnp.ggr           : Core genome alignment (GGR format)"
echo "  - parsnp.vcf           : Core genome SNPs (VCF format)"
echo "  - parsnp.xmfa          : Multi-alignment (XMFA format)"
echo "  - assemblies/          : Collected assembly files"
echo ""

# Display tree if it exists
if [ -f "${OUTPUT_DIR}/parsnp.tree" ]; then
    echo "Phylogenetic tree:"
    echo "=================="
    cat "${OUTPUT_DIR}/parsnp.tree"
    echo ""
fi

# Count SNPs if VCF exists
if [ -f "${OUTPUT_DIR}/parsnp.vcf" ]; then
    SNP_COUNT=$(grep -v "^#" "${OUTPUT_DIR}/parsnp.vcf" | wc -l)
    echo "Core genome SNPs detected: $SNP_COUNT"
    echo ""
fi

echo "Visualization options:"
echo "  1. View tree in FigTree: figtree ${OUTPUT_DIR}/parsnp.tree"
echo "  2. Upload tree to iTOL: https://itol.embl.de/"
echo "  3. Use Harvest Gingr: gingr ${OUTPUT_DIR}/parsnp.ggr"
echo ""

echo "Next steps:"
echo "  - Review phylogenetic relationships"
echo "  - Identify strain-specific SNPs"
echo "  - Analyze core genome variants"
echo ""
