#!/usr/bin/bash

# --- run_parsnp_analysis.sh ---
# Parsnp analysis for core genome alignment and phylogenetic tree construction
# Performs rapid core genome multi-alignment for closely related genomes

set -e

# --- CONDA SETUP ---
source ~/miniforge3/etc/profile.d/conda.sh

# --- USER-DEFINED VARIABLES ---

# Reference genome (will be used as anchor) - using Nissle RefSeq genome
REFERENCE_GENOME="/bulkpool/reference_data/Nissel1917/ncbi_dataset/data/GCF_000714595.1/GCF_000714595.1_ASM71459v1_genomic.fna"

# Number of threads (default: 16)
THREADS=${1:-16}

# Output directory
OUTPUT_DIR="parsnp_results_nissle_refseq"

# Minimum length for core genome alignment (default: 1000bp)
MIN_LENGTH=${2:-1000}

echo "============================================"
echo "   Parsnp Core Genome Analysis"
echo "============================================"
echo "Reference: Nissle RefSeq (GCF_000714595.1)"
echo "Threads: $THREADS"
echo "Minimum alignment length: $MIN_LENGTH bp"
echo "Output: $OUTPUT_DIR"
echo "============================================"
echo ""

# Activate assembly-qc environment (contains parsnp)
conda activate assembly-qc

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Create a directory for assemblies
ASSEMBLY_DIR="${OUTPUT_DIR}/assemblies"
mkdir -p "$ASSEMBLY_DIR"

# Find and copy all completed assemblies
echo "Collecting completed assemblies..."
ASSEMBLY_COUNT=0

for asm in $(find . -path "*/unicycler_output/assembly.fasta" -type f | sort); do
    sample=$(echo "$asm" | sed 's|./\(.*\)_hybrid_assembly/.*|\1|')
    echo "  ✓ $sample"
    cp "$asm" "${ASSEMBLY_DIR}/${sample}.fasta"
    ASSEMBLY_COUNT=$((ASSEMBLY_COUNT + 1))
done

if [ $ASSEMBLY_COUNT -eq 0 ]; then
    echo "ERROR: No assembly.fasta files found"
    exit 1
fi

echo ""
echo "Found $ASSEMBLY_COUNT assemblies"
echo ""

# Check reference
if [ ! -f "$REFERENCE_GENOME" ]; then
    echo "ERROR: Reference genome not found: $REFERENCE_GENOME"
    exit 1
fi

# Copy reference to assembly directory
echo "Using Nissle RefSeq genome as anchor..."
cp "$REFERENCE_GENOME" "${ASSEMBLY_DIR}/reference_Nissle_RefSeq.fasta"

echo ""
echo "Running Parsnp..."
echo "This will:"
echo "  1. Identify core genome regions conserved across all samples"
echo "  2. Construct whole-genome alignment"
echo "  3. Build maximum likelihood phylogenetic tree"
echo "  4. Identify SNPs and structural variants"
echo ""
echo "This may take 10-30 minutes depending on number of samples..."
echo ""

# Run parsnp
# -r: reference genome (anchor)
# -d: directory containing genomes
# -p: number of threads
# -c: force inclusion of all genomes (no filtering)
# -o: output directory
# -v: verbose output

parsnp \
    -r "${ASSEMBLY_DIR}/reference_Nissle_RefSeq.fasta" \
    -d "$ASSEMBLY_DIR" \
    -p $THREADS \
    -c \
    -o "$OUTPUT_DIR" \
    -v

echo ""
echo "============================================"
echo "   Parsnp Analysis Complete!"
echo "============================================"
echo ""
echo "Results saved to: $OUTPUT_DIR/"
echo ""
echo "Key output files:"
echo "  - parsnp.tree          : Newick format phylogenetic tree"
echo "  - parsnp.ggr           : Harvest GGR format (core genome alignment)"
echo "  - parsnp.xmfa          : XMFA format multi-alignment"
echo "  - parsnp.vcf           : Variant calls (SNPs, indels)"
echo ""

# Convert tree to better format if possible
if command -v harvesttools &> /dev/null; then
    echo "Generating additional outputs with harvesttools..."

    # Generate SNP matrix
    harvesttools -i "${OUTPUT_DIR}/parsnp.ggr" -S "${OUTPUT_DIR}/snp_matrix.tsv"
    echo "  ✓ SNP matrix: ${OUTPUT_DIR}/snp_matrix.tsv"

    # Generate phylip alignment
    harvesttools -i "${OUTPUT_DIR}/parsnp.ggr" -M "${OUTPUT_DIR}/alignment.phylip"
    echo "  ✓ Phylip alignment: ${OUTPUT_DIR}/alignment.phylip"

    echo ""
fi

echo "Visualize the tree:"
echo "  1. Use FigTree, iTOL, or online tools"
echo "  2. Tree file: ${OUTPUT_DIR}/parsnp.tree"
echo ""

# Parse VCF for SNP summary
if [ -f "${OUTPUT_DIR}/parsnp.vcf" ]; then
    echo "SNP Summary:"
    echo "============"

    total_snps=$(grep -v "^#" "${OUTPUT_DIR}/parsnp.vcf" | wc -l)
    echo "Total core genome SNPs: $total_snps"
    echo ""

    # Count SNPs per sample (from VCF)
    echo "SNPs per sample (relative to reference):"
    grep -v "^#" "${OUTPUT_DIR}/parsnp.vcf" | head -1 | cut -f10- | tr '\t' '\n' | nl -v 0 | while read num sample; do
        if [ ! -z "$sample" ]; then
            snp_count=$(grep -v "^#" "${OUTPUT_DIR}/parsnp.vcf" | cut -f$((10+num)) | grep -c "1")
            echo "  Sample $num: $snp_count SNPs"
        fi
    done 2>/dev/null || echo "  (Run harvesttools for detailed SNP matrix)"
fi

echo ""
echo "For downstream analysis:"
echo "  - View tree: cat ${OUTPUT_DIR}/parsnp.tree"
echo "  - SNP variants: ${OUTPUT_DIR}/parsnp.vcf"
echo "  - Core genome alignment: ${OUTPUT_DIR}/parsnp.xmfa"
echo ""

conda deactivate

echo "============================================"
