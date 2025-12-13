#!/usr/bin/bash

# --- run_snippy_core_analysis.sh ---
# Snippy variant calling and core SNP analysis
# Runs snippy on reads against reference, then snippy-core for core genome SNPs

set -e

# --- CONDA SETUP ---
source ~/miniforge3/etc/profile.d/conda.sh

# --- USER-DEFINED VARIABLES ---

# Reference genome (using N2SKTQ_1_1 hybrid assembly)
REFERENCE_GENOME="N2SKTQ_1_1_hybrid_assembly/unicycler_output/assembly.fasta"

# Number of threads per snippy job (default: 8)
THREADS_PER_JOB=${1:-8}

# Output directory
OUTPUT_DIR="snippy_results"

# Samples to process (excluding reference sample 1_1)
SAMPLES=(
    "N2SKTQ_2_2"
    "N2SKTQ_3_3"
    "N2SKTQ_4_4"
    "N2SKTQ_5_5"
    "N2SKTQ_6_6"
    "N2SKTQ_7_7"
    "N2SKTQ_8_8"
    "N2SKTQ_9_9"
    "N2SKTQ_10_10"
    "N2SKTQ_11_11"
    "N2SKTQ_12_12"
)

echo "============================================"
echo "   Snippy Core SNP Analysis"
echo "============================================"
echo "Reference: N2SKTQ_1_1 hybrid assembly"
echo "Threads per job: $THREADS_PER_JOB"
echo "Samples to process: ${#SAMPLES[@]}"
echo "Output: $OUTPUT_DIR"
echo "============================================"
echo ""

# Activate snippy environment
conda activate snippy

# Check reference exists
if [ ! -f "$REFERENCE_GENOME" ]; then
    echo "ERROR: Reference genome not found: $REFERENCE_GENOME"
    exit 1
fi

# Create output directory
mkdir -p "$OUTPUT_DIR"

echo "Step 1: Running snippy on each sample against reference..."
echo ""

# Run snippy for each sample
for sample in "${SAMPLES[@]}"; do
    echo "Processing $sample..."

    # Find Illumina reads
    R1="/fastpool/active_data/ecoli_genomics/N2SKTQ_results/${sample}_illumina_R1.fastq.gz"
    R2="/fastpool/active_data/ecoli_genomics/N2SKTQ_results/${sample}_illumina_R2.fastq.gz"

    if [ ! -f "$R1" ] || [ ! -f "$R2" ]; then
        echo "  WARNING: Reads not found for $sample, skipping..."
        continue
    fi

    # Output directory for this sample
    OUTDIR="${OUTPUT_DIR}/${sample}"

    if [ -d "$OUTDIR" ]; then
        echo "  ✓ $sample already processed, skipping..."
        continue
    fi

    echo "  Running snippy for $sample..."
    snippy \
        --cpus $THREADS_PER_JOB \
        --outdir "$OUTDIR" \
        --ref "$REFERENCE_GENOME" \
        --R1 "$R1" \
        --R2 "$R2" \
        --cleanup \
        2>&1 | tee "${OUTPUT_DIR}/${sample}_snippy.log"

    echo "  ✓ $sample complete"
    echo ""
done

echo ""
echo "Step 2: Running snippy-core to generate core SNP alignment..."
echo ""

# Collect all snippy output directories
SNIPPY_DIRS=()
for sample in "${SAMPLES[@]}"; do
    OUTDIR="${OUTPUT_DIR}/${sample}"
    if [ -d "$OUTDIR" ]; then
        SNIPPY_DIRS+=("$OUTDIR")
    fi
done

if [ ${#SNIPPY_DIRS[@]} -eq 0 ]; then
    echo "ERROR: No snippy results found"
    exit 1
fi

echo "Found ${#SNIPPY_DIRS[@]} snippy results"
echo ""

# Run snippy-core
echo "Running snippy-core..."
snippy-core \
    --ref "$REFERENCE_GENOME" \
    --prefix "${OUTPUT_DIR}/core" \
    "${SNIPPY_DIRS[@]}"

echo ""
echo "============================================"
echo "   Snippy Core Analysis Complete!"
echo "============================================"
echo ""
echo "Results saved to: $OUTPUT_DIR/"
echo ""
echo "Key output files:"
echo "  - core.aln           : Core SNP alignment (FASTA)"
echo "  - core.tab           : Core SNP table"
echo "  - core.vcf           : Core SNPs in VCF format"
echo "  - core.txt           : Core SNP summary"
echo ""

# Build phylogenetic tree if snp-sites and iqtree available
if command -v snp-sites &> /dev/null && command -v iqtree &> /dev/null; then
    echo "Building phylogenetic tree..."

    # Extract variable sites
    snp-sites -o "${OUTPUT_DIR}/core.snps.aln" "${OUTPUT_DIR}/core.aln"

    # Build tree with IQ-TREE
    iqtree \
        -s "${OUTPUT_DIR}/core.snps.aln" \
        -m GTR+G \
        -nt AUTO \
        -pre "${OUTPUT_DIR}/core.tree"

    echo ""
    echo "Phylogenetic tree: ${OUTPUT_DIR}/core.tree.treefile"
fi

# Generate SNP matrix if snippy-core produced one
if [ -f "${OUTPUT_DIR}/core.tab" ]; then
    echo ""
    echo "SNP Summary:"
    echo "============"

    # Count total core SNPs
    total_snps=$(tail -n +2 "${OUTPUT_DIR}/core.tab" | wc -l)
    echo "Total core SNPs: $total_snps"

    # SNP counts per sample
    echo ""
    echo "SNPs per sample:"
    for sample in "${SAMPLES[@]}"; do
        if [ -f "${OUTPUT_DIR}/${sample}/snps.txt" ]; then
            snp_count=$(tail -n +2 "${OUTPUT_DIR}/${sample}/snps.txt" | wc -l)
            echo "  $sample: $snp_count SNPs"
        fi
    done
fi

echo ""
echo "For visualization:"
echo "  - View alignment: less ${OUTPUT_DIR}/core.aln"
echo "  - View SNP table: less ${OUTPUT_DIR}/core.tab"
echo "  - View tree: cat ${OUTPUT_DIR}/core.tree.treefile"
echo ""

conda deactivate

echo "============================================"
