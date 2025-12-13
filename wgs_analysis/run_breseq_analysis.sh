#!/usr/bin/bash

# --- run_breseq_analysis.sh ---
# breseq analysis for bacterial resequencing
# Identifies mutations in evolved/variant strains compared to reference

set -e

# --- CONDA SETUP ---
source ~/miniforge3/etc/profile.d/conda.sh

# --- USER-DEFINED VARIABLES ---

# Sample prefix (e.g., N2SKTQ_1_1)
SAMPLE_PREFIX=$1

# Reference genome - using Nissle 1917 RefSeq
REFERENCE_GENOME="/bulkpool/reference_data/Nissel1917/ncbi_dataset/data/GCF_000714595.1/GCF_000714595.1_ASM71459v1_genomic.fna"
REFERENCE_GBK="/bulkpool/reference_data/Nissel1917/ncbi_dataset/data/GCF_000714595.1/genomic.gbff"

# Optional: Number of threads (default: 16)
THREADS=${2:-16}

# --- INPUT VALIDATION ---
if [ -z "$SAMPLE_PREFIX" ]; then
    echo "Usage: ./run_breseq_analysis.sh <sample_prefix> [threads]"
    echo ""
    echo "This script can analyze:"
    echo "  1. Illumina reads (paired-end)"
    echo "  2. Nanopore reads"
    echo "  3. Both (for consensus analysis)"
    echo ""
    echo "Example: ./run_breseq_analysis.sh N2SKTQ_1_1 16"
    exit 1
fi

# Check for reference
if [ ! -f "$REFERENCE_GENOME" ]; then
    echo "ERROR: Reference genome not found: $REFERENCE_GENOME"
    exit 1
fi

echo "============================================"
echo "   breseq Mutation Analysis"
echo "============================================"
echo "Sample: $SAMPLE_PREFIX"
echo "Reference: Nissle 1917 RefSeq (GCF_000714595.1)"
echo "Threads: $THREADS"
echo "============================================"
echo ""

# Activate breseq environment
conda activate breseq

# Check for input files
ILLUMINA_R1="${SAMPLE_PREFIX}_illumina_R1.fastq.gz"
ILLUMINA_R2="${SAMPLE_PREFIX}_illumina_R2.fastq.gz"
NANOPORE="${SAMPLE_PREFIX}_nanopore.fastq.gz"

HAS_ILLUMINA=false
HAS_NANOPORE=false

if [ -f "$ILLUMINA_R1" ] && [ -f "$ILLUMINA_R2" ]; then
    HAS_ILLUMINA=true
    echo "✓ Found Illumina reads"
fi

if [ -f "$NANOPORE" ]; then
    HAS_NANOPORE=true
    echo "✓ Found Nanopore reads"
fi

if [ "$HAS_ILLUMINA" = false ] && [ "$HAS_NANOPORE" = false ]; then
    echo "ERROR: No sequencing data found for $SAMPLE_PREFIX"
    echo "Expected:"
    echo "  - ${SAMPLE_PREFIX}_illumina_R1.fastq.gz & R2"
    echo "  - ${SAMPLE_PREFIX}_nanopore.fastq.gz"
    exit 1
fi

echo ""

# ============================================
# CONSENSUS ANALYSIS
# ============================================
# Use all available data for the most accurate mutation calls
OUTPUT_DIR="${SAMPLE_PREFIX}_breseq"

# Build read list based on what's available
READ_FILES=()

if [ "$HAS_ILLUMINA" = true ]; then
    READ_FILES+=("$ILLUMINA_R1" "$ILLUMINA_R2")
fi

if [ "$HAS_NANOPORE" = true ]; then
    READ_FILES+=("$NANOPORE")
fi

# Determine analysis type message
if [ "$HAS_ILLUMINA" = true ] && [ "$HAS_NANOPORE" = true ]; then
    echo "=== Running breseq consensus analysis (Illumina + Nanopore) ==="
    echo "-> Using both Illumina and Nanopore reads for most accurate results..."
    echo "   This may take 60-90 minutes..."
elif [ "$HAS_ILLUMINA" = true ]; then
    echo "=== Running breseq on Illumina reads ==="
    echo "-> Using Illumina paired-end reads..."
    echo "   This may take 30-60 minutes..."
elif [ "$HAS_NANOPORE" = true ]; then
    echo "=== Running breseq on Nanopore reads ==="
    echo "-> Using Nanopore reads..."
    echo "   This may take 30-60 minutes..."
fi

breseq -j $THREADS \
       -o "$OUTPUT_DIR" \
       -r "$REFERENCE_GBK" \
       "${READ_FILES[@]}"

echo "-> Analysis complete!"
echo "   HTML report: ${OUTPUT_DIR}/output/index.html"
echo "   Summary: ${OUTPUT_DIR}/output/summary.html"
echo ""

# ============================================
# SUMMARY
# ============================================
echo "============================================"
echo "   breseq Analysis Complete!"
echo "============================================"
echo ""
echo "Output directory:"
echo "  - ${SAMPLE_PREFIX}_breseq/"
echo ""

echo "Key output files:"
echo "  - output/index.html        : Main results page"
echo "  - output/summary.html      : Mutation summary"
echo "  - output/evidence/         : Detailed evidence for each variant"
echo "  - data/output.gd           : Genome diff file (for gdtools)"
echo "  - data/output.vcf          : Variants in VCF format"
echo ""

echo "To compare multiple samples:"
echo "  1. Use gdtools to compare .gd files"
echo "  2. Example: gdtools COMPARE -o comparison.html -r $REFERENCE_GBK sample1/data/output.gd sample2/data/output.gd"
echo ""

conda deactivate

echo "============================================"
