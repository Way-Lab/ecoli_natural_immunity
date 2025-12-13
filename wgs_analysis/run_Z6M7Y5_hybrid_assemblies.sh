#!/bin/bash
# Hybrid assembly script for Z6M7Y5 samples using Unicycler
# All samples now have paired-end Illumina data (R1 + R2)
# Creates assemblies with simplified names: Z6M7Y5_1, Z6M7Y5_2, etc.

set -e

# --- CONDA SETUP ---
source ~/miniforge3/etc/profile.d/conda.sh

# Configuration
THREADS=${1:-32}  # Default to 32 threads, or use command line argument
ILLUMINA_DIR="Z6M7Y5_illumina_prepared"
NANOPORE_DIR="Z6M7Y5_nanopore_prepared"
OUTPUT_BASE="Z6M7Y5_assemblies"

echo "=========================================="
echo "Z6M7Y5 Hybrid Assembly Pipeline"
echo "=========================================="
echo "Using ${THREADS} threads"
echo "Assembler: Unicycler (hybrid mode)"
echo ""

# Check if data is prepared
if [ ! -d "${ILLUMINA_DIR}" ] || [ ! -d "${NANOPORE_DIR}" ]; then
    echo "ERROR: Data not prepared!"
    echo "Please run: ./prepare_Z6M7Y5_data.sh first"
    exit 1
fi

# Create output directory
mkdir -p ${OUTPUT_BASE}

# Function to run assembly for a single sample
run_assembly() {
    # Source conda in the function for parallel execution
    source ~/miniforge3/etc/profile.d/conda.sh

    local SAMPLE=$1
    local R1="${ILLUMINA_DIR}/Z6M7Y5_${SAMPLE}_R1.fastq.gz"
    local R2="${ILLUMINA_DIR}/Z6M7Y5_${SAMPLE}_R2.fastq.gz"
    local NANOPORE="${NANOPORE_DIR}/Z6M7Y5_${SAMPLE}_nanopore.fastq.gz"
    local OUT_DIR="${OUTPUT_BASE}/Z6M7Y5_${SAMPLE}"

    echo ""
    echo "=========================================="
    echo "Processing: Z6M7Y5_${SAMPLE}"
    echo "=========================================="
    echo "Start time: $(date)"

    # Activate unicycler environment
    conda activate unicycler-env

    # Check if R2 exists (sample 5 is single-end)
    if [ -f "${R2}" ]; then
        echo "Mode: Paired-end Illumina + Nanopore"
        unicycler -t ${THREADS} \
                  -1 ${R1} \
                  -2 ${R2} \
                  -l ${NANOPORE} \
                  -o ${OUT_DIR}
    else
        echo "Mode: Single-end Illumina + Nanopore (sample 5)"
        unicycler -t ${THREADS} \
                  -s ${R1} \
                  -l ${NANOPORE} \
                  -o ${OUT_DIR}
    fi

    local EXIT_CODE=$?

    conda deactivate

    # Create simplified output filename
    if [ ${EXIT_CODE} -eq 0 ] && [ -f "${OUT_DIR}/assembly.fasta" ]; then
        cp "${OUT_DIR}/assembly.fasta" "${OUT_DIR}/Z6M7Y5_${SAMPLE}.fasta"
        echo ""
        echo "✓ Assembly complete: ${OUT_DIR}/Z6M7Y5_${SAMPLE}.fasta"
    else
        echo ""
        echo "✗ ERROR: Assembly failed for Z6M7Y5_${SAMPLE}"
        return 1
    fi

    echo "End time: $(date)"
}

# Export function for parallel execution
export -f run_assembly
export THREADS ILLUMINA_DIR NANOPORE_DIR OUTPUT_BASE

# Run assemblies in parallel (all 6 samples)
echo "Starting parallel assemblies for 6 samples..."
echo "Each assembly will use ${THREADS} threads"
echo ""

# Use GNU parallel to run 1 assembly at a time (since each uses many threads)
# Adjust -j parameter based on your system resources
parallel -j 1 run_assembly ::: {1..6}

echo ""
echo "=========================================="
echo "All Assemblies Complete!"
echo "=========================================="
echo ""
echo "Output directory: ${OUTPUT_BASE}/"
echo ""
echo "Assembly files:"
for i in {1..6}; do
    ASSEMBLY="${OUTPUT_BASE}/Z6M7Y5_${i}/Z6M7Y5_${i}.fasta"
    if [ -f "${ASSEMBLY}" ]; then
        SIZE=$(du -h "${ASSEMBLY}" | cut -f1)
        CONTIGS=$(grep -c "^>" "${ASSEMBLY}")
        echo "  ✓ Z6M7Y5_${i}.fasta - ${SIZE}, ${CONTIGS} contigs"
    else
        echo "  ✗ Z6M7Y5_${i}.fasta - FAILED"
    fi
done

echo ""
echo "Next steps:"
echo "  1. Quality check: ./run_Z6M7Y5_assembly_qc.sh ${THREADS}"
echo "  2. Comparative analysis: ./run_Z6M7Y5_comparative_analysis.sh ${THREADS}"
