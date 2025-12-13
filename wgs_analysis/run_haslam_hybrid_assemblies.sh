#!/usr/bin/bash

# --- run_haslam_hybrid_assemblies.sh ---
# Run hybrid assemblies for all 6 Haslam samples in parallel
# Uses GNU parallel for efficient job management

set -e

# --- CONFIGURATION ---

# Number of parallel jobs (samples to process simultaneously)
PARALLEL_JOBS=2

# Threads per assembly job
THREADS_PER_JOB=32

# Samples to process
SAMPLES=(
    "Control-1"
    "Control-2"
    "Control-3"
    "Control-4"
    "NIL-1"
    "NIL-2"
)

# Log file for batch processing
BATCH_LOG="haslam_batch_assemblies.log"

echo "============================================"
echo "   Haslam Hybrid Assembly - Batch Mode"
echo "============================================"
echo "Samples: ${#SAMPLES[@]}"
echo "Parallel jobs: $PARALLEL_JOBS"
echo "Threads per job: $THREADS_PER_JOB"
echo "Log file: $BATCH_LOG"
echo "============================================"
echo ""

# Initialize log file
echo "Batch hybrid assembly started at $(date)" > "$BATCH_LOG"
echo "Processing ${#SAMPLES[@]} samples" >> "$BATCH_LOG"
echo "" >> "$BATCH_LOG"

# Function to check if sample is ready for assembly
check_sample_ready() {
    local sample=$1
    local r1="${sample}_illumina_R1.fastq.gz"
    local r2="${sample}_illumina_R2.fastq.gz"
    local nanopore="${sample}_nanopore.fastq.gz"

    if [ -f "$r1" ] && [ -f "$r2" ] && [ -f "$nanopore" ]; then
        return 0  # Ready
    else
        return 1  # Not ready
    fi
}

# Check which samples are ready
echo "Checking sample readiness..."
READY_SAMPLES=()
for sample in "${SAMPLES[@]}"; do
    if check_sample_ready "$sample"; then
        echo "  ✓ $sample - ready"
        READY_SAMPLES+=("$sample")
    else
        echo "  ✗ $sample - missing input files"
        echo "    Expected:"
        echo "      - ${sample}_illumina_R1.fastq.gz"
        echo "      - ${sample}_illumina_R2.fastq.gz"
        echo "      - ${sample}_nanopore.fastq.gz"
    fi
done

echo ""
echo "Ready to assemble: ${#READY_SAMPLES[@]} samples"
echo ""

if [ ${#READY_SAMPLES[@]} -eq 0 ]; then
    echo "ERROR: No samples ready for assembly"
    echo "Please run data preparation scripts first:"
    echo "  ./prepare_haslam_illumina.sh"
    echo "  ./prepare_haslam_nanopore.sh"
    exit 1
fi

# Check if GNU parallel is available
if ! command -v parallel &> /dev/null; then
    echo "WARNING: GNU parallel not found. Running sequentially..."
    echo ""

    # Sequential processing
    for sample in "${READY_SAMPLES[@]}"; do
        echo "=========================================="
        echo "Processing $sample"
        echo "=========================================="
        echo ""

        ./process_hybrid_assembly.sh "$sample" $THREADS_PER_JOB 2>&1 | tee -a "$BATCH_LOG"

        echo "" >> "$BATCH_LOG"
        echo "Completed $sample at $(date)" >> "$BATCH_LOG"
        echo "" >> "$BATCH_LOG"
    done
else
    # Parallel processing with GNU parallel
    echo "Using GNU parallel for job management..."
    echo ""

    # Export function and variables for parallel
    export -f check_sample_ready
    export THREADS_PER_JOB
    export BATCH_LOG

    # Run assemblies in parallel
    printf "%s\n" "${READY_SAMPLES[@]}" | \
        parallel -j $PARALLEL_JOBS --joblog parallel_assembly.log \
        "./process_hybrid_assembly.sh {} $THREADS_PER_JOB 2>&1 | tee -a haslam_{}_assembly.log"

    echo ""
    echo "All parallel jobs completed"
fi

echo ""
echo "============================================"
echo "   Batch Assembly Complete!"
echo "============================================"
echo "Completed at $(date)" | tee -a "$BATCH_LOG"
echo ""

# Summary of results
echo "Assembly Summary:"
echo "================"
for sample in "${READY_SAMPLES[@]}"; do
    assembly_dir="${sample}_hybrid_assembly"
    if [ -d "$assembly_dir" ]; then
        assembly_file="${assembly_dir}/${sample}_hybrid_assembly.fasta"
        if [ -f "$assembly_file" ]; then
            size=$(du -h "$assembly_file" | cut -f1)
            echo "  ✓ $sample - Assembly complete ($size)"
        else
            echo "  ⚠ $sample - Assembly directory exists but assembly file not found"
        fi
    else
        echo "  ✗ $sample - Assembly failed (directory not found)"
    fi
done

echo ""
echo "Individual logs available:"
for sample in "${READY_SAMPLES[@]}"; do
    log_file="${sample}_hybrid_assembly/${sample}_hybrid_processing.log"
    if [ -f "$log_file" ]; then
        echo "  - $log_file"
    fi
done

echo ""
echo "Next steps:"
echo "  1. Review assembly quality: ./run_assembly_qc.sh"
echo "  2. Run comparative analysis: ./run_haslam_comparative_analysis.sh"
echo ""
