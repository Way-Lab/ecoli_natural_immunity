#!/usr/bin/bash

# --- run_parallel_assemblies.sh ---
# Parallel hybrid assembly for multiple samples
# Optimized for 128-core system

set -e

# Configuration
TOTAL_CORES=128
SAMPLES_TO_PROCESS=(2 3 4 5 6 7 8 9 10 11 12)  # Skip 1 (running), 13, 14 (bad quality)

# Calculate optimal parallelization
NUM_SAMPLES=${#SAMPLES_TO_PROCESS[@]}
THREADS_PER_SAMPLE=$((TOTAL_CORES / 4))  # 32 threads per sample, run 4 at a time

echo "============================================"
echo "   Parallel Hybrid Assembly Pipeline"
echo "============================================"
echo "Total cores: $TOTAL_CORES"
echo "Samples to process: ${SAMPLES_TO_PROCESS[@]}"
echo "Strategy: 4 samples at a time, $THREADS_PER_SAMPLE threads each"
echo "============================================"
echo ""

# Check for required files
echo "Checking for input files..."
VALID_SAMPLES=()
for i in "${SAMPLES_TO_PROCESS[@]}"; do
    SAMPLE="N2SKTQ_${i}_${i}"
    R1="${SAMPLE}_illumina_R1.fastq.gz"
    R2="${SAMPLE}_illumina_R2.fastq.gz"
    ONT="${SAMPLE}_nanopore.fastq.gz"

    if [ -f "$R1" ] && [ -f "$R2" ] && [ -f "$ONT" ]; then
        echo "  ✓ $SAMPLE - All files present"
        VALID_SAMPLES+=($i)
    else
        echo "  ✗ $SAMPLE - Missing files, skipping"
    fi
done

echo ""
echo "Valid samples to process: ${VALID_SAMPLES[@]}"
echo "Total: ${#VALID_SAMPLES[@]} samples"
echo ""

if [ ${#VALID_SAMPLES[@]} -eq 0 ]; then
    echo "ERROR: No valid samples found!"
    exit 1
fi

# Create log directory
mkdir -p batch_assembly_logs

# Function to process a single sample
process_sample() {
    local sample_num=$1
    local threads=$2
    local sample="N2SKTQ_${sample_num}_${sample_num}"
    local log_file="batch_assembly_logs/${sample}_batch.log"

    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Starting $sample with $threads threads" | tee -a "$log_file"

    if ./process_hybrid_assembly.sh "$sample" "$threads" >> "$log_file" 2>&1; then
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] ✓ SUCCESS: $sample" | tee -a "$log_file"
        return 0
    else
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] ✗ FAILED: $sample" | tee -a "$log_file"
        return 1
    fi
}

# Export function for parallel
export -f process_sample
export THREADS_PER_SAMPLE

# Start time
START_TIME=$(date +%s)

echo "Starting parallel assembly at $(date)"
echo "============================================"
echo ""

# Check if GNU parallel is available
if command -v parallel &> /dev/null; then
    echo "Using GNU parallel for job management..."
    printf '%s\n' "${VALID_SAMPLES[@]}" | \
        parallel -j 4 --eta --bar \
        "process_sample {} $THREADS_PER_SAMPLE"
else
    echo "GNU parallel not found, using background jobs..."

    # Process samples in batches of 4
    for sample_num in "${VALID_SAMPLES[@]}"; do
        # Wait if we have 4 jobs running
        while [ $(jobs -r | wc -l) -ge 4 ]; do
            sleep 10
        done

        # Start new job
        process_sample "$sample_num" "$THREADS_PER_SAMPLE" &
        echo "Started background job for N2SKTQ_${sample_num}_${sample_num}"
    done

    echo ""
    echo "All jobs started. Waiting for completion..."
    wait
fi

# Calculate runtime
END_TIME=$(date +%s)
RUNTIME=$((END_TIME - START_TIME))
RUNTIME_MIN=$((RUNTIME / 60))
RUNTIME_SEC=$((RUNTIME % 60))

echo ""
echo "============================================"
echo "   Batch Assembly Complete!"
echo "============================================"
echo "Runtime: ${RUNTIME_MIN}m ${RUNTIME_SEC}s"
echo "Completion time: $(date)"
echo ""

# Summary
echo "=== RESULTS SUMMARY ==="
SUCCESS_COUNT=0
FAILED_SAMPLES=()

for sample_num in "${VALID_SAMPLES[@]}"; do
    sample="N2SKTQ_${sample_num}_${sample_num}"
    if [ -f "${sample}_hybrid_assembly/${sample}_hybrid_assembly.fasta" ]; then
        SUCCESS_COUNT=$((SUCCESS_COUNT + 1))
        echo "  ✓ $sample"
    else
        FAILED_SAMPLES+=($sample)
        echo "  ✗ $sample"
    fi
done

echo ""
echo "Success: $SUCCESS_COUNT / ${#VALID_SAMPLES[@]} samples"

if [ ${#FAILED_SAMPLES[@]} -gt 0 ]; then
    echo ""
    echo "Failed samples:"
    for sample in "${FAILED_SAMPLES[@]}"; do
        echo "  - $sample"
        echo "    Log: batch_assembly_logs/${sample}_batch.log"
    done
fi

echo ""
echo "All results in: *_hybrid_assembly/"
echo "All logs in: batch_assembly_logs/"
echo ""
