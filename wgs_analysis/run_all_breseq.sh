#!/usr/bin/bash

# --- run_all_breseq.sh ---
# Run breseq analysis on all samples in order
# Runs samples 2-9, 11-12 first (completed assemblies)
# Then runs samples 1 and 10 last (still assembling)

set -e

# --- CONDA SETUP ---
source ~/miniforge3/etc/profile.d/conda.sh

# Number of threads per sample
THREADS=32

# Log directory
LOG_DIR="breseq_logs"
mkdir -p "$LOG_DIR"

echo "============================================"
echo "   Running breseq on all N2SKTQ samples"
echo "============================================"
echo "Threads per sample: $THREADS"
echo "Logs will be saved to: $LOG_DIR/"
echo "============================================"
echo ""

# Function to run breseq for a sample
run_breseq_sample() {
    local sample=$1
    local timestamp=$(date '+%Y-%m-%d %H:%M:%S')

    echo "[$timestamp] Starting breseq for $sample"

    # Run breseq and log output
    if ./run_breseq_analysis.sh "$sample" "$THREADS" > "${LOG_DIR}/${sample}_breseq.log" 2>&1; then
        timestamp=$(date '+%Y-%m-%d %H:%M:%S')
        echo "[$timestamp] ✓ COMPLETED: $sample"
    else
        timestamp=$(date '+%Y-%m-%d %H:%M:%S')
        echo "[$timestamp] ✗ FAILED: $sample (check ${LOG_DIR}/${sample}_breseq.log)"
    fi
    echo ""
}

# All samples (1-12)
SAMPLES=(
    "N2SKTQ_1_1"
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

echo "=== Running breseq on all 12 samples ==="
echo "Samples: ${SAMPLES[@]}"
echo ""

for sample in "${SAMPLES[@]}"; do
    run_breseq_sample "$sample"
done

echo "============================================"
echo "   All breseq analyses complete!"
echo "============================================"
echo ""
echo "Summary:"
echo "--------"
completed_count=$(ls -d N2SKTQ_*_breseq_* 2>/dev/null | wc -l)
echo "Total analyses completed: $completed_count"
echo ""
echo "Output directories:"
ls -d N2SKTQ_*_breseq_* 2>/dev/null || echo "  None found"
echo ""
echo "To view results:"
echo "  firefox N2SKTQ_X_X_breseq_illumina/output/index.html"
echo ""
