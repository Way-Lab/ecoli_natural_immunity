#!/bin/bash
#
# Run Bakta annotation on all identified plasmids in parallel
#

set -u

# Configuration
WORKING_DIR="/fastpool/active_data/ecoli_genomics"
OUTPUT_DIR="${WORKING_DIR}/plasmid_analysis_results"
THREADS=32
JOBS=6

# Bakta configuration
BAKTA_ENV="bakta"
BAKTA_DB="/bulkpool/reference_data/bakta_db/db"

echo "===== Running Bakta annotation on all plasmids ====="
echo "Threads: ${THREADS}"
echo "Parallel jobs: ${JOBS}"
echo ""

# Function to annotate a single plasmid
annotate_plasmid() {
    local plasmid_fasta=$1
    local sample_dir=$(dirname "${plasmid_fasta}")
    local sample_name=$(basename "${sample_dir}")
    local plasmid_id=$(basename "${plasmid_fasta}" .fasta | sed 's/plasmid_//g')
    local output_dir="${OUTPUT_DIR}/bakta_plasmid_annotations/${sample_name}_${plasmid_id}"

    echo "[$(date '+%H:%M:%S')] Annotating ${sample_name} plasmid ${plasmid_id}"

    # Activate Bakta environment
    source /home/david/miniforge3/bin/activate ${BAKTA_ENV} 2>/dev/null || true

    bakta \
        --db "${BAKTA_DB}" \
        --output "${output_dir}" \
        --prefix "${sample_name}_${plasmid_id}" \
        --threads $((THREADS / JOBS)) \
        --compliant \
        --force \
        "${plasmid_fasta}" \
        &> "${output_dir}.log"

    conda deactivate 2>/dev/null || true

    echo "[$(date '+%H:%M:%S')] Completed ${sample_name} plasmid ${plasmid_id}"
}

export -f annotate_plasmid
export OUTPUT_DIR BAKTA_ENV BAKTA_DB THREADS JOBS

# Find all plasmid FASTA files and annotate in parallel
find "${OUTPUT_DIR}/mobsuite_results" -name "plasmid_*.fasta" | \
    parallel -j ${JOBS} --progress annotate_plasmid {}

echo ""
echo "===== Bakta annotation complete! ====="
echo "Results in: ${OUTPUT_DIR}/bakta_plasmid_annotations/"
