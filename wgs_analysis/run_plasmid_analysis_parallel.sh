#!/bin/bash
#
# Parallel Plasmid Analysis Pipeline
# Uses GNU parallel to analyze all 18 samples simultaneously
#

set -u

# Configuration
WORKING_DIR="/fastpool/active_data/ecoli_genomics"
OUTPUT_DIR="${WORKING_DIR}/plasmid_analysis_results"
THREADS=32
JOBS=6  # Number of parallel jobs (adjust based on available resources)

# Conda environments
MOBSUITE_ENV="mobsuite"
BAKTA_ENV="bakta"

# Bakta database location
BAKTA_DB="/bulkpool/reference_data/bakta_db/db"

# Assembly locations
N2SKTQ_ASSEMBLIES=(
    "${WORKING_DIR}/N2SKTQ_results/N2SKTQ_1_1_hybrid_assembly/N2SKTQ_1_1_hybrid_assembly.fasta"
    "${WORKING_DIR}/N2SKTQ_results/N2SKTQ_2_2_hybrid_assembly/N2SKTQ_2_2_hybrid_assembly.fasta"
    "${WORKING_DIR}/N2SKTQ_results/N2SKTQ_3_3_hybrid_assembly/N2SKTQ_3_3_hybrid_assembly.fasta"
    "${WORKING_DIR}/N2SKTQ_results/N2SKTQ_4_4_hybrid_assembly/N2SKTQ_4_4_hybrid_assembly.fasta"
    "${WORKING_DIR}/N2SKTQ_results/N2SKTQ_5_5_hybrid_assembly/N2SKTQ_5_5_hybrid_assembly.fasta"
    "${WORKING_DIR}/N2SKTQ_results/N2SKTQ_6_6_hybrid_assembly/N2SKTQ_6_6_hybrid_assembly.fasta"
    "${WORKING_DIR}/N2SKTQ_results/N2SKTQ_7_7_hybrid_assembly/N2SKTQ_7_7_hybrid_assembly.fasta"
    "${WORKING_DIR}/N2SKTQ_results/N2SKTQ_8_8_hybrid_assembly/N2SKTQ_8_8_hybrid_assembly.fasta"
    "${WORKING_DIR}/N2SKTQ_results/N2SKTQ_9_9_hybrid_assembly/N2SKTQ_9_9_hybrid_assembly.fasta"
    "${WORKING_DIR}/N2SKTQ_results/N2SKTQ_10_10_hybrid_assembly/N2SKTQ_10_10_hybrid_assembly.fasta"
    "${WORKING_DIR}/N2SKTQ_results/N2SKTQ_11_11_hybrid_assembly/N2SKTQ_11_11_hybrid_assembly.fasta"
    "${WORKING_DIR}/N2SKTQ_results/N2SKTQ_12_12_hybrid_assembly/N2SKTQ_12_12_hybrid_assembly.fasta"
)

Z6M7Y5_ASSEMBLIES=(
    "${WORKING_DIR}/Z6M7Y5_assemblies/Z6M7Y5_1/Z6M7Y5_1.fasta"
    "${WORKING_DIR}/Z6M7Y5_assemblies/Z6M7Y5_2/Z6M7Y5_2.fasta"
    "${WORKING_DIR}/Z6M7Y5_assemblies/Z6M7Y5_3/Z6M7Y5_3.fasta"
    "${WORKING_DIR}/Z6M7Y5_assemblies/Z6M7Y5_4/Z6M7Y5_4.fasta"
    "${WORKING_DIR}/Z6M7Y5_assemblies/Z6M7Y5_5/Z6M7Y5_5.fasta"
    "${WORKING_DIR}/Z6M7Y5_assemblies/Z6M7Y5_6/Z6M7Y5_6.fasta"
)

# Color output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Logging functions
log_info() {
    echo -e "${BLUE}[INFO]${NC} $(date '+%Y-%m-%d %H:%M:%S') - $1"
}

log_success() {
    echo -e "${GREEN}[SUCCESS]${NC} $(date '+%Y-%m-%d %H:%M:%S') - $1"
}

log_error() {
    echo -e "${RED}[ERROR]${NC} $(date '+%Y-%m-%d %H:%M:%S') - $1"
}

# Create output directory structure
mkdir -p "${OUTPUT_DIR}"/{mobsuite_results,bakta_plasmid_annotations,summary_reports,plasmid_fastas}

# Function to process a single sample with MOB-suite
process_sample_mobsuite() {
    local assembly_file=$1
    local sample_name=$(basename "${assembly_file}" | sed 's/_hybrid_assembly\.fasta//g' | sed 's/\.fasta//g')
    local output_dir="${OUTPUT_DIR}/mobsuite_results/${sample_name}"

    echo "[$(date '+%H:%M:%S')] Starting MOB-suite for ${sample_name}"

    # Activate MOB-suite environment
    source /home/david/miniforge3/bin/activate ${MOBSUITE_ENV} 2>/dev/null || true

    # Run mob_recon
    mob_recon \
        --infile "${assembly_file}" \
        --outdir "${output_dir}" \
        --num_threads $((THREADS / JOBS)) \
        --force \
        --keep_tmp \
        &> "${output_dir}_mob_recon.log"

    # Run mob_typer on identified plasmids
    for plasmid_fasta in "${output_dir}"/plasmid_*.fasta; do
        if [ -f "${plasmid_fasta}" ]; then
            local plasmid_name=$(basename "${plasmid_fasta}" .fasta)
            mob_typer \
                --infile "${plasmid_fasta}" \
                --out_file "${output_dir}/${plasmid_name}_mobtyper.txt" \
                --num_threads $((THREADS / JOBS)) \
                &>> "${output_dir}_mob_typer.log"
        fi
    done

    conda deactivate 2>/dev/null || true

    echo "[$(date '+%H:%M:%S')] Completed MOB-suite for ${sample_name}"
}

# Export function and variables for parallel
export -f process_sample_mobsuite
export OUTPUT_DIR MOBSUITE_ENV THREADS JOBS

log_info "===== Starting Parallel Plasmid Analysis Pipeline ====="
log_info "Working directory: ${WORKING_DIR}"
log_info "Output directory: ${OUTPUT_DIR}"
log_info "Total threads: ${THREADS}"
log_info "Parallel jobs: ${JOBS}"
log_info "Threads per job: $((THREADS / JOBS))"

# Combine all assemblies
ALL_ASSEMBLIES=("${N2SKTQ_ASSEMBLIES[@]}" "${Z6M7Y5_ASSEMBLIES[@]}")

log_info "===== PHASE 1: Running MOB-suite on all ${#ALL_ASSEMBLIES[@]} assemblies in parallel ====="

# Run MOB-suite in parallel
printf '%s\n' "${ALL_ASSEMBLIES[@]}" | parallel -j ${JOBS} --progress process_sample_mobsuite {}

log_success "MOB-suite analysis complete for all samples"

log_info "===== PHASE 2: Annotating all plasmids with Bakta ====="

# Function to annotate a single plasmid
annotate_plasmid() {
    local plasmid_fasta=$1
    local sample_dir=$(dirname "${plasmid_fasta}")
    local sample_name=$(basename "${sample_dir}")
    local plasmid_id=$(basename "${plasmid_fasta}" .fasta | sed 's/plasmid_//g')
    local output_dir="${OUTPUT_DIR}/bakta_plasmid_annotations/${sample_name}_${plasmid_id}"

    echo "[$(date '+%H:%M:%S')] Annotating ${sample_name} plasmid ${plasmid_id}"

    mkdir -p "${output_dir}"

    # Activate Bakta environment
    source /home/david/miniforge3/bin/activate ${BAKTA_ENV} 2>/dev/null || true

    bakta \
        --db "${BAKTA_DB}" \
        --output "${output_dir}" \
        --prefix "${sample_name}_${plasmid_id}" \
        --threads $((THREADS / JOBS)) \
        --compliant \
        "${plasmid_fasta}" \
        &> "${output_dir}/bakta.log"

    # Copy plasmid to plasmid_fastas directory
    cp "${plasmid_fasta}" "${OUTPUT_DIR}/plasmid_fastas/${sample_name}_${plasmid_id}.fasta"

    conda deactivate 2>/dev/null || true

    echo "[$(date '+%H:%M:%S')] Completed ${sample_name} plasmid ${plasmid_id}"
}

export -f annotate_plasmid
export BAKTA_ENV BAKTA_DB

# Find all plasmid FASTA files and annotate in parallel
find "${OUTPUT_DIR}/mobsuite_results" -name "plasmid_*.fasta" | \
    parallel -j ${JOBS} --progress annotate_plasmid {}

log_success "Bakta annotation complete for all plasmids"

log_info "===== PHASE 3: Generating summary reports ====="

# Generate summary reports
summary_file="${OUTPUT_DIR}/summary_reports/PLASMID_SUMMARY.txt"
csv_file="${OUTPUT_DIR}/summary_reports/plasmid_summary.csv"

cat > "${summary_file}" << EOF
================================================================================
                    PLASMID ANALYSIS SUMMARY REPORT
================================================================================
Analysis Date: $(date '+%Y-%m-%d %H:%M:%S')
Total Strains Analyzed: 18
Pipeline: MOB-suite v3.1.9 + Bakta
Processing: Parallel (${JOBS} concurrent jobs)

================================================================================
                           PLASMID INVENTORY
================================================================================

EOF

# CSV header
echo "Strain,Plasmid_ID,Size_bp,Circularity,Predicted_Mobility,Replicon_Type,Primary_Cluster_ID" > "${csv_file}"

# Parse MOB-suite results
for mobsuite_dir in "${OUTPUT_DIR}"/mobsuite_results/*/; do
    sample_name=$(basename "${mobsuite_dir}")

    echo "Strain: ${sample_name}" >> "${summary_file}"
    echo "----------------------------------------" >> "${summary_file}"

    if [ -f "${mobsuite_dir}/contig_report.txt" ]; then
        grep -v "^contig_id" "${mobsuite_dir}/contig_report.txt" | while read line; do
            echo "  ${line}" >> "${summary_file}"

            # Parse for CSV (simplified - adjust field numbers as needed)
            contig_id=$(echo "$line" | awk '{print $1}')
            size=$(echo "$line" | awk '{print $2}')
            circular=$(echo "$line" | awk '{print $3}')
            mob_type=$(echo "$line" | awk '{print $4}')
            rep_type=$(echo "$line" | awk '{print $5}')
            cluster=$(echo "$line" | awk '{print $6}')

            echo "${sample_name},${contig_id},${size},${circular},${mob_type},${rep_type},${cluster}" >> "${csv_file}"
        done
    else
        echo "  No plasmids detected" >> "${summary_file}"
    fi

    echo "" >> "${summary_file}"
done

log_success "Summary reports generated"

log_success "===== Plasmid Analysis Pipeline Complete! ====="
log_info "Results are in: ${OUTPUT_DIR}"
log_info ""
log_info "Summary files:"
log_info "  - ${summary_file}"
log_info "  - ${csv_file}"
log_info ""
log_info "Plasmid count:"
find "${OUTPUT_DIR}/plasmid_fastas" -name "*.fasta" | wc -l | xargs echo "  Total plasmids identified:"
