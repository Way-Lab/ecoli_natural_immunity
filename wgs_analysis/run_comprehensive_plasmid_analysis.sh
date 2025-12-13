#!/bin/bash
#
# Comprehensive Plasmid Analysis Pipeline
# Combines MOB-suite + Bakta for detailed plasmid characterization
# across all 14 E. coli strain assemblies
#

set -u  # Exit on undefined variable
# NOTE: Removed set -e to prevent script from exiting on conda activation issues

# Configuration
WORKING_DIR="/fastpool/active_data/ecoli_genomics"
OUTPUT_DIR="${WORKING_DIR}/plasmid_analysis_results"
THREADS=32

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

log_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $(date '+%Y-%m-%d %H:%M:%S') - $1"
}

log_error() {
    echo -e "${RED}[ERROR]${NC} $(date '+%Y-%m-%d %H:%M:%S') - $1"
}

# Create output directory structure
create_output_dirs() {
    log_info "Creating output directory structure..."
    mkdir -p "${OUTPUT_DIR}"/{mobsuite_results,bakta_plasmid_annotations,summary_reports,plasmid_fastas}
    log_success "Output directories created"
}

# Run MOB-suite on a single assembly
run_mobsuite() {
    local assembly_file=$1
    local sample_name=$(basename "${assembly_file}" | sed 's/_hybrid_assembly\.fasta//g' | sed 's/\.fasta//g')
    local output_dir="${OUTPUT_DIR}/mobsuite_results/${sample_name}"

    log_info "Running MOB-suite on ${sample_name}..."

    # Activate MOB-suite environment and run (don't exit on activation failure)
    set +e
    source /home/david/miniforge3/bin/activate ${MOBSUITE_ENV} 2>/dev/null || conda activate ${MOBSUITE_ENV}
    set -e

    # Run mob_recon to identify and reconstruct plasmids
    mob_recon \
        --infile "${assembly_file}" \
        --outdir "${output_dir}" \
        --num_threads ${THREADS} \
        --force \
        --keep_tmp \
        2>&1 | tee "${output_dir}/${sample_name}_mob_recon.log"

    if [ $? -ne 0 ]; then
        log_error "mob_recon failed for ${sample_name}"
        conda deactivate || true
        return 1
    fi

    # Run mob_typer on identified plasmids
    for plasmid_fasta in "${output_dir}"/plasmid_*.fasta; do
        if [ -f "${plasmid_fasta}" ]; then
            local plasmid_name=$(basename "${plasmid_fasta}" .fasta)
            log_info "Running mob_typer on ${plasmid_name}..."

            mob_typer \
                --infile "${plasmid_fasta}" \
                --out_file "${output_dir}/${plasmid_name}_mobtyper.txt" \
                --num_threads ${THREADS} \
                2>&1 | tee -a "${output_dir}/${sample_name}_mob_typer.log"
        fi
    done

    log_success "MOB-suite completed for ${sample_name}"
    conda deactivate || true
}

# Annotate plasmid with Bakta
annotate_plasmid_with_bakta() {
    local plasmid_fasta=$1
    local sample_name=$2
    local plasmid_id=$3
    local output_dir="${OUTPUT_DIR}/bakta_plasmid_annotations/${sample_name}_${plasmid_id}"

    log_info "Annotating ${sample_name} plasmid ${plasmid_id} with Bakta..."

    # Create output directory
    mkdir -p "${output_dir}"

    # Activate Bakta environment (don't exit on activation failure)
    set +e
    source /home/david/miniforge3/bin/activate ${BAKTA_ENV} 2>/dev/null || conda activate ${BAKTA_ENV}
    set -e

    bakta \
        --db "${BAKTA_DB}" \
        --output "${output_dir}" \
        --prefix "${sample_name}_${plasmid_id}" \
        --threads ${THREADS} \
        --compliant \
        "${plasmid_fasta}" \
        2>&1 | tee "${output_dir}/${sample_name}_${plasmid_id}_bakta.log"

    if [ $? -eq 0 ]; then
        log_success "Bakta annotation completed for ${sample_name} plasmid ${plasmid_id}"
    else
        log_warning "Bakta annotation failed for ${sample_name} plasmid ${plasmid_id}"
    fi

    conda deactivate || true
}

# Process all assemblies
process_all_assemblies() {
    log_info "===== PHASE 1: Running MOB-suite on all assemblies ====="

    local total_assemblies=$((${#N2SKTQ_ASSEMBLIES[@]} + ${#Z6M7Y5_ASSEMBLIES[@]}))
    local current=0

    # Process N2SKTQ assemblies
    for assembly in "${N2SKTQ_ASSEMBLIES[@]}"; do
        current=$((current + 1))
        log_info "Processing ${current}/${total_assemblies}: $(basename ${assembly})"
        run_mobsuite "${assembly}"
    done

    # Process Z6M7Y5 assemblies
    for assembly in "${Z6M7Y5_ASSEMBLIES[@]}"; do
        current=$((current + 1))
        log_info "Processing ${current}/${total_assemblies}: $(basename ${assembly})"
        run_mobsuite "${assembly}"
    done

    log_success "MOB-suite analysis complete for all ${total_assemblies} assemblies"
}

# Annotate all identified plasmids
annotate_all_plasmids() {
    log_info "===== PHASE 2: Annotating all plasmids with Bakta ====="

    local plasmid_count=0

    # Find all plasmid FASTA files output by MOB-suite
    for mobsuite_dir in "${OUTPUT_DIR}"/mobsuite_results/*/; do
        sample_name=$(basename "${mobsuite_dir}")

        # Look for plasmid FASTAs (MOB-suite names them as plasmid_*.fasta)
        for plasmid_file in "${mobsuite_dir}"plasmid_*.fasta; do
            if [ -f "${plasmid_file}" ]; then
                plasmid_id=$(basename "${plasmid_file}" .fasta | sed 's/plasmid_//g')
                plasmid_count=$((plasmid_count + 1))

                # Copy to plasmid_fastas directory for easy access
                cp "${plasmid_file}" "${OUTPUT_DIR}/plasmid_fastas/${sample_name}_${plasmid_id}.fasta"

                # Annotate with Bakta
                annotate_plasmid_with_bakta "${plasmid_file}" "${sample_name}" "${plasmid_id}"
            fi
        done
    done

    log_success "Annotated ${plasmid_count} plasmids with Bakta"
}

# Generate summary reports
generate_summary_reports() {
    log_info "===== PHASE 3: Generating summary reports ====="

    local summary_file="${OUTPUT_DIR}/summary_reports/PLASMID_SUMMARY.txt"
    local csv_file="${OUTPUT_DIR}/summary_reports/plasmid_summary.csv"

    # Create summary header
    cat > "${summary_file}" << 'EOF'
================================================================================
                    PLASMID ANALYSIS SUMMARY REPORT
================================================================================
Analysis Date: $(date '+%Y-%m-%d %H:%M:%S')
Total Strains Analyzed: 18
Pipeline: MOB-suite v3.1.9 + Bakta

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

        # Check for contig report
        if [ -f "${mobsuite_dir}/contig_report.txt" ]; then
            grep -v "^contig_id" "${mobsuite_dir}/contig_report.txt" | while read line; do
                echo "  ${line}" >> "${summary_file}"

                # Parse for CSV
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
    log_info "Summary text file: ${summary_file}"
    log_info "Summary CSV file: ${csv_file}"
}

# Main execution
main() {
    log_info "===== Starting Comprehensive Plasmid Analysis Pipeline ====="
    log_info "Working directory: ${WORKING_DIR}"
    log_info "Output directory: ${OUTPUT_DIR}"
    log_info "Threads: ${THREADS}"

    create_output_dirs
    process_all_assemblies
    annotate_all_plasmids
    generate_summary_reports

    log_success "===== Plasmid Analysis Pipeline Complete! ====="
    log_info "Results are in: ${OUTPUT_DIR}"
    log_info ""
    log_info "Next steps:"
    log_info "  1. Review summary reports in: ${OUTPUT_DIR}/summary_reports/"
    log_info "  2. Examine individual plasmid annotations in: ${OUTPUT_DIR}/bakta_plasmid_annotations/"
    log_info "  3. Check plasmid typing results in: ${OUTPUT_DIR}/mobsuite_results/"
}

# Run the pipeline
main "$@"
