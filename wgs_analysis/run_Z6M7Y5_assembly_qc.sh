#!/bin/bash
# Assembly QC for Z6M7Y5 samples using QUAST and CheckM

set -e

# --- CONDA SETUP ---
source ~/miniforge3/etc/profile.d/conda.sh

# Configuration
THREADS=${1:-32}
ASSEMBLY_DIR="Z6M7Y5_assemblies"
QC_OUTPUT="Z6M7Y5_assembly_qc"

echo "=========================================="
echo "Z6M7Y5 Assembly Quality Control"
echo "=========================================="
echo "Using ${THREADS} threads"
echo ""

# Check if assemblies exist
if [ ! -d "${ASSEMBLY_DIR}" ]; then
    echo "ERROR: Assembly directory not found!"
    echo "Please run: ./run_Z6M7Y5_hybrid_assemblies.sh first"
    exit 1
fi

# Create QC output directory
mkdir -p ${QC_OUTPUT}

# Collect all assembly files
ASSEMBLIES=""
for i in {1..6}; do
    ASSEMBLY="${ASSEMBLY_DIR}/Z6M7Y5_${i}/Z6M7Y5_${i}.fasta"
    if [ -f "${ASSEMBLY}" ]; then
        ASSEMBLIES="${ASSEMBLIES} ${ASSEMBLY}"
    else
        echo "WARNING: Assembly not found: ${ASSEMBLY}"
    fi
done

if [ -z "${ASSEMBLIES}" ]; then
    echo "ERROR: No assemblies found!"
    exit 1
fi

echo "Found assemblies:"
echo "${ASSEMBLIES}" | tr ' ' '\n' | grep -v "^$"
echo ""

# Run QUAST for assembly statistics
echo "=========================================="
echo "Running QUAST..."
echo "=========================================="
conda activate assembly-qc
quast.py ${ASSEMBLIES} \
    -o ${QC_OUTPUT}/quast_results \
    -t ${THREADS} \
    --min-contig 500 \
    --glimmer \
    --circos

echo ""
echo "✓ QUAST complete"
conda deactivate
echo ""

# Run BUSCO for completeness assessment (Enterobacterales lineage)
echo "=========================================="
echo "Running BUSCO..."
echo "=========================================="

# Check if BUSCO is available
if command -v busco &> /dev/null; then
    mkdir -p ${QC_OUTPUT}/busco_results

    for i in {1..6}; do
        ASSEMBLY="${ASSEMBLY_DIR}/Z6M7Y5_${i}/Z6M7Y5_${i}.fasta"
        if [ ! -f "${ASSEMBLY}" ]; then
            continue
        fi

        echo ""
        echo "Processing Z6M7Y5_${i}..."

        busco -i ${ASSEMBLY} \
            -o Z6M7Y5_${i}_busco \
            -m genome \
            -l enterobacterales_odb10 \
            --cpu ${THREADS} \
            --out_path ${QC_OUTPUT}/busco_results \
            -f

        echo "✓ BUSCO complete for Z6M7Y5_${i}"
    done
else
    echo "WARNING: BUSCO not available, skipping completeness assessment"
fi

echo ""
echo "=========================================="
echo "Quality Control Complete!"
echo "=========================================="
echo ""
echo "Results directory: ${QC_OUTPUT}/"
echo ""
echo "Key output files:"
echo "  - ${QC_OUTPUT}/quast_results/report.html"
echo "  - ${QC_OUTPUT}/quast_results/report.tsv"
echo "  - ${QC_OUTPUT}/busco_results/*/short_summary.*.txt"
echo ""

# Generate summary report
echo "=========================================="
echo "Assembly Quality Summary"
echo "=========================================="
echo ""

if [ -f "${QC_OUTPUT}/quast_results/report.tsv" ]; then
    echo "QUAST Statistics:"
    echo "----------------"
    # Extract key metrics from QUAST report
    head -30 "${QC_OUTPUT}/quast_results/report.tsv" | awk 'BEGIN{FS="\t"; OFS="\t"} {print}' | sed 's/\t/  |  /g'
    echo ""
fi

# BUSCO summary
echo "BUSCO Completeness:"
echo "-------------------"
for i in {1..6}; do
    BUSCO_SUMMARY="${QC_OUTPUT}/busco_results/Z6M7Y5_${i}_busco/short_summary.*.txt"
    if ls ${BUSCO_SUMMARY} 1> /dev/null 2>&1; then
        SUMMARY_FILE=$(ls ${BUSCO_SUMMARY})
        COMPLETE=$(grep "Complete BUSCOs" ${SUMMARY_FILE} | head -1 | awk '{print $1}')
        echo "  Z6M7Y5_${i}: ${COMPLETE}% complete BUSCOs"
    fi
done

echo ""
echo "Next step: Comparative analysis with reference"
echo "  ./run_Z6M7Y5_comparative_analysis.sh ${THREADS}"
