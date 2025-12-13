#!/bin/bash
# Complete workflow for Z6M7Y5 hybrid assembly and comparative analysis
# This script runs the entire pipeline in a screen session to avoid losing work due to disconnection
# Compares 6 new Z6M7Y5 isolates with 12 N2SKTQ strains

set -e

THREADS=${1:-32}
SCREEN_NAME="Z6M7Y5_workflow"

echo "=========================================="
echo "Z6M7Y5 Complete Workflow"
echo "=========================================="
echo ""
echo "This workflow will:"
echo "  1. Prepare raw sequencing data (symlinks)"
echo "  2. Run hybrid assemblies (Unicycler) - ~6-12 hours"
echo "  3. Quality control (QUAST, BUSCO) - ~30 min"
echo "  4. Comparative analysis (Parsnp with N2SKTQ strains) - ~1 hour"
echo ""
echo "Total estimated time: 8-14 hours"
echo "Using ${THREADS} threads per assembly"
echo ""
echo "All work will run in a screen session: ${SCREEN_NAME}"
echo ""
read -p "Press Enter to start, or Ctrl+C to cancel..."
echo ""

# Check if screen is available
if ! command -v screen &> /dev/null; then
    echo "ERROR: screen is not installed!"
    echo "Please install: sudo apt-get install screen"
    exit 1
fi

# Check if screen session already exists
if screen -list | grep -q "${SCREEN_NAME}"; then
    echo "WARNING: Screen session '${SCREEN_NAME}' already exists!"
    echo ""
    echo "Options:"
    echo "  1. Resume existing session: screen -r ${SCREEN_NAME}"
    echo "  2. Kill existing session: screen -X -S ${SCREEN_NAME} quit"
    echo "  3. Use a different name"
    exit 1
fi

# Create a script that will run inside screen
SCREEN_SCRIPT="/tmp/Z6M7Y5_workflow_$(date +%Y%m%d_%H%M%S).sh"

cat > ${SCREEN_SCRIPT} << 'WORKFLOW_EOF'
#!/bin/bash
set -e

THREADS=THREADS_PLACEHOLDER
WORKFLOW_LOG="Z6M7Y5_workflow_$(date +%Y%m%d_%H%M%S).log"

# Function to log with timestamp
log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" | tee -a ${WORKFLOW_LOG}
}

log "=========================================="
log "Z6M7Y5 Workflow Started"
log "=========================================="
log "Log file: ${WORKFLOW_LOG}"
log "Threads: ${THREADS}"
log ""

# Step 1: Data preparation
log "=========================================="
log "STEP 1: Data Preparation"
log "=========================================="
if [ ! -d "Z6M7Y5_illumina_prepared" ] || [ ! -d "Z6M7Y5_nanopore_prepared" ]; then
    log "Running: ./prepare_Z6M7Y5_data.sh"
    ./prepare_Z6M7Y5_data.sh 2>&1 | tee -a ${WORKFLOW_LOG}
    log "✓ Data preparation complete"
else
    log "⚠ Data already prepared, skipping..."
fi
log ""

# Step 2: Hybrid assemblies
log "=========================================="
log "STEP 2: Hybrid Assemblies (Unicycler)"
log "=========================================="
log "This will take 6-12 hours..."
log "Start time: $(date)"

if [ ! -d "Z6M7Y5_assemblies" ]; then
    log "Running: ./run_Z6M7Y5_hybrid_assemblies.sh ${THREADS}"
    ./run_Z6M7Y5_hybrid_assemblies.sh ${THREADS} 2>&1 | tee -a ${WORKFLOW_LOG}
    log "✓ Hybrid assemblies complete"
else
    log "⚠ Assemblies directory exists, checking for completeness..."
    COMPLETE=true
    for i in {1..6}; do
        if [ ! -f "Z6M7Y5_assemblies/Z6M7Y5_${i}/Z6M7Y5_${i}.fasta" ]; then
            COMPLETE=false
            break
        fi
    done

    if [ "$COMPLETE" = true ]; then
        log "✓ All assemblies already complete, skipping..."
    else
        log "⚠ Some assemblies missing, re-running..."
        ./run_Z6M7Y5_hybrid_assemblies.sh ${THREADS} 2>&1 | tee -a ${WORKFLOW_LOG}
    fi
fi

log "End time: $(date)"
log ""

# Step 3: Assembly QC
log "=========================================="
log "STEP 3: Assembly Quality Control"
log "=========================================="
log "Running QUAST and BUSCO..."

if [ ! -d "Z6M7Y5_assembly_qc" ]; then
    log "Running: ./run_Z6M7Y5_assembly_qc.sh ${THREADS}"
    ./run_Z6M7Y5_assembly_qc.sh ${THREADS} 2>&1 | tee -a ${WORKFLOW_LOG}
    log "✓ Quality control complete"
else
    log "⚠ QC directory exists, skipping... (delete Z6M7Y5_assembly_qc to re-run)"
fi
log ""

# Step 4: Comparative analysis
log "=========================================="
log "STEP 4: Comparative Analysis (Parsnp)"
log "=========================================="
log "Comparing Z6M7Y5 isolates with N2SKTQ strains..."
log "Reference: N2SKTQ_1_1"
log "Genomes: 6 Z6M7Y5 + 12 N2SKTQ = 18 total"

if [ ! -d "Z6M7Y5_comparative_analysis" ]; then
    log "Running: ./run_Z6M7Y5_comparative_analysis.sh ${THREADS}"
    ./run_Z6M7Y5_comparative_analysis.sh ${THREADS} 2>&1 | tee -a ${WORKFLOW_LOG}
    log "✓ Comparative analysis complete"
else
    log "⚠ Comparative analysis directory exists, skipping... (delete Z6M7Y5_comparative_analysis to re-run)"
fi
log ""

# Summary
log "=========================================="
log "WORKFLOW COMPLETE!"
log "=========================================="
log "End time: $(date)"
log ""
log "Output directories:"
log "  - Z6M7Y5_assemblies/ - Hybrid assembly results"
log "  - Z6M7Y5_assembly_qc/ - QUAST and BUSCO results"
log "  - Z6M7Y5_comparative_analysis/ - Parsnp phylogenetic analysis"
log ""
log "Key files:"
log "  - Z6M7Y5_assembly_qc/quast_results/report.html"
log "  - Z6M7Y5_comparative_analysis/parsnp_results/parsnp.tree"
log "  - Z6M7Y5_comparative_analysis/parsnp_results/snp_distance_matrix.tsv"
log ""
log "Full log saved to: ${WORKFLOW_LOG}"
log ""
log "Press Ctrl+A, D to detach from screen session"
log "To reattach: screen -r Z6M7Y5_workflow"

# Keep screen session open
echo ""
echo "Workflow complete! Press Enter to close this screen session..."
read
WORKFLOW_EOF

# Replace placeholder with actual thread count
sed -i "s/THREADS_PLACEHOLDER/${THREADS}/g" ${SCREEN_SCRIPT}
chmod +x ${SCREEN_SCRIPT}

# Start screen session
echo "=========================================="
echo "Starting Screen Session"
echo "=========================================="
echo ""
echo "Screen session name: ${SCREEN_NAME}"
echo "Workflow script: ${SCREEN_SCRIPT}"
echo ""
echo "To detach from screen: Press Ctrl+A, then D"
echo "To reattach: screen -r ${SCREEN_NAME}"
echo "To view all sessions: screen -ls"
echo ""
echo "Starting in 3 seconds..."
sleep 3

# Execute workflow in screen
screen -dmS ${SCREEN_NAME} bash -c "cd $(pwd) && ${SCREEN_SCRIPT}; exec bash"

echo ""
echo "✓ Screen session started!"
echo ""
echo "To attach to the session and monitor progress:"
echo "  screen -r ${SCREEN_NAME}"
echo ""
echo "The workflow will run in the background even if you disconnect."
echo ""
