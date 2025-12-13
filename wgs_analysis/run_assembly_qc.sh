#!/usr/bin/bash

# --- run_assembly_qc.sh ---
# Comprehensive assembly quality assessment using QUAST
# Evaluates assembly completeness, contiguity, and accuracy

set -e

# --- CONDA SETUP ---
source ~/miniforge3/etc/profile.d/conda.sh

# --- USER-DEFINED VARIABLES ---

# Reference genome (optional but recommended)
REFERENCE_GENOME="/bulkpool/reference_data/Nissel1917/ncbi_dataset/data/GCF_000714595.1/GCF_000714595.1_ASM71459v1_genomic.fna"
REFERENCE_GFF="/bulkpool/reference_data/Nissel1917/ncbi_dataset/data/GCF_000714595.1/genomic.gff"

# Number of threads (default: 16)
THREADS=${1:-16}

# Output directory
OUTPUT_DIR="assembly_qc_results"

echo "============================================"
echo "   Assembly Quality Control with QUAST"
echo "============================================"
echo "Reference: E. coli Nissle 1917"
echo "Threads: $THREADS"
echo "Output: $OUTPUT_DIR"
echo "============================================"
echo ""

# Activate assembly-qc environment
conda activate assembly-qc

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Find all completed assemblies
echo "Searching for completed assemblies..."
ASSEMBLIES=$(find . -path "*/unicycler_output/assembly.fasta" -type f | sort)

if [ -z "$ASSEMBLIES" ]; then
    echo "ERROR: No assembly.fasta files found"
    echo "Expected pattern: */unicycler_output/assembly.fasta"
    exit 1
fi

echo "Found assemblies:"
echo "$ASSEMBLIES" | while read asm; do
    sample=$(echo "$asm" | sed 's|./\(.*\)_hybrid_assembly/.*|\1|')
    echo "  ✓ $sample"
done
echo ""

# Convert to comma-separated list for QUAST
ASSEMBLY_LIST=$(echo "$ASSEMBLIES" | tr '\n' ' ')

# Run QUAST with reference
if [ -f "$REFERENCE_GENOME" ]; then
    echo "Running QUAST with reference genome..."
    echo "This will provide:"
    echo "  - Assembly statistics (N50, L50, contigs)"
    echo "  - Reference-based metrics (genome fraction, misassemblies)"
    echo "  - Gene annotation completeness"
    echo ""

    quast.py $ASSEMBLY_LIST \
        -r "$REFERENCE_GENOME" \
        --features "$REFERENCE_GFF" \
        -o "$OUTPUT_DIR" \
        -t $THREADS \
        --min-contig 500 \
        --circos \
        --glimmer \
        --rna-finding
else
    echo "Running QUAST without reference (de novo mode)..."

    quast.py $ASSEMBLY_LIST \
        -o "$OUTPUT_DIR" \
        -t $THREADS \
        --min-contig 500 \
        --glimmer \
        --rna-finding
fi

echo ""
echo "============================================"
echo "   QUAST Analysis Complete!"
echo "============================================"
echo ""
echo "Results saved to: $OUTPUT_DIR/"
echo ""
echo "Key output files:"
echo "  - report.html          : Interactive HTML report"
echo "  - report.txt           : Text summary"
echo "  - report.tsv           : Tab-separated metrics"
echo "  - transposed_report.tsv: Transposed metrics (samples as columns)"
echo ""
echo "View the report:"
echo "  firefox $OUTPUT_DIR/report.html"
echo ""

# Generate summary statistics
echo "Assembly Summary:"
echo "================"
if [ -f "$OUTPUT_DIR/transposed_report.tsv" ]; then
    # Extract key metrics
    echo ""
    echo "Samples and key metrics:"
    awk -F'\t' 'NR==1 {for(i=2;i<=NF;i++) samples[i]=$i}
                /^# contigs \(>= 0 bp\)/ {printf "Total contigs: "; for(i=2;i<=NF;i++) printf "%s=%s ", samples[i], $i; print ""}
                /^Total length \(>= 0 bp\)/ {printf "Total length: "; for(i=2;i<=NF;i++) printf "%s=%s ", samples[i], $i; print ""}
                /^N50/ {printf "N50: "; for(i=2;i<=NF;i++) printf "%s=%s ", samples[i], $i; print ""}
                /^Genome fraction/ {printf "Genome fraction: "; for(i=2;i<=NF;i++) printf "%s=%s%% ", samples[i], $i; print ""}' \
        "$OUTPUT_DIR/transposed_report.tsv"
fi

echo ""
conda deactivate

echo "============================================"
