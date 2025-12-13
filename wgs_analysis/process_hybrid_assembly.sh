#!/usr/bin/bash

# --- process_hybrid_assembly.sh ---
# Hybrid assembly pipeline using Unicycler for bacterial genomes
# Combines Illumina short reads and Nanopore long reads
# Includes assembly, polishing, QC, annotation, and comparison to reference

set -e # Exit immediately if a command exits with a non-zero status

# --- CONDA SETUP ---
source ~/miniforge3/etc/profile.d/conda.sh

# --- USER-DEFINED VARIABLES ---

# Sample prefix (e.g., N2SKTQ_1_1)
SAMPLE_PREFIX=$1

# Optional: Number of threads (default: 16)
THREADS=${2:-16}

# Reference genome paths
REFERENCE_GENOME="/bulkpool/reference_data/Nissel1917/ncbi_dataset/data/GCF_000714595.1/GCF_000714595.1_ASM71459v1_genomic.fna"
REFERENCE_GBK="/bulkpool/reference_data/Nissel1917/ncbi_dataset/data/GCF_000714595.1/genomic.gbff"
REFERENCE_GFF="/bulkpool/reference_data/Nissel1917/ncbi_dataset/data/GCF_000714595.1/genomic.gff"

# --- INPUT VALIDATION ---
if [ -z "$SAMPLE_PREFIX" ]; then
    echo "Usage: ./process_hybrid_assembly.sh <sample_prefix> [threads]"
    echo "Example: ./process_hybrid_assembly.sh N2SKTQ_1_1 16"
    exit 1
fi

# Check for input files
ILLUMINA_R1="${SAMPLE_PREFIX}_illumina_R1.fastq.gz"
ILLUMINA_R2="${SAMPLE_PREFIX}_illumina_R2.fastq.gz"
NANOPORE="${SAMPLE_PREFIX}_nanopore.fastq.gz"

if [ ! -f "$ILLUMINA_R1" ]; then
    echo "ERROR: Illumina R1 file not found: $ILLUMINA_R1"
    exit 1
fi

if [ ! -f "$ILLUMINA_R2" ]; then
    echo "ERROR: Illumina R2 file not found: $ILLUMINA_R2"
    exit 1
fi

if [ ! -f "$NANOPORE" ]; then
    echo "ERROR: Nanopore file not found: $NANOPORE"
    exit 1
fi

if [ ! -f "$REFERENCE_GENOME" ]; then
    echo "ERROR: Reference genome not found: $REFERENCE_GENOME"
    exit 1
fi

# --- SETUP OUTPUT DIRECTORY ---
OUTPUT_DIR="${SAMPLE_PREFIX}_hybrid_assembly"
mkdir -p "$OUTPUT_DIR"
LOG_FILE="${OUTPUT_DIR}/${SAMPLE_PREFIX}_hybrid_processing.log"

echo "============================================"
echo "   Hybrid Assembly Pipeline - Unicycler"
echo "============================================"
echo "Sample: $SAMPLE_PREFIX"
echo "Threads: $THREADS"
echo "Output: $OUTPUT_DIR"
echo "Reference: E. coli Nissle 1917"
echo "============================================"
echo ""
echo "Processing started at $(date)" > "$LOG_FILE"

# Get file sizes
R1_SIZE=$(du -h "$ILLUMINA_R1" | cut -f1)
R2_SIZE=$(du -h "$ILLUMINA_R2" | cut -f1)
ONT_SIZE=$(du -h "$NANOPORE" | cut -f1)

echo "Input files:"
echo "  Illumina R1: $ILLUMINA_R1 ($R1_SIZE)"
echo "  Illumina R2: $ILLUMINA_R2 ($R2_SIZE)"
echo "  Nanopore: $NANOPORE ($ONT_SIZE)"
echo ""

# ============================================
# STEP 1: QC on Raw Reads
# ============================================
echo "=== STEP 1: Quality Control on Raw Reads ==="

# Activate nanopore-qc for QC tools
if conda env list | grep -q "nanopore-qc"; then
    echo "-> Activating nanopore-qc environment..."
    conda activate nanopore-qc

    # QC on Nanopore reads
    echo "-> Running NanoPlot on Nanopore reads..."
    NanoPlot --fastq "$NANOPORE" \
             --outdir "${OUTPUT_DIR}/qc_nanopore_raw" \
             --threads $THREADS \
             --plots kde hex \
             --N50 \
             --title "Raw Nanopore - $SAMPLE_PREFIX" \
             2>&1 | tee -a "$LOG_FILE"

    echo "-> QC complete. Check ${OUTPUT_DIR}/qc_nanopore_raw/ for reports"

    conda deactivate
else
    echo "WARNING: nanopore-qc environment not found. Skipping Nanopore QC."
fi

echo "Step 1 completed at $(date)" >> "$LOG_FILE"
echo ""

# ============================================
# STEP 2: Hybrid Assembly with Unicycler
# ============================================
echo "=== STEP 2: Hybrid Assembly with Unicycler ==="
echo "-> Activating unicycler-env..."
conda activate unicycler-env

ASSEMBLY_DIR="${OUTPUT_DIR}/unicycler_output"

echo "-> Running Unicycler hybrid assembly..."
echo "   This may take 30-60 minutes depending on data size..."

unicycler -1 "$ILLUMINA_R1" \
          -2 "$ILLUMINA_R2" \
          -l "$NANOPORE" \
          -o "$ASSEMBLY_DIR" \
          -t $THREADS \
          --verbosity 2 \
          2>&1 | tee -a "$LOG_FILE"

# Check if assembly succeeded
if [ ! -f "${ASSEMBLY_DIR}/assembly.fasta" ]; then
    echo "ERROR: Unicycler assembly failed. Check log file: $LOG_FILE"
    exit 1
fi

# Copy main assembly file
ASSEMBLY_FASTA="${OUTPUT_DIR}/${SAMPLE_PREFIX}_hybrid_assembly.fasta"
cp "${ASSEMBLY_DIR}/assembly.fasta" "$ASSEMBLY_FASTA"

echo "-> Assembly complete: $ASSEMBLY_FASTA"
echo "Step 2 completed at $(date)" >> "$LOG_FILE"
echo ""

# ============================================
# STEP 3: Assembly Statistics
# ============================================
echo "=== STEP 3: Assembly Statistics ==="
echo "-> Calculating assembly statistics..."

# Basic stats using seqkit or custom script
python3 << EOF
import gzip

def parse_fasta(filename):
    sequences = {}
    current_id = None
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                current_id = line[1:].split()[0]
                sequences[current_id] = ''
            else:
                sequences[current_id] += line
    return sequences

seqs = parse_fasta("${ASSEMBLY_FASTA}")
contig_lengths = sorted([len(seq) for seq in seqs.values()], reverse=True)

total_length = sum(contig_lengths)
num_contigs = len(contig_lengths)
longest = contig_lengths[0] if contig_lengths else 0

# Calculate N50
cumsum = 0
n50 = 0
for length in contig_lengths:
    cumsum += length
    if cumsum >= total_length / 2:
        n50 = length
        break

print(f"Number of contigs: {num_contigs}")
print(f"Total assembly length: {total_length:,} bp")
print(f"Longest contig: {longest:,} bp")
print(f"N50: {n50:,} bp")
print(f"Mean contig length: {int(total_length/num_contigs):,} bp" if num_contigs > 0 else "N/A")

# Save to file
with open("${OUTPUT_DIR}/${SAMPLE_PREFIX}_assembly_stats.txt", 'w') as out:
    out.write(f"Assembly Statistics for {SAMPLE_PREFIX}\n")
    out.write("="*50 + "\n")
    out.write(f"Number of contigs: {num_contigs}\n")
    out.write(f"Total assembly length: {total_length:,} bp\n")
    out.write(f"Longest contig: {longest:,} bp\n")
    out.write(f"N50: {n50:,} bp\n")
    out.write(f"Mean contig length: {int(total_length/num_contigs):,} bp\n" if num_contigs > 0 else "Mean contig length: N/A\n")
    out.write(f"\nContig lengths:\n")
    for i, length in enumerate(contig_lengths, 1):
        out.write(f"  Contig {i}: {length:,} bp\n")
EOF

cat "${OUTPUT_DIR}/${SAMPLE_PREFIX}_assembly_stats.txt"

echo "Step 3 completed at $(date)" >> "$LOG_FILE"
echo ""

# ============================================
# STEP 4: Compare Assembly to Reference
# ============================================
echo "=== STEP 4: Compare Assembly to Reference ==="
echo "-> Switching to nanopore-wgs environment for alignment..."
conda deactivate
conda activate nanopore-wgs

# Align assembly to reference using minimap2
ALIGNED_PAF="${OUTPUT_DIR}/${SAMPLE_PREFIX}_vs_reference.paf"
ALIGNED_SAM="${OUTPUT_DIR}/${SAMPLE_PREFIX}_vs_reference.sam"

echo "-> Aligning assembly to reference genome..."
minimap2 -ax asm5 \
         -t $THREADS \
         "$REFERENCE_GENOME" \
         "$ASSEMBLY_FASTA" \
         > "$ALIGNED_SAM" 2>&1 | tee -a "$LOG_FILE"

# Convert to sorted BAM
ALIGNED_BAM="${OUTPUT_DIR}/${SAMPLE_PREFIX}_vs_reference.sorted.bam"
samtools view -bS -@ $THREADS "$ALIGNED_SAM" | \
    samtools sort -@ $THREADS -o "$ALIGNED_BAM"
samtools index -@ $THREADS "$ALIGNED_BAM"

# Generate alignment statistics
echo "-> Generating alignment statistics..."
samtools flagstat "$ALIGNED_BAM" > "${OUTPUT_DIR}/${SAMPLE_PREFIX}_alignment_stats.txt"
samtools coverage "$ALIGNED_BAM" > "${OUTPUT_DIR}/${SAMPLE_PREFIX}_reference_coverage.txt"

echo "-> Alignment statistics:"
cat "${OUTPUT_DIR}/${SAMPLE_PREFIX}_alignment_stats.txt"
echo ""
echo "-> Reference coverage:"
cat "${OUTPUT_DIR}/${SAMPLE_PREFIX}_reference_coverage.txt"

# Use MUMmer for detailed comparison
if command -v nucmer &> /dev/null; then
    echo "-> Running MUMmer alignment for detailed comparison..."
    MUMMER_PREFIX="${OUTPUT_DIR}/${SAMPLE_PREFIX}_mummer"

    nucmer --prefix="$MUMMER_PREFIX" \
           -t $THREADS \
           "$REFERENCE_GENOME" \
           "$ASSEMBLY_FASTA" \
           2>&1 | tee -a "$LOG_FILE"

    # Generate alignment statistics
    dnadiff -d "${MUMMER_PREFIX}.delta" -p "$MUMMER_PREFIX" 2>&1 | tee -a "$LOG_FILE"

    echo "-> MUMmer report generated: ${MUMMER_PREFIX}.report"

    if [ -f "${MUMMER_PREFIX}.report" ]; then
        echo ""
        echo "=== MUMmer Alignment Summary ==="
        grep -A 20 "TotalBases" "${MUMMER_PREFIX}.report" || echo "Summary not found in expected format"
    fi
fi

echo "Step 4 completed at $(date)" >> "$LOG_FILE"
echo ""

# ============================================
# STEP 5: Call Variants Against Reference
# ============================================
echo "=== STEP 5: Variant Calling Against Reference ==="
echo "-> Calling variants using bcftools..."

# Call variants on the assembly vs reference alignment
RAW_VCF="${OUTPUT_DIR}/${SAMPLE_PREFIX}_variants.raw.vcf.gz"

bcftools mpileup -f "$REFERENCE_GENOME" \
                 "$ALIGNED_BAM" \
                 -Ou \
                 --threads $THREADS | \
    bcftools call -mv -Oz -o "$RAW_VCF" --threads $THREADS

tabix -p vcf "$RAW_VCF"

# Filter variants
FILTERED_VCF="${OUTPUT_DIR}/${SAMPLE_PREFIX}_variants.filtered.vcf.gz"
bcftools view -i 'QUAL>20 && DP>5' "$RAW_VCF" -Oz -o "$FILTERED_VCF"
tabix -p vcf "$FILTERED_VCF"

# Count variants
RAW_COUNT=$(bcftools view -H "$RAW_VCF" 2>/dev/null | wc -l)
FILTERED_COUNT=$(bcftools view -H "$FILTERED_VCF" 2>/dev/null | wc -l)

echo "-> Variants detected:"
echo "   Raw variants: $RAW_COUNT"
echo "   Filtered variants (QUAL>20, DP>5): $FILTERED_COUNT"

# Generate variant statistics
bcftools stats "$FILTERED_VCF" > "${OUTPUT_DIR}/${SAMPLE_PREFIX}_variant_stats.txt"

echo "Step 5 completed at $(date)" >> "$LOG_FILE"
echo ""

# ============================================
# STEP 6: Annotation (if Bakta is available)
# ============================================
echo "=== STEP 6: Genome Annotation ==="

# Check if bakta is available in any environment
if command -v bakta &> /dev/null || conda run -n base bakta --version &> /dev/null; then
    echo "-> Running Bakta annotation..."

    ANNOTATION_DIR="${OUTPUT_DIR}/annotation"

    bakta --outdir "$ANNOTATION_DIR" \
           --prefix "$SAMPLE_PREFIX" \
           --genus Escherichia \
           --species "coli" \
           --strain "Nissle_1917_variant" \
           --cpus $THREADS \
           --force \
           "$ASSEMBLY_FASTA" \
           2>&1 | tee -a "$LOG_FILE"

    if [ -f "${ANNOTATION_DIR}/${SAMPLE_PREFIX}.gff" ]; then
        echo "-> Annotation complete!"
        echo "   GFF: ${ANNOTATION_DIR}/${SAMPLE_PREFIX}.gff"
        echo "   GBK: ${ANNOTATION_DIR}/${SAMPLE_PREFIX}.gbk"

        # Count features
        GENES=$(grep -c "gene" "${ANNOTATION_DIR}/${SAMPLE_PREFIX}.gff" || echo 0)
        CDS=$(grep -c "CDS" "${ANNOTATION_DIR}/${SAMPLE_PREFIX}.gff" || echo 0)
        echo "   Predicted genes: $GENES"
        echo "   Predicted CDS: $CDS"
    fi
else
    echo "WARNING: Bakta not available. Skipping annotation."
    echo "To enable annotation, install Bakta:"
    echo "  conda install -c bioconda bakta"
fi

echo "Step 6 completed at $(date)" >> "$LOG_FILE"
echo ""

# ============================================
# FINAL SUMMARY
# ============================================
echo "============================================"
echo "       Processing Complete!"
echo "============================================"
echo "Processing completed at $(date)" >> "$LOG_FILE"

SUMMARY_FILE="${OUTPUT_DIR}/${SAMPLE_PREFIX}_SUMMARY.txt"

cat > "$SUMMARY_FILE" << SUMMARY_EOF
=================================================
  HYBRID ASSEMBLY SUMMARY - ${SAMPLE_PREFIX}
=================================================
Date: $(date)
Reference: E. coli Nissle 1917

INPUT DATA:
-----------
Illumina R1: $ILLUMINA_R1 ($R1_SIZE)
Illumina R2: $ILLUMINA_R2 ($R2_SIZE)
Nanopore:    $NANOPORE ($ONT_SIZE)

ASSEMBLY STATISTICS:
--------------------
SUMMARY_EOF

cat "${OUTPUT_DIR}/${SAMPLE_PREFIX}_assembly_stats.txt" >> "$SUMMARY_FILE"

cat >> "$SUMMARY_FILE" << SUMMARY_EOF

ALIGNMENT TO REFERENCE:
-----------------------
SUMMARY_EOF

cat "${OUTPUT_DIR}/${SAMPLE_PREFIX}_alignment_stats.txt" >> "$SUMMARY_FILE"

cat >> "$SUMMARY_FILE" << SUMMARY_EOF

VARIANT CALLING:
----------------
Raw variants: $RAW_COUNT
Filtered variants: $FILTERED_COUNT

OUTPUT FILES:
-------------
Assembly:          $ASSEMBLY_FASTA
Assembly graph:    ${ASSEMBLY_DIR}/assembly.gfa
Unicycler log:     ${ASSEMBLY_DIR}/unicycler.log
Alignment (BAM):   $ALIGNED_BAM
Variants (VCF):    $FILTERED_VCF
Assembly stats:    ${OUTPUT_DIR}/${SAMPLE_PREFIX}_assembly_stats.txt
QC reports:        ${OUTPUT_DIR}/qc_*/
Summary:           $SUMMARY_FILE
Full log:          $LOG_FILE

=================================================
SUMMARY_EOF

echo ""
cat "$SUMMARY_FILE"
echo ""
echo "Full summary saved to: $SUMMARY_FILE"
echo "Full log saved to: $LOG_FILE"
echo ""
echo "To view assembly graph, open in Bandage:"
echo "  ${ASSEMBLY_DIR}/assembly.gfa"
echo ""

conda deactivate

echo "============================================"
echo "          Pipeline Complete!"
echo "============================================"
