#!/usr/bin/bash

# --- prepare_haslam_illumina.sh ---
# Concatenate multi-lane Illumina reads for Haslam samples
# Each sample has reads split across 6 lanes (L003-L008)

set -e

ILLUMINA_DIR="HaslamIlluminaSeq"
OUTPUT_DIR="."

# Sample names
SAMPLES=(
    "Control-1"
    "Control-2"
    "Control-3"
    "Control-4"
    "NIL-1"
    "NIL-2"
)

echo "============================================"
echo "  Concatenating Illumina Reads (Multi-Lane)"
echo "============================================"
echo ""

for sample in "${SAMPLES[@]}"; do
    echo "Processing $sample..."

    # Output files
    R1_OUT="${OUTPUT_DIR}/${sample}_illumina_R1.fastq.gz"
    R2_OUT="${OUTPUT_DIR}/${sample}_illumina_R2.fastq.gz"

    # Check if already processed
    if [ -f "$R1_OUT" ] && [ -f "$R2_OUT" ]; then
        echo "  ✓ Already processed, skipping..."
        continue
    fi

    # Find all R1 files for this sample (sorted by lane)
    R1_FILES=$(ls ${ILLUMINA_DIR}/${sample}_*_R1_*.fastq.gz 2>/dev/null | sort)
    R2_FILES=$(ls ${ILLUMINA_DIR}/${sample}_*_R2_*.fastq.gz 2>/dev/null | sort)

    if [ -z "$R1_FILES" ]; then
        echo "  ⚠ WARNING: No Illumina files found for $sample"
        continue
    fi

    # Count files
    R1_COUNT=$(echo "$R1_FILES" | wc -l)
    R2_COUNT=$(echo "$R2_FILES" | wc -l)

    echo "  Found $R1_COUNT R1 files and $R2_COUNT R2 files"

    # Concatenate R1 files
    echo "  Concatenating R1 files..."
    cat $R1_FILES > "$R1_OUT"

    # Concatenate R2 files
    echo "  Concatenating R2 files..."
    cat $R2_FILES > "$R2_OUT"

    # Get file sizes
    R1_SIZE=$(du -h "$R1_OUT" | cut -f1)
    R2_SIZE=$(du -h "$R2_OUT" | cut -f1)

    echo "  ✓ Complete!"
    echo "    R1: $R1_OUT ($R1_SIZE)"
    echo "    R2: $R2_OUT ($R2_SIZE)"
    echo ""
done

echo "============================================"
echo "  Illumina Preparation Complete!"
echo "============================================"
echo ""
echo "Output files created:"
for sample in "${SAMPLES[@]}"; do
    R1_OUT="${OUTPUT_DIR}/${sample}_illumina_R1.fastq.gz"
    R2_OUT="${OUTPUT_DIR}/${sample}_illumina_R2.fastq.gz"
    if [ -f "$R1_OUT" ] && [ -f "$R2_OUT" ]; then
        echo "  - ${sample}_illumina_R1.fastq.gz"
        echo "  - ${sample}_illumina_R2.fastq.gz"
    fi
done
echo ""
