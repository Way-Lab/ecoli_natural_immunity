#!/usr/bin/bash

# --- prepare_haslam_nanopore.sh ---
# Concatenate multi-file Nanopore reads for Haslam samples
# Each sample has 52 individual fastq files that need to be combined

set -e

NANOPORE_DIR="HaslamNanoporeSeq"
OUTPUT_DIR="."

# Sample names (with Haslam- prefix for directory names)
SAMPLES=(
    "Haslam-Control-1"
    "Haslam-Control-2"
    "Haslam-Control-3"
    "Haslam-Control-4"
    "Haslam-NIL-1"
    "Haslam-NIL-2"
)

# Output names (without Haslam- prefix)
OUTPUT_NAMES=(
    "Control-1"
    "Control-2"
    "Control-3"
    "Control-4"
    "NIL-1"
    "NIL-2"
)

echo "============================================"
echo "  Concatenating Nanopore Reads (Multi-File)"
echo "============================================"
echo ""

for i in "${!SAMPLES[@]}"; do
    sample="${SAMPLES[$i]}"
    output_name="${OUTPUT_NAMES[$i]}"

    echo "Processing $sample..."

    # Output file
    NANOPORE_OUT="${OUTPUT_DIR}/${output_name}_nanopore.fastq.gz"

    # Check if already processed
    if [ -f "$NANOPORE_OUT" ]; then
        echo "  ✓ Already processed, skipping..."
        continue
    fi

    # Sample directory
    SAMPLE_DIR="${NANOPORE_DIR}/${sample}"

    if [ ! -d "$SAMPLE_DIR" ]; then
        echo "  ⚠ WARNING: Directory not found: $SAMPLE_DIR"
        continue
    fi

    # Find all fastq files (not in subdirectories)
    FASTQ_FILES=$(find "$SAMPLE_DIR" -maxdepth 1 -name "*.fastq" -o -maxdepth 1 -name "*.fastq.gz" | sort)

    if [ -z "$FASTQ_FILES" ]; then
        echo "  ⚠ WARNING: No fastq files found in $SAMPLE_DIR"
        continue
    fi

    # Count files
    FILE_COUNT=$(echo "$FASTQ_FILES" | wc -l)
    echo "  Found $FILE_COUNT fastq files"

    # Concatenate files
    # If files are already gzipped, cat them; if not, gzip while concatenating
    echo "  Concatenating and compressing..."

    # Create temporary combined file
    TMP_FILE="${OUTPUT_DIR}/${output_name}_nanopore_tmp.fastq"

    for file in $FASTQ_FILES; do
        if [[ $file == *.gz ]]; then
            gunzip -c "$file" >> "$TMP_FILE"
        else
            cat "$file" >> "$TMP_FILE"
        fi
    done

    # Compress the final file
    echo "  Compressing..."
    gzip -c "$TMP_FILE" > "$NANOPORE_OUT"
    rm "$TMP_FILE"

    # Get file size
    OUT_SIZE=$(du -h "$NANOPORE_OUT" | cut -f1)

    echo "  ✓ Complete!"
    echo "    Output: $NANOPORE_OUT ($OUT_SIZE)"
    echo ""
done

echo "============================================"
echo "  Nanopore Preparation Complete!"
echo "============================================"
echo ""
echo "Output files created:"
for output_name in "${OUTPUT_NAMES[@]}"; do
    NANOPORE_OUT="${OUTPUT_DIR}/${output_name}_nanopore.fastq.gz"
    if [ -f "$NANOPORE_OUT" ]; then
        echo "  - ${output_name}_nanopore.fastq.gz"
    fi
done
echo ""
