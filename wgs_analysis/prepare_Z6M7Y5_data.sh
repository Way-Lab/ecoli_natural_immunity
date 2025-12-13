#!/bin/bash
# Data preparation script for Z6M7Y5 samples
# No concatenation needed - data is already consolidated
# Creates symlinks with simplified names

set -e

echo "=========================================="
echo "Z6M7Y5 Data Preparation"
echo "=========================================="

# Create working directories
mkdir -p Z6M7Y5_illumina_prepared
mkdir -p Z6M7Y5_nanopore_prepared

cd Z6M7Y5_results

# Process each of the 6 samples
for i in {1..6}; do
    echo ""
    echo "Processing sample ${i}..."

    # Nanopore data - create symlink with simplified name
    nanopore_file="Z6M7Y5_${i}_sample_${i}__IS-*_nanopore.fastq.gz"
    if ls ${nanopore_file} 1> /dev/null 2>&1; then
        ln -sf "$(pwd)/$(ls ${nanopore_file})" ../Z6M7Y5_nanopore_prepared/Z6M7Y5_${i}_nanopore.fastq.gz
        echo "  ✓ Nanopore: Z6M7Y5_${i}_nanopore.fastq.gz"
    else
        echo "  ✗ ERROR: Nanopore file not found for sample ${i}"
        exit 1
    fi

    # Illumina R1 - always present
    r1_file="Z6M7Y5_${i}_sample_${i}__IS-*_illumina_R1.fastq.gz"
    if ls ${r1_file} 1> /dev/null 2>&1; then
        ln -sf "$(pwd)/$(ls ${r1_file})" ../Z6M7Y5_illumina_prepared/Z6M7Y5_${i}_R1.fastq.gz
        echo "  ✓ Illumina R1: Z6M7Y5_${i}_R1.fastq.gz"
    else
        echo "  ✗ ERROR: Illumina R1 file not found for sample ${i}"
        exit 1
    fi

    # Illumina R2 - special handling for sample 5 (only has R1)
    if [ ${i} -eq 5 ]; then
        echo "  ⚠ Sample 5: R2 not available (single-end data)"
    else
        r2_file="Z6M7Y5_${i}_sample_${i}__IS-*_illumina_R2.fastq.gz"
        if ls ${r2_file} 1> /dev/null 2>&1; then
            ln -sf "$(pwd)/$(ls ${r2_file})" ../Z6M7Y5_illumina_prepared/Z6M7Y5_${i}_R2.fastq.gz
            echo "  ✓ Illumina R2: Z6M7Y5_${i}_R2.fastq.gz"
        else
            echo "  ✗ ERROR: Illumina R2 file not found for sample ${i}"
            exit 1
        fi
    fi
done

cd ..

echo ""
echo "=========================================="
echo "Data Preparation Complete!"
echo "=========================================="
echo ""
echo "Summary:"
echo "  Illumina data → Z6M7Y5_illumina_prepared/"
echo "  Nanopore data → Z6M7Y5_nanopore_prepared/"
echo ""
echo "Note: Sample 5 only has R1 (single-end Illumina data)"
echo ""
echo "Next step: Run hybrid assemblies"
echo "  ./run_Z6M7Y5_hybrid_assemblies.sh"
