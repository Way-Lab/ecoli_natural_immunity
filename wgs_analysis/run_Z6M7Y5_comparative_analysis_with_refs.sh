#!/bin/bash
# Comparative genomic analysis of Z6M7Y5 samples with N2SKTQ strains using Parsnp
# UPDATED: Now includes reference E. coli genomes for broader phylogenetic context
# Uses N2SKTQ_1_1 as the reference genome for alignment
# Includes all 12 N2SKTQ strains + 6 Z6M7Y5 isolates + reference genomes

set -e

# --- CONDA SETUP ---
source ~/miniforge3/etc/profile.d/conda.sh

# Configuration
THREADS=${1:-32}
Z6M7Y5_ASSEMBLY_DIR="Z6M7Y5_assemblies"
N2SKTQ_ASSEMBLY_DIR="N2SKTQ_results"
REFERENCE="N2SKTQ_results/N2SKTQ_1_1_hybrid_assembly/N2SKTQ_1_1_hybrid_assembly.fasta"
REFERENCE_GENOMES_DIR="HaslamNanoporeSeq/phylogenetic_analysis/reference_genomes"
OUTPUT_DIR="Z6M7Y5_comparative_analysis_with_refs"

echo "=========================================="
echo "Z6M7Y5 + N2SKTQ + Reference Genomes"
echo "Comparative Analysis"
echo "=========================================="
echo "Tool: Parsnp (core genome alignment)"
echo "Reference: N2SKTQ_1_1"
echo "Using ${THREADS} threads"
echo ""

# Check if reference exists
if [ ! -f "${REFERENCE}" ]; then
    echo "ERROR: Reference genome not found!"
    echo "Expected: ${REFERENCE}"
    exit 1
fi

# Check if Z6M7Y5 assemblies exist
if [ ! -d "${Z6M7Y5_ASSEMBLY_DIR}" ]; then
    echo "ERROR: Z6M7Y5 assembly directory not found!"
    echo "Please run: ./run_Z6M7Y5_hybrid_assemblies.sh first"
    exit 1
fi

# Create temporary directory for Parsnp input
PARSNP_INPUT="${OUTPUT_DIR}/genomes"
mkdir -p ${PARSNP_INPUT}

echo "Preparing genomes for Parsnp..."
echo ""

# Reference genome - will be specified with -r flag only (NOT copied to genomes directory)
echo "Reference:"
echo "  ✓ N2SKTQ_1_1 (reference - will be used with -r flag)"
echo ""

# Copy all N2SKTQ assemblies (2-12, since 1_1 is the reference)
echo "N2SKTQ strains (2-12):"
N2SKTQ_COUNT=0
for i in {2..12}; do
    ASSEMBLY="${N2SKTQ_ASSEMBLY_DIR}/N2SKTQ_${i}_${i}_hybrid_assembly/N2SKTQ_${i}_${i}_hybrid_assembly.fasta"
    if [ -f "${ASSEMBLY}" ]; then
        cp ${ASSEMBLY} ${PARSNP_INPUT}/N2SKTQ_${i}_${i}.fasta
        echo "  ✓ N2SKTQ_${i}_${i}"
        N2SKTQ_COUNT=$((N2SKTQ_COUNT + 1))
    else
        echo "  ✗ WARNING: N2SKTQ_${i}_${i} assembly not found"
    fi
done

echo ""
echo "Z6M7Y5 strains (1-6):"
# Copy Z6M7Y5 assemblies
Z6M7Y5_COUNT=0
for i in {1..6}; do
    ASSEMBLY="${Z6M7Y5_ASSEMBLY_DIR}/Z6M7Y5_${i}/Z6M7Y5_${i}.fasta"
    if [ -f "${ASSEMBLY}" ]; then
        cp ${ASSEMBLY} ${PARSNP_INPUT}/Z6M7Y5_${i}.fasta
        echo "  ✓ Z6M7Y5_${i}"
        Z6M7Y5_COUNT=$((Z6M7Y5_COUNT + 1))
    else
        echo "  ✗ WARNING: Z6M7Y5_${i} assembly not found"
    fi
done

echo ""
echo "Reference E. coli genomes:"
# Copy reference genomes
REF_COUNT=0
declare -a REFERENCE_STRAINS=(
    "ecoli_K12_MG1655"
    "ecoli_CFT073"
    "ecoli_UTI89"
    "ecoli_EC958"
    "ecoli_MS7163"
    "ecoli_Nissle1917"
)

for ref_strain in "${REFERENCE_STRAINS[@]}"; do
    REF_FILE="${REFERENCE_GENOMES_DIR}/${ref_strain}.fna"
    if [ -f "${REF_FILE}" ]; then
        cp ${REF_FILE} ${PARSNP_INPUT}/${ref_strain}.fasta
        echo "  ✓ ${ref_strain}"
        REF_COUNT=$((REF_COUNT + 1))
    else
        echo "  ✗ WARNING: ${ref_strain} not found at ${REF_FILE}"
    fi
done

TOTAL_GENOMES=$((N2SKTQ_COUNT + Z6M7Y5_COUNT + REF_COUNT + 1))
echo ""
echo "Total genomes for analysis: ${TOTAL_GENOMES}"
echo "  - N2SKTQ_1_1 (reference): 1"
echo "  - N2SKTQ strains (2-12): ${N2SKTQ_COUNT}"
echo "  - Z6M7Y5 strains (1-6): ${Z6M7Y5_COUNT}"
echo "  - Reference genomes: ${REF_COUNT}"
echo ""

# Run Parsnp
echo "=========================================="
echo "Running Parsnp..."
echo "=========================================="
echo "Start time: $(date)"
echo ""

conda activate assembly-qc
parsnp \
    -r ${REFERENCE} \
    -d ${PARSNP_INPUT} \
    -o ${OUTPUT_DIR}/parsnp_results \
    -p ${THREADS} \
    -c -v

echo ""
echo "End time: $(date)"
echo "✓ Parsnp complete"
conda deactivate
echo ""

# Generate SNP distance matrix
echo "=========================================="
echo "Generating SNP distance matrix..."
echo "=========================================="

cd ${OUTPUT_DIR}/parsnp_results

# Use harvesttools to convert parsnp output to various formats
if command -v harvesttools &> /dev/null; then
    # Extract core SNP alignment
    harvesttools -i parsnp.ggr -V parsnp_core.vcf
    echo "  ✓ Generated: parsnp_core.vcf"

    # Generate alignment in fasta format
    harvesttools -i parsnp.ggr -M parsnp_core.fasta
    echo "  ✓ Generated: parsnp_core.fasta"
fi

cd ../..

# Create SNP distance matrix from VCF
if [ -f "${OUTPUT_DIR}/parsnp_results/parsnp.vcf" ]; then
    echo ""
    echo "Creating SNP distance matrix..."

    # Python script to generate SNP matrix from VCF
    python3 << 'EOF'
import sys
import os
from collections import defaultdict

vcf_file = "Z6M7Y5_comparative_analysis_with_refs/parsnp_results/parsnp.vcf"
output_file = "Z6M7Y5_comparative_analysis_with_refs/parsnp_results/snp_distance_matrix.tsv"

if not os.path.exists(vcf_file):
    print(f"ERROR: VCF file not found: {vcf_file}")
    sys.exit(1)

# Parse VCF to count SNPs between samples
samples = []
snp_data = []

with open(vcf_file, 'r') as f:
    for line in f:
        if line.startswith('#CHROM'):
            samples = line.strip().split('\t')[9:]
            break

# Parse genotype data properly
with open(vcf_file, 'r') as f:
    for line in f:
        if line.startswith('#'):
            continue
        fields = line.strip().split('\t')
        if len(fields) > 9:
            # Extract genotypes - typically in format "GT:..." where GT is "0" or "1"
            genotypes = []
            for gt_field in fields[9:]:
                # Split by ':' and take first field (GT)
                gt = gt_field.split(':')[0]
                genotypes.append(gt)
            snp_data.append(genotypes)

# Calculate pairwise SNP distances
n_samples = len(samples)
snp_matrix = [[0] * n_samples for _ in range(n_samples)]

for i in range(n_samples):
    for j in range(i+1, n_samples):
        diff_count = 0
        for snp in snp_data:
            if i < len(snp) and j < len(snp):
                # Compare actual genotypes
                if snp[i] != snp[j]:
                    diff_count += 1
        snp_matrix[i][j] = diff_count
        snp_matrix[j][i] = diff_count

# Write matrix to file
with open(output_file, 'w') as f:
    # Header
    f.write("Sample\t" + "\t".join(samples) + "\n")
    # Data rows
    for i, sample in enumerate(samples):
        f.write(sample + "\t")
        f.write("\t".join(str(snp_matrix[i][j]) for j in range(n_samples)))
        f.write("\n")

print(f"  ✓ SNP distance matrix: {output_file}")
print(f"\nSample count: {len(samples)}")
print(f"Total SNPs analyzed: {len(snp_data)}")
EOF

fi

echo ""
echo "=========================================="
echo "Comparative Analysis Complete!"
echo "=========================================="
echo ""
echo "Output directory: ${OUTPUT_DIR}/"
echo ""
echo "Key output files:"
echo "  - ${OUTPUT_DIR}/parsnp_results/parsnp.tree"
echo "  - ${OUTPUT_DIR}/parsnp_results/parsnp.vcf"
echo "  - ${OUTPUT_DIR}/parsnp_results/parsnp_core.fasta"
echo "  - ${OUTPUT_DIR}/parsnp_results/snp_distance_matrix.tsv"
echo ""
echo "To visualize the tree, you can use:"
echo "  - FigTree: Load parsnp.tree"
echo "  - Online: https://itol.embl.de/ (upload parsnp.tree)"
echo ""

# Display tree if available
if [ -f "${OUTPUT_DIR}/parsnp_results/parsnp.tree" ]; then
    echo "Phylogenetic tree (Newick format):"
    echo "-----------------------------------"
    cat ${OUTPUT_DIR}/parsnp_results/parsnp.tree
    echo ""
fi

# Display SNP matrix if available
if [ -f "${OUTPUT_DIR}/parsnp_results/snp_distance_matrix.tsv" ]; then
    echo ""
    echo "SNP Distance Matrix (first 10x10):"
    echo "------------------------------------"
    head -11 ${OUTPUT_DIR}/parsnp_results/snp_distance_matrix.tsv | column -t -s $'\t'
    echo ""
    echo "(Full matrix available in: ${OUTPUT_DIR}/parsnp_results/snp_distance_matrix.tsv)"
fi

echo ""
echo "=========================================="
echo "Analysis pipeline complete!"
echo "=========================================="
