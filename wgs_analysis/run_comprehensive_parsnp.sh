#!/bin/bash

# Activate conda environment
source /home/david/miniforge3/bin/activate assembly-qc

# Run parsnp with Nissle Stock-1 as reference
echo "Starting parsnp analysis at $(date)"
echo "Reference: N2SKTQ_1_1.fasta (Nissle Stock-1)"
echo "Samples: 13 Nissle isolates"
echo "Threads: 32"
echo ""

parsnp -r comprehensive_parsnp_analysis/assemblies/N2SKTQ_1_1.fasta \
       -d comprehensive_parsnp_analysis/assemblies \
       -o comprehensive_parsnp_analysis/parsnp_output \
       -c -v -p 32

echo ""
echo "Parsnp analysis completed at $(date)"
echo ""

# Run snp-dists to generate SNP distance matrix
if [ -f comprehensive_parsnp_analysis/parsnp_output/parsnp.xmfa ]; then
    echo "Generating SNP distance matrix..."
    snp-dists comprehensive_parsnp_analysis/parsnp_output/parsnp.xmfa > comprehensive_parsnp_analysis/parsnp_output/parsnp_snp_distances.tsv
    echo "SNP distance matrix saved to: comprehensive_parsnp_analysis/parsnp_output/parsnp_snp_distances.tsv"
fi

echo ""
echo "All done at $(date)"
