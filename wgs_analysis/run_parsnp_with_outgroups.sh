#!/bin/bash

# Activate conda environment
source /home/david/miniforge3/bin/activate assembly-qc

# Run parsnp with Nissle Stock-1 as reference, including outgroups
echo "Starting parsnp analysis with outgroups at $(date)"
echo "Reference: N2SKTQ_1_1.fasta (Nissle Stock-1)"
echo "Samples: 13 Nissle isolates + 5 outgroups + 1 EcN RefSeq"
echo "Threads: 32"
echo ""

parsnp -r parsnp_with_outgroups/assemblies/N2SKTQ_1_1.fasta \
       -d parsnp_with_outgroups/assemblies \
       -o parsnp_with_outgroups/parsnp_output \
       -c -v -p 32

echo ""
echo "Parsnp analysis completed at $(date)"
echo ""

# Convert to VCF
if [ -f parsnp_with_outgroups/parsnp_output/parsnp.ggr ]; then
    echo "Converting to VCF format..."
    harvesttools -i parsnp_with_outgroups/parsnp_output/parsnp.ggr \
                 -V parsnp_with_outgroups/parsnp_output/parsnp.vcf
    echo "VCF saved to: parsnp_with_outgroups/parsnp_output/parsnp.vcf"
fi

echo ""
echo "All done at $(date)"
