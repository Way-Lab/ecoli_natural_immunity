#!/usr/bin/env python3
"""
Analyze SNPs from parsnp analysis with outgroups (more restrictive core genome).
Calculate distances and list SNPs by sample.
"""

import pandas as pd
import numpy as np
from collections import defaultdict
import re

# Mapping from file names to display names
FILE_TO_NAME = {
    'N2SKTQ_1_1.fasta.ref': 'Nissle Stock-1 (ref)',
    'N2SKTQ_1_1.fasta': 'Nissle Stock-1',
    'N2SKTQ_2_2.fasta': 'Nissle Stock-2',
    'N2SKTQ_3_3.fasta': 'Nissle Control-1',
    'N2SKTQ_4_4.fasta': 'Nissle Control-2',
    'N2SKTQ_5_5.fasta': 'Nissle Control-3',
    'N2SKTQ_6_6.fasta': 'Nissle Control-4',
    'N2SKTQ_7_7.fasta': 'Nissle Mark stock',
    'N2SKTQ_8_8.fasta': 'CFT073/rpoS',
    'N2SKTQ_9_9.fasta': 'Nissle Preg-1',
    'N2SKTQ_10_10.fasta': 'Nissle Preg-2',
    'N2SKTQ_11_11.fasta': 'MS8707',
    'N2SKTQ_12_12.fasta': 'MS25509',
    'Z6M7Y5_1.fasta': 'Nissle Preg-5',
    'Z6M7Y5_2.fasta': 'Nissle Preg-6',
    'Z6M7Y5_3.fasta': 'Nissle Preg-13',
    'Z6M7Y5_4.fasta': 'Nissle Preg-14',
    'Z6M7Y5_5.fasta': 'Nissle Preg-16',
    'Z6M7Y5_6.fasta': 'Nissle Preg-18',
    'reference_Nissle_RefSeq.fasta': 'EcN RefSeq',
}

# Exclude outgroups AND the reference (we want to show SNPs in query samples, not the ref)
EXCLUDE_FROM_ANALYSIS = ['CFT073/rpoS', 'MS8707', 'MS25509', 'Nissle Preg-16', 'Nissle Preg-18', 'Nissle Stock-1 (ref)']

def parse_vcf(vcf_file):
    """Parse VCF file and extract SNP data."""

    # Read header to get sample names
    with open(vcf_file, 'r') as f:
        for line in f:
            if line.startswith('#CHROM'):
                header = line.strip().split('\t')
                sample_names = header[9:]  # Sample names start after FORMAT column
                break

    # Map to display names
    display_names = [FILE_TO_NAME.get(name, name) for name in sample_names]

    # Read variant data
    variants = []
    with open(vcf_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            chrom = fields[0]
            pos = int(fields[1])
            ref = fields[3]
            alt = fields[4]
            genotypes = fields[9:]

            # Convert genotypes to integers (0=ref, 1=alt, .=missing)
            gt_values = []
            for gt in genotypes:
                if gt == '0':
                    gt_values.append(0)
                elif gt == '1':
                    gt_values.append(1)
                else:
                    gt_values.append(-1)  # Missing

            variants.append({
                'chrom': chrom,
                'pos': pos,
                'ref': ref,
                'alt': alt,
                'genotypes': gt_values
            })

    return display_names, variants

def list_snps_by_sample(sample_names, variants, reference_idx=0, exclude_samples=None):
    """List all SNPs for each sample compared to reference."""

    if exclude_samples is None:
        exclude_samples = []

    snps_by_sample = defaultdict(list)

    for var in variants:
        ref_gt = var['genotypes'][reference_idx]

        for i, sample in enumerate(sample_names):
            if i == reference_idx:
                continue  # Skip reference itself

            if sample in exclude_samples:
                continue  # Skip outgroups

            gt = var['genotypes'][i]

            # If this sample has alternate allele
            if gt == 1:
                snps_by_sample[sample].append({
                    'position': var['pos'],
                    'ref': var['ref'],
                    'alt': var['alt'],
                    'change': f"{var['ref']}>{var['alt']}"
                })

    return snps_by_sample

def main():
    print("="*80)
    print("SNP Analysis from Parsnp with Outgroups (Restrictive Core Genome)")
    print("="*80)
    print()

    vcf_file = 'parsnp_with_outgroups/parsnp_output/parsnp.vcf'

    # Parse VCF
    print("Parsing VCF file...")
    sample_names, variants = parse_vcf(vcf_file)
    print(f"Found {len(sample_names)} samples and {len(variants)} variant sites\n")

    # Get SNPs by sample (using first sample as reference, exclude outgroups)
    print("="*80)
    print(f"SNPs by Sample (compared to {sample_names[0]})")
    print("Excluding outgroups from detailed analysis")
    print("="*80)

    snps_by_sample = list_snps_by_sample(sample_names, variants, reference_idx=0,
                                         exclude_samples=EXCLUDE_FROM_ANALYSIS)

    # Save detailed SNP list
    snp_details_file = 'parsnp_with_outgroups/parsnp_output/snp_details_by_sample_filtered.txt'
    with open(snp_details_file, 'w') as f:
        f.write(f"SNP Details by Sample (compared to {sample_names[0]})\n")
        f.write("Core genome defined with outgroups (CFT073/rpoS, MS8707, MS25509, Preg-16, Preg-18)\n")
        f.write("="*80 + "\n\n")

        for sample in sorted(snps_by_sample.keys()):
            snps = snps_by_sample[sample]
            f.write(f"{sample}: {len(snps)} SNPs\n")
            f.write("-" * 80 + "\n")

            if len(snps) > 0:
                for snp in snps:
                    f.write(f"  Position {snp['position']:>7}: {snp['change']:>5} "
                           f"(Ref: {snp['ref']}, Alt: {snp['alt']})\n")
            f.write("\n")

    print(f"\nDetailed SNP list saved to: {snp_details_file}\n")

    # Print summary
    print("Summary:")
    print("-" * 80)
    for sample in sorted(snps_by_sample.keys()):
        n_snps = len(snps_by_sample[sample])
        print(f"{sample:30s}: {n_snps:3d} SNPs")

    print()
    print("="*80)
    print("Complete!")
    print("="*80)

if __name__ == "__main__":
    main()
