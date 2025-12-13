#!/usr/bin/env python3
"""
Create a comprehensive summary of SNPs by sample from filtered analysis.
"""

import pandas as pd

def main():
    # Load annotated SNPs
    df = pd.read_csv('parsnp_with_outgroups/parsnp_output/snps_annotated_filtered.tsv', sep='\t')

    # Create output file
    output_file = 'parsnp_with_outgroups/parsnp_output/snp_summary_by_sample_filtered.txt'

    with open(output_file, 'w') as f:
        f.write('='*100 + '\n')
        f.write('SNP ANALYSIS - SNPs by Sample with Gene Annotations\n')
        f.write('Reference: Nissle Stock-1 (N2SKTQ_1_1)\n')
        f.write('Core genome defined with outgroups (CFT073/rpoS, MS8707, MS25509, Preg-16, Preg-18)\n')
        f.write('='*100 + '\n\n')

        # Get unique samples
        all_samples = set()
        for samples_str in df['samples']:
            if pd.notna(samples_str) and samples_str != 'reference':
                all_samples.update([s.strip() for s in samples_str.split(',')])

        # Print for each sample
        for sample in sorted(all_samples):
            sample_snps = df[df['samples'].str.contains(sample, na=False, regex=False)]

            f.write(f'\n{sample}: {len(sample_snps)} SNPs\n')
            f.write('-' * 100 + '\n')

            for idx, row in sample_snps.iterrows():
                if row['gene'] != 'intergenic':
                    gene_info = f"{row['gene']} - {row['product']}"
                else:
                    gene_info = 'intergenic region'

                f.write(f"  Position {row['position']:>7}: {row['change']:>5}  in  {gene_info}\n")

        f.write('\n' + '='*100 + '\n')
        f.write('Complete!\n')
        f.write('='*100 + '\n')

    print(f"SNP summary by sample saved to: {output_file}")

    # Also print to console
    print()
    print('='*100)
    print('SNP ANALYSIS - SNPs by Sample with Gene Annotations')
    print('Reference: Nissle Stock-1 (N2SKTQ_1_1)')
    print('Core genome defined with outgroups')
    print('='*100)
    print()

    # Get unique samples
    all_samples = set()
    for samples_str in df['samples']:
        if pd.notna(samples_str) and samples_str != 'reference':
            all_samples.update([s.strip() for s in samples_str.split(',')])

    # Print summary counts
    print("Summary:")
    print("-" * 80)
    for sample in sorted(all_samples):
        sample_snps = df[df['samples'].str.contains(sample, na=False, regex=False)]
        print(f'{sample:30s}: {len(sample_snps):3d} SNPs')

    print()
    print("="*80)
    print("Note: The restrictive core genome excludes regions like SfaX")
    print("      where Stock-2 had 13 unique SNPs in the comprehensive analysis.")
    print("="*80)

if __name__ == "__main__":
    main()
