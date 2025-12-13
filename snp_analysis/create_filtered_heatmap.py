#!/usr/bin/env python3
"""
Create SNP distance heatmap excluding outgroups.
Analysis includes outgroups to define core genome,
but heatmap only shows target Nissle strains + EcN RefSeq.
"""

import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns

# Set font types for editable text in Illustrator
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

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

# Target strains for heatmap (exclude outgroups)
# Note: Stock-1 is the reference, so it will show as 0 SNPs from itself
TARGET_STRAINS = [
    'Nissle Stock-1 (ref)',  # This is the reference
    'Nissle Stock-2',
    'Nissle Mark stock',
    'Nissle Control-1',
    'Nissle Control-2',
    'Nissle Control-3',
    'Nissle Control-4',
    'Nissle Preg-1',
    'Nissle Preg-2',
    'Nissle Preg-5',
    'Nissle Preg-6',
    'Nissle Preg-13',
    'Nissle Preg-14',
    'EcN RefSeq'
]

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

    # Add reference genome (Stock-1) as first sample
    # It will have genotype 0 (reference allele) at all positions
    display_names.insert(0, 'Nissle Stock-1 (ref)')

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

            # Insert reference genome genotype (always 0) at the beginning
            gt_values.insert(0, 0)

            variants.append({
                'chrom': chrom,
                'pos': pos,
                'ref': ref,
                'alt': alt,
                'genotypes': gt_values
            })

    return display_names, variants

def calculate_snp_distances(sample_names, variants):
    """Calculate pairwise SNP distances."""

    n_samples = len(sample_names)
    distance_matrix = np.zeros((n_samples, n_samples), dtype=int)

    # Calculate pairwise distances
    for i in range(n_samples):
        for j in range(i+1, n_samples):
            dist = 0
            for var in variants:
                gt_i = var['genotypes'][i]
                gt_j = var['genotypes'][j]
                # Only count if both are not missing and different
                if gt_i != -1 and gt_j != -1 and gt_i != gt_j:
                    dist += 1
            distance_matrix[i, j] = dist
            distance_matrix[j, i] = dist

    # Create DataFrame
    df = pd.DataFrame(distance_matrix, index=sample_names, columns=sample_names)
    return df

def create_heatmap(dist_df, output_file='filtered_parsnp_snp_heatmap.pdf'):
    """Create a publication-quality heatmap of SNP distances."""

    # Set up the plot
    fig, ax = plt.subplots(figsize=(14, 10))

    # Create heatmap with actual values annotated
    sns.heatmap(dist_df,
                annot=True,           # Show numbers
                fmt='d',              # Integer format
                cmap='YlOrRd',        # Yellow to red colormap
                square=True,          # Square cells
                linewidths=0.5,       # Grid lines
                cbar_kws={'label': 'Core Genome SNP Distance'},
                ax=ax,
                vmin=0,
                vmax=80)              # Cap color scale

    # Customize labels
    ax.set_xlabel('', fontsize=12)
    ax.set_ylabel('', fontsize=12)
    plt.title('Pairwise SNP Distances - Nissle Isolates\n(Reference: Nissle Stock-1, Core genome defined with outgroups)',
              fontsize=14, fontweight='bold', pad=20)

    # Rotate x-axis labels
    plt.xticks(rotation=45, ha='right')
    plt.yticks(rotation=0)

    # Adjust layout
    plt.tight_layout()

    # Save
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"✓ Heatmap saved to: {output_file}")

    # Also save as PNG
    png_file = output_file.replace('.pdf', '.png')
    plt.savefig(png_file, dpi=300, bbox_inches='tight')
    print(f"✓ Heatmap saved to: {png_file}")

    plt.close()

def main():
    print("="*80)
    print("SNP Distance Heatmap (Filtered, Core Genome from Analysis with Outgroups)")
    print("="*80)
    print()

    # Parse VCF
    vcf_file = 'parsnp_with_outgroups/parsnp_output/parsnp.vcf'
    print(f"Parsing VCF file: {vcf_file}")
    sample_names, variants = parse_vcf(vcf_file)
    print(f"Found {len(sample_names)} samples and {len(variants)} variant sites\n")

    # Calculate distances for all samples
    print("Calculating pairwise SNP distances for all samples...")
    full_dist_df = calculate_snp_distances(sample_names, variants)

    # Save full distance matrix
    full_dist_df.to_csv('parsnp_with_outgroups/parsnp_output/snp_distance_matrix_all.tsv', sep='\t')
    print(f"✓ Full distance matrix saved\n")

    # Filter to target strains only
    available_targets = [s for s in TARGET_STRAINS if s in full_dist_df.index]
    print(f"Filtering to {len(available_targets)} target strains for heatmap:")
    for strain in available_targets:
        print(f"  ✓ {strain}")
    print()

    filtered_dist_df = full_dist_df.loc[available_targets, available_targets]

    # Save filtered matrix
    output_tsv = 'filtered_parsnp_snp_distances.tsv'
    filtered_dist_df.to_csv(output_tsv, sep='\t')
    print(f"✓ Filtered distance matrix saved to: {output_tsv}\n")

    # Print summary statistics
    print("="*80)
    print("SNP Distance Summary (Filtered)")
    print("="*80)

    # Get upper triangle (exclude diagonal)
    mask = np.triu(np.ones_like(filtered_dist_df, dtype=bool), k=1)
    distances = filtered_dist_df.values[mask]

    print(f"Number of strain pairs: {len(distances)}")
    print(f"Minimum SNP distance: {distances.min()}")
    print(f"Maximum SNP distance: {distances.max()}")
    print(f"Mean SNP distance: {distances.mean():.2f}")
    print(f"Median SNP distance: {np.median(distances):.1f}")

    # Genome size
    genome_size = 5_000_000
    max_divergence = (distances.max() / genome_size) * 100
    print(f"\nMaximum divergence: {distances.max()} SNPs / ~{genome_size/1e6:.1f} Mb = {max_divergence:.4f}%")

    # Create heatmap
    print("\n" + "="*80)
    print("Creating heatmap visualization...")
    print("="*80)
    create_heatmap(filtered_dist_df, 'filtered_parsnp_snp_heatmap.pdf')

    print("\n" + "="*80)
    print("Complete!")
    print("="*80)
    print("\nNote: Outgroups (CFT073/rpoS, MS8707, MS25509, Preg-16, Preg-18)")
    print("      were included in parsnp analysis to define core genome")
    print("      but are excluded from this heatmap visualization.")

if __name__ == "__main__":
    main()
