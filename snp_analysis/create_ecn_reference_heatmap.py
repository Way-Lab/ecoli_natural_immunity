#!/usr/bin/env python3
"""
Create SNP distance heatmap using the EcN reference genome parsnp analysis.
Filters to only the E. coli strains from the most recent analysis.
"""

import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns

# Set font types for editable text in Illustrator
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

# Target strains for the heatmap (from most recent analysis)
# Include EcN Reference to see SNP differences between isolates and reference genome
TARGET_STRAINS = [
    'EcN Reference',
    'Nissle Stock-1',
    'Nissle Stock-2',
    'Nissle Control-1',
    'Nissle Control-2',
    'Nissle Control-3',
    'Nissle Control-4',
    'Nissle Preg-1',
    'Nissle Preg-2',
    'Nissle Preg-5',
    'Nissle Preg-6',
    'Nissle Preg-13',
    'Nissle Preg-14'
]

# Mapping from parsnp file names to our display names
# Based on FileNameConversion.csv
FILE_TO_NAME = {
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
    'N2SKTQ_13_13.fasta': 'Nissle Preg-9',
    'N2SKTQ_14_14.fasta': 'Nissle Preg-11',
    'Z6M7Y5_1.fasta': 'Nissle Preg-5',
    'Z6M7Y5_2.fasta': 'Nissle Preg-6',
    'Z6M7Y5_3.fasta': 'Nissle Preg-13',
    'Z6M7Y5_4.fasta': 'Nissle Preg-14',
    'reference_Nissle_RefSeq.fasta.ref': 'EcN Reference',
    'reference_Nissle_RefSeq.fasta': 'EcN Reference'
}

def load_and_filter_matrix(file_path, target_strains, file_to_name):
    """Load parsnp SNP matrix and filter to target strains."""

    # Load the full matrix (skip version line if present)
    df = pd.read_csv(file_path, sep='\t', index_col=0)

    print(f"Loaded SNP distance matrix: {df.shape}")
    print(f"Samples in matrix: {df.shape[0]}\n")

    print("Original sample names in matrix:")
    for idx in df.index:
        print(f"  - {idx}")
    print()

    # Remove the duplicate reference entry (keep only the .ref version)
    if 'reference_Nissle_RefSeq.fasta' in df.index:
        print("Removing duplicate reference entry: reference_Nissle_RefSeq.fasta")
        df = df.drop('reference_Nissle_RefSeq.fasta', axis=0)
        df = df.drop('reference_Nissle_RefSeq.fasta', axis=1)
        print(f"After removing duplicate: {df.shape}\n")

    # Map file names to our names
    df.index = [file_to_name.get(idx, idx) for idx in df.index]
    df.columns = [file_to_name.get(col, col) for col in df.columns]

    print("After name mapping:")
    for idx in df.index:
        print(f"  - {idx}")
    print()

    # Filter to target strains (only include those that exist in matrix)
    available_targets = [s for s in target_strains if s in df.index]

    if len(available_targets) != len(target_strains):
        missing = set(target_strains) - set(available_targets)
        print(f"WARNING: Missing strains from EcN reference analysis:")
        for strain in missing:
            print(f"  ✗ {strain}")
        print()

    print(f"Extracting {len(available_targets)} available target strains:")
    for strain in available_targets:
        print(f"  ✓ {strain}")
    print()

    # Extract submatrix
    filtered_df = df.loc[available_targets, available_targets]

    return filtered_df

def create_heatmap(dist_df, output_file='ecn_reference_snp_heatmap.pdf'):
    """Create a publication-quality heatmap of SNP distances."""

    # Set up the plot - wider figure with shorter height
    fig, ax = plt.subplots(figsize=(13, 7))

    # Create heatmap with actual values annotated
    sns.heatmap(dist_df,
                annot=True,           # Show numbers
                fmt='d',              # Integer format
                cmap='YlOrRd',        # Yellow to red colormap
                square=False,         # Allow rectangular cells
                linewidths=0.5,       # Grid lines
                cbar_kws={'label': 'Core Genome SNP Distance'},
                ax=ax,
                vmin=0,
                vmax=70)              # Cap color scale to emphasize differences

    # Customize labels
    ax.set_xlabel('', fontsize=12)
    ax.set_ylabel('', fontsize=12)
    plt.title('Pairwise SNP Distances - Nissle Isolates\n(Aligned to EcN Reference Genome via Parsnp)',
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
    print("SNP Distance Heatmap from EcN Reference Parsnp Analysis")
    print("="*80)
    print()

    # Load and filter matrix
    input_file = 'N2SKTQ_results/parsnp_results_nissle_refseq/parsnp_snp_distances.tsv'
    dist_df = load_and_filter_matrix(input_file, TARGET_STRAINS, FILE_TO_NAME)

    # Save filtered matrix
    output_tsv = 'ecn_reference_snp_distances.tsv'
    dist_df.to_csv(output_tsv, sep='\t')
    print(f"✓ Filtered distance matrix saved to: {output_tsv}\n")

    # Print summary statistics
    print("="*80)
    print("SNP Distance Summary")
    print("="*80)

    # Get upper triangle (exclude diagonal)
    mask = np.triu(np.ones_like(dist_df, dtype=bool), k=1)
    distances = dist_df.values[mask]

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
    create_heatmap(dist_df, 'ecn_reference_snp_heatmap.pdf')

    print("\n" + "="*80)
    print("Complete!")
    print("="*80)

if __name__ == "__main__":
    main()
