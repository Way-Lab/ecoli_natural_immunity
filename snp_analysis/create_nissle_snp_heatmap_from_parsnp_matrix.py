#!/usr/bin/env python3
"""
Create SNP distance heatmap for Nissle isolates using the actual parsnp SNP matrix.
"""

import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns

# Set font types for editable text in Illustrator
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

# Target strains for the heatmap (user requested)
TARGET_STRAINS = [
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

# Mapping from file names to our names
FILE_TO_NAME = {
    'N2SKTQ_1_1_hybrid_assembly.fasta.ref': 'Nissle Stock-1',
    'Nissle Stock-2.fasta': 'Nissle Stock-2',
    'Nissle Control-1.fasta': 'Nissle Control-1',
    'Nissle Control-2.fasta': 'Nissle Control-2',
    'Nissle Control-3.fasta': 'Nissle Control-3',
    'Nissle Control-4.fasta': 'Nissle Control-4',
    'Nissle Preg-1.fasta': 'Nissle Preg-1',
    'Nissle Preg-2.fasta': 'Nissle Preg-2',
    'Nissle Preg-5.fasta': 'Nissle Preg-5',
    'Nissle Preg-6.fasta': 'Nissle Preg-6',
    'Nissle Preg-13.fasta': 'Nissle Preg-13',
    'Nissle Preg-14.fasta': 'Nissle Preg-14'
}

def load_and_filter_matrix(file_path, target_strains, file_to_name):
    """Load parsnp SNP matrix and filter to target strains."""

    # Load the full matrix
    df = pd.read_csv(file_path, sep='\t', index_col=0)

    print(f"Loaded SNP distance matrix: {df.shape}")
    print(f"Samples in matrix: {df.shape[0]}\n")

    # Map file names to our names
    df.index = [file_to_name.get(idx, idx) for idx in df.index]
    df.columns = [file_to_name.get(col, col) for col in df.columns]

    # Filter to target strains
    available_targets = [s for s in target_strains if s in df.index]

    if len(available_targets) != len(target_strains):
        missing = set(target_strains) - set(available_targets)
        print(f"WARNING: Missing strains: {missing}")

    print(f"Extracting {len(available_targets)} target strains:")
    for strain in available_targets:
        print(f"  ✓ {strain}")
    print()

    # Extract submatrix
    filtered_df = df.loc[available_targets, available_targets]

    return filtered_df

def create_heatmap(dist_df, output_file='nissle_snp_heatmap_from_parsnp.pdf'):
    """Create a publication-quality heatmap of SNP distances."""

    # Set up the plot - wider figure with shorter height
    fig, ax = plt.subplots(figsize=(14, 7))

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
                vmax=50)              # Cap color scale at 50 to emphasize low divergence

    # Customize labels
    ax.set_xlabel('', fontsize=12)
    ax.set_ylabel('', fontsize=12)
    plt.title('Pairwise SNP Distances - Nissle Isolates\n(Core Genome Alignment via Parsnp)',
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
    print("SNP Distance Heatmap from Parsnp Matrix")
    print("="*80)
    print()

    # Load and filter matrix
    input_file = 'Z6M7Y5_comparative_analysis/renamed_parsnp_output/snp_distance_matrix.tsv'
    dist_df = load_and_filter_matrix(input_file, TARGET_STRAINS, FILE_TO_NAME)

    # Save filtered matrix
    output_tsv = 'nissle_snp_distances_from_parsnp.tsv'
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
    create_heatmap(dist_df, 'nissle_snp_heatmap_from_parsnp.pdf')

    print("\n" + "="*80)
    print("Complete!")
    print("="*80)

if __name__ == "__main__":
    main()
