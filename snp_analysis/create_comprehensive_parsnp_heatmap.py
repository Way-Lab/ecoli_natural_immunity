#!/usr/bin/env python3
"""
Create SNP distance heatmap from comprehensive parsnp analysis.
"""

import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns

# Set font types for editable text in Illustrator
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

# Target strains for the heatmap (exclude reference samples, include EcN RefSeq for comparison)
TARGET_STRAINS = [
    'Nissle Stock-1',
    'Nissle Stock-2',
    'Nissle Control-1',
    'Nissle Control-2',
    'Nissle Control-3',
    'Nissle Control-4',
    'Nissle Mark stock',
    'Nissle Preg-1',
    'Nissle Preg-2',
    'Nissle Preg-5',
    'Nissle Preg-6',
    'Nissle Preg-13',
    'Nissle Preg-14',
    'EcN RefSeq'
]

def load_and_filter_matrix(file_path, target_strains):
    """Load SNP distance matrix and filter to target strains."""

    # Load the full matrix
    df = pd.read_csv(file_path, sep='\t', index_col=0)

    print(f"Loaded SNP distance matrix: {df.shape}")
    print(f"Samples in matrix: {df.shape[0]}\n")

    # Filter to target strains (keep order)
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

def create_heatmap(dist_df, output_file='comprehensive_parsnp_snp_heatmap.pdf'):
    """Create a publication-quality heatmap of SNP distances."""

    # Set up the plot - wider figure with good height
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
                vmax=80)              # Cap color scale to show full range

    # Customize labels
    ax.set_xlabel('', fontsize=12)
    ax.set_ylabel('', fontsize=12)
    plt.title('Pairwise SNP Distances - Comprehensive Nissle Analysis\n(Reference: Nissle Stock-1 via Parsnp)',
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
    print("SNP Distance Heatmap from Comprehensive Parsnp Analysis")
    print("="*80)
    print()

    # Load and filter matrix
    input_file = 'comprehensive_parsnp_analysis/parsnp_output/snp_distance_matrix.tsv'
    dist_df = load_and_filter_matrix(input_file, TARGET_STRAINS)

    # Save filtered matrix
    output_tsv = 'comprehensive_parsnp_snp_distances.tsv'
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
    create_heatmap(dist_df, 'comprehensive_parsnp_snp_heatmap.pdf')

    print("\n" + "="*80)
    print("Complete!")
    print("="*80)

if __name__ == "__main__":
    main()
