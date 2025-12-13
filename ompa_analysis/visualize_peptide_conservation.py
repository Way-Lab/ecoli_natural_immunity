#!/usr/bin/env python3
"""
Create visualization of peptide conservation across Enterobacteriaceae.
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import Rectangle
import numpy as np

def create_conservation_heatmap(csv_file, output_file):
    """
    Create a heatmap showing conservation at each position for each peptide.
    """
    df = pd.read_csv(csv_file)

    # Get unique peptides
    peptides = df['Peptide'].unique()

    # Create figure with subplots for each peptide
    fig, axes = plt.subplots(len(peptides), 1, figsize=(14, 2*len(peptides)))

    if len(peptides) == 1:
        axes = [axes]

    for idx, peptide in enumerate(peptides):
        peptide_data = df[df['Peptide'] == peptide]

        # Get all positions
        positions = sorted(peptide_data['Position'].unique())

        # Create conservation matrix
        conservation_values = []
        for pos in positions:
            pos_data = peptide_data[peptide_data['Position'] == pos]
            conservation = pos_data['Conservation_Percent'].iloc[0]
            conservation_values.append(conservation)

        # Plot
        ax = axes[idx]
        bars = ax.bar(positions, conservation_values, color='steelblue', edgecolor='black', linewidth=0.5)

        # Color bars based on conservation
        for i, (bar, cons) in enumerate(zip(bars, conservation_values)):
            if cons >= 90:
                bar.set_color('darkgreen')
            elif cons >= 75:
                bar.set_color('orange')
            else:
                bar.set_color('red')

        ax.set_ylabel('Conservation (%)', fontsize=10)
        ax.set_xlabel('Position', fontsize=10)
        ax.set_title(f'{peptide}', fontsize=11, fontweight='bold')
        ax.set_ylim(0, 105)
        ax.axhline(y=90, color='gray', linestyle='--', linewidth=0.8, alpha=0.5)
        ax.axhline(y=75, color='gray', linestyle='--', linewidth=0.8, alpha=0.5)
        ax.grid(axis='y', alpha=0.3)
        ax.set_xticks(positions)

    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Conservation heatmap saved to: {output_file}")


def create_variability_matrix(csv_file, output_file):
    """
    Create a matrix showing amino acid variability at each position.
    """
    df = pd.read_csv(csv_file)

    # Get peptides with found sequences
    peptides_with_data = []
    for peptide in df['Peptide'].unique():
        peptide_data = df[df['Peptide'] == peptide]
        if len(peptide_data) > 0:
            peptides_with_data.append(peptide)

    if not peptides_with_data:
        print("No peptides with data to visualize")
        return

    # Filter to most interesting peptides (those with some variability)
    interesting_peptides = []
    for peptide in peptides_with_data[:8]:  # Limit to first 8 for readability
        peptide_data = df[df['Peptide'] == peptide]
        min_conservation = peptide_data['Conservation_Percent'].min()
        if min_conservation < 100:
            interesting_peptides.append(peptide)

    if not interesting_peptides:
        interesting_peptides = peptides_with_data[:4]

    # Create figure
    n_peptides = len(interesting_peptides)
    fig, axes = plt.subplots(n_peptides, 1, figsize=(16, 3*n_peptides))

    if n_peptides == 1:
        axes = [axes]

    for idx, peptide in enumerate(interesting_peptides):
        peptide_data = df[df['Peptide'] == peptide]

        # Get positions and create matrix
        positions = sorted(peptide_data['Position'].unique())
        max_pos = len(positions)

        # Get all amino acids observed
        all_aas = sorted(peptide_data['Observed_AA'].unique())

        # Create matrix
        matrix = np.zeros((len(all_aas), max_pos))
        aa_to_idx = {aa: i for i, aa in enumerate(all_aas)}

        for _, row in peptide_data.iterrows():
            pos_idx = row['Position'] - 1
            aa_idx = aa_to_idx[row['Observed_AA']]
            matrix[aa_idx, pos_idx] = row['Frequency_Percent']

        # Plot
        ax = axes[idx]
        im = ax.imshow(matrix, aspect='auto', cmap='YlOrRd', interpolation='nearest')

        # Set ticks
        ax.set_xticks(range(max_pos))
        ax.set_xticklabels(positions, fontsize=9)
        ax.set_yticks(range(len(all_aas)))
        ax.set_yticklabels(all_aas, fontsize=9)

        ax.set_xlabel('Position', fontsize=10)
        ax.set_ylabel('Amino Acid', fontsize=10)
        ax.set_title(f'{peptide} - Amino Acid Frequency', fontsize=11, fontweight='bold')

        # Add colorbar
        cbar = plt.colorbar(im, ax=ax)
        cbar.set_label('Frequency (%)', fontsize=9)

        # Add text annotations for values > 5%
        for i in range(len(all_aas)):
            for j in range(max_pos):
                if matrix[i, j] > 5:
                    text = ax.text(j, i, f'{matrix[i, j]:.0f}',
                                   ha="center", va="center", color="black", fontsize=7)

    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Variability matrix saved to: {output_file}")


def generate_summary_statistics(csv_file, output_file):
    """
    Generate summary statistics table.
    """
    df = pd.read_csv(csv_file)

    summary = []
    for peptide in df['Peptide'].unique():
        peptide_data = df[df['Peptide'] == peptide]

        if len(peptide_data) == 0:
            continue

        positions = peptide_data['Position'].unique()
        avg_conservation = peptide_data.groupby('Position')['Conservation_Percent'].first().mean()
        min_conservation = peptide_data.groupby('Position')['Conservation_Percent'].first().min()
        max_conservation = peptide_data.groupby('Position')['Conservation_Percent'].first().max()

        variable_positions = len(peptide_data[peptide_data['Conservation_Percent'] < 90]['Position'].unique())

        summary.append({
            'Peptide': peptide,
            'Length': len(positions),
            'Avg_Conservation_%': round(avg_conservation, 2),
            'Min_Conservation_%': round(min_conservation, 2),
            'Max_Conservation_%': round(max_conservation, 2),
            'Variable_Positions': variable_positions,
        })

    summary_df = pd.DataFrame(summary)
    summary_df.to_csv(output_file, index=False)
    print(f"Summary statistics saved to: {output_file}")

    # Print to console
    print("\n" + "="*80)
    print("PEPTIDE CONSERVATION SUMMARY")
    print("="*80)
    print(summary_df.to_string(index=False))
    print("="*80)


def main():
    import sys

    if len(sys.argv) < 2:
        print("Usage: python visualize_peptide_conservation.py <csv_file> [output_prefix]")
        sys.exit(1)

    csv_file = sys.argv[1]
    output_prefix = sys.argv[2] if len(sys.argv) > 2 else 'peptide_conservation'

    # Generate visualizations
    create_conservation_heatmap(csv_file, f'{output_prefix}_conservation.png')
    create_variability_matrix(csv_file, f'{output_prefix}_variability.png')
    generate_summary_statistics(csv_file, f'{output_prefix}_summary.csv')


if __name__ == "__main__":
    main()
