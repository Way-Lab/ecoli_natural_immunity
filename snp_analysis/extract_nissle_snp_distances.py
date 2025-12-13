#!/usr/bin/env python3
"""
Extract SNP distances from parsnp alignment and create heatmap
for closely related Nissle isolates.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

# Parsnp alignment file
XMFA_FILE = "Z6M7Y5_comparative_analysis_with_refs/renamed_parsnp_output/parsnp.xmfa"

# Strains to include
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

def clean_strain_name(record_id):
    """Extract clean strain name from record ID."""
    # Remove .fasta extension
    name = record_id.replace('.fasta', '')
    # Handle reference format
    if '.ref' in name:
        name = name.replace('.ref', '')
    if '_hybrid_assembly' in name:
        name = name.replace('_hybrid_assembly', '')

    # Keep the full Nissle strain names (Preg-X, Stock-X, Control-X, etc.)
    # Don't simplify further
    return name

def parse_parsnp_xmfa(xmfa_file):
    """Parse parsnp XMFA file and extract alignment blocks."""
    print(f"Parsing parsnp XMFA file: {xmfa_file}")

    # First, read the header to get sequence mapping
    seq_map = {}  # index -> filename
    with open(xmfa_file, 'r') as f:
        for line in f:
            if line.startswith('##SequenceIndex'):
                idx = int(line.split()[1])
            elif line.startswith('##SequenceFile'):
                # Get everything after "##SequenceFile " (handles spaces in filenames)
                filename = line.split(maxsplit=1)[1].strip()
                seq_map[idx] = clean_strain_name(filename)
            elif line.startswith('>'):
                # End of header, start of alignment blocks
                break

    print(f"Found {len(seq_map)} sequences in header")
    for idx, name in sorted(seq_map.items()):
        print(f"  {idx}: {name}")

    # Now parse alignment blocks
    alignments = {idx: [] for idx in seq_map.keys()}

    with open(xmfa_file, 'r') as f:
        # Skip header
        for line in f:
            if line.startswith('>'):
                break

        current_seq = None
        current_data = []

        for line in f:
            line = line.strip()

            if not line or line.startswith('#'):
                continue

            if line.startswith('>'):
                # Save previous sequence
                if current_seq is not None and current_data:
                    alignments[current_seq].append(''.join(current_data))
                    current_data = []

                # Parse new sequence header
                # Format: >index:start-end +/- ...
                parts = line[1:].split()
                seq_idx = int(parts[0].split(':')[0])
                current_seq = seq_idx

            elif line == '=':
                # End of alignment block
                if current_seq is not None and current_data:
                    alignments[current_seq].append(''.join(current_data))
                    current_data = []
                current_seq = None

            else:
                # Sequence data
                current_data.append(line)

    # Concatenate all blocks for each sequence
    concatenated = {}
    for idx, blocks in alignments.items():
        concatenated[idx] = ''.join(blocks)

    return seq_map, concatenated


def calculate_snp_distances(alignment_file, target_strains):
    """Calculate pairwise SNP distances from XMFA alignment."""
    print(f"Reading alignment from {alignment_file}...")

    # Parse the XMFA file
    seq_map, alignments = parse_parsnp_xmfa(alignment_file)

    # Find target strains
    target_indices = []
    target_names_clean = []

    print("\nMatching target strains:")
    for target in target_strains:
        for idx, name in seq_map.items():
            if target.lower() in name.lower() or name.lower() in target.lower():
                if idx in alignments:
                    target_indices.append(idx)
                    target_names_clean.append(name)
                    print(f"  ✓ {target} -> {name}")
                    break

    if len(target_names_clean) == 0:
        print("\nERROR: No target strains found!")
        print("\nAvailable strains:")
        for idx, name in seq_map.items():
            print(f"  - {name}")
        raise ValueError("No target strains found in alignment")

    print(f"\nFound {len(target_indices)} target strains")

    # Calculate pairwise SNP distances
    n = len(target_indices)
    distances = np.zeros((n, n), dtype=int)

    print("\nCalculating pairwise SNP distances...")
    for i in range(n):
        for j in range(i+1, n):
            idx_i = target_indices[i]
            idx_j = target_indices[j]

            seq_i = alignments[idx_i]
            seq_j = alignments[idx_j]

            # Make sure sequences are same length
            min_len = min(len(seq_i), len(seq_j))

            # Count differences (ignoring gaps)
            diffs = sum(1 for k in range(min_len)
                       if seq_i[k] != seq_j[k] and seq_i[k] != '-' and seq_j[k] != '-')

            distances[i, j] = diffs
            distances[j, i] = diffs

        if (i + 1) % 5 == 0:
            print(f"  Processed {i+1}/{n} strains...")

    # Create DataFrame
    dist_df = pd.DataFrame(distances,
                          index=target_names_clean,
                          columns=target_names_clean)

    return dist_df

def create_heatmap(dist_df, output_file='nissle_snp_heatmap.pdf'):
    """Create a publication-quality heatmap of SNP distances."""

    # Set up the plot
    fig, ax = plt.subplots(figsize=(12, 10))

    # Create heatmap with actual values annotated
    sns.heatmap(dist_df,
                annot=True,           # Show numbers
                fmt='d',              # Integer format
                cmap='YlOrRd',        # Yellow to red colormap
                square=True,          # Square cells
                linewidths=0.5,       # Grid lines
                cbar_kws={'label': 'SNP Distance'},
                ax=ax,
                vmin=0,
                vmax=dist_df.max().max() if dist_df.max().max() > 0 else 10)

    # Customize labels
    ax.set_xlabel('', fontsize=12)
    ax.set_ylabel('', fontsize=12)
    plt.title('Pairwise SNP Distances - Nissle Isolates\n(Core Genome Alignment)',
              fontsize=14, fontweight='bold', pad=20)

    # Rotate x-axis labels
    plt.xticks(rotation=45, ha='right')
    plt.yticks(rotation=0)

    # Adjust layout
    plt.tight_layout()

    # Save
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"\n✓ Heatmap saved to: {output_file}")

    # Also save as PNG
    png_file = output_file.replace('.pdf', '.png')
    plt.savefig(png_file, dpi=300, bbox_inches='tight')
    print(f"✓ Heatmap saved to: {png_file}")

    plt.close()

def main():
    """Main function."""
    print("="*80)
    print("SNP Distance Extraction for Nissle Isolates")
    print("="*80)
    print()

    # Calculate distances
    dist_df = calculate_snp_distances(XMFA_FILE, TARGET_STRAINS)

    # Save distance matrix
    output_tsv = 'nissle_snp_distances.tsv'
    dist_df.to_csv(output_tsv, sep='\t')
    print(f"\n✓ Distance matrix saved to: {output_tsv}")

    # Print summary statistics
    print("\n" + "="*80)
    print("SNP Distance Summary")
    print("="*80)

    # Get upper triangle (exclude diagonal)
    mask = np.triu(np.ones_like(dist_df, dtype=bool), k=1)
    distances = dist_df.values[mask]

    print(f"Number of strain pairs: {len(distances)}")
    print(f"Minimum SNP distance: {distances.min()}")
    print(f"Maximum SNP distance: {distances.max()}")
    print(f"Mean SNP distance: {distances.mean():.2f}")
    print(f"Median SNP distance: {np.median(distances):.2f}")

    # Genome size estimate (typical E. coli)
    genome_size = 5_000_000
    max_divergence = (distances.max() / genome_size) * 100
    print(f"\nMaximum divergence: {distances.max()} SNPs / ~{genome_size/1e6:.1f} Mb = {max_divergence:.4f}%")

    # Create heatmap
    print("\n" + "="*80)
    print("Creating heatmap visualization...")
    print("="*80)
    create_heatmap(dist_df)

    print("\n" + "="*80)
    print("Complete!")
    print("="*80)

if __name__ == "__main__":
    main()
