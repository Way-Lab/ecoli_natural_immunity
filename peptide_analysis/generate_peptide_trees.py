#!/usr/bin/env python3
"""
Generate phylogenetic trees (dendrograms) for each peptide alignment.
Uses distance-based methods suitable for short sequences.
"""

from Bio import SeqIO, Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
import os

# Define peptides
PEPTIDES = {
    'Peptide_1': 'WSQYHDTGFINNNGPTHEN',
    'Peptide_2': 'GRMPYKGSVENGAYKA',
    'Peptide_3a': 'KSNVYGKN',
    'Peptide_3b': 'KANVPGGASFKD',
    'Peptide_4': 'TNNIGDAHTIGTRPDNGM',
    'Peptide_2b': 'GRMPYKGDNINGAYKA',
}


def create_alignment_from_fasta(fasta_file):
    """
    Read FASTA and create a MultipleSeqAlignment object.
    Make sure all sequence IDs are unique.
    """
    sequences = list(SeqIO.parse(fasta_file, 'fasta'))

    if len(sequences) < 2:
        return None, 0

    # Make IDs unique by adding index
    seen_ids = {}
    for seq in sequences:
        original_id = seq.id
        if original_id in seen_ids:
            seen_ids[original_id] += 1
            seq.id = f"{original_id}_{seen_ids[original_id]}"
            seq.name = seq.id
        else:
            seen_ids[original_id] = 0

    # Create alignment (sequences are already aligned - same length)
    alignment = MultipleSeqAlignment(sequences)

    return alignment, len(sequences)


def simplify_name(name):
    """
    Simplify sequence names for better tree readability.
    """
    # Remove REFERENCE_ prefix
    if name.startswith('REFERENCE_'):
        return 'REFERENCE'

    # Extract main species name
    parts = name.split('|')
    if len(parts) > 0:
        species_part = parts[0].strip()
        # Shorten long names
        if len(species_part) > 30:
            # Try to get species and accession
            if '_' in species_part:
                parts = species_part.split('_')
                if len(parts) >= 2:
                    return f"{parts[0]}_{parts[1]}"[:25]
        return species_part[:25]

    return name[:25]


def create_tree_plot(tree, peptide_name, output_file, title=None):
    """
    Create a publication-quality tree plot.
    """
    # Determine figure size based on number of terminals
    num_terminals = len(tree.get_terminals())

    # Calculate appropriate figure height (0.3 inches per terminal, min 6, max 20)
    fig_height = max(6, min(20, num_terminals * 0.3))
    fig_width = 12

    fig, ax = plt.subplots(figsize=(fig_width, fig_height))

    # Simplify terminal names
    for terminal in tree.get_terminals():
        terminal.name = simplify_name(terminal.name)

    # Draw the tree
    Phylo.draw(tree, axes=ax, do_show=False, show_confidence=False)

    # Add title
    if title:
        ax.set_title(title, fontsize=14, fontweight='bold', pad=20)
    else:
        ax.set_title(f'Phylogenetic Tree: {peptide_name}', fontsize=14, fontweight='bold', pad=20)

    # Improve layout
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)

    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()

    print(f"    Tree plot saved: {output_file}")


def create_circular_tree_plot(tree, peptide_name, output_file):
    """
    Create a circular/radial tree plot.
    """
    num_terminals = len(tree.get_terminals())

    # For circular, make it square
    fig_size = max(8, min(16, num_terminals * 0.2))

    fig, ax = plt.subplots(figsize=(fig_size, fig_size))

    # Simplify terminal names
    for terminal in tree.get_terminals():
        terminal.name = simplify_name(terminal.name)

    # Draw circular tree
    Phylo.draw(tree, axes=ax, do_show=False, show_confidence=False)

    # Convert to radial/circular
    ax.set_title(f'{peptide_name} - Circular Tree', fontsize=14, fontweight='bold', pad=20)

    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()

    print(f"    Circular tree saved: {output_file}")


def generate_tree_for_peptide(fasta_file, peptide_name, output_dir):
    """
    Generate phylogenetic tree for a single peptide.
    """
    print(f"\n{'='*80}")
    print(f"Generating tree for {peptide_name}")
    print(f"{'='*80}")

    # Read alignment
    alignment, num_seqs = create_alignment_from_fasta(fasta_file)

    if alignment is None:
        print(f"  WARNING: Not enough sequences (<2) to build tree")
        return False

    print(f"  Sequences: {num_seqs}")
    print(f"  Alignment length: {alignment.get_alignment_length()} aa")

    # Calculate distance matrix
    print(f"  Calculating distance matrix...")
    calculator = DistanceCalculator('identity')
    distance_matrix = calculator.get_distance(alignment)

    # Build tree using UPGMA (good for closely related sequences)
    print(f"  Building UPGMA tree...")
    constructor = DistanceTreeConstructor()
    upgma_tree = constructor.upgma(distance_matrix)

    # Build tree using Neighbor Joining (alternative method)
    print(f"  Building Neighbor-Joining tree...")
    nj_tree = constructor.nj(distance_matrix)

    # Save trees in Newick format
    newick_upgma = f"{output_dir}/{peptide_name}_upgma.nwk"
    newick_nj = f"{output_dir}/{peptide_name}_nj.nwk"

    Phylo.write(upgma_tree, newick_upgma, 'newick')
    Phylo.write(nj_tree, newick_nj, 'newick')

    print(f"  Newick files saved:")
    print(f"    UPGMA: {newick_upgma}")
    print(f"    NJ: {newick_nj}")

    # Create visualizations
    print(f"  Creating visualizations...")

    # UPGMA rectangular tree
    create_tree_plot(
        upgma_tree,
        peptide_name,
        f"{output_dir}/{peptide_name}_upgma_tree.png",
        title=f"{peptide_name} - UPGMA Tree ({num_seqs} sequences)"
    )

    # NJ rectangular tree
    create_tree_plot(
        nj_tree,
        peptide_name,
        f"{output_dir}/{peptide_name}_nj_tree.png",
        title=f"{peptide_name} - Neighbor-Joining Tree ({num_seqs} sequences)"
    )

    # Also save as PDF (better for publications)
    try:
        # UPGMA PDF
        fig, ax = plt.subplots(figsize=(12, max(6, num_seqs * 0.3)))
        for terminal in upgma_tree.get_terminals():
            terminal.name = simplify_name(terminal.name)
        Phylo.draw(upgma_tree, axes=ax, do_show=False, show_confidence=False)
        ax.set_title(f'{peptide_name} - UPGMA Tree', fontsize=14, fontweight='bold', pad=20)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_visible(False)
        plt.tight_layout()
        plt.savefig(f"{output_dir}/{peptide_name}_upgma_tree.pdf", bbox_inches='tight')
        plt.close()
        print(f"    PDF saved: {output_dir}/{peptide_name}_upgma_tree.pdf")
    except Exception as e:
        print(f"    PDF save failed: {e}")

    # Create distance matrix heatmap
    create_distance_heatmap(distance_matrix, peptide_name, output_dir, num_seqs)

    return True


def create_distance_heatmap(distance_matrix, peptide_name, output_dir, num_seqs):
    """
    Create a heatmap of the distance matrix.
    """
    import numpy as np
    import seaborn as sns

    # Extract matrix data
    names = [simplify_name(name) for name in distance_matrix.names]

    # Convert to numpy array
    matrix_data = np.zeros((len(names), len(names)))
    for i, name1 in enumerate(distance_matrix.names):
        for j, name2 in enumerate(distance_matrix.names):
            if i == j:
                matrix_data[i][j] = 0
            else:
                matrix_data[i][j] = distance_matrix[name1, name2]

    # Create heatmap
    fig_size = max(8, min(20, num_seqs * 0.3))
    fig, ax = plt.subplots(figsize=(fig_size, fig_size))

    sns.heatmap(matrix_data,
                xticklabels=names,
                yticklabels=names,
                cmap='YlOrRd',
                square=True,
                cbar_kws={'label': 'Distance'},
                ax=ax,
                annot=False if num_seqs > 30 else True,
                fmt='.2f' if num_seqs <= 30 else None,
                linewidths=0.5 if num_seqs <= 30 else 0)

    ax.set_title(f'{peptide_name} - Pairwise Distance Matrix',
                 fontsize=14, fontweight='bold', pad=20)

    # Rotate labels for readability
    plt.xticks(rotation=90, ha='right', fontsize=8)
    plt.yticks(rotation=0, fontsize=8)

    plt.tight_layout()
    plt.savefig(f"{output_dir}/{peptide_name}_distance_matrix.png", dpi=300, bbox_inches='tight')
    plt.close()

    print(f"    Distance heatmap saved: {output_dir}/{peptide_name}_distance_matrix.png")


def create_summary_document(output_dir, results):
    """
    Create a summary document about the trees.
    """
    with open(f"{output_dir}/TREES_README.txt", 'w') as f:
        f.write("PHYLOGENETIC TREES FOR PEPTIDE SEQUENCES\n")
        f.write("="*80 + "\n\n")

        f.write("This directory contains phylogenetic trees for each peptide sequence.\n\n")

        f.write("FILES GENERATED FOR EACH PEPTIDE:\n")
        f.write("-"*80 + "\n")
        f.write("1. {peptide}_upgma_tree.png - UPGMA tree visualization (PNG)\n")
        f.write("2. {peptide}_upgma_tree.pdf - UPGMA tree visualization (PDF)\n")
        f.write("3. {peptide}_nj_tree.png - Neighbor-Joining tree (PNG)\n")
        f.write("4. {peptide}_upgma.nwk - UPGMA tree in Newick format\n")
        f.write("5. {peptide}_nj.nwk - NJ tree in Newick format\n")
        f.write("6. {peptide}_distance_matrix.png - Pairwise distance heatmap\n")
        f.write("\n")

        f.write("TREE CONSTRUCTION METHODS:\n")
        f.write("-"*80 + "\n")
        f.write("UPGMA (Unweighted Pair Group Method with Arithmetic Mean):\n")
        f.write("  - Assumes constant evolution rate (molecular clock)\n")
        f.write("  - Good for closely related sequences\n")
        f.write("  - Creates ultrametric tree (all leaves same distance from root)\n")
        f.write("\n")
        f.write("Neighbor-Joining (NJ):\n")
        f.write("  - Does NOT assume molecular clock\n")
        f.write("  - More robust to rate variation\n")
        f.write("  - Generally preferred for phylogenetic analysis\n")
        f.write("\n")

        f.write("DISTANCE METRIC:\n")
        f.write("-"*80 + "\n")
        f.write("Identity-based distance:\n")
        f.write("  - Measures proportion of different amino acids\n")
        f.write("  - Distance = (# differences) / (alignment length)\n")
        f.write("  - Range: 0.0 (identical) to 1.0 (completely different)\n")
        f.write("\n")

        f.write("HOW TO USE THESE FILES:\n")
        f.write("-"*80 + "\n")
        f.write("Newick files (.nwk) can be opened in:\n")
        f.write("  - FigTree: http://tree.bio.ed.ac.uk/software/figtree/\n")
        f.write("  - iTOL: https://itol.embl.de/\n")
        f.write("  - MEGA: https://www.megasoftware.net/\n")
        f.write("  - Dendroscope\n")
        f.write("  - R (using ape, ggtree packages)\n")
        f.write("\n")
        f.write("Image files (.png, .pdf) can be:\n")
        f.write("  - Used directly in presentations/papers\n")
        f.write("  - Edited in Illustrator/Inkscape\n")
        f.write("\n")

        f.write("PEPTIDES ANALYZED:\n")
        f.write("-"*80 + "\n")
        for peptide_name, success, num_seqs in results:
            if success:
                f.write(f"{peptide_name:15} - {num_seqs:3} sequences - ✓ Trees generated\n")
            else:
                f.write(f"{peptide_name:15} - Not enough sequences for tree\n")
        f.write("\n")

        f.write("INTERPRETATION TIPS:\n")
        f.write("-"*80 + "\n")
        f.write("- Branch length = evolutionary distance (amino acid changes)\n")
        f.write("- Sequences that cluster together are more similar\n")
        f.write("- Look for species-specific clades (groups)\n")
        f.write("- Compare UPGMA vs NJ trees for consistency\n")
        f.write("- Reference sequence helps anchor the comparison\n")
        f.write("\n")

        f.write("NOTES:\n")
        f.write("-"*80 + "\n")
        f.write("- Trees are unrooted (no ancestral direction implied)\n")
        f.write("- Short sequences may have limited phylogenetic signal\n")
        f.write("- Distance-based methods used (appropriate for peptides)\n")
        f.write("- All sequences included (not just unique variants)\n")
        f.write("\n")


def main():
    import sys

    if len(sys.argv) < 2:
        print("Usage: python generate_peptide_trees.py <alignment_dir> [output_dir]")
        sys.exit(1)

    alignment_dir = sys.argv[1]
    output_dir = sys.argv[2] if len(sys.argv) > 2 else "peptide_trees"

    # Create output directory
    os.makedirs(output_dir, exist_ok=True)

    print(f"\n{'='*80}")
    print("GENERATING PHYLOGENETIC TREES FOR PEPTIDE SEQUENCES")
    print(f"{'='*80}\n")
    print(f"Input directory: {alignment_dir}/")
    print(f"Output directory: {output_dir}/\n")

    results = []

    for peptide_name in PEPTIDES.keys():
        fasta_file = f"{alignment_dir}/{peptide_name}_all_sequences.fasta"

        if not os.path.exists(fasta_file):
            print(f"WARNING: {fasta_file} not found, skipping...")
            results.append((peptide_name, False, 0))
            continue

        # Count sequences
        seqs = list(SeqIO.parse(fasta_file, 'fasta'))
        num_seqs = len(seqs)

        success = generate_tree_for_peptide(fasta_file, peptide_name, output_dir)
        results.append((peptide_name, success, num_seqs))

    # Create summary
    create_summary_document(output_dir, results)

    print(f"\n{'='*80}")
    print("COMPLETE!")
    print(f"{'='*80}\n")
    print(f"All tree files written to: {output_dir}/\n")
    print("Summary of results:")
    for peptide_name, success, num_seqs in results:
        status = "✓" if success else "✗"
        print(f"  {status} {peptide_name:15} ({num_seqs:3} sequences)")
    print()


if __name__ == "__main__":
    main()
