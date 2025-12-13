#!/usr/bin/env python3

import sys
import matplotlib.pyplot as plt
from Bio import Phylo
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend

def plot_cladogram(tree_file, output_file="OmpA_cladogram.pdf", title="OmpA Phylogenetic Tree"):
    """
    Create a cladogram visualization from a Newick tree file.
    """
    # Read the tree
    tree = Phylo.read(tree_file, "newick")

    # Create figure with appropriate size based on number of taxa
    num_taxa = len(tree.get_terminals())
    fig_height = max(10, num_taxa * 0.3)  # Scale height with number of sequences
    fig_width = 12

    fig = plt.figure(figsize=(fig_width, fig_height))
    ax = fig.add_subplot(111)

    # Draw the tree as a cladogram
    Phylo.draw(tree, axes=ax, do_show=False, show_confidence=False)

    # Customize the plot
    ax.set_title(title, fontsize=14, fontweight='bold', pad=20)
    ax.set_xlabel("", fontsize=0)  # Remove x-label for cleaner look

    # Adjust layout to prevent label cutoff
    plt.tight_layout()

    # Save the figure
    plt.savefig(output_file, format='pdf', dpi=300, bbox_inches='tight')
    print(f"Cladogram saved to {output_file}")

    # Also save as PNG for easier viewing
    png_file = output_file.replace('.pdf', '.png')
    plt.savefig(png_file, format='png', dpi=150, bbox_inches='tight')
    print(f"Cladogram also saved to {png_file}")

    # Print some tree statistics
    print(f"\nTree statistics:")
    print(f"Number of taxa: {num_taxa}")
    print(f"Total branch length: {tree.total_branch_length():.6f}")

    # Find and highlight our strains of interest
    strains_of_interest = ['RS218', 'SCB60', 'SCB61', 'EcN']
    print(f"\nStrains of interest:")
    for terminal in tree.get_terminals():
        if any(strain in terminal.name for strain in strains_of_interest):
            print(f"  - {terminal.name}")

    plt.close()

if __name__ == "__main__":
    tree_file = "OmpA_tree.nwk"
    output_file = "OmpA_cladogram.pdf"

    if len(sys.argv) > 1:
        tree_file = sys.argv[1]
    if len(sys.argv) > 2:
        output_file = sys.argv[2]

    plot_cladogram(tree_file, output_file)