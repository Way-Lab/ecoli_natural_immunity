
from Bio import Phylo

tree = Phylo.read('ompA_tree_strains.nwk', 'newick')
Phylo.draw_ascii(tree)
