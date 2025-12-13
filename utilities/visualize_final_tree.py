
from Bio import Phylo

tree = Phylo.read('ompA_tree_final.nwk', 'newick')
Phylo.draw_ascii(tree)
