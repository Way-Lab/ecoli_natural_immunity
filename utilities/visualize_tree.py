
from Bio import Phylo

tree = Phylo.read('ompA_tree.nwk', 'newick')
Phylo.draw_ascii(tree)
