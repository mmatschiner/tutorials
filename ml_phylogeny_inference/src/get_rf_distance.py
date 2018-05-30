# m_matschiner Mon May 28 18:53:32 CEST 2018

# Import libraries.
import sys
from ete3 import Tree

# Get the command line arguments.
tree1_file = open(sys.argv[1])
tree2_file = open(sys.argv[2])

t1_str = tree1_file.read()
t2_str = tree2_file.read()

t1 = Tree(t1_str)
t2 = Tree(t2_str)
rf = t1.robinson_foulds(t2, unrooted_trees=True)[0]
print rf
