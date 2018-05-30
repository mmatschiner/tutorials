# m_matschiner Tue May 29 14:24:10 CEST 2018

# Load the ape library.
library(ape)

# Define a function to print the node ages.
print_mrca_table <- function(tree){
	ancestors = mrca(tree)
	cat("\n")
	cat("#ancestors\n")
	print(ancestors)
	node_ages = branching.times(tree)
	cat("\n")
	cat("#node_ages\n")
	print(node_ages)	
}

# Get the command-line arguments.
args <- commandArgs(trailingOnly = TRUE)
gene_tree_file_name <- args[1]

# Read the trees from the gene-tree file.
trees <- read.nexus(gene_tree_file_name)

# Set the size of the output wide enough for the table.
options(width=1000)

# Apply the function to print the node ages of all trees.
lapply(trees, print_mrca_table)
