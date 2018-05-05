# m_matschiner Sat May 5 17:37:00 CEST 2018

# Load the ape and phangorn libraries.
library("ape")
library("phangorn")

# Get the command-line arguments.
args <- commandArgs(trailingOnly = TRUE)
tree1_file_name <- args[1]
tree2_file_name <- args[2]

# Read the two trees.
tree1 <- read.tree(tree1_file_name)
tree2 <- read.tree(tree2_file_name)

# Calculate the Robinson-Foulds distance.
cat("\nRobinson-Foulds distance:", RF.dist(tree1, tree2), "\n\n")