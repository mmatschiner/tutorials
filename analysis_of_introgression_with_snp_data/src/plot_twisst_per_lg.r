# m_matschiner Thu May 24 00:32:49 CEST 2018

# Include simon martin's plot_twisst.R script.
source ("twisst/plot_twisst.R")

# Define a colour scheme (using http://ethanschoonover.com/solarized).
cols = c("#2aa198", "#d33682", "#b58900") 

# Get the command-line arguments.
args <- commandArgs(trailingOnly = TRUE)
weights_file_name <- args[1]
window_data_file_name <- args[2]
lg_name <- args[3]
lg_length <- as.numeric(args[4])
rect_plot_file_name <- args[5]
smooth_plot_file_name <- args[6]

# Read the input files.
weights = read.table(weights_file_name, header = T)
window_data = read.table(window_data_file_name, header = T)

# Normalize weights.
weights = weights / apply(weights, 1, sum)

# Get the subset of sites without na values.
good_rows = which(is.na(apply(weights, 1, sum)) == F)
weights <- weights[good_rows,]
window_data <- window_data[good_rows,]

# Smooth the weights.
weights_smooth = smooth.weights(window_position=window_data$mid, weights_dataframe=weights, span=0.03, window_sites=window_data$sites)

# Plot the smoothed weights.
pdf(smooth_plot_file_name, height=7, width=7)
plot.weights(weights_dataframe=weights_smooth, positions=window_data$mid, line_cols=cols, fill_cols=cols, main=lg_name, xlim =c(0,lg_length), stacked=TRUE)
dev.off()

# Plot the weights as rectangles.
pdf(rect_plot_file_name, height=7, width=7)
plot(c(0,lg_length), c(0,1), type="n", main=lg_name, xlab="Position", ylab="Weights")
rect(window_data$start, rep(0,length(window_data$start)), window_data$end, weights$topo1, col="#2aa198", border="NA")
rect(window_data$start, weights$topo1, window_data$end, weights$topo1+weights$topo2, col="#d33682", border="NA")
rect(window_data$start, weights$topo1+weights$topo2, window_data$end, rep(1,length(window_data$start)), col="#b58900", border="NA")
dev.off()