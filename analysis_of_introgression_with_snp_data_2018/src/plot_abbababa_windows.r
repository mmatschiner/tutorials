# m_matschiner Wed May 23 19:16:56 CEST 2018

# Get the command-line arguments.
args <- commandArgs(trailingOnly = TRUE)
window_data_file_name <- args[1]
lg_lengths_file_name <- args[2]
plot_file_name <- args[3]

# Read the input files.
window_data = read.table(window_data_file_name, sep=",", header = T)
lg_length_data = read.table(lg_lengths_file_name, header = F)
lgs <- lg_length_data$V1
lg_lengths <- lg_length_data$V2

# Plot the weights as rectangles.
pdf(plot_file_name, height=7, width=7)
for(x in 1:length(lgs)) {
	lg <- as.character(lgs[x])
	lg_length <- as.numeric(lg_lengths[x])
	lg_data <- window_data[ which(window_data$scaffold==lg), ]
	if( nrow(lg_data) > 0 ){
		plot(lg_data$mid, lg_data$fd, type="l", ylim=c(0,1), xlab="Position", ylab="fd", main=lg)
	}
}
dev.off()