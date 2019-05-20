# m_matschiner Mon Apr 16 12:47:43 CEST 2018

# Define functions to calculate the numbers of abba and baba patterns.
bbaa = function(p1, p2, p3, p4) p1 * p2 * (1 - p3) * (1 - p4)
abba = function(p1, p2, p3, p4) (1 - p1) * p2 * p3 * (1 - p4)
baba = function(p1, p2, p3, p4) p1 * (1 - p2) * p3 * (1 - p4)
D.stat = function(dataframe) (sum(dataframe$ABBA) - sum(dataframe$BABA)) / (sum(dataframe$ABBA) + sum(dataframe$BABA))
fd.stat = function(p1, p2, p3, p4) {
    pd = pmax(p2, p3)
    (sum(abba(p1, p2, p3, p4)) - sum(baba(p1, p2, p3, p4))) / (sum(abba(p1, pd, pd, p4)) - sum(baba(p1, pd, pd, p4)))
}

# Define functions for jackkniving.
get_genome_blocks <- function(block_size, lg_lengths) {
    block_starts <- sapply(lg_lengths, function(l) seq(1, l, block_size))
    data.frame(start = unlist(block_starts),
               end = unlist(block_starts) + block_size - 1,
               lg = rep(names(block_starts), sapply(block_starts, length)))
    }
get_genome_jackknife_indices <- function(lg, position, block_info){
    lapply(1:nrow(block_info), function(x) !(lg == block_info$lg[x] &
                                             position >= block_info$start[x] &
                                             position <= block_info$end[x]))
    }
get_jackknife_sd <- function(FUN, input_dataframe, jackknife_indices){
    n_blocks <- length(jackknife_indices)
    overall_mean <- FUN(input_dataframe)
    sd(sapply(1:n_blocks, function(i) overall_mean*n_blocks - FUN(input_dataframe[jackknife_indices[[i]],])*(n_blocks-1)))
    }

# Get the command-line arguments.
args <- commandArgs(trailingOnly = TRUE)
allele_freqs_file_name <- args[1]
output_file_name <- args[2]
spc_p1 <- args[3]
spc_p2 <- args[4]
spc_p3 <- args[5]
spc_o <- args[6]
lg_length_table_file_name <- args[7]

# Read the allele-frequencies table.
freq_table = read.table(allele_freqs_file_name, header=T, as.is=T)

# Open the output file.
output_file <- file(output_file_name, "w")

# Output.
write(paste("Species 1: ", spc_p1, sep=""), output_file, append=T)
write(paste("Species 2: ", spc_p2, sep=""), output_file, append=T)
write(paste("Species 3: ", spc_p3, sep=""), output_file, append=T)
write(paste("Species O: ", spc_o, sep=""), output_file, append=T)
write("", output_file, append=T)
write(paste("Number of sites: ", nrow(freq_table), sep=""), output_file, append=T)

# Get the allele frequencies.
p1 = freq_table[,spc_p1]
p2 = freq_table[,spc_p2]
p3 = freq_table[,spc_p3]
p4 = freq_table[,spc_o]

# Calculate the number of abba and baba patterns.
BBAA = bbaa(p1, p2, p3, p4)
ABBA = abba(p1,	p2, p3,	p4)
BABA = baba(p1,	p2, p3,	p4)

# Calculate the d statistic.
ABBA_BABA_df = as.data.frame(cbind(ABBA,BABA))
D = D.stat(ABBA_BABA_df)

# Calculate Simon Martin's fd statistic.
fd = fd.stat(p1, p2, p3, p4)

# Output.
write(paste("Number of BBAA sites: ", sum(BBAA), sep=""), output_file, append=T)
write(paste("Number of ABBA sites: ", sum(ABBA), sep=""), output_file, append=T)
write(paste("Number of BABA sites: ", sum(BABA), sep=""), output_file, append=T)
write(paste("D statistic: ", D, sep=""), output_file, append=T)
write(paste("fd statistic: ", fd, sep=""), output_file, append=T)
if( sum(ABBA) > sum(BBAA)){
    cat(paste("\nWARNING: The number of ABBA sites (", sum(ABBA) , ") is greater than the number of BBAA sites (", sum(BBAA) , "), indicating that ", spc_p2, " and ", spc_p3, " are more closely related than ", spc_p1, " and ", spc_p2, ". You should swap ", spc_p1, " and ", spc_p3, " to get the correct D-statistic.\n\n", sep=""))
}

# Read the lg length table
lg_table = read.table(lg_length_table_file_name)
lg_lengths = lg_table[,2]
names(lg_lengths) = lg_table[,1]

# Run jackknife iterations to estimate the standard deviation in d.
blocks = get_genome_blocks(block_size=1e6, lg_lengths=lg_lengths)
n_blocks = nrow(blocks)
indices = get_genome_jackknife_indices(lg=freq_table$scaffold,position=freq_table$position,block_info=blocks)
D_sd = get_jackknife_sd(FUN=D.stat, input_dataframe=as.data.frame(cbind(ABBA,BABA)),jackknife_indices=indices)

# Calculate the p-value for d.
D_err <- D_sd/sqrt(n_blocks)
D_Z <- D / D_err
D_p <- 2*pnorm(-abs(D_Z))

# Feedback.
write(paste("p-value: ", D_p, sep=""), output_file, append=T)

# Close the output file.
close(output_file)

# Feedback.
cat(paste("\nWrote results to file ", output_file_name, ".\n\n", sep=""))
