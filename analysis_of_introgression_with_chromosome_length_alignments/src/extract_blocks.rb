# m_matschiner Thu May 16 09:44:25 CEST 2019

# This script reads a chromosome-length alignment in fasta format
# and extracts non-overlapping alignment blocks of a specified size.
# These alignment blocks will be written to a specified directory.
# Optionally, a maximum proportion of missing data can be specified,
# alignments with a greater proportion of missing data will then 
# not be written to the output directory.

# This script should be run e.g. with
# ruby extract_blocks.rb align.fasta chromosome_alignment_dir 50000 0.8

# Include library FileUtils.
require 'fileutils'

# Get the command line arguments.
fasta_file_name = ARGV[0]
block_dir = ARGV[1]
block_size = ARGV[2].to_i
max_missing_proportion_raw = ARGV[3]
if ARGV[3] == nil
	max_missing_proportion = 1.0
else
	max_missing_proportion = max_missing_proportion_raw.to_f
end
FileUtils::mkdir_p block_dir

# Read the input file.
fasta_file = File.open(fasta_file_name)
print "Reading file #{fasta_file_name}..."
fasta_lines = fasta_file.readlines
puts " done."
fasta_ids = []
fasta_seqs = []
fasta_lines.each do |l|
	if l[0] == ">"
		fasta_ids << l[1..-1].strip
	elsif l.strip != ""
		fasta_seqs << l.strip
	end
end

# Determine the maximum id size.
max_id_size = 0
fasta_ids.each {|i| max_id_size = i.size if i.size > max_id_size}

# Make sure that all sequences have the same length.
fasta_seqs.each do |s|
	if s.size != fasta_seqs[0].size
		puts "ERROR: Not all sequences have the same size!"
		exit 1
	end
end

# Split the chromosome-length alignment it into blocks.
block_start = 0
block_end = block_start + block_size - 1
blocks_written = 0
blocks_excluded = 0
while block_end < fasta_seqs[0].size
	block_seqs = []
	fasta_seqs.each do |s|
		block_seqs << s[block_start..block_end]
	end
	block_start = block_start += block_size
	block_end = block_start + block_size - 1

	# Determine the proportion of missing data per block.
	write_block = true
	if max_missing_proportion < 1.0
		n_missing = 0
		block_seqs.each do |s|
			n_missing += s.count("N") + s.count("n") + s.count("?") + s.count("-")
		end
		write_block = false if n_missing > max_missing_proportion * block_seqs.size * block_size
	end

	# Prepare and write the block alignment in nexus format unless there is too much missing data.
	if write_block
		block_alignment_string = "#nexus\n"
		block_alignment_string << "begin data;\n"
		block_alignment_string << "dimensions  ntax=#{block_seqs.size} nchar=#{block_seqs[0].size};\n"
		block_alignment_string << "format datatype=DNA gap=- missing=?;\n"
		block_alignment_string << "matrix\n"
		fasta_ids.size.times do |x|
			block_alignment_string << "  #{fasta_ids[x].ljust(max_id_size+2)}"
			block_alignment_string << "#{block_seqs[x]}\n"
		end
		block_alignment_string << ";\n"
		block_alignment_string << "end;\n"

		# Write local alignment.
		n_digits = fasta_seqs[0].size.to_s.size
		block_start_string = (block_start+1).to_s.rjust(n_digits).gsub(" ","0")
		block_end_string = (block_end).to_s.rjust(n_digits).gsub(" ","0")
		block_alignment_file_name = "#{block_dir}/#{fasta_file_name.chomp(".fasta").chomp(".fa")}_#{block_start_string}_#{block_end_string}.nex"
		block_alignment_file = File.open(block_alignment_file_name,"w")
		block_alignment_file.write(block_alignment_string)
		blocks_written += 1
	else
		blocks_excluded += 1
	end
end

# Feedback.
puts "Wrote #{blocks_written} alignment blocks of #{block_size} bp."
if blocks_excluded > 0
	puts "Excluded #{blocks_excluded} blocks due to a proportion of missing data above #{max_missing_proportion}."
end
