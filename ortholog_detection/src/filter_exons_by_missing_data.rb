# m_matschiner Tue May 15 15:52:44 CEST 2018

# Get the command line arguments.
alignment_directory_in = ARGV[0].chomp("/")
alignment_directory_out = ARGV[1].chomp("/")
maximum_number_of_missing_sequences = ARGV[2].to_i
minimum_alignment_length = ARGV[3].to_i

# Create the output directory if it does not exist yet.
unless Dir.exists?(alignment_directory_out)
	Dir.mkdir(alignment_directory_out)
end

# Collect names of nucleotide fasta files in the input directory.
dir_entries_in = Dir.entries(alignment_directory_in)
filenames_in = []
dir_entries_in.each {|e| filenames_in << e if e.match(/.*\.fasta/)}

# Do for each fasta file in the input directory.
n_removed_exons = 0
filenames_in.each do |f|

	# Read the fasta file.
	fasta_file = File.open("#{alignment_directory_in}/#{f}")
	fasta_lines = fasta_file.readlines
	fasta_ids = []
	fasta_bitscores = []
	fasta_hits = []
	fasta_seqs = []
	fasta_lines.each do |l|
		if l[0] == ">"
			fasta_ids << l[1..-1].strip
			fasta_seqs << ""
		else
			fasta_seqs.last << l.strip
		end
	end

	# Count the number of missing sequences.
	n_missing_seqs = 0
	fasta_ids.size.times do |x|
		if fasta_seqs[x].match(/^-+$/)
			n_missing_seqs += 1
		end
	end

	# See whether the alignment matches criteria for completeness and length.
	if n_missing_seqs <= maximum_number_of_missing_sequences and fasta_seqs[0].size >= minimum_alignment_length
	
		# Prepare the string for a new fasta file.
		new_fasta_string = ""
		fasta_ids.size.times do |x|
			new_fasta_string << ">#{fasta_ids[x]}\n"
			new_fasta_string << "#{fasta_seqs[x]}\n"
		end

		# Write the new fasta file.
		new_fasta_file = File.open("#{alignment_directory_out}/#{f}","w")
		new_fasta_file.write(new_fasta_string)
		
	else

		n_removed_exons += 1

	end

end

# Feedback.
puts "Removed #{n_removed_exons} exons."
