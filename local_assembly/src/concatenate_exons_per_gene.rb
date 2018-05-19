# m_matschiner Tue May 15 18:16:32 CEST 2018

# Get the command line arguments.
alignment_directory_in = ARGV[0].chomp("/")
alignment_directory_out = ARGV[1].chomp("/")
exon_info_file_name = ARGV[2]

# Create the output directory if it does not exist yet.
unless Dir.exists?(alignment_directory_out)
	Dir.mkdir(alignment_directory_out)
end

# Collect names of nucleotide fasta files in the input directory.
dir_entries_in = Dir.entries(alignment_directory_in).sort
filenames_in = []
dir_entries_in.each {|e| filenames_in << e if e.match(/.*_nucl.fasta/)}

# Get exon ids of alignments from the filenames.
exon_ids_with_alignments = []
filenames_in.each do |f|
	exon_ids_with_alignments << f.chomp("_nucl.fasta")
end

# Read the exons info file.
exons_info_exon_ids = []
exons_info_gene_ids = []
exon_info_file = File.open(exon_info_file_name)
exon_info_lines = exon_info_file.readlines
exon_info_lines.each do |l|
	line_ary = l.split
	exons_info_exon_ids << line_ary[0]
	exons_info_gene_ids << line_ary[1]
end
gene_ids = exons_info_gene_ids.uniq

# For each gene id, read all corresponding exon alignments.
gene_ids.each do |g|
	gene_ids = []
	gene_seqs = []
	exons_info_exon_ids.size.times do |x|
		if exons_info_gene_ids[x] == g
			exon_ids = []
			exon_seqs = []
			# Read this exon alignment file.
			exon_alignment_file_name = "#{alignment_directory_in}/#{exons_info_exon_ids[x]}_nucl.fasta"
			if File.exists?(exon_alignment_file_name)
				exon_alignment_file = File.open(exon_alignment_file_name)
				exon_alignment_lines = exon_alignment_file.readlines
				exon_alignment_lines.each do |l|
					if l[0] == ">"
						exon_ids << l[1..-1].gsub(/\[.+\]/,"").strip
						exon_seqs << ""
					elsif l.strip != ""
						exon_seqs.last << l.strip
					end
				end
				# Make sure that exon ids are identical among all files.
				if gene_ids == []
					gene_ids = exon_ids
				elsif gene_ids != exon_ids
					puts "ERROR: Exon IDs appear to differ between alignment files!"
					exit 1
				end
				# Concatenate exon sequences.
				if gene_seqs == []
					gene_seqs = exon_seqs
				else
					exon_seqs.size.times do |y|
						gene_seqs[y] << exon_seqs[y]
					end
				end
			end
		end
	end

	unless gene_ids == []
		# Prepare an output string for the concatenated alignment for this gene.
		gene_alignment_string = ""
		gene_ids.size.times do |x|
			gene_alignment_string << ">#{gene_ids[x]}\n"
			gene_alignment_string << "#{gene_seqs[x]}\n"
		end

		# Write the output string.
		gene_alignment_file_name = "#{alignment_directory_out}/#{g}.fasta"
		gene_alignment_file = File.open(gene_alignment_file_name,"w")
		gene_alignment_file.write(gene_alignment_string)

		# Feedback.
		puts "Wrote file #{gene_alignment_file_name}."
	end
end
