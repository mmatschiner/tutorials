# m_matschiner Thu May 24 15:34:45 CEST 2018

# Get the command-line arguments.
vcf_input_file_name = ARGV[0]
ref_seq_file_name = ARGV[1]
mask_file_name = ARGV[2]
outfile_name = ARGV[3]

# Read the reference sequence file.
ref_seq_file = File.open(ref_seq_file_name)
ref_lines = ref_seq_file.readlines
ref_ids = []
ref_seqs = []
ref_lines.each do |l|
	if l[0] == ">"
		ref_ids << l[1..-1].strip.split(".")[0]
		ref_seqs << ""
	elsif l.strip != ""
		ref_seqs.last << l.strip
	end
end
if ref_ids.size > 1
	puts "ERROR: Expected a single sequence in file #{ref_seq_file_name} but found #{ref_ids.size}!"
	exit 1
end
ref_id = ref_ids[0]
ref_seq = ref_seqs[0]

# Start an array for the output sequence, and fill it with the reference sequence for positions that are not masked.
out_seq_ary = []
ref_seq.size.times do |x|
	out_seq_ary << ref_seq[x]
end

# Mask the output sequence using the input mask file.
mask_file = File.open(mask_file_name)
mask_lines = mask_file.readlines
first_mask_line_ary = mask_lines[0].split
if first_mask_line_ary.size != 3
	puts "ERROR: Expected three columns in file #{mask_file_name} but found #{line_ary.size} columns!"
	exit 1
elsif first_mask_line_ary[1].match(/[a-zA-Z]+/)
	mask_lines = mask_lines[1..-1]
end
mask_lines.each do |l|
	line_ary = l.split
	from = line_ary[1].to_i
	to = line_ary[2].to_i
	(from-1).upto(to-1) do |x|
		out_seq_ary[x] = "N"
	end
end

# Duplicate the outseq array.
out_seq_ary1 = out_seq_ary
out_seq_ary2 = out_seq_ary.dup

# Read the vcf input file and update the output sequence according to SNP information.
vcf_input_file = File.open(vcf_input_file_name)
vcf_lines = vcf_input_file.readlines
if vcf_lines[0][0..15] != "##fileformat=VCF"
	puts "ERROR: The input file #{vcf_input_file_name} does not seem to be in VCF file format!"
	exit 1
end
sample_id = nil
vcf_lines.each do |l|
	if l[0..1] == "##"
		next
	elsif l[0] == "#"
		sample_ids = l.split[9..-1]
		if sample_ids.size > 1
			puts "ERROR: Expected a single sample ID in file #{vcf_input_file_name} but found #{l[9..-1].size} sample IDs in the header line!"
			exit 1
		end
		sample_id = sample_ids[0]
	else
		chr = l.split[0]
		if chr != ref_id
			puts "ERROR: Expected only records for chromosome #{ref_id} in file #{vcf_input_file_name}, but found #{chr}!"
			exit 1
		end
		pos = l.split[1].to_i
		ref = l.split[3]
		alts = l.split[4].split(",")
		gt = l.split[9]
		alleles = []
		if gt.include?("|")
			alleles = gt.split("|")
		elsif gt.include?("/")
			puts "ERROR: Expected phased genotypes but found genotype #{gt}!"
			exit 1
		else
			alleles = [".", "."]
		end
		if alleles.size != 2
			puts "ERROR: Expected two alleles per genotype but found #{alleles.size} alleles!"
			exit 1
		end
		if alleles[0] != "0"
			if alleles[0] == "."
				out_seq_ary1[pos-1] = "N"
			else
				if alleles[0].to_i > alts.size
					puts "ERROR: Expected #{alts.size} alternate alleles at pos #{pos} but found #{alleles[0].to_i}!"
					exit 1
				end
				out_seq_ary1[pos-1] = alts[alleles[0].to_i-1]
			end
		end
		if alleles[1] != "0"
			if alleles[1] == "."
				out_seq_ary2[pos-1] = "N"
			else
				if alleles[1].to_i > alts.size
					puts "ERROR: Expected #{alts.size} alternate alleles at pos #{pos} but found #{alleles[1].to_i}!"
					exit 1
				end
				out_seq_ary2[pos-1] = alts[alleles[1].to_i-1]
			end
		end
	end
end

# Prepare the output string in fasta format.
outstring = ">#{sample_id}_A\n"
out_seq_ary1.each do |i|
	outstring << i
end
outstring << "\n"
outstring << ">#{sample_id}_B\n"
out_seq_ary2.each do |i|
	outstring << i
end
outstring << "\n"

# Write the output string to a fasta file.
outfile = File.open(outfile_name,"w")
outfile.write(outstring)
