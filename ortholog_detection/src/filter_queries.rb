# m_matschiner Mon May 14 01:43:51 CEST 2018

# Get the command-line arguments.
input_file_name = ARGV[0]
output_file_name = ARGV[1]
minimum_sequence_length = ARGV[2].to_i
minimum_bitscore_threshold = ARGV[3].to_i

# Read the input file.
input_file = File.open(input_file_name)
lines = input_file.readlines
ids = []
seqs = []
bitscores = []
lines.each do |l|
	if l[0] == ">"
		id = l[1..-1].strip
		bitscore = 0
		if id.match(/\[&bitscore=([\d\.]+)\]/)
			bitscore = $1.to_f
		end
		ids << id
		bitscores << bitscore
		seqs << ""
	elsif l.strip != ""
		seqs.last << l.strip
	end
end

# Prepare an output string for a filtered set of sequences.
outstring = ""
seq_count = 0
ids.size.times do |x|
	if seqs[x].size >= minimum_sequence_length and bitscores[x] >= minimum_bitscore_threshold
		outstring << ">#{ids[x]}\n"
		outstring << "#{seqs[x]}\n"
		outstring << "\n"
		seq_count += 1
	end
end

# Write the output string to the output file.
output_file = File.open(output_file_name,"w")
output_file.write(outstring)

# Feedback.
puts "Removed #{ids.size-seq_count} out of #{ids.size} sequences."