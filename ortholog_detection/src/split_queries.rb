# m_matschiner Mon May 14 11:53:12 CEST 2018

# Get the command-line arguments.
input_file_name = ARGV[0]
n_parts = ARGV[1].to_i

# Read the input file.
input_file = File.open(input_file_name)
lines = input_file.readlines
ids = []
seqs = []
lines.each do |l|
	if l[0] == ">"
		id = l[1..-1].strip
		ids << id
		seqs << ""
	elsif l.strip != ""
		seqs.last << l.strip
	end
end

# Prepare output strings.
outstrings = []
n_parts.times {outstrings << ""}
n_seqs_per_set = (ids.size/n_parts.to_f).ceil
ids.size.times do |x|
	outstring_index = x/n_seqs_per_set
	outstrings[outstring_index] << ">#{ids[x]}\n"
	outstrings[outstring_index] << "#{seqs[x]}\n"
	outstrings[outstring_index] << "\n"
end

# Write the output strings to output files.
output_file_names = []
n_parts.times do |x|
	output_file_names << input_file_name.chomp(".fasta").chomp(".fa").chomp(".fn") + ".#{(x+1).to_s.rjust(2).gsub(" ","0")}.fasta"
end
n_parts.times do |x|
	output_file = File.open(output_file_names[x],"w")
	output_file.write(outstrings[x])
end

# Feedback.
feedback_string = "Wrote #{n_parts} files named"
(output_file_names.size-1).times do |x|
	feedback_string << " #{output_file_names[x]},"
end
feedback_string << " and #{output_file_names.last}."
puts feedback_string