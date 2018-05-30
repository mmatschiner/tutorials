# Michael Matschiner, 2016-08-05

# This script calcualates the proportion of missing
# data in an alignment in nexus format.
#
# It should be called as
# ruby get_proportion_of_missing_data.rb alignment.nex

input_file_name = ARGV[0]
input_file = File.open(input_file_name)
input_lines = input_file.readlines

# Read the sequences and ids.
ids = []
seqs = []
in_matrix = false
input_lines.each do |l|
	if l.strip.downcase == "matrix"
		in_matrix = true
	elsif l.strip.downcase == "end;"
		in_matrix = false
	elsif in_matrix and l.strip != ";"
		line_ary = l.strip.split
		raise "ERROR: File could not be read properly!" if line_ary.size != 2
		ids << line_ary[0]
		seqs << line_ary[1]
	end
end

# Make sure all sequences have the same length.
seqs[1..-1].each do |s|
	raise "ERROR: Sequences have different lengths!" if s.size != seqs[0].size
end

# Get the number of parsimony-informative sites.
number_of_missing_bases = 0
seqs.each do |s|
	number_of_missing_bases += s.count("N")
	number_of_missing_bases += s.count("n")
	number_of_missing_bases += s.count("-")
	number_of_missing_bases += s.count("?")
end

# Output the number of parsimony-informative sites.
puts number_of_missing_bases/(seqs.size*seqs[0].size).to_f