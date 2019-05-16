# Michael Matschiner, 2016-07-07

# This script calcualates the number of variable
# sites in an alignment in nexus format.
#
# It should be called as
# ruby get_number_of_variable_sites.rb alignment.nex

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
	elsif in_matrix and l.strip != ";" and l.strip != ""
		line_ary = l.strip.split
		if line_ary.size != 2
			puts "ERROR: File could not be read properly!"
			puts line_ary
			exit 1
		end
		ids << line_ary[0]
		seqs << line_ary[1]
	end
end

# Make sure all sequences have the same length.
seqs[1..-1].each do |s|
	raise "ERROR: Sequences have different lengths!" if s.size != seqs[0].size
end

# Get the number of variable sites.
variable_sites = 0
seqs[0].size.times do |x|
	alleles_at_this_site = []
	seqs.each do |s|
		alleles_at_this_site << s[x].upcase if ["A","C","G","T"].include?(s[x].upcase)
	end
	uniq_alleles_at_this_site = alleles_at_this_site.uniq
	if uniq_alleles_at_this_site.size > 1
		variable_sites += 1
	end
end

# Output the number of parsimony-informative sites.
puts variable_sites