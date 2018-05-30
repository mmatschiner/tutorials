# m_matschiner Sun May 27 23:46:11 CEST 2018

# This script uses input and output of the software saguaro to
# generate local alignments according to the segments inferred
# by saguaro. These alignments will be written in nexus format.

# This script should be run e.g. with
# ruby generate_alignments.rb LocalTrees.out chromosome_alignment_dir local_alignment_dir
# where 'chromosome_alignment_dir' should be replaced with the
# path to a directory in which the chromosome-length alignments
# that were used as input for saguaro are stored, and
# segment_alignment_dir should be replaced with the path to a
# directory, in which output alignments should be saved.
# Note that the chromosome alignment directory must contain an
# alignment file in fasta format for each chromosome that was
# used with saguaro.

# Include library FileUtils.
require 'fileutils'

# Check if command line arguments are provided, and print help
# text if they aren't.
if ARGV.size < 3 or ARGV.include?("-h") or ARGV.include?("--help")
	puts
	puts "generate_alignments.rb"
	puts
	puts "This script uses input and output of the software saguaro to"
	puts "generate local alignments according to the segments inferred"
	puts "by saguaro. These alignments will be written in nexus format."
	puts
	puts "This script should be run e.g. with"
	puts "ruby generate_alignments.rb LocalTrees.out chromosome_alignment_dir local_alignment_dir local_alignment_length"
	puts "where 'chromosome_alignment_dir' should be replaced with the"
	puts "path to a directory in which the chromosome-length alignments"
	puts "that were used as input for saguaro are stored, and"
	puts "local_alignment_dir should be replaced with the path to a"
	puts "directory, in which output alignments should be saved."
	puts "'local_alignment_length' should be replaced with the length"
	puts "(in bp) that all alignments should have."
	puts "Note that the chromosome alignment directory must contain an"
	puts "alignment file in fasta format for each chromosome that was"
	puts "used with saguaro."
	exit
end

# Get the command line arguments.
input_file_name = ARGV[0]
chromosome_alignment_dir = ARGV[1]
local_alignment_dir = ARGV[2]
local_alignment_length = ARGV[3].to_i
FileUtils::mkdir_p local_alignment_dir

# Read the input file.
input_file = File.open(input_file_name)
input_lines = input_file.readlines
filtered_lines = []
input_lines.each do |l|
	if l.include?("cactus") and l.include?("length") and l.include?("score")
		filtered_lines << l
	end
end

# Analyze the input.
lg_ids = []
local_alignment_ids = []
local_alignment_starts = []
local_alignment_ends = []
filtered_lines.each do |l|
	line_ary = l.split
	lg_id = line_ary[1].chomp(":")
	lg_ids << lg_id
	segment_start = line_ary[2].to_i
	segment_end = line_ary[4].to_i
	segment_center = (segment_start+segment_end)/2
	number_of_local_alignments_within_this_segment = (segment_end-segment_start+1)/local_alignment_length
	if number_of_local_alignments_within_this_segment > 0
		first_local_alignment_start = (segment_center - 0.5 * number_of_local_alignments_within_this_segment * local_alignment_length).to_i
		number_of_local_alignments_within_this_segment.times do |x|
			local_alignment_start = first_local_alignment_start + x * local_alignment_length
			local_alignment_end = first_local_alignment_start + (x+1) * local_alignment_length - 1
			local_alignent_id = "#{lg_id}_#{local_alignment_start.to_s.rjust(8).gsub(" ","0")}_#{local_alignment_end.to_s.rjust(8).gsub(" ","0")}"
			local_alignment_ids << local_alignent_id
			local_alignment_starts << local_alignment_start
			local_alignment_ends << local_alignment_end
		end
	end
end
uniq_lg_ids = lg_ids.uniq.sort

# For each of the unique linkage groups, read the chromosome-length
# alignment, and split it into smaller local alignments.
uniq_lg_ids.each do |lg|

	# Read the chromosome-length alignment file.
	chromosome_alignment_file_name = "#{chromosome_alignment_dir}/#{lg}.fasta"
	chromosome_alignment_file = File.open(chromosome_alignment_file_name,"r")
	chromosome_alignment_lines = chromosome_alignment_file.readlines
	chromosome_alignment_seq_ids = []
	chromosome_alignment_seq_ids_max_length = 0
	chromosome_alignment_seqs = []
	chromosome_alignment_lines.each do |l|
		if l[0] == ">"
			chromosome_alignment_seq_id = l[1..-1].strip
			chromosome_alignment_seq_ids << chromosome_alignment_seq_id
			chromosome_alignment_seq_ids_max_length = chromosome_alignment_seq_id.size if chromosome_alignment_seq_id.size > chromosome_alignment_seq_ids_max_length
		elsif l.strip != ""
			chromosome_alignment_seqs << l.strip
		end
	end

	# Produce and write local alignments.
	local_alignment_ids.size.times do |x|
		# Produce local alignment.
		local_alignment_seqs = []
		chromosome_alignment_seqs.each do |s|
			local_alignment_seqs << s[(local_alignment_starts[x]-1)..(local_alignment_ends[x]-1)]
		end

		# Prepare local alignment.
		local_alignment_string = "#nexus\n"
		local_alignment_string << "begin data;\n"
		local_alignment_string << "dimensions  ntax=#{local_alignment_seqs.size} nchar=#{local_alignment_seqs[0].size};\n"
		local_alignment_string << "format datatype=DNA gap=- missing=?;\n"
		local_alignment_string << "matrix\n"
		chromosome_alignment_seq_ids.size.times do |y|
			local_alignment_string << "  #{chromosome_alignment_seq_ids[y].ljust(chromosome_alignment_seq_ids_max_length+2)}"
			local_alignment_string << "#{local_alignment_seqs[y]}\n"
		end
		local_alignment_string << ";\n"
		local_alignment_string << "end;\n"

		# Write local alignment.
		local_alignment_file_name = "#{local_alignment_dir}/#{local_alignment_ids[x]}.nex"
		local_alignment_file = File.open(local_alignment_file_name,"w")
		local_alignment_file.write(local_alignment_string)
	end
end
