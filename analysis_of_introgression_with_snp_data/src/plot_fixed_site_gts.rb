# Michael Matschiner, 2016-03-09

# This script reads the output of script 
# get_fixed_site_gts.rb
# and produces plots for the ancestry of hybrid
# individuals.

# Read the input file.
infile_name = ARGV[0]
infile = File.open(infile_name)
lines = infile.readlines
outfile_name = ARGV[1]
required_completeness = ARGV[2].to_f
thinning_distance = ARGV[3].to_i

# Define the ids of all specimens.
specimen_ids = lines[0].split[4..-1]
n_specimens = specimen_ids.size

# Prepare an empty matrix as the basis for the plots.
specimen_allele1 = []
specimen_allele2 = []
n_specimens.times do |x|
	specimen_allele1[x] = []
	specimen_allele2[x] = []
end

# Analyse the presence of parent alleles in all specimens.
scaffolds_for_all_selected_sites = []
positions_for_all_selected_sites = []
new_scaffold = []
previous_scaffold = nil
previous_position = nil
lines[1..-1].each do |l|
	line_ary = l.split
	scaffold = line_ary[0]
	position = line_ary[1].to_i
	parent1_allele = line_ary[2]
	parent2_allele = line_ary[3]
	alleles_this_pos = line_ary[4..-1]
	# Make sure the number of alleles at this position is identical to that of specimen ids.
	unless alleles_this_pos.size == n_specimens
		puts "ERROR: The number of alleles found at position #{position} does not match the number of specimens #{specimens}!"
		exit
	end
	# Check if the completeness at this position is sufficient.
	if alleles_this_pos.count("./.") + alleles_this_pos.count(".|.") <= alleles_this_pos.size*(1-required_completeness)
		# Check if this position is sufficiently distant to the previous position.
		if previous_position == nil or position < previous_position or position > previous_position + thinning_distance
			alleles_this_pos.size.times do |x|
				if alleles_this_pos[x].include?("/")
					alleles_this_pos_this_specimen = alleles_this_pos[x].split("/")
				elsif alleles_this_pos[x].include?("|")
					alleles_this_pos_this_specimen = alleles_this_pos[x].split("|")
				else
					puts "ERROR: Expected alleles to be separated eiter by '/' or '|' but found #{alleles_this_pos[x]}!"
					exit
				end
				if alleles_this_pos_this_specimen[0] == parent1_allele and alleles_this_pos_this_specimen[1] == parent1_allele
					specimen_allele1[x] << "p1"
					specimen_allele2[x] << "p1"
				elsif alleles_this_pos_this_specimen[0] == parent1_allele and alleles_this_pos_this_specimen[1] == parent2_allele
					specimen_allele1[x] << "p2"
					specimen_allele2[x] << "p1"
				elsif alleles_this_pos_this_specimen[0] == parent1_allele and alleles_this_pos_this_specimen[1] == "."
					specimen_allele1[x] << "p1"
					specimen_allele2[x] << "m"
				elsif alleles_this_pos_this_specimen[0] == parent2_allele and alleles_this_pos_this_specimen[1] == parent1_allele
					specimen_allele1[x] << "p2"
					specimen_allele2[x] << "p1"
				elsif alleles_this_pos_this_specimen[0] == parent2_allele and alleles_this_pos_this_specimen[1] == parent2_allele
					specimen_allele1[x] << "p2"
					specimen_allele2[x] << "p2"
				elsif alleles_this_pos_this_specimen[0] == parent2_allele and alleles_this_pos_this_specimen[1] == "."
					specimen_allele1[x] << "p2"
					specimen_allele2[x] << "m"
				elsif alleles_this_pos_this_specimen[0] == "." and alleles_this_pos_this_specimen[1] == parent1_allele
					specimen_allele1[x] << "p1"
					specimen_allele2[x] << "m"
				elsif alleles_this_pos_this_specimen[0] == "." and alleles_this_pos_this_specimen[1] == parent2_allele
					specimen_allele1[x] << "p2"
					specimen_allele2[x] << "m"
				elsif alleles_this_pos_this_specimen[0] == "." and alleles_this_pos_this_specimen[1] == "."
					specimen_allele1[x] << "m"
					specimen_allele2[x] << "m"
				else
					puts "WARNING: Expected either one of the parental alleles #{parent1_allele} and #{parent2_allele} or missing data coded with '.' but found #{alleles_this_pos_this_specimen[0]} and #{alleles_this_pos_this_specimen[1]}."
					specimen_allele1[x] << "m"
					specimen_allele2[x] << "m"
				end
			end
			if previous_scaffold == nil
				new_scaffold << false
			else
				if scaffold == previous_scaffold
					new_scaffold << false
				else
					new_scaffold << true
				end
			end
			previous_scaffold = scaffold
			previous_position = position
			scaffolds_for_all_selected_sites << scaffold
			positions_for_all_selected_sites << position
		end
	end
end

# Feedback, report number of homo- and heterozygous sites per specimen.
output_string = "\n"
output_string << "specimen".ljust(20)
output_string << "n_hom(p1)".rjust(12)
output_string << "n_het".rjust(12)
output_string << "n_hom(p2)".rjust(12)
output_string << "\n"
n_specimens.times do |x|
	n_hom_p1 = 0
	n_het = 0
	n_hom_p2 = 0
	specimen_allele1[x].size.times do |pos|
		if specimen_allele1[x][pos] == "p1" and specimen_allele2[x][pos] == "p1"
			n_hom_p1 += 1
		elsif specimen_allele1[x][pos] == "p2" and specimen_allele2[x][pos] == "p2"
			n_hom_p2 += 1
		elsif specimen_allele1[x][pos] == "p1" and specimen_allele2[x][pos] == "p2"
			n_het += 1
		elsif specimen_allele1[x][pos] == "p2" and specimen_allele2[x][pos] == "p1"
			n_het += 1
		end
	end
	output_string << "#{specimen_ids[x].ljust(20)}#{n_hom_p1.to_s.rjust(12)}#{n_het.to_s.rjust(12)}#{n_hom_p2.to_s.rjust(12)}\n"
end
output_string << "\n"
puts output_string

# Get the number of sites included in the plot.
n_sites = specimen_allele1[0].size

# Feedback.
puts "Found #{n_sites} positions with the required completeness."

# Some specifications for the SVG string.
top_margin = 10
bottom_margin = 10
left_margin = 10
left_text_block_width = 50
right_margin = 10
dim_x = 800
dim_y = 600
line_color = "5b5b5b"
color1 = "ef2746"
color2 = "27abd0"
font_size = 6
font_correction_y = 1.5
height = ((dim_y - top_margin - bottom_margin) / n_specimens.to_f) / 2.0
width = (dim_x - left_margin - left_text_block_width - right_margin) / n_sites.to_f

# Build the SVG header.
svgstring = "<?xml version=\"1.0\" standalone=\"no\"?>\n"
svgstring << "<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.0//EN\" \"http://www.w3.org/TR/2001/REC-SVG-20010904/DTD/svg10.dtd\">\n"
svgstring << "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"#{dim_x}\" height=\"#{dim_y}\" viewBox=\"0 0 #{dim_x} #{dim_y}\" xmlns:xlink=\"htttp://www.w3.org/1999/xlink\">\n\n"
svgstring << "\n"

# Draw each individual block.
n_specimens.times do |y|
	top_left_y = top_margin + y * (2 * height)
	# Draw blocks for the first alleles of this specimen.
	this_allele_first_site = 0
	top_left_x = left_margin + left_text_block_width
	n_sites.times do |x|
		if x < 1
			prev_allele = nil
		else
			prev_allele = specimen_allele1[y][x-1]
		end
		this_allele = specimen_allele1[y][x]
		if x > n_sites-2
			next_allele = nil
		else
			next_allele = specimen_allele1[y][x+1]
		end
		this_allele_is_a_first = false
		this_allele_is_a_last = false
		this_allele_is_a_first = true if prev_allele == nil or prev_allele != this_allele
		this_allele_is_a_last = true if next_allele == nil or this_allele != next_allele
		if this_allele_is_a_first and this_allele_is_a_last
			this_allele_first_site = x
			top_left_x = left_margin + left_text_block_width + this_allele_first_site * width
			this_allele_last_site = x
			top_right_x = left_margin + left_text_block_width + (this_allele_last_site+1) * width
			if this_allele == "p1"
				svgstring << "<rect fill=\"##{color1}\" x=\"#{top_left_x.round(3)}\" y=\"#{top_left_y.round(4)}\" width=\"#{(top_right_x-top_left_x).round(4)}\" height=\"#{height.round(3)}\" />\n"
			elsif this_allele == "p2"
				svgstring << "<rect fill=\"##{color2}\" x=\"#{top_left_x.round(3)}\" y=\"#{top_left_y.round(4)}\" width=\"#{(top_right_x-top_left_x).round(4)}\" height=\"#{height.round(3)}\" />\n"
			end
		elsif this_allele_is_a_first
			this_allele_first_site = x
			top_left_x = left_margin + left_text_block_width + this_allele_first_site * width
		elsif this_allele_is_a_last
			this_allele_last_site = x
			top_right_x = left_margin + left_text_block_width + (this_allele_last_site+1) * width
			if this_allele == "p1"
				svgstring << "<rect fill=\"##{color1}\" x=\"#{top_left_x.round(3)}\" y=\"#{top_left_y.round(4)}\" width=\"#{(top_right_x-top_left_x).round(4)}\" height=\"#{height.round(3)}\" />\n"
			elsif this_allele == "p2"
				svgstring << "<rect fill=\"##{color2}\" x=\"#{top_left_x.round(3)}\" y=\"#{top_left_y.round(4)}\" width=\"#{(top_right_x-top_left_x).round(4)}\" height=\"#{height.round(3)}\" />\n"
			end
		end
	end
	# Draw blocks for the second alleles of this specimen.
	this_allele_first_site = 0
	top_left_x = left_margin + left_text_block_width
	n_sites.times do |x|
		if x < 1
			prev_allele = nil
		else
			prev_allele = specimen_allele2[y][x-1]
		end
		this_allele = specimen_allele2[y][x]
		if x > n_sites-2
			next_allele = nil
		else
			next_allele = specimen_allele2[y][x+1]
		end
		this_allele_is_a_first = false
		this_allele_is_a_last = false
		this_allele_is_a_first = true if prev_allele == nil or prev_allele != this_allele
		this_allele_is_a_last = true if next_allele == nil or this_allele != next_allele
		if this_allele_is_a_first and this_allele_is_a_last
			this_allele_first_site = x
			top_left_x = left_margin + left_text_block_width + this_allele_first_site * width
			this_allele_last_site = x
			top_right_x = left_margin + left_text_block_width + (this_allele_last_site+1) * width
			if this_allele == "p1"
				svgstring << "<rect fill=\"##{color1}\" x=\"#{top_left_x.round(3)}\" y=\"#{(top_left_y+height).round(4)}\" width=\"#{(top_right_x-top_left_x).round(4)}\" height=\"#{height.round(3)}\" />\n"
			elsif this_allele == "p2"
				svgstring << "<rect fill=\"##{color2}\" x=\"#{top_left_x.round(3)}\" y=\"#{(top_left_y+height).round(4)}\" width=\"#{(top_right_x-top_left_x).round(4)}\" height=\"#{height.round(3)}\" />\n"
			end
		elsif this_allele_is_a_first
			this_allele_first_site = x
			top_left_x = left_margin + left_text_block_width + this_allele_first_site * width
		elsif this_allele_is_a_last
			this_allele_last_site = x
			top_right_x = left_margin + left_text_block_width + (this_allele_last_site+1) * width
			if this_allele == "p1"
				svgstring << "<rect fill=\"##{color1}\" x=\"#{top_left_x.round(3)}\" y=\"#{(top_left_y+height).round(4)}\" width=\"#{(top_right_x-top_left_x).round(4)}\" height=\"#{height.round(3)}\" />\n"
			elsif this_allele == "p2"
				svgstring << "<rect fill=\"##{color2}\" x=\"#{top_left_x.round(3)}\" y=\"#{(top_left_y+height).round(4)}\" width=\"#{(top_right_x-top_left_x).round(4)}\" height=\"#{height.round(3)}\" />\n"
			end
		end
	end

		# if specimen_presence_for_parent1_allele[y][x]
		# 	if specimen_presence_for_parent2_allele[y][x]
		# 		svgstring << "<rect fill=\"##{color1}\" x=\"#{top_left_x.round(3)}\" y=\"#{top_left_y.round(4)}\" width=\"#{width.round(4)}\" height=\"#{height.round(3)}\" />\n"
		# 		svgstring << "<rect fill=\"##{color2}\" x=\"#{top_left_x.round(3)}\" y=\"#{(top_left_y+height).round(4)}\" width=\"#{width.round(4)}\" height=\"#{height.round(3)}\" />\n"
		# 	else
		# 		svgstring << "<rect fill=\"##{color1}\" x=\"#{top_left_x.round(3)}\" y=\"#{top_left_y.round(4)}\" width=\"#{width.round(4)}\" height=\"#{(height*2).round(3)}\" />\n"
		# 	end
		# else
		# 	if specimen_presence_for_parent2_allele[y][x]
		# 		svgstring << "<rect fill=\"##{color2}\" x=\"#{top_left_x.round(3)}\" y=\"#{(top_left_y).round(4)}\" width=\"#{width.round(4)}\" height=\"#{height*2.round(3)}\" />\n"
		# 	end
		# end

end

# Draw a horizontal line above the blocks for this specimen.
n_specimens.times do |y|
	top_left_y = top_margin + y * (2 * height)
	unless y == 0
		svgstring << "<line stroke=\"##{line_color}\" x1=\"#{left_margin + left_text_block_width}\" y1=\"#{top_left_y}\" x2=\"#{dim_x - right_margin}\" y2=\"#{top_left_y}\" />\n"
	end
end

# Add the specimen name at the left of the plot.
n_specimens.times do |y|
	top_left_y = top_margin + y * (2 * height)
	svgstring << "<text x=\"#{left_margin}\" y=\"#{top_left_y+height+font_correction_y}\" fill=\"black\" font-family=\"Helvetica\" font-size=\"#{font_size}\">#{specimen_ids[y]}</text>\n"
end

# Draw the scaffold boundaries.
height = (dim_y - top_margin - bottom_margin)
n_sites.times do |x|
	if new_scaffold[x]
		top_left_x = left_margin + left_text_block_width + x * width
		top_left_y = top_margin
		svgstring << "<line stroke=\"##{line_color}\" x1=\"#{top_left_x}\" y1=\"#{top_left_y}\" x2=\"#{top_left_x}\" y2=\"#{top_left_y+height}\" />\n"
	end
end

# Draw the frame.
svgstring << "<rect stroke=\"##{line_color}\" fill=\"none\" x=\"#{left_margin+left_text_block_width}\" y=\"#{top_margin}\" width=\"#{dim_x-left_margin-left_text_block_width-right_margin}\" height=\"#{dim_y-top_margin-bottom_margin}\" />\n"

# Finalize the svg string.
svgstring << "\n"
svgstring << "</svg>\n"

# Write the svg string to file.
outfile = File.open(outfile_name,"w")
outfile.write(svgstring)
outfile.close

# Feedback.
puts "Wrote file #{outfile_name}."
