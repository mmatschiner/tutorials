# Michael Matschiner, 2016-02-26

# This script uses output from the software Saguaro to paint
# chromosomes according to Saguaro cacti. The output will be
# written in svg format.

# This script should be run e.g. with
# ruby paint_chromosomes.rb LocalTrees.out LocalTrees.svg
# where 'LocalTrees.out' should be replaced with the actual path
# to the Saguaro output file.

# Check if command line arguments are provided, and print help
# text if they aren't.
if ARGV == [] or ARGV.include?("-h") or ARGV.include?("--help")
	puts
	puts "paint_chromosomes.rb"
	puts
	puts "This script uses output from the software Saguaro to paint"
	puts "chromosomes according to Saguaro cacti. The output will be"
	puts "written in svg format."
	puts
	puts "This script should be run e.g. with"
	puts "ruby paint_chromosomes.rb LocalTrees.out LocalTrees.svg"
	puts "where 'LocalTrees.out' should be replaced with the actual path"
	puts "to the Saguaro output file."
	exit
end

# Get the command line arguments.
input_file_name = ARGV[0]
if ARGV[1] == nil
	output_file_name = ARGV[0].chomp(".out").chomp(".txt") + ".svg"
else
	output_file_name = ARGV[1]
end

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
cactus_ids = []
lg_ids = []
segment_starts = []
segment_ends = []
filtered_lines.each do |l|
	line_ary = l.split
	cactus_ids << line_ary[0]
	lg_ids << line_ary[1].chomp(":").split(".")[0]
	segment_starts << line_ary[2].to_i
	segment_ends << line_ary[4].to_i
end
uniq_lg_ids = lg_ids.uniq.sort
uniq_cactus_ids = cactus_ids.uniq.sort

# Prepare the svg string.
top_margin = 10
bottom_margin = 10
scale_bar_height = 40
scale_bar_font_shift = 10
legend_font_shift = 30
left_margin = 10
text_width = 70
right_margin = 10
height_per_chromosome = 100
space_between_chromosomes = 10
dim_x = 800
stroke_width = 1
font_size = 12
font_correction_y = 4
tick_length = 2
minor_tick_distance = 5000000
major_tick_distance = 10000000
colors = ["5b5b5b","ef2746","f59031","27abd0","bad531","606cb3"]
rest_color = "bfbfbf"
dim_y = top_margin + bottom_margin + scale_bar_height + uniq_lg_ids.size * height_per_chromosome + (uniq_lg_ids.size-1) * space_between_chromosomes
svgstring = "<?xml version=\"1.0\" standalone=\"no\"?>\n"
svgstring << "<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.0//EN\" \"http://www.w3.org/TR/2001/REC-SVG-20010904/DTD/svg10.dtd\">\n"
svgstring << "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"#{dim_x}\" height=\"#{dim_y}\" viewBox=\"0 0 #{dim_x} #{dim_y}\" xmlns:xlink=\"htttp://www.w3.org/1999/xlink\">\n\n"
svgstring << "  <defs>\n"
svgstring << "    <style type=\"text/css\">\n"
svgstring << "      <![CDATA[\n"
svgstring << "        line {fill:none; stroke:#93a1a1; stroke-width:1.0px;}\n"
svgstring << "      ]]>\n"
svgstring << "    </style>\n"
svgstring << "  </defs>\n"
svgstring << "\n"

# Get the longest chromosome length.
max_chromosome_length = 0
segment_ends.each do |e|
	max_chromosome_length = e if e > max_chromosome_length
end

# Get the total length for each cactus.
uniq_cactus_total_lengths = []
uniq_cactus_ids.each do |c|
	total_length_for_this_cactus = 0
	cactus_ids.size.times do |x|
		if cactus_ids[x] == c
			total_length_for_this_cactus += segment_ends[x] - segment_starts[x]
		end
	end
	uniq_cactus_total_lengths << total_length_for_this_cactus
end

# Sort unique cactus ids according to their total length.
sorted = false
until sorted
	sorted = true
	0.upto(uniq_cactus_ids.size-2) do |x|
		(x+1).upto(uniq_cactus_ids.size-1) do |y|
			if uniq_cactus_total_lengths[y] > uniq_cactus_total_lengths[x]
				uniq_cactus_total_lengths[x],uniq_cactus_total_lengths[y] = uniq_cactus_total_lengths[y],uniq_cactus_total_lengths[x]
				uniq_cactus_ids[x],uniq_cactus_ids[y] = uniq_cactus_ids[y],uniq_cactus_ids[x]
				sorted = false
			end
		end
	end
end

# Write the painted chromosomes part of the svg.
top_left_y = top_margin
height = height_per_chromosome
uniq_lg_ids.each do |lg|
	svgstring << "<text x=\"#{left_margin}\" y=\"#{top_left_y+(0.5*height_per_chromosome)+font_correction_y}\" fill=\"black\" font-family=\"Helvetica\" font-size=\"#{font_size}\" font-weight=\"bold\">#{lg}</text>\n"
	segment_starts_for_this_lg = []
	segment_ends_for_this_lg = []
	lg_ids.size.times do |x|
		if lg_ids[x] == lg
			segment_starts_for_this_lg << segment_starts[x]
			segment_ends_for_this_lg << segment_ends[x]
			top_left_x = left_margin + text_width + segment_starts[x] * ((dim_x - left_margin - text_width - right_margin)/max_chromosome_length.to_f)
			width = (segment_ends[x]-segment_starts[x]) * ((dim_x - left_margin - text_width - right_margin)/max_chromosome_length.to_f)
			if uniq_cactus_ids.index(cactus_ids[x]) < colors.size
				color = colors[uniq_cactus_ids.index(cactus_ids[x])]
			else
				color = rest_color
			end
			svgstring << "<rect style=\"stroke:none; fill:##{color}\" x=\"#{top_left_x}\" y=\"#{top_left_y}\" width=\"#{width}\" height=\"#{height}\" />\n"
		end
	end
	frame_top_left_x = left_margin + text_width + segment_starts_for_this_lg.min * ((dim_x - left_margin - text_width - right_margin)/max_chromosome_length.to_f)
	frame_width = (segment_ends_for_this_lg.max-segment_starts_for_this_lg.min) * ((dim_x - left_margin - text_width - right_margin)/max_chromosome_length.to_f)
	svgstring << "<rect style=\"stroke:black; stroke-width:#{stroke_width}px; fill:none\" x=\"#{frame_top_left_x}\" y=\"#{top_left_y}\" width=\"#{frame_width}\" height=\"#{height}\" />\n"
	top_left_y += height_per_chromosome
	top_left_y += space_between_chromosomes
end

# Write the scale bar part of the svg.
svgstring << "<line style=\"stroke:black; stroke-width:#{stroke_width}px; fill:none\" x1=\"#{left_margin+text_width}\" y1=\"#{top_left_y}\" x2=\"#{dim_x-right_margin}\" y2=\"#{top_left_y}\"/>\n"
((max_chromosome_length/minor_tick_distance)+1).times do |x|
	svgstring << "<line style=\"stroke:black; stroke-width:#{stroke_width}px; fill:none\" x1=\"#{left_margin+text_width+((x*minor_tick_distance)/max_chromosome_length.to_f)*(dim_x-left_margin-text_width-right_margin)}\" y1=\"#{top_left_y}\" x2=\"#{left_margin+text_width+((x*minor_tick_distance)/max_chromosome_length.to_f)*(dim_x-left_margin-text_width-right_margin)}\" y2=\"#{top_left_y+tick_length}\"/>\n"
end
((max_chromosome_length/major_tick_distance)+1).times do |x|
	svgstring << "<text text-anchor=\"middle\" x=\"#{left_margin+text_width+((x*major_tick_distance)/max_chromosome_length.to_f)*(dim_x-left_margin-text_width-right_margin)}\" y=\"#{top_left_y+tick_length+scale_bar_font_shift+font_correction_y}\" fill=\"black\" font-family=\"Helvetica\" font-size=\"#{font_size}\">#{x*major_tick_distance}</text>\n"
end
svgstring << "<text text-anchor=\"middle\" x=\"#{left_margin+text_width+((dim_x-left_margin-text_width-right_margin)/2.0)}\" y=\"#{top_left_y+legend_font_shift+font_correction_y}\" fill=\"black\" font-family=\"Helvetica\" font-size=\"#{font_size}\" font-weight=\"bold\">Chromosome position (bp)</text>\n"

# Finalize the svg string.
svgstring << "\n"
svgstring << "</svg>\n"

# Write the svg string to file.
output_file = File.open(output_file_name,"w")
output_file.write(svgstring)
output_file.close

# Feedback.
puts "Wrote file #{output_file_name}."
