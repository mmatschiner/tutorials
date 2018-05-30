# m_matschiner Tue May 29 17:20:38 CEST 2018

# Get the command-line argument.
inmatrix_file_name = ARGV[0]
assignment_file_name = ARGV[1]
outmatrix_file_name = ARGV[2]

# Read the input matrix.
inmatrix_file = File.open(inmatrix_file_name)
inmatrix_lines = inmatrix_file.readlines
inmatrix_ids = inmatrix_lines[0].split
inmatrix_cells = []
inmatrix_vertical_ids = []
inmatrix_lines[1..-1].each do |l|
	line_ary = l.split
	inmatrix_vertical_ids << line_ary[0]
	inmatrix_row = []
	line_ary[1..-1].each do |i|
		inmatrix_row << i.to_f
	end
	inmatrix_cells << inmatrix_row
end

# Make sure that the matrix is ordered in the same way horizontally and vertically.
if inmatrix_ids != inmatrix_vertical_ids
	puts "ERROR: The matrix is not ordered identically horizontally and vertically!"
	exit 1
end

# Read the assignment file.
assignment_file = File.open(assignment_file_name)
assignment_lines = assignment_file.readlines
smaller_units = []
larger_units = []
assignment_lines.each do |l|
	line_ary = l.split
	smaller_units << line_ary[0]
	larger_units << line_ary[1]
end

# Make sure that more smaller units than larger units are in the assignment file.
if smaller_units.uniq.size < larger_units.uniq.size
	puts "WARNING: The number of unique smaller units is smaller than the number of unique larger units. These will be swapped."
	smaller_units, larger_units = larger_units, smaller_units
end

# Prepare the matrix for the larger units.
outmatrix_ids = larger_units.sort.uniq
outmatrix_cells = []
outmatrix_ids.size.times do |x|
	row_ary = []
	outmatrix_ids.size.times do |y|
		row_ary << 0
	end
	outmatrix_cells << row_ary
end

# Make the matrix for the larger units.
outmatrix_ids.size.times do |a|
	u = outmatrix_ids[a]
	smaller_units_for_the_first_larger_unit = []
	smaller_units.size.times do |x|
		smaller_units_for_the_first_larger_unit << smaller_units[x] if larger_units[x] == u
	end
	outmatrix_ids.size.times do |b|
		v = outmatrix_ids[b]
		smaller_units_for_the_second_larger_unit = []
		smaller_units.size.times do |x|
			smaller_units_for_the_second_larger_unit << smaller_units[x] if larger_units[x] == v
		end
		# Get the mean value from pairwise comparisons of the smaller units of the first and second larter unit.
		comparison_values = []
		inmatrix_ids.size.times do |x|
			inmatrix_ids.size.times do |y|
				if smaller_units_for_the_first_larger_unit.include?(inmatrix_ids[x]) and smaller_units_for_the_second_larger_unit.include?(inmatrix_ids[y])
					comparison_values << inmatrix_cells[x][y]
				end
			end
		end
		# Make sure that pairwise comparisons are found.
		if comparison_values.size == 0
			puts "ERROR: No pairwise comparisons of the larger units #{u} and #{v} were found!"
			exit 1
		end
		comparison_values_sum = 0
		comparison_values.each{|i| comparison_values_sum += i}
		outmatrix_cells[a][b] = ((comparison_values_sum)/comparison_values.size.to_f).round(5)
	end
end

# Prepare the string for the output matrix.
outstring = ""
cell_width = 0
outmatrix_ids.each{|i| cell_width = i.size if i.size > cell_width}
outmatrix_cells.each do |row|
	row.each{|i| cell_width = i.to_s.size if i.to_s.size > cell_width}
end
cell_width += 2
outstring << "".ljust(cell_width)
outmatrix_ids.each {|i| outstring << "#{i.ljust(cell_width)}"}
outstring << "\n"
outmatrix_ids.size.times do |x|
	outstring << "#{outmatrix_ids[x].ljust(cell_width)}".ljust(cell_width)
	outmatrix_ids.size.times do |y|
		outstring << "#{outmatrix_cells[x][y].to_s.ljust(cell_width)}".ljust(cell_width)
	end
	outstring << "\n"
end

# Write the output file.
outmatrix_file = File.open(outmatrix_file_name,"w")
outmatrix_file.write(outstring)
