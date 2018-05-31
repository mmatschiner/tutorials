# m_matschiner Tue May 29 21:35:40 CEST 2018

# Get the command-line argument.
inmatrix_file_name = ARGV[0]
outmatrix_file_name = ARGV[1]
unless ARGV[2] == nil
	exclude_ids = ARGV[2].split(",")
end

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

if exclude_ids == nil
	species_ids = inmatrix_ids
else
	species_ids = []
	inmatrix_ids.each {|i| species_ids << i unless exclude_ids.include?(i)}
	inmatrix_cells.size.times do |x|
		inmatrix_cells[x] = nil if exclude_ids.include?(inmatrix_ids[x])
	end
	inmatrix_cells.compact!
	inmatrix_cells.size.times do |z|
		row = inmatrix_cells[z]
		row.size.times do |x|
			row[x] = nil if exclude_ids.include?(inmatrix_ids[x])
		end
		row.compact!
	end
end

# Prepare the matrix for age reductions.
outmatrix_cells = []
species_ids.size.times do |x|
	row_ary = []
	species_ids.size.times do |y|
		row_ary << 0
	end
	outmatrix_cells << row_ary
end

# Make the matrix for age reductions.
species_ids.size.times do |a|
	species_ids.size.times do |c|
		mrca_reduction_pairwise_ac = 0
		species_ids.size.times do |b|
			mrca_ab_mean_age = inmatrix_cells[a][b]
			mrca_ac_mean_age = inmatrix_cells[a][c]
			mrca_bc_mean_age = inmatrix_cells[b][c]
			if mrca_ab_mean_age < [mrca_ac_mean_age,mrca_bc_mean_age].min
				mrca_reduction = [mrca_bc_mean_age-mrca_ac_mean_age,0].max
			else
				mrca_reduction = 0
			end
			mrca_reduction_pairwise_ac = mrca_reduction if mrca_reduction > mrca_reduction_pairwise_ac
		end
		outmatrix_cells[a][c] = mrca_reduction_pairwise_ac.round(5)
	end
end

# Determine the cell width.
cell_width = 0
species_ids.each{|i| cell_width = i.size if i.size > cell_width}
outmatrix_cells.each do |row|
	row.each{|i| cell_width = i.to_s.size if i.to_s.size > cell_width}
end
cell_width += 2

# Produce a matrix with either the maximum or the sum of mrca reductions for any combination of species A and C.
outstring = "".ljust(cell_width)
species_ids.each {|i| outstring << "#{i.ljust(cell_width)}"}
outstring << "\n"
species_ids.size.times do |x|
	outstring << "#{species_ids[x].ljust(cell_width)}".ljust(cell_width)
	species_ids.size.times do |y|
		outstring << "#{outmatrix_cells[x][y].to_s.ljust(cell_width)}".ljust(cell_width)
	end
	outstring << "\n"
end

# Write the output matrix.
outmatrix_file = File.open(outmatrix_file_name,"w")
outmatrix_file.write(outstring)
