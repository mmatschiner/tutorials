# m_matschiner Mon Dec 3 15:59:17 CET 2018

# Define a class for lines of the SVG graph.
class Line
	def initialize(x_start,x_end,y_start,y_end,color,stroke,opacity,stroke_dasharray)
		@x_start = x_start
		@x_end = x_end
		@y_start = y_start
		@y_end = y_end
		@color = color
		@stroke = stroke
		@opacity = opacity
		@stroke_dasharray = stroke_dasharray
	end
	def to_svg
		if @stroke_dasharray == "none"
			svg = "<line x1=\"#{@x_start.round(3)}\" y1=\"#{@y_start.round(3)}\" x2=\"#{@x_end.round(3)}\" y2=\"#{@y_end.round(3)}\" stroke=\"#{@color}\" stroke-width=\"#{@stroke}\" stroke-opacity=\"#{@opacity}\" stroke-linecap=\"round\" />"
		else
			svg = "<line x1=\"#{@x_start.round(3)}\" y1=\"#{@y_start.round(3)}\" x2=\"#{@x_end.round(3)}\" y2=\"#{@y_end.round(3)}\" stroke=\"#{@color}\" stroke-width=\"#{@stroke}\" stroke-opacity=\"#{@opacity}\" stroke-linecap=\"round\" stroke-dasharray=\"#{@stroke_dasharray}\" />"
		end
		svg
	end
end

# Define class for rectangles of the SVG graph.
class Rectangle
	def initialize(x,y,width,height,fill,color,stroke,opacity)
		@x = x
		@y = y
		@width = width
		@height = height
		@fill = fill
		@color = color
		@stroke = stroke
		@opacity = opacity
	end
	def to_svg
		svg = "<rect x=\"#{@x}\" y=\"#{@y}\" width=\"#{@width}\" height=\"#{@height}\" fill=\"#{@fill}\" stroke=\"#{@color}\" stroke-width=\"#{@stroke}\" fill-opacity=\"#{@opacity}\" stroke-linejoin=\"round\" />"
		svg
	end
end

# Define a class for text of the SVG graph.
class Text
	def initialize(x,y,font,font_size,weight,color,string,anchor,transform)
		@x = x
		@y = y
		@font = font
		@font_size = font_size
		@weight = weight
		@color = color
		@string = string
		@anchor = anchor
		@transform = transform
	end
	def to_svg
		if @transform == "none"
			if @anchor == "none"
				svg = "<text font-family=\"#{@font}\" font-size=\"#{@font_size}\" fill=\"#{@color}\" font-weight=\"#{@weight}\" x=\"#{@x}\" y=\"#{@y}\">#{@string}</text>"
			else
				svg = "<text text-anchor=\"#{@anchor}\" font-family=\"#{@font}\" font-size=\"#{@font_size}\" fill=\"#{@color}\" font-weight=\"#{@weight}\" x=\"#{@x}\" y=\"#{@y}\">#{@string}</text>"
			end
		else
			if @anchor == "none"
				svg = "<text font-family=\"#{@font}\" font-size=\"#{@font_size}\" fill=\"#{@color}\" font-weight=\"#{@weight}\" x=\"#{@x}\" y=\"#{@y}\" transform=\"#{@transform}\">#{@string}</text>"
			else
				svg = "<text text-anchor=\"#{@anchor}\" font-family=\"#{@font}\" font-size=\"#{@font_size}\" fill=\"#{@color}\" font-weight=\"#{@weight}\" x=\"#{@x}\" y=\"#{@y}\" transform=\"#{@transform}\">#{@string}</text>"
			end
		end
		svg
	end
end

# Get the command-line arguments.
asymmetry_table_file_name = ARGV[0]
heatmap_name_order_file_name = ARGV[1]
max_d_value = ARGV[2].to_f
plot_file_name = ARGV[3]

# Read the order file.
order_file = File.open(heatmap_name_order_file_name)
order_lines = order_file.readlines
names_ordered = []
order_lines.each do |l|
	names_ordered << l.strip unless l[0] == "#"
end

# Read the asymmetry table file.
asymmetry_table_file = File.open(asymmetry_table_file_name)
asymmetry_table_lines = asymmetry_table_file.readlines

# Calculate d and p for each pairwise comparison of two taxa.
d_matrix = []
p_matrix = []
names_ordered.size.times do
	d_matrix << []
	p_matrix << []
end
names_ordered.size.times do |x|
	names_ordered.size.times do |y|
		selected_d_value = 0.0
		selected_p_value = 1.0
		asymmetry_table_lines.each do |l|
			asymmetry_table_line_ary = l.split
			p2 = asymmetry_table_line_ary[1]
			p3 = asymmetry_table_line_ary[2]
			d_value = asymmetry_table_line_ary[3].to_f
                        if asymmetry_table_line_ary[5] == "nan"
				p_value = 1.0
			else
				p_value = asymmetry_table_line_ary[5].to_f
			end
			if [p2,p3] == [names_ordered[x],names_ordered[y]] or [p3,p2] == [names_ordered[x],names_ordered[y]]
				if p_value < selected_p_value
					selected_d_value = d_value
					selected_p_value = p_value
				elsif p_value == selected_p_value and d_value > selected_d_value
					selected_d_value = d_value
				end
			end
		end
		d_matrix[x][y] = selected_d_value
		d_matrix[y][x] = selected_d_value
		p_matrix[x][y] = selected_p_value
		p_matrix[y][x] = selected_p_value
	end
end

# Set up variables for svg elements.
rectangles = []
texts = []
lines = []
svg_width = 170
svg_height = svg_width
svg_margin = 2
text_width = 50
heatmap_left_margin = svg_margin
heatmap_right_margin = svg_margin + text_width
heatmap_top_margin = svg_margin
heatmap_bottom_margin = svg_margin + text_width
font_size = 2.8225806452
legend_font_size = 2.4697580646
font = "Helvetica"
heatmap_height = svg_height - heatmap_top_margin - heatmap_bottom_margin
heatmap_width = svg_width - heatmap_left_margin - heatmap_right_margin
names_horizontal_shift = 2
names_vertical_shift = 1
legend_width = text_width-10
legend_height = legend_width
legend_x = svg_width-svg_margin-legend_width
legend_y = svg_height-svg_margin-legend_height

# Prepare the elements of the svg.
heatmap_frame = Rectangle.new(heatmap_left_margin,heatmap_top_margin,heatmap_width,heatmap_height,"none","black",1,1)
legend_frame = Rectangle.new(legend_x,legend_y,legend_width,legend_height,"none","black",1,1)
(d_matrix.size-1).times do |z|
	# Draw vertical lines.
	x1 = heatmap_left_margin + (z+1) * heatmap_width/d_matrix.size.to_f
	x2 = x1
	y1 = heatmap_top_margin
	y2 = heatmap_top_margin + heatmap_height
	lines << Line.new(x1,x2,y1,y2,"black",0.5,0.5,"none")
	# Draw horizontal lines.
	x1 = heatmap_left_margin
	x2 = heatmap_left_margin + heatmap_width
	y1 = heatmap_top_margin + (z+1) * heatmap_height/d_matrix.size.to_f
	y2 = y1
	lines << Line.new(x1,x2,y1,y2,"black",0.5,0.5,"none")
end
# Write the names.
names_ordered.size.times do |z|
	x = heatmap_left_margin + heatmap_width + names_horizontal_shift
	y = heatmap_top_margin + (z+0.5) * heatmap_height/d_matrix.size.to_f + names_vertical_shift
	texts << Text.new(x,y,font,font_size,"normal","black",names_ordered[z],"start","none")
	x = heatmap_left_margin + (z+0.5) * heatmap_width/d_matrix.size.to_f + names_vertical_shift
	y = heatmap_top_margin + heatmap_height + names_horizontal_shift
	texts << Text.new(x,y,font,font_size,"normal","black",names_ordered[z],"end","rotate(-90,#{x},#{y})")
end

# Fill the heatmap cells with color and opacity according to d and p values.
r_cold = 2       # 38
g_cold = 202     # 139
b_cold = 238     # 210
r_warm = 236     # 220
g_warm = 4   # 102     # 50
b_warm = 31  # 4       # 47
cell_width = heatmap_width/d_matrix.size.to_f
cell_height = heatmap_height/d_matrix.size.to_f
(d_matrix.size).times do |z1|
	x = heatmap_left_margin + z1 * cell_width
	(d_matrix.size).times do |z2|
		y = heatmap_left_margin + z2 * cell_height
		r_cell = (r_cold + (r_warm-r_cold) * [d_matrix[z1][z2]/max_d_value,1].min).to_i
		g_cell = (g_cold + (g_warm-g_cold) * [d_matrix[z1][z2]/max_d_value,1].min).to_i
		b_cell = (b_cold + (b_warm-b_cold) * [d_matrix[z1][z2]/max_d_value,1].min).to_i
		color = "##{r_cell.to_s(16).rjust(2).gsub(" ","0")}#{g_cell.to_s(16).rjust(2).gsub(" ","0")}#{b_cell.to_s(16).rjust(2).gsub(" ","0")}"
		opacity = 0
		if p_matrix[z1][z2] == 0.00000000
			opacity = 1.0
		else
			neg_log_p = -Math.log(p_matrix[z1][z2],10)
			opacity = [(neg_log_p-1)/7.0,0].max
		end
		if opacity > 0
			rectangles << Rectangle.new(x,y,cell_width,cell_height,"#{color}","none",0,opacity)
		end
	end
end

# Fill the legend cells.
n_legend_cells = 20
legend_cell_width = legend_width/n_legend_cells.to_f
legend_cell_height = legend_height/n_legend_cells.to_f
n_legend_cells.times do |z1|
	x = legend_x + z1 * legend_cell_width
	r_legend_cell = (r_cold + (r_warm-r_cold) * z1/(n_legend_cells-1).to_f).to_i
	g_legend_cell = (g_cold + (g_warm-g_cold) * z1/(n_legend_cells-1).to_f).to_i
	b_legend_cell = (b_cold + (b_warm-b_cold) * z1/(n_legend_cells-1).to_f).to_i
	color = "##{r_legend_cell.to_s(16).rjust(2).gsub(" ","0")}#{g_legend_cell.to_s(16).rjust(2).gsub(" ","0")}#{b_legend_cell.to_s(16).rjust(2).gsub(" ","0")}"
	n_legend_cells.times do |z2|
		y = legend_x + z2 * legend_cell_height
		opacity = z2/(n_legend_cells-1).to_f
		rectangles << Rectangle.new(x,y,legend_cell_width,legend_cell_height,"#{color}","none",0,opacity)
	end
end

# Write the horizontal legend labels.
legend_label_shift1 = 4
legend_label_shift2 = 1
legend_label_shift3 = 2
legend_label_shift4 = 2
x = legend_x + 0.5 * legend_width
y = legend_y - legend_label_shift1
texts << Text.new(x,y,font,legend_font_size,"normal","black","D","middle","none")
y = legend_y - legend_label_shift2
((max_d_value*10).to_i+1).times do |z|
	label = z/10.0
	x = legend_x + legend_label_shift2 + (z/(max_d_value*10).to_i.to_f)*(legend_width - 2*legend_label_shift2)
	texts << Text.new(x,y,font,legend_font_size,"normal","black",label,"middle","none")
end
# Write the vertical legend labels.
x = legend_x - legend_label_shift2
1.upto(8) do |z|
	y = legend_y + legend_label_shift3 + ((z-1)/7.0)*(legend_height-legend_label_shift3)
	texts << Text.new(x,y,font,legend_font_size,"normal","black","-#{z}","end","none")
end
x = legend_x - legend_label_shift1 - legend_label_shift4
y = legend_y + 0.5 * legend_height
texts << Text.new(x,y,font,legend_font_size,"normal","black","log(p)","middle","rotate(-90,#{x},#{y})")

# Prepare the svg string.
svg_string = ""
svg_string << "<?xml version=\"1.0\" standalone=\"no\"?>\n"
svg_string << "<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\" \"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">\n"
svg_string << "<svg width=\"#{svg_width}mm\" height=\"#{svg_height}mm\" viewBox=\"0 0 #{svg_width} #{svg_height}\" xmlns=\"http://www.w3.org/2000/svg\" version=\"1.1\">\n"
rectangles.each {|r| svg_string << "    #{r.to_svg}\n"}
texts.each {|t| svg_string << "    #{t.to_svg}\n"}
lines.each {|l| svg_string << "    #{l.to_svg}\n"}
svg_string << "    #{heatmap_frame.to_svg}\n"
svg_string << "    #{legend_frame.to_svg}\n"
svg_string << "</svg>\n"

# Write the svg string to file.
svg_file = File.new(plot_file_name,"w")
svg_file.write(svg_string)
