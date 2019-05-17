# m_matschiner Sun Dec 2 15:28:49 CET 2018

# Define a class for lines of the SVG graph.
class Line
	def initialize(x_start,x_end,y_start,y_end,color,stroke,opacity,stroke_dasharray,transform)
		@x_start = x_start
		@x_end = x_end
		@y_start = y_start
		@y_end = y_end
		@color = color
		@stroke = stroke
		@opacity = opacity
		@stroke_dasharray = stroke_dasharray
		@transform = transform
	end
	def to_s
		string = ""
		string << "Start:  #{@x_start.round(3)},#{@y_start.round(3)}\n"
		string << "End:    #{@x_end.round(3)},#{@y_end.round(3)}\n"
		string << "Color:  #{@color}\n"
		string << "Stroke: #{@stroke}\n"
		string << "\n"
		string
	end
	def to_svg
		if @stroke_dasharray == "none"
			if @transform == "none"
				svg = "<line x1=\"#{@x_start.round(3)}\" y1=\"#{@y_start.round(3)}\" x2=\"#{@x_end.round(3)}\" y2=\"#{@y_end.round(3)}\" stroke=\"#{@color}\" stroke-width=\"#{@stroke}\" stroke-opacity=\"#{@opacity}\" stroke-linecap=\"round\" />"
			else
				svg = "<line x1=\"#{@x_start.round(3)}\" y1=\"#{@y_start.round(3)}\" x2=\"#{@x_end.round(3)}\" y2=\"#{@y_end.round(3)}\" stroke=\"#{@color}\" stroke-width=\"#{@stroke}\" stroke-opacity=\"#{@opacity}\" stroke-linecap=\"round\" transform=\"#{@transform}\" />"
			end
		else
			if @transform == "none"
				svg = "<line x1=\"#{@x_start.round(3)}\" y1=\"#{@y_start.round(3)}\" x2=\"#{@x_end.round(3)}\" y2=\"#{@y_end.round(3)}\" stroke=\"#{@color}\" stroke-width=\"#{@stroke}\" stroke-opacity=\"#{@opacity}\" stroke-linecap=\"round\" stroke-dasharray=\"#{@stroke_dasharray}\" />"
			else
				svg = "<line x1=\"#{@x_start.round(3)}\" y1=\"#{@y_start.round(3)}\" x2=\"#{@x_end.round(3)}\" y2=\"#{@y_end.round(3)}\" stroke=\"#{@color}\" stroke-width=\"#{@stroke}\" stroke-opacity=\"#{@opacity}\" stroke-linecap=\"round\" stroke-dasharray=\"#{@stroke_dasharray}\" transform=\"#{@transform}\" />"
			end
		end
		svg
	end
	def max_y
		[@y_start,@y_end].max
	end
end

# Define a class for text of the SVG graph.
class Text
	def initialize(x,y,font,font_size,color,string,transform)
		@x = x
		@y = y
		@font = font
		@font_size = font_size
		@color = color
		@string = string
		@transform = transform
	end
	def to_svg
		if @transform == "none"
			svg = "<text text-anchor=\"middle\" font-family=\"#{@font}\" font-size=\"#{@font_size}\" fill=\"#{@color}\" font-weight=\"bold\" x=\"#{@x}\" y=\"#{@y}\">#{@string}</text>"
		else
			svg = "<text text-anchor=\"middle\" font-family=\"#{@font}\" font-size=\"#{@font_size}\" fill=\"#{@color}\" font-weight=\"bold\" x=\"#{@x}\" y=\"#{@y}\" transform=\"#{@transform}\">#{@string}</text>"
		end
		svg
	end
end

# Define class for paths of the SVG graph.
class Path
	def initialize(x,y,fill,color,stroke,opacity)
		@x = [x]
		@y = [y]
		@fill = fill
		@color = color
		@stroke = stroke
		@opacity = opacity
	end
	def add_point(x,y)
		@x << x
		@y << y
	end
	def to_svg
		svg = "<path d=\"M #{@x[0]} #{@y[0]} "
		if @x.size > 1
			1.upto(@x.size-1) do |z|
				svg << "L #{@x[z]} #{@y[z]} "
			end
		end
		svg << "z\" fill=\"#{@fill}\" stroke=\"#{@color}\" stroke-width=\"#{@stroke}\" fill-opacity=\"#{@opacity}\" stroke-linejoin=\"round\" />"
		svg
	end
end

# Define a class for circles of the SVG graph.
class Circle
	def initialize(cx,cy,r,fill,color,stroke,opacity) # (cx,cy,cr,"gray","black",1,0.5)
		@cx = cx
		@cy = cy
		@r = r
		@fill = fill
		@color = color
		@stroke = stroke # Stroke color.
		@opacity = opacity
	end
	def to_svg
		if @color == "none"
			svg = "<circle cx=\"#{@cx}\" cy=\"#{@cy}\" r=\"#{@r}\" fill=\"#{@fill}\" stroke=\"#{@color}\" opacity=\"#{@opacity}\" />"
		else
			svg = "<circle cx=\"#{@cx}\" cy=\"#{@cy}\" r=\"#{@r}\" fill=\"#{@fill}\" stroke=\"#{@color}\" stroke-width=\"#{@stroke}\" opacity=\"#{@opacity}\" />"
		end
		svg
	end
end

# Add methods for all enumerables.
module Enumerable
	def sum
		self.inject(0){|accum, i| accum + i }
	end
	def mean
		if self.length == 0
			nil
		else
			self.sum/self.length.to_f
		end
	end
end

# Get the command-line arguments.
likelihood_table_file_name = ARGV[0]
svg_file_name = ARGV[1]

# Read the likelihood table.
likelihood_table_file = File.open(likelihood_table_file_name)
likelihood_table_lines = likelihood_table_file.readlines

# Make sure the header is as expected.
header_ary = likelihood_table_lines[0].split
unless header_ary[1] == "best" and header_ary[2] == "delta_lik" and header_ary.size == 6
	puts "ERROR: Unexpected header line: #{likelihood_table_lines[0]}!"
	exit 1
end

# Get the likelihoods for three hypotheses.
t1_lik_diffs = []
t2_lik_diffs = []
t3_lik_diffs = []
t1_over_t2 = 0
t2_over_t1 = 0
t1_over_t3 = 0
t3_over_t1 = 0
t2_over_t3 = 0
t3_over_t2 = 0
max_lik_diff = 0
likelihood_table_lines[1..-1].each do |l|
	line_ary = l.split
	t1_lik = line_ary[3].to_f
	t2_lik = line_ary[4].to_f
	t3_lik = line_ary[5].to_f
	min_lik = [t1_lik,t2_lik,t3_lik].min
	max_lik = [t1_lik,t2_lik,t3_lik].max
	lik_diff = max_lik-min_lik
	max_lik_diff = lik_diff if lik_diff > max_lik_diff
	t1_lik_diffs << t1_lik-min_lik
	t2_lik_diffs << t2_lik-min_lik
	t3_lik_diffs << t3_lik-min_lik
	t1_over_t2 += 1 if t1_lik > t2_lik
	t2_over_t1 += 1 if t2_lik > t1_lik
	t1_over_t3 += 1 if t1_lik > t3_lik
	t3_over_t1 += 1 if t3_lik > t1_lik
	t2_over_t3 += 1 if t2_lik > t3_lik
	t3_over_t2 += 1 if t3_lik > t2_lik
end

# Get the mean values.
t1_lik_diff_mean = t1_lik_diffs.mean
t2_lik_diff_mean = t2_lik_diffs.mean
t3_lik_diff_mean = t3_lik_diffs.mean

# Prepare the elements of a svg.
texts = []
paths = []
lines = []
circles = []
svg_width = 60
svg_height = 60
svg_margin = 5
font_size_l = 2.8218694885
font_size_s = 2.1164021164
cr = 0.5
lower_label_vertical_shift = 2
upper_label_vertical_shift = -2
central_x = svg_width/2.0
max_lik_diff_for_plot = (max_lik_diff*1.3).ceil
triangle_width = svg_width - 2*svg_margin
triangle_height = Math.sqrt(3) * triangle_width / 2.0
triangle_top_bottom_space = (svg_height - triangle_height)/2.0
if triangle_top_bottom_space < svg_margin
	puts "ERROR: Horizontal dimensions are not large enough to fit the triangle plot!"
	exit 1
end
central_y = triangle_top_bottom_space+(2*triangle_height/3.0)
t1_red = 42
t1_green = 161
t1_blue = 152
t2_red = 108
t2_green = 113
t2_blue = 196
t3_red = 203
t3_green = 75
t3_blue = 22
center_red = 147
center_green = 161
center_blue = 161
scale_bar_shift = 4
scale_bar_legend_shift = 2.5
legend_shift = 1
arrow_shift = 0.75

# Add the triangle.
triangle = Path.new(svg_margin,triangle_top_bottom_space+triangle_height,"none","black",1,1)
triangle.add_point(svg_margin+triangle_width,triangle_top_bottom_space+triangle_height)
triangle.add_point(central_x,triangle_top_bottom_space)
paths << triangle

# Add labels.
texts << Text.new(0.5*svg_margin,triangle_top_bottom_space+triangle_height+lower_label_vertical_shift,"Helvetica",font_size_l,"##{t1_red.to_s(16)}#{t1_green.to_s(16)}#{t1_blue.to_s(16)}","T1","none")
texts << Text.new(1.5*svg_margin+triangle_width,triangle_top_bottom_space+triangle_height+lower_label_vertical_shift,"Helvetica",font_size_l,"##{t2_red.to_s(16)}#{t2_green.to_s(16)}#{t2_blue.to_s(16)}","T2","none")
texts << Text.new(central_x,triangle_top_bottom_space+upper_label_vertical_shift,"Helvetica",font_size_l,"##{t3_red.to_s(16)}#{t3_green.to_s(16)}#{t3_blue.to_s(16)}","T3","none")

# Add lines from the corners to the center.
lines << Line.new(svg_margin,central_x,triangle_top_bottom_space+triangle_height,central_y,"black",0.5,1,1,"none")
lines << Line.new(svg_margin+triangle_width,central_x,triangle_top_bottom_space+triangle_height,central_y,"black",0.5,1,1,"none")
lines << Line.new(central_x,central_x,triangle_top_bottom_space,central_y,"black",0.5,1,1,"none")

# Add circles for each likelihood comparison.
t1_lik_diffs.size.times do |x|
	cx = central_x - (t1_lik_diffs[x]/max_lik_diff_for_plot.to_f)*(triangle_width/2.0) + (t2_lik_diffs[x]/max_lik_diff_for_plot.to_f)*(triangle_width/2.0)
	cy = central_y + (t1_lik_diffs[x]/max_lik_diff_for_plot.to_f)*(triangle_height/3.0) + (t2_lik_diffs[x]/max_lik_diff_for_plot.to_f)*(triangle_height/3.0) - (t3_lik_diffs[x]/max_lik_diff_for_plot.to_f)*(2*triangle_height/3.0)
	c_red = (center_red + (t1_lik_diffs[x]/max_lik_diff_for_plot.to_f) * (t1_red-center_red) + (t2_lik_diffs[x]/max_lik_diff_for_plot.to_f) * (t2_red-center_red) + (t3_lik_diffs[x]/max_lik_diff_for_plot.to_f) * (t3_red-center_red)).to_i
	c_green = (center_green + (t1_lik_diffs[x]/max_lik_diff_for_plot.to_f) * (t1_green-center_green) + (t2_lik_diffs[x]/max_lik_diff_for_plot.to_f) * (t2_green-center_green) + (t3_lik_diffs[x]/max_lik_diff_for_plot.to_f) * (t3_green-center_green)).to_i
	c_blue = (center_blue + (t1_lik_diffs[x]/max_lik_diff_for_plot.to_f) * (t1_blue-center_blue) + (t2_lik_diffs[x]/max_lik_diff_for_plot.to_f) * (t2_blue-center_blue) + (t3_lik_diffs[x]/max_lik_diff_for_plot.to_f) * (t3_blue-center_blue)).to_i
	if [c_red,c_green,c_blue].min >= 0 and [c_red,c_green,c_blue].max <= 255
		c_red_hex = c_red.to_s(16)
		c_red_hex = "0#{c_red_hex}" if c_red_hex.size == 1
		c_green_hex = c_green.to_s(16)
		c_green_hex = "0#{c_green_hex}" if c_green_hex.size == 1
		c_blue_hex = c_blue.to_s(16)
		c_blue_hex = "0#{c_blue_hex}" if c_blue_hex.size == 1
		circles << Circle.new(cx,cy,cr,"##{c_red_hex}#{c_green_hex}#{c_blue_hex}","none",0,0.5)
	end
end

# Add a triangle for the mean.
mx = central_x - (t1_lik_diff_mean/max_lik_diff_for_plot.to_f)*(triangle_width/2.0)
my = central_y + (t1_lik_diff_mean/max_lik_diff_for_plot.to_f)*(triangle_height/3.0)
mean_triangle = Path.new(mx,my,"none","black",1,1)
mx = central_x + (t2_lik_diff_mean/max_lik_diff_for_plot)*(triangle_width/2.0)
my = central_y + (t2_lik_diff_mean/max_lik_diff_for_plot)*(triangle_height/3.0)
mean_triangle.add_point(mx,my)
mx = central_x
my = central_y - (t3_lik_diff_mean/max_lik_diff_for_plot)*(2*triangle_height/3.0)
mean_triangle.add_point(mx,my)
paths << mean_triangle

# Add a small circle for the center.
cx = central_x
cy = central_y
circles << Circle.new(cx,cy,0.1,"black","none",0,1)

# Add the scale bar.
lik_diff_for_scale_bar = 10
scale_bar_width = (lik_diff_for_scale_bar/max_lik_diff_for_plot.to_f) * (2*triangle_height/3.0)
lines << Line.new(central_x-(scale_bar_width/2.0),central_x+(scale_bar_width/2.0),triangle_top_bottom_space+triangle_height+scale_bar_shift,triangle_top_bottom_space+triangle_height+scale_bar_shift,"black",1,1,"none","none")
texts << Text.new(central_x,triangle_top_bottom_space+triangle_height+scale_bar_shift+scale_bar_legend_shift,"Helvetica",font_size_s,"black","&#916;L = #{lik_diff_for_scale_bar}","none")
texts << Text.new(central_x,triangle_top_bottom_space+triangle_height-legend_shift,"Helvetica",font_size_s,"black","n = #{t1_lik_diffs.size}","none")

# Add the arrows.
x1 = svg_margin
x2 = central_x-0.5
y1 = triangle_top_bottom_space+triangle_height+arrow_shift
y2 = triangle_top_bottom_space+triangle_height+arrow_shift
lines << Line.new(x1,x2,y1,y2,"##{t1_red.to_s(16)}#{t1_green.to_s(16)}#{t1_blue.to_s(16)}",1,1,"none","none")
texts << Text.new((x1+x2)/2.0,((y1+y2)/2.0)+3,"Helvetica",font_size_s,"##{t1_red.to_s(16)}#{t1_green.to_s(16)}#{t1_blue.to_s(16)}","n = #{t1_over_t2}","none")
x1 = central_x+0.5
x2 = svg_margin+triangle_width
y1 = triangle_top_bottom_space+triangle_height+arrow_shift
y2 = triangle_top_bottom_space+triangle_height+arrow_shift
lines << Line.new(x1,x2,y1,y2,"##{t2_red.to_s(16)}#{t2_green.to_s(16)}#{t2_blue.to_s(16)}",1,1,"none","none")
texts << Text.new((x1+x2)/2.0,((y1+y2)/2.0)+3,"Helvetica",font_size_s,"##{t2_red.to_s(16)}#{t2_green.to_s(16)}#{t2_blue.to_s(16)}","n = #{t2_over_t1}","none")
x1 = svg_margin
x2 = central_x-0.5
y1 = triangle_top_bottom_space+triangle_height-arrow_shift
y2 = triangle_top_bottom_space+triangle_height-arrow_shift
lines << Line.new(x1,x2,y1,y2,"##{t1_red.to_s(16)}#{t1_green.to_s(16)}#{t1_blue.to_s(16)}",1,1,"none","rotate(-60,#{svg_margin},#{triangle_top_bottom_space+triangle_height})")
texts << Text.new((x1+x2)/2.0,((y1+y2)/2.0)-1.3,"Helvetica",font_size_s,"##{t1_red.to_s(16)}#{t1_green.to_s(16)}#{t1_blue.to_s(16)}","n = #{t1_over_t3}","rotate(-60,#{svg_margin},#{triangle_top_bottom_space+triangle_height})")
x1 = central_x+0.5
x2 = svg_margin+triangle_width
y1 = triangle_top_bottom_space+triangle_height-arrow_shift
y2 = triangle_top_bottom_space+triangle_height-arrow_shift
lines << Line.new(x1,x2,y1,y2,"##{t3_red.to_s(16)}#{t3_green.to_s(16)}#{t3_blue.to_s(16)}",1,1,"none","rotate(-60,#{svg_margin},#{triangle_top_bottom_space+triangle_height})")
texts << Text.new((x1+x2)/2.0,((y1+y2)/2.0)-1.3,"Helvetica",font_size_s,"##{t3_red.to_s(16)}#{t3_green.to_s(16)}#{t3_blue.to_s(16)}","n = #{t3_over_t1}","rotate(-60,#{svg_margin},#{triangle_top_bottom_space+triangle_height})")
x1 = svg_margin
x2 = central_x-0.5
y1 = triangle_top_bottom_space+triangle_height-arrow_shift
y2 = triangle_top_bottom_space+triangle_height-arrow_shift
lines << Line.new(x1,x2,y1,y2,"##{t3_red.to_s(16)}#{t3_green.to_s(16)}#{t3_blue.to_s(16)}",1,1,"none","rotate(60,#{svg_margin+triangle_width},#{triangle_top_bottom_space+triangle_height})")
texts << Text.new((x1+x2)/2.0,((y1+y2)/2.0)-1.3,"Helvetica",font_size_s,"##{t3_red.to_s(16)}#{t3_green.to_s(16)}#{t3_blue.to_s(16)}","n = #{t3_over_t2}","rotate(60,#{svg_margin+triangle_width},#{triangle_top_bottom_space+triangle_height})")
x1 = central_x+0.5
x2 = svg_margin+triangle_width
y1 = triangle_top_bottom_space+triangle_height-arrow_shift
y2 = triangle_top_bottom_space+triangle_height-arrow_shift
lines << Line.new(x1,x2,y1,y2,"##{t2_red.to_s(16)}#{t2_green.to_s(16)}#{t2_blue.to_s(16)}",1,1,"none","rotate(60,#{svg_margin+triangle_width},#{triangle_top_bottom_space+triangle_height})")
texts << Text.new((x1+x2)/2.0,((y1+y2)/2.0)-1.3,"Helvetica",font_size_s,"##{t2_red.to_s(16)}#{t2_green.to_s(16)}#{t2_blue.to_s(16)}","n = #{t2_over_t3}","rotate(60,#{svg_margin+triangle_width},#{triangle_top_bottom_space+triangle_height})")

# Prepare the svg string.
svg_string = ""
svg_string << "<?xml version=\"1.0\" standalone=\"no\"?>\n"
svg_string << "<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\" \"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">\n"
svg_string << "<svg width=\"#{svg_width}mm\" height=\"#{svg_height}mm\" viewBox=\"0 0 #{svg_width} #{svg_height}\" xmlns=\"http://www.w3.org/2000/svg\" version=\"1.1\">\n"
texts.each {|t| svg_string << "    #{t.to_svg}\n"}
lines.each {|l| svg_string << "    #{l.to_svg}\n"}
circles.each {|c| svg_string << "    #{c.to_svg}\n"}
paths.each {|p| svg_string << "    #{p.to_svg}\n"}
svg_string << "</svg>\n"

# Write the svg string to file.
svg_file = File.new(svg_file_name,"w")
svg_file.write(svg_string)

# Feedback.
puts "Wrote file #{svg_file_name}."
