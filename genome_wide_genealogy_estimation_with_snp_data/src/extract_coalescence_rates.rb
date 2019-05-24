# m_matschiner Sun May 19 14:18:27 CEST 2019

# Get the command-line arguments.
coalescence_rates_file_name = ARGV[0]
generation_time = ARGV[1].to_f
species = ARGV[2]
table_file_name = ARGV[3]

# Read the coalescence rates file.
coalescence_rates_file = File.open(coalescence_rates_file_name)
coalescence_rates_lines = coalescence_rates_file.readlines

# Get the populations.
pops = coalescence_rates_lines[0].split

# Get the time bins.
time_breaks = coalescence_rates_lines[1].split.map {|i| i.to_f * generation_time}

# Get the coalescence times among populations.
coalescence_rates_per_pair = []
pops.size.times do |x|
	pops.size.times do |y|
		if y > x
			if pops[x] == species or pops[y] == species
				coalescence_rates_this_pair = nil
				coalescence_rates_lines[3..-1].each do |l|
					if l.split[0].to_i == x and l.split[1].to_i == y
						coalescence_rates_this_pair = l.split[2..-1]
					end
				end
				if coalescence_rates_this_pair == nil
					puts "ERROR: Couldn't find cross-coalescence rates for populations #{pops[x]} and #{pops[y]}!"
					exit 1
				end
				coalescence_rates_per_pair << coalescence_rates_this_pair
			end
		end
	end
end

# Prepare the output table.
outstring = ""
(time_breaks.size-1).times do |x|
	bin_center = 0.5 * (time_breaks[x+1] + time_breaks[x])
	outstring << "#{bin_center}\t"
	coalescence_rates_per_pair.size.times do |y|
		outstring << "#{coalescence_rates_per_pair[y][x]}\t"
	end
	outstring.chomp!("\t")
	outstring << "\n"
end

# Write the output table.
table_file = File.open(table_file_name, "w")
table_file.write(outstring)
