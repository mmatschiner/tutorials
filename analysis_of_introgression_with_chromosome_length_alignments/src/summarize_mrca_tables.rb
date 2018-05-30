# m_matschiner Tue May 29 15:24:01 CEST 2018

# Get the command-line arguments.
table_dir = ARGV[0]
summary_table_means_file_name = ARGV[1]

# Read the tables with mrca times.
table_file_names = []
Dir.entries(table_dir).each {|e| table_file_names << e unless e[0] == "."}

# Read the table file until taxon ids can be read from a line.
taxon_ids = []
table_file = File.open("#{table_dir}/#{table_file_names[0]}")
table_lines = table_file.readlines
table_lines.size.times do |x|
	if table_lines[x-1].strip == "#ancestors"
		taxon_ids = table_lines[x].split
		break
	end
end

# Prepare a summary table with samples as units.
summary_table = []
taxon_ids.size.times do |x|
	summary_table << []
	taxon_ids.size.times do |y|
		summary_table.last << []
	end
end
table_file_names.each do |t|
	print "Analyzing #{t}..."
	# Read the table file and split it into 100 parts for the individual gene tree's mrca tables.
	table_file = File.open("#{table_dir}/#{t}")
	table_lines = table_file.readlines
	table_file.close
	sets_of_gene_tree_table_lines = []
	gene_tree_table_lines = []
	valid_line_count = 0
	table_lines.size.times do |x|
		if table_lines[x].strip != "" and table_lines[x].strip != "WARNING: ignoring environment value of R_HOME"
			valid_line_count += 1
			if table_lines[x].strip == "#ancestors" and valid_line_count > 1
				sets_of_gene_tree_table_lines << gene_tree_table_lines
				gene_tree_table_lines = [table_lines[x]]
			elsif x == table_lines.size-1
				gene_tree_table_lines << table_lines[x]
				sets_of_gene_tree_table_lines << gene_tree_table_lines
			else
				gene_tree_table_lines << table_lines[x]
			end
		end
	end
	# Analyse each mrca table, and add the mrca times to the summary table.
	sets_of_gene_tree_table_lines.each do |set_lines|
		ancestor_table_lines = []
		node_ages_table_lines = []
		in_ancestors_table = false
		in_node_ages_table = false
		set_lines.each do |l|
			if l.strip == "#ancestors"
				in_ancestors_table = true
			elsif l.strip == "#node_ages"
				in_ancestors_table = false
				in_node_ages_table = true
			elsif in_ancestors_table
				ancestor_table_lines << l
			elsif in_node_ages_table
				node_ages_table_lines << l
			end
		end
		node_ages_table_keys = node_ages_table_lines[0].split
		node_ages_table_values = node_ages_table_lines[1].split
		# column names and row names are always the same, therefore only column
		# names are considered.
		ancestor_table_taxon_ids = ancestor_table_lines[0].split
		# Make sure the taxon ids in this ancestor table are identical to that in the first table of the first file.
		if ancestor_table_taxon_ids != taxon_ids
			raise "ERROR! Taxon Ids differ between tables!"
		end
		if ancestor_table_lines[1..-1].size != ancestor_table_taxon_ids.size
			raise "ERROR! The number of lines in an ancestor table does not match the number of taxa for that table!"
		end
		1.upto(ancestor_table_lines.size-1) do |y|
			line_ary = ancestor_table_lines[y].split[1..-1]
			row_num = y-1
			if line_ary.size != ancestor_table_taxon_ids.size
				raise "ERROR! The number of items on a table line does not match that of the species ids of that table!"
			end
			line_ary.size.times do |x|
				col_num = x
				unless col_num == row_num
					descendant_taxon_ids = [ancestor_table_taxon_ids[col_num],ancestor_table_taxon_ids[row_num]]
					if descendant_taxon_ids[0] == descendant_taxon_ids[1]
						raise "ERROR: Two identical taxon ids found as descendants!"
					end
					mrca_key = line_ary[x]
					unless node_ages_table_keys.count(mrca_key) == 1
						raise "ERROR: Ambiguous node age table key: #{mrca_key}!"
					end
					node_age = node_ages_table_values[node_ages_table_keys.index(mrca_key)].to_f
					taxon_ids_index1 = taxon_ids.index(descendant_taxon_ids[0])
					taxon_ids_index2 = taxon_ids.index(descendant_taxon_ids[1])
					summary_table[taxon_ids_index1][taxon_ids_index2] << node_age
					summary_table[taxon_ids_index2][taxon_ids_index1] << node_age
				end
			end
		end
	end
	puts " done."
end

# Make a reduced table containing only the mean node ages.
summary_table_means = []
taxon_ids.size.times do |y|
	summary_table_means << []
	taxon_ids.size.times do |x|
		if x == y
			summary_table_means[y][x] = 0
		else
			node_ages = summary_table[y][x]
			node_ages_sum = 0
			node_ages.each {|i| node_ages_sum += i}
			node_ages_mean = node_ages_sum/node_ages.size.to_f
			summary_table_means[y][x] = node_ages_mean
		end
	end
end			

# Print the summary table.
cell_width = 10
outstring = "".rjust(cell_width)
taxon_ids.size.times do |x|
	outstring << "#{taxon_ids[x].rjust(cell_width)}"
end
outstring << "\n"
taxon_ids.size.times do |y|
	outstring << "#{taxon_ids[y].ljust(cell_width)}"
	taxon_ids.size.times do |x|
		node_age_string = "%.2f" % summary_table_means[y][x]
		outstring << "#{node_age_string.rjust(cell_width)}"
	end
	outstring << "\n"
end
summary_table_means_file = File.open(summary_table_means_file_name,"w")
summary_table_means_file.write(outstring)
summary_table_means_file.close

# Feedback.
puts "Wrote summary table to file #{summary_table_means_file_name}."
