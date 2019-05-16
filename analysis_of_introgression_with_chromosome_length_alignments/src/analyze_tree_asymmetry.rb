# m_matschiner Fri Nov 30 14:56:54 CET 2018

# See https://en.wikipedia.org/wiki/Binomial_coefficient
def n_choose_k(n, k)
	if k > n
		puts "ERROR: k should never be larger than n!"
		exit 1
	elsif k == 0
		return 1
	elsif k > n/2
		return n_choose_k(n, n-k)
	end
	r = 1
	1.upto(k) do |j|
		r = r * (n + 1 - j) / j
	end
	r
end

# Get the command-line arguments.
tree_file_name = ARGV[0]
table_file_name = ARGV[1]

# Read the tree file.
tree_file = File.open(tree_file_name)
trees = []
tree_file.readlines.each { |l| trees << l.strip.gsub(/:\d+\.\d+/,"").chomp(";") }

# Get the complete set of taxa.
taxa = []
trees.each do |t|
	taxa_this_tree = t.scan(/[A-Za-z][A-Za-z]+/)
	taxa_this_tree.each do |ttt|
		taxa << ttt
	end
end
sorted_taxa = taxa.sort.uniq

# Open the output file.
table_file = File.open(table_file_name, "w")

# Prepare the output string.
outstring = "p1        p2        p3         n(p1,p2)  n(p2,p3)  n(p1,p3)    d         p\n"
table_file.write(outstring)

# Pick a trio for the analysis.
sorted_taxa.size.times do |x|
	taxon1 = sorted_taxa[x]
	(x+1).upto(sorted_taxa.size-1) do |y|
		taxon2 = sorted_taxa[y]
		(y+1).upto(sorted_taxa.size-1) do |z|
			taxon3 = sorted_taxa[z]

			# Analyze all trees.
			n_taxa12_monophyletic = 0
			n_taxa13_monophyletic = 0
			n_taxa23_monophyletic = 0
			trees.each do |t|
				# Get the taxon positions in the tree and convert the tree to a series of levels.
				levels = []
				taxon1_index = nil
				taxon2_index = nil
				taxon3_index = nil
				level = 0
				t.size.times do |pos|
					if taxon1_index == nil
						taxon1_index = pos if t[pos..pos+taxon1.size-1] == taxon1
					end
					if taxon2_index == nil
						taxon2_index = pos if t[pos..pos+taxon2.size-1] == taxon2
					end
					if taxon3_index == nil
						taxon3_index = pos if t[pos..pos+taxon3.size-1] == taxon3
					end
					if t[pos] == "("
						level += 1
					elsif t[pos] == ")"
						level -= 1
					end
					levels << level
				end

				# Make sure that all taxa are included in this tree, if not, move to the next tree.
				if taxon1_index != nil and taxon2_index != nil and taxon3_index != nil

					# Make sure the first and last levels are as expected.
					if levels[0] != 1
						puts "ERROR: Expected the first level to be 1 but found #{levels[-1]}!"
						exit 1
					end
					if levels[-1] != 0
						puts "ERROR: Expected the last level to be 0 but found #{levels[-1]}!"
						exit 1
					end

					# Determine which of the two taxa are monophyletic.
					taxa12_monophyletic = false
					taxa13_monophyletic = false
					taxa23_monophyletic = false
					sorted_taxon_indices = [taxon1_index,taxon2_index,taxon3_index].sort
					if sorted_taxon_indices == [taxon1_index,taxon2_index,taxon3_index]
						if levels[taxon1_index..taxon2_index].min > levels[taxon2_index..taxon3_index].min
							taxa12_monophyletic = true
						elsif levels[taxon1_index..taxon2_index].min < levels[taxon2_index..taxon3_index].min
							taxa23_monophyletic = true
						end
					elsif sorted_taxon_indices == [taxon1_index,taxon3_index,taxon2_index]
						if levels[taxon1_index..taxon3_index].min > levels[taxon3_index..taxon2_index].min
							taxa13_monophyletic = true
						elsif levels[taxon1_index..taxon3_index].min < levels[taxon3_index..taxon2_index].min
							taxa23_monophyletic = true
						end
					elsif sorted_taxon_indices == [taxon2_index,taxon1_index,taxon3_index]
						if levels[taxon2_index..taxon1_index].min > levels[taxon1_index..taxon3_index].min
							taxa12_monophyletic = true
						elsif levels[taxon2_index..taxon1_index].min < levels[taxon1_index..taxon3_index].min
							taxa13_monophyletic = true
						end
					elsif sorted_taxon_indices == [taxon2_index,taxon3_index,taxon1_index]
						if levels[taxon2_index..taxon3_index].min > levels[taxon3_index..taxon1_index].min
							taxa23_monophyletic = true
						elsif levels[taxon2_index..taxon3_index].min < levels[taxon3_index..taxon1_index].min
							taxa13_monophyletic = true
						end
					elsif sorted_taxon_indices == [taxon3_index,taxon1_index,taxon2_index]
						if levels[taxon3_index..taxon1_index].min > levels[taxon1_index..taxon2_index].min
							taxa13_monophyletic = true
						elsif levels[taxon3_index..taxon1_index].min < levels[taxon1_index..taxon2_index].min
							taxa12_monophyletic = true
						end
					elsif sorted_taxon_indices == [taxon3_index,taxon2_index,taxon1_index]
						if levels[taxon3_index..taxon2_index].min > levels[taxon2_index..taxon1_index].min
							taxa23_monophyletic = true
						elsif levels[taxon3_index..taxon2_index].min < levels[taxon2_index..taxon1_index].min
							taxa12_monophyletic = true
						end
					else
						puts "ERROR: Unexpected sequence of sorted indices: #{sorted_taxon_indices[0]}, #{sorted_taxon_indices[2]}, #{sorted_taxon_indices[3]}!"
						exit 1
					end

					# Add to overall numbers of monophylies.
					n_taxa12_monophyletic += 1 if taxa12_monophyletic
					n_taxa13_monophyletic += 1 if taxa13_monophyletic
					n_taxa23_monophyletic += 1 if taxa23_monophyletic
				end
			end

			# Sort the results according to numbers of monophylies.
			sorted_n_monophyletic = [n_taxa12_monophyletic,n_taxa13_monophyletic,n_taxa23_monophyletic].sort.reverse

			# Get the numbers of concordant and discordant trees.
			n_obs_concordant = sorted_n_monophyletic[0]
			n_obs_discordant = sorted_n_monophyletic[1] + sorted_n_monophyletic[2]

			# Calculate the d statistic.
			if sorted_n_monophyletic[1] + sorted_n_monophyletic[2] == 0
				d = 0
			else
				d = (sorted_n_monophyletic[1] - sorted_n_monophyletic[2])/(sorted_n_monophyletic[1] + sorted_n_monophyletic[2]).to_f
			end

			# Test if the imbalance between discordant trees is greater than expected by chance.
			pr = 0
			n = n_obs_discordant
			k = (sorted_n_monophyletic[2]+1)
			k.times do |i|
				pr += Math.exp(Math.log(n_choose_k(n,i)) + Math.log(0.5**i) + Math.log((1-0.5)**(n-i)))
			end
			pr = pr * 2 # Make the test two-tailed instead of one-tailed.
			pr = 1.0 if pr > 1.0 # This is possible when both alternative topologies are equally frequent.

			# Prepare the output string for this trio.
			outstring = ""
			if sorted_n_monophyletic == [n_taxa12_monophyletic,n_taxa13_monophyletic,n_taxa23_monophyletic]
				outstring << "#{taxon2.ljust(10)}#{taxon1.ljust(10)}#{taxon3.ljust(10)} #{sorted_n_monophyletic[0].to_s.rjust(8)}  #{sorted_n_monophyletic[1].to_s.rjust(8)}  #{sorted_n_monophyletic[2].to_s.rjust(8)} #{('%.3f' % d).rjust(8)}  #{('%.8f' % pr).rjust(12)}\n"
			elsif sorted_n_monophyletic == [n_taxa12_monophyletic,n_taxa23_monophyletic,n_taxa13_monophyletic]
				outstring << "#{taxon1.ljust(10)}#{taxon2.ljust(10)}#{taxon3.ljust(10)} #{sorted_n_monophyletic[0].to_s.rjust(8)}  #{sorted_n_monophyletic[1].to_s.rjust(8)}  #{sorted_n_monophyletic[2].to_s.rjust(8)} #{('%.3f' % d).rjust(8)}  #{('%.8f' % pr).rjust(12)}\n"
			elsif sorted_n_monophyletic == [n_taxa13_monophyletic,n_taxa12_monophyletic,n_taxa23_monophyletic]
				outstring << "#{taxon3.ljust(10)}#{taxon1.ljust(10)}#{taxon2.ljust(10)} #{sorted_n_monophyletic[0].to_s.rjust(8)}  #{sorted_n_monophyletic[1].to_s.rjust(8)}  #{sorted_n_monophyletic[2].to_s.rjust(8)} #{('%.3f' % d).rjust(8)}  #{('%.8f' % pr).rjust(12)}\n"
			elsif sorted_n_monophyletic == [n_taxa13_monophyletic,n_taxa23_monophyletic,n_taxa12_monophyletic]
				outstring << "#{taxon1.ljust(10)}#{taxon3.ljust(10)}#{taxon2.ljust(10)} #{sorted_n_monophyletic[0].to_s.rjust(8)}  #{sorted_n_monophyletic[1].to_s.rjust(8)}  #{sorted_n_monophyletic[2].to_s.rjust(8)} #{('%.3f' % d).rjust(8)}  #{('%.8f' % pr).rjust(12)}\n"
			elsif sorted_n_monophyletic == [n_taxa23_monophyletic,n_taxa12_monophyletic,n_taxa13_monophyletic]
				outstring << "#{taxon3.ljust(10)}#{taxon2.ljust(10)}#{taxon1.ljust(10)} #{sorted_n_monophyletic[0].to_s.rjust(8)}  #{sorted_n_monophyletic[1].to_s.rjust(8)}  #{sorted_n_monophyletic[2].to_s.rjust(8)} #{('%.3f' % d).rjust(8)}  #{('%.8f' % pr).rjust(12)}\n"
			elsif sorted_n_monophyletic == [n_taxa23_monophyletic,n_taxa13_monophyletic,n_taxa12_monophyletic]
				outstring << "#{taxon2.ljust(10)}#{taxon3.ljust(10)}#{taxon1.ljust(10)} #{sorted_n_monophyletic[0].to_s.rjust(8)}  #{sorted_n_monophyletic[1].to_s.rjust(8)}  #{sorted_n_monophyletic[2].to_s.rjust(8)} #{('%.3f' % d).rjust(8)}  #{('%.8f' % pr).rjust(12)}\n"
			end

			# Write the output.
			table_file.write(outstring)
		end
	end
end
