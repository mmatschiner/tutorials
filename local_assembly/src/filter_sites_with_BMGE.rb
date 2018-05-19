# m_matschiner Tue May 15 13:35:40 CEST 2018

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
	def sample_variance
		if self.length == 0
			nil
		else
			m = self.mean
			sum = self.inject(0){|accum, i| accum +(i-m)**2 }
			sum/(self.length - 1).to_f
		end
	end
	def standard_deviation
		if self.length == 0
			nil
		else
			return Math.sqrt(self.sample_variance)
		end
	end
	def median
		if self.length == 0
			nil
		else
			sorted_array = self.sort
			if self.size.modulo(2) == 1 
				sorted_array[self.size/2]
			else
				(sorted_array[(self.size/2)-1]+sorted_array[self.size/2])/2.0
			end
		end
	end
	def rough_percentile(p)
		if self.length == 0
			nil
		else
			raise "p should be between 0 and 1!" if p < 0.0 or p > 1.0
			sorted_array = self.sort
			index = (sorted_array.size*p).floor
			if index > sorted_array.size-1
				index = sorted_array.size-1
			elsif index < 0
				index = 0
			end
			sorted_array[index]
		end
	end
	def lower_quartile(method = 1)
		# The first and second quartile computing methods on http://en.wikipedia.org/wiki/Quartile are implemented.
		sorted_array = self.sort
		lower_half = []
		if method == 1
			sorted_array.each {|i| lower_half << i if i < self.median}
		elsif method == 2
			sorted_array.each {|i| lower_half << i if i <= self.median}
		else
			raise "Unknown quartile method!"
		end
		lower_half.median
	end
	def upper_quartile(method = 1)
		# The first and second quartile computing methods on http://en.wikipedia.org/wiki/Quartile are implemented.
		sorted_array = self.sort
		upper_half = []
		if method == 1
			sorted_array.each {|i| upper_half << i if i > self.median}
		elsif method == 2
			sorted_array.each {|i| upper_half << i if i >= self.median}
		else
			raise "Unknown quartile method!"
		end
		upper_half.median
	end
	def inter_quartile_range(method = 1)
		# The first and second quartile computing methods on http://en.wikipedia.org/wiki/Quartile are implemented.
		if method == 1
			self.upper_quartile(method = 1) - self.lower_quartile(method = 1)
		elsif method == 2
			self.upper_quartile(method = 2) - self.lower_quartile(method = 2)
		else
			raise "Unknown quartile method!"
		end
	end
	def lower_whisker(method = 1)
		# The first and second quartile computing methods on http://en.wikipedia.org/wiki/Quartile are implemented.
		sorted_array = self.sort
		lw = 0
		if method == 1
			(sorted_array.size-1).downto(0) {|x| lw = sorted_array[x] if sorted_array[x] >= self.lower_quartile(method = 1) - 1.5*self.inter_quartile_range(method = 1)}
		elsif method == 2
			(sorted_array.size-1).downto(0) {|x| lw = sorted_array[x] if sorted_array[x] >= self.lower_quartile(method = 2) - 1.5*self.inter_quartile_range(method = 2)}
		else
			raise "Unknown quartile method!"
		end
		return lw
	end
	def upper_whisker(method = 1)
		# The first and second quartile computing methods on http://en.wikipedia.org/wiki/Quartile are implemented.
		sorted_array = self.sort
		uw = 0
		if method == 1
			0.upto(sorted_array.size-1) {|x| uw = sorted_array[x] if sorted_array[x] <= self.upper_quartile(method = 1) + 1.5*self.inter_quartile_range(method = 1)}
		elsif method == 2
			0.upto(sorted_array.size-1) {|x| uw = sorted_array[x] if sorted_array[x] <= self.upper_quartile(method = 2) + 1.5*self.inter_quartile_range(method = 2)}
		else
			raise "Unknown quartile method!"
		end
		return uw
	end
	def hpd_lower(proportion)
		raise "The interval should be between 0 and 1!" if proportion >= 1 or proportion <= 0
		sorted_array = self.sort
		hpd_index = 0
		min_range = sorted_array[-1]
		diff = (proportion*self.size).round
		(self.size-diff).times do |i|
			min_value = sorted_array[i]
			max_value = sorted_array[i+diff-1]
			range = max_value - min_value
			if range < min_range
				min_range = range
				hpd_index = i
			end
		end
		sorted_array[hpd_index]
	end
	def hpd_upper(proportion)
		raise "The interval should be between 0 and 1!" if proportion >= 1 or proportion <= 0
		sorted_array = self.sort
		hpd_index = 0
		min_range = sorted_array[-1]
		diff = (proportion*self.size).round
		(self.size-diff).times do |i|
			min_value = sorted_array[i]
			max_value = sorted_array[i+diff-1]
			range = max_value - min_value
			if range < min_range
				min_range = range
				hpd_index = i
			end
		end
		sorted_array[hpd_index+diff-1]
	end
end

# Get the command line arguments.
alignment_directory_in = ARGV[0].chomp("/")
alignment_directory_out = ARGV[1].chomp("/")
bmge_path = ARGV[2]
gap_rate_cut_off = ARGV[3].to_f
entropy_like_score = ARGV[4].to_f

# Create the output directory if it does not exist yet.
unless Dir.exists?(alignment_directory_out)
	Dir.mkdir(alignment_directory_out)
end

# Collect names of nucleotide fasta files in the input directory.
dir_entries_in = Dir.entries(alignment_directory_in).sort
filenames_in = []
dir_entries_in.each do |e|
	if e.match(/.*_nucl.fasta/)
		filenames_in << e
	end
end

# Do for each fasta file in the input directory.
n_discarded_sites_overall = 0
filenames_in.each do |f|

	# Feedback.
	print "Analysing file #{f}..."

	# Read the fasta file.
	fasta_file = File.open("#{alignment_directory_in}/#{f}")
	fasta_lines = fasta_file.readlines
	fasta_ids = []
	fasta_seqs = []
	fasta_lines.each do |l|
		if l[0] == ">"
			fasta_ids << l[1..-1].strip.gsub(/\[.+\]/,"")
			fasta_seqs << ""
		else
			fasta_seqs.last << l.strip
		end
	end

	# Make sure that all sequences are the same length.
	fasta_seqs.each do |s|
		if s.size != fasta_seqs[0].size
			puts "ERROR: Sequences have different lengths!"
			exit 1
		end
	end

	# Make sure that the file actually contains non-zero length sequences.
	discarded_sites = []
	if fasta_seqs[0].size > 0

		# Run BMGE.
		system("java -jar #{bmge_path} -i #{alignment_directory_in}/#{f} -t DNA -g 0.99 -oh tmp.html > /dev/null")
		
		# Read the BMGE HTML file, and then delete it.
		html_file = File.open("tmp.html")
		html_lines = html_file.readlines
		File.delete("tmp.html")

		# Parse the HTML file.
		smoothed_entropies = []
		gap_rates = []
		in_table = false
		html_lines.each do |l|
			if l.match(/ch\.\s+entropy\s+smooth\. entr.\s+gap rate/)
				in_table = true
			elsif l.match(/<\/span>/)
				in_table = false
			elsif in_table
				smoothed_entropies << l.split[2].to_f
				gap_rates << l.split[3].to_f
			end
		end

		# Make sure that entropy scores and gap rates are found for each site.
		if fasta_seqs[0].size != smoothed_entropies.size or fasta_seqs[0].size != gap_rates.size
			raise "Alignment scores were not found for all positions!"
		end
		
		# Determine the sites to be discarded.
		(fasta_seqs[0].size/3).times do |codon_pos|
			keep_codon = true
			if smoothed_entropies[3*codon_pos] > entropy_like_score or gap_rates[3*codon_pos] > gap_rate_cut_off
				keep_codon = false
			elsif smoothed_entropies[3*codon_pos+1] > entropy_like_score or gap_rates[3*codon_pos+1] > gap_rate_cut_off
				keep_codon = false
			elsif smoothed_entropies[3*codon_pos+2] > entropy_like_score or gap_rates[3*codon_pos+2] > gap_rate_cut_off
				keep_codon = false
			end
			if keep_codon == false
				discarded_sites << 3*codon_pos
				discarded_sites << 3*codon_pos+1
				discarded_sites << 3*codon_pos+2
			end
		end

	end
	
	# Prepare the string for a new fasta file.
	new_fasta_string = ""
	fasta_ids.size.times do |x|
		new_fasta_string << ">#{fasta_ids[x]}\n"
		fasta_seqs[x].size.times do |pos|
			unless discarded_sites.include?(pos)
				new_fasta_string << fasta_seqs[x][pos]
			end
		end
		new_fasta_string << "\n"
	end

	# Write the corrected fasta file.
	new_fasta_file = File.open("#{alignment_directory_out}/#{f}","w")
	new_fasta_file.write(new_fasta_string)

	# Feedback.
	if discarded_sites.size > 0
		puts " done. Removed #{discarded_sites.size} sites."
		n_discarded_sites_overall += discarded_sites.size
	else
		puts " done."
	end

end

# Feedback.
puts "Removed #{n_discarded_sites_overall} sites overall."
