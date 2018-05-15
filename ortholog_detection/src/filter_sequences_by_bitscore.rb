# m_matschiner Mon May 14 17:09:17 CEST 2018

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
relative_bitscore_threshold = ARGV[2].to_f

# Create the output directory if it does not exist yet.
unless Dir.exists?(alignment_directory_out)
	Dir.mkdir(alignment_directory_out)
end

# Collect names of nucleotide fasta files in the input directory.
dir_entries_in = Dir.entries(alignment_directory_in).sort
filenames_in = []
dir_entries_in.each {|e| filenames_in << e if e.match(/.*_nucl.fasta/)}

# Do for each fasta file in the input directory.
n_excluded_overall = 0
filenames_in.each do |f|

	# Feedback.
	print "Analysing file #{f}..."

	# Read the fasta file.
	n_excluded_this_alignment = 0
	fasta_file = File.open("#{alignment_directory_in}/#{f}")
	fasta_lines = fasta_file.readlines
	fasta_ids = []
	fasta_bitscores = []
	fasta_hits = []
	fasta_seqs = []
	fasta_lines.each do |l|
		if l[0] == ">"
			fasta_header = l[1..-1].strip
			fasta_ids << fasta_header.sub(/\[.*\]/,"")
			fasta_header.match(/bitscore=(.+?)[,\]]/)
			fasta_bitscores << $1
			fasta_header.match(/nhits=(.+?)[,\]]/)
			fasta_hits << $1.to_i
			fasta_seqs << ""
		else
			fasta_seqs.last << l.strip
		end
	end

	# Determine the minimum bitscore needed to qualify as ortholog.
	# It is calculates as relative_bitscore_threshold * maximum bitscore (excluding the reference).
	maximum_ingroup_bitscore = 0
	n_known_bases_in_best_seq = 0
	1.upto(fasta_bitscores.size-1) do |x|
		bitscore = fasta_bitscores[x]
		unless bitscore == "None"
			if maximum_ingroup_bitscore < bitscore.to_f
				maximum_ingroup_bitscore = bitscore.to_f
				n_known_bases_in_best_seq = fasta_seqs[x].upcase.count("A") + fasta_seqs[x].upcase.count("C") + fasta_seqs[x].upcase.count("G") + fasta_seqs[x].upcase.count("T")
			end
		end
	end

	# Prepare the string for a new fasta file.
	new_fasta_string = ""
	fasta_ids.size.times do |x|
		new_fasta_string << ">#{fasta_ids[x]}\n"
		if fasta_bitscores[x] == "None"
			fasta_seqs[0].size.times {new_fasta_string << "-"}
			new_fasta_string << "\n"
		else
			n_known_bases_in_this_seq = fasta_seqs[x].upcase.count("A") + fasta_seqs[x].upcase.count("C") + fasta_seqs[x].upcase.count("G") + fasta_seqs[x].upcase.count("T")
			if n_known_bases_in_this_seq < n_known_bases_in_best_seq
				# This correction attempts to penalize sequences with more missing data to a lesser degree.
				# It reduces the required bitscore according to sequence completeness, but sequence
				# completeness is only accounted for by 50%. Thus, if for example a sequence is 60%
				# complete, the required bitscore is 80% or the original bitscore.
				correction = ((n_known_bases_in_this_seq/n_known_bases_in_best_seq.to_f)+1.0)/2.0
			else
				correction = 1.0
			end
			if fasta_bitscores[x].to_f < correction * relative_bitscore_threshold * maximum_ingroup_bitscore
				fasta_seqs[0].size.times {new_fasta_string << "-"}
				new_fasta_string << "\n"
				n_excluded_this_alignment += 1
				n_excluded_overall += 1
			else
				new_fasta_string << "#{fasta_seqs[x]}\n"
			end
		end
	end

	# Write the corrected fasta file.
	new_fasta_file = File.open("#{alignment_directory_out}/#{f}","w")
	new_fasta_file.write(new_fasta_string)

	# Feedback.
	print " done."
	if n_excluded_this_alignment == 1
		puts " Removed 1 sequence."
	elsif n_excluded_this_alignment > 1
		puts " Removed #{n_excluded_this_alignment} sequences."
	else
		puts
	end

end

# Feedback.
puts "Removed #{n_excluded_overall} sequences overall."
