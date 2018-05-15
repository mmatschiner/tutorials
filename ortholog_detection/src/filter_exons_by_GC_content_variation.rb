# m_matschiner Tue May 15 16:26:05 CEST 2018

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
maximum_gc_content_variation = ARGV[2].to_f

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

	# Read the fasta file.
	fasta_file = File.open("#{alignment_directory_in}/#{f}")
	fasta_lines = fasta_file.readlines
	fasta_ids = []
	new_fasta_ids = []
	fasta_seqs = []
	fasta_lines.each do |l|
		if l[0] == ">"
			fasta_ids << l[1..-1].strip
			new_fasta_ids << l[1..-1].strip
			fasta_seqs << ""
		else
			fasta_seqs.last << l.strip
		end
	end

	# Get the variation in GC content per codon position.
	gc_contents = []
	fasta_ids.size.times do |x|
		unless fasta_seqs[x].match(/^-+$/)
			bases = []
			fasta_seqs[x].size.times do |pos|
				bases << fasta_seqs[x][pos]
			end
			base_count = bases.count("A") + bases.count("C") + bases.count("G") + bases.count("T")
			gc_content = (bases.count("C") + bases.count("G"))/base_count.to_f
			gc_contents << gc_content
			new_fasta_ids[x] = "#{new_fasta_ids[x]}[&gc=#{gc_content.round(3)}]"
		end
	end

	if gc_contents.standard_deviation <= maximum_gc_content_variation

		# Prepare the string for a new fasta file.
		new_fasta_string = ""
		fasta_ids.size.times do |x|
			new_fasta_string << ">#{new_fasta_ids[x]}\n"
			new_fasta_string << "#{fasta_seqs[x]}\n"
		end

		# Write the new fasta file.
		new_fasta_file = File.open("#{alignment_directory_out}/#{f}","w")
		new_fasta_file.write(new_fasta_string)

	else

		puts "Removed file #{f}."
		n_excluded_overall += 1

	end

end

# Feedback.
if n_excluded_overall == 1
	puts "Removed #{n_excluded_overall} exon overall."
else
	puts "Removed #{n_excluded_overall} exons overall."
end
