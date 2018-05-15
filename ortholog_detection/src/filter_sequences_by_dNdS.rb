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
omega_threshold = ARGV[2].to_f

# Create the output directory if it does not exist yet.
unless Dir.exists?(alignment_directory_out)
	Dir.mkdir(alignment_directory_out)
end

# Prepare a codeml.ctl file.
control_string = ""
control_string << "      seqfile = tmp.phy   * sequence data file name\n"
control_string << "      outfile = results.txt   * main result file name\n"
control_string << "\n"
control_string << "        noisy = 0      * 0,1,2,3,9: how much rubbish on the screen\n"
control_string << "      verbose = 0      * 1:detailed output\n"
control_string << "      runmode = -2     * -2:pairwise\n"
control_string << "\n"
control_string << "      seqtype = 1      * 1:codons\n"
control_string << "    CodonFreq = 1      * 0:equal, 1:F1X4, 2:F3X4, 3:F61\n"
control_string << "        model = 0      *\n"
control_string << "      NSsites = 0      *\n"
control_string << "        icode = 0      * 0:universal code\n"
control_string << "\n"
control_string << "    fix_kappa = 0      * 1:kappa fixed, 0:kappa to be estimated\n"
control_string << "        kappa = 1      * initial or fixed kappa\n"
control_string << "\n"
control_string << "    fix_omega = 0      * 1:omega fixed, 0:omega to be estimated\n"
control_string << "        omega = 0.5    * initial omega value\n"

# Write the control file to the current directory.
control_file = File.open("codeml.ctl","w")
control_file.write(control_string)
control_file.close

# Collect names of nucleotide fasta files in the input directory.
dir_entries_in = Dir.entries("#{alignment_directory_in}").sort
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

	# For each non-empty sequence except the reference, write a reduced fasta 
	# file with only that sequence and the reference, then run codeml to estimate
	# omega for this 2-sequence comparison.
	reference_empty = false
	species_ids_compared_with_reference = []
	omega_estimates = []
	reference_phylip_string = "2 #{fasta_seqs[0].size}\n"
	reference_empty = true if fasta_seqs[0].match(/^-+$/) or fasta_seqs[0].strip == ""
	reference_phylip_string << "#{fasta_ids[0].ljust(12)}#{fasta_seqs[0]}\n"
	unless reference_empty
		1.upto(fasta_ids.size-1) do |x|
			unless fasta_seqs[x].match(/^-+$/)

				# Memorize that this species has been used for dNdS estimation.
				species_ids_compared_with_reference << fasta_ids[x]

				# Add the sequence of this taxon to the fasta string that already contains
				# the reference sequence.
				reduced_phylip_string = Marshal.load(Marshal.dump(reference_phylip_string))
				reduced_phylip_string << "#{fasta_ids[x].ljust(12)}#{fasta_seqs[x]}\n"

				# Write the reduced fasta file with these two sequences.
				reduced_phylip_file = File.open("tmp.phy","w")
				reduced_phylip_file.write(reduced_phylip_string)
				reduced_phylip_file.close

				# Run codeml with codeml.ctl and file tmp.fasta
				system("codeml > /dev/null")
				
				# Read the codeml results file.
				result_file = File.open("results.txt")
				result_lines = result_file.readlines

				# Delete codeml output files.
				File.delete("results.txt")
				File.delete("2ML.dN")
				File.delete("2ML.dS")
				File.delete("2ML.t")
				File.delete("2NG.dN")
				File.delete("2NG.dS")
				File.delete("2NG.t")
				File.delete("rst")
				File.delete("rub")
				File.delete("rst1")

				# Delete the reduced phylip file.
				File.delete("tmp.phy")

				# Make sure codeml has run without errors.
				unless result_lines.size == 89
					raise "Unexpected length of codeml result file!"
				end

				# Parse the result file and store the omega estimate.
				result_lines[-1].match(/dN\/dS=\s+(\d+\.\d+)/)
				if $1 == nil
					puts "The omega estimate could not be parsed!"
					puts result_lines
					exit
				else
					omega_estimate = $1.to_f
					omega_estimates << omega_estimate
				end
				new_fasta_ids[x] = "#{fasta_ids[x]}[&omega=#{omega_estimate}]"

			end
		end
	end

	# Prepare the string for a new fasta file.
	new_fasta_string = ""
	fasta_ids.size.times do |x|
		new_fasta_string << ">#{new_fasta_ids[x]}\n"
		if species_ids_compared_with_reference.include?(fasta_ids[x]) and omega_estimates[species_ids_compared_with_reference.index(fasta_ids[x])] > omega_threshold
			fasta_seqs[0].size.times {new_fasta_string << "-"}
			new_fasta_string << "\n"
			n_excluded_this_alignment += 1
			n_excluded_overall += 1
		else
			new_fasta_string << "#{fasta_seqs[x]}\n"
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

# Delete the codeml.ctl file.
File.delete("codeml.ctl")

# Feedback.
puts "Removed #{n_excluded_overall} sequences overall."
