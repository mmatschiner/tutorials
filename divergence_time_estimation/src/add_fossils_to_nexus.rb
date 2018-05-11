# m_matschiner Wed May 9 01:28:20 CEST 2018

# Get the command-line arguments.
nexus_infile_name = ARGV[0]
fossil_list_file_name = ARGV[1]
nexus_out_file_name = ARGV[2]

# Read the file with the list of fossil ids and ages.
fossil_ids = []
fossil_ages = []
fossil_list_file = File.open(fossil_list_file_name)
fossil_list_lines = fossil_list_file.readlines
fossil_list_lines[1..-1].each do |l|
	line_ary = l.split
	fossil_ids << line_ary[0]
	fossil_ages << (line_ary[1].to_f + rand * (line_ary[2].to_f-line_ary[1].to_f)).round(2)
end

# Read the nexus input file and compose an output string.
nexus_infile = File.open(nexus_infile_name)
nexus_lines = nexus_infile.readlines
nexus_outstring = ""
in_marker = false
extant_seq_length = 0
nexus_lines.each do |l|
	if l.match(/([Nn][Tt][Aa][Xx]\s*=\s*)(\d+)(\s+)/)
		nexus_outstring << l.sub("#{$1}#{$2}#{$3}","#{$1}#{$2.to_i + fossil_ids.size}#{$3}")
	elsif l[0] == "["
		in_marker = true
		nexus_outstring << l
	elsif in_marker
		if l.strip == "" or l.strip == ";"
			in_marker = false
			fossil_ids.size.times do |x|
				fossil_id_with_age = "#{fossil_ids[x]}_#{fossil_ages[x]}"
				nexus_outstring << "#{fossil_id_with_age.ljust(40)}"
				extant_seq_length.times {nexus_outstring << "?"}
				nexus_outstring << "\n"
			end
			nexus_outstring << l
		else
			line_ary = l.split
			extant_id = line_ary[0].strip
			extant_seq = line_ary[1].strip
			extant_seq_length = extant_seq.size
			extant_id_with_age = "#{extant_id}_0.00"
			nexus_outstring << "#{extant_id_with_age.ljust(40)}#{extant_seq}\n"
		end	
	else
		nexus_outstring << l
	end
end

# Write the output nexus file.
nexus_outfile = File.open(nexus_out_file_name,"w")
nexus_outfile.write(nexus_outstring)
