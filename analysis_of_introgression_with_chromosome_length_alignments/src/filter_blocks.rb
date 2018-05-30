# m_matschiner Mon May 28 10:44:03 CEST 2018

# Get the command line arguments.
input_file_name = ARGV[0]
max_proportion_of_missing_data = ARGV[1].to_f
min_number_of_pi_sites = ARGV[2].to_i
min_bootstrap_support_value = ARGV[3].to_f
min_p_value_for_recombination = ARGV[4].to_f

# Load the input file.
input_file = File.open(input_file_name)
input_lines = input_file.readlines

# Filter the lines.
count = 0
puts input_lines[0]
input_lines[1..-1].each do |l|
	count += 1
	line_ary = l.split
	block_id = line_ary[0]
	proportion_of_missing_data = line_ary[1].to_f
	number_of_pi_sites = line_ary[2].to_i
	bootstrap_support_value = line_ary[3].to_f
	p_value_for_recombination = line_ary[4].to_f

	good_line = true
	good_line = false if proportion_of_missing_data > max_proportion_of_missing_data
	good_line = false if number_of_pi_sites < min_number_of_pi_sites
	good_line = false if bootstrap_support_value < min_bootstrap_support_value
	good_line = false if p_value_for_recombination < min_p_value_for_recombination
	puts l if good_line
end