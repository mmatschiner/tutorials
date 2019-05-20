# m_matschiner Wed May 23 22:36:31 CEST 2018

# Get the command-line arguments.
unimputed_vcf_file_name = ARGV[0]
imputed_vcf_file_name = ARGV[1]
masked_vcf_file_name = ARGV[2]

# Open the output file.
masked_vcf_file = File.open(masked_vcf_file_name,"w")

# Open the input files.
unimputed_vcf_file = File.open(unimputed_vcf_file_name)
imputed_vcf_file = File.open(imputed_vcf_file_name)

#unimputed_gzvcf_file.each.zip(imputed_gzvcf_file.each).each do |l_unimp, l_imp|
until unimputed_vcf_file.eof do
  current_lines = []
  [unimputed_vcf_file,imputed_vcf_file].each do |f|
    current_lines << f.readline
  end
  l_unimp = current_lines[0]
  l_imp = current_lines[1]
  ary_unimp = l_unimp.split
  ary_imp = l_imp.split
  unless ary_unimp.size == ary_imp.size
    puts "ERROR: The two input files have different numbers of columns!"
    puts "The number of columns in the unimputed file: #{ary_unimp.size}."
    puts "The number of columns in the imputed file: #{ary_imp.size}"
    puts
    puts "The column values in the unimputed file:"
    puts ary_unimp
    puts
    puts "The column values in the imputed file:"
    puts ary_imp
    puts
    puts "The line in the unimputed file:"
    puts l_unimp
    puts
    puts "The line in the imputed file:"
    puts l_imp
    exit 1
  end
  ary_unimp.size.times do |x|
    ary_imp[x] = ".|." if ary_unimp[x].split(":")[0] == "./."
  end
  str_masked = ""
  ary_imp.each do |i|
    str_masked << "#{i}\t"
  end
  str_masked.strip!
  str_masked << "\n"
  masked_vcf_file.write(str_masked) 
end
