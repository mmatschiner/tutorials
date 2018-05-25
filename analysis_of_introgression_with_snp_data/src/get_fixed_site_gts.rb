# m_matschiner Wed May 2 12:30:55 CEST 2018

# This script reads a vcf file and determines the
# sites that represent bi-allelic SNS that are fixed
# in the parental populations.

# Get the command-line arguments.
vcf_file_name = ARGV[0]
output_table_file_name = ARGV[1]
parent1_ids = ARGV[2].split(",")
parent2_ids = ARGV[3].split(",")
hybrid_ids = ARGV[4].split(",")
required_completeness_in_parents = ARGV[5].to_f

# Initiate the output.
output_table_file = File.open(output_table_file_name,"w")
string = "scaffold\tpos\tp1\tp2"
parent1_ids.each do |i|
  string << "\t#{i}"
end
hybrid_ids.each do |i|
  string << "\t#{i}"
end
parent2_ids.each do |i|
  string << "\t#{i}"
end
string << "\n"
output_table_file.write(string)
output_table_file.close

# Read the vcf file and analyse the alleles at sites fixed in the two parental species.
output_table_file = File.open(output_table_file_name,"a")
parent1_indices = []
parent2_indices = []
hybrid_indices = []
in_data = false
File.open(vcf_file_name) do |f|
  f.each_line do |l|
    if l[0..5] == "#CHROM"
      in_data = true
      header = l
      header_ary = header.split
      header_ary.size.times do |x|
        if parent1_ids.include?(header_ary[x])
          parent1_indices << x
        elsif parent2_ids.include?(header_ary[x])
          parent2_indices << x
        elsif hybrid_ids.include?(header_ary[x])
          hybrid_indices << x
        end
      end
    elsif in_data
      line_ary = l.split
      chromosome = line_ary[0]
      position = line_ary[1]
      parent1_alleles = []
      parent2_alleles = []
      parent1_indices.each do |x|
        gt = line_ary[x].split(":")[0]
        if gt.include?("/")
          gt1 = gt.split("/")[0]
          gt2 = gt.split("/")[1]
        elsif gt.include?("|")
          gt1 = gt.split("|")[0]
          gt2 = gt.split("|")[1]
        else
          puts "ERROR: Expected genotypes to be separated with either '/' or '|' but found '#{gt}'!"
          exit 1
        end
        parent1_alleles << gt1 unless gt1 == "."
        parent1_alleles << gt2 unless gt2 == "."
      end
      parent2_indices.each do |x|
        gt = line_ary[x].split(":")[0]
        if gt.include?("/")
          gt1 = gt.split("/")[0]
          gt2 = gt.split("/")[1]
        elsif gt.include?("|")
          gt1 = gt.split("|")[0]
          gt2 =gt.split("|")[1]
        else
          puts "ERROR: Expected genotypes to be separated with either '/' or '|' but found '#{gt}'!"
          exit 1
        end
        parent2_alleles << gt1 unless gt1 == "."
        parent2_alleles << gt2 unless gt2 == "."
      end
      site_is_complete_enough = true
      site_is_complete_enough = false if parent1_alleles.size < 2*parent1_indices.size*required_completeness_in_parents
      site_is_complete_enough = false if parent2_alleles.size < 2*parent2_indices.size*required_completeness_in_parents
      if site_is_complete_enough
        if parent1_alleles.uniq.size == 1 and parent2_alleles.uniq.size == 1
          if parent1_alleles[0] != parent2_alleles[0]
            string = "#{chromosome}\t#{position}\t#{parent1_alleles[0]}\t#{parent2_alleles[0]}"
            parent1_indices.each do |x|
              gt = line_ary[x].split(":")[0]
              string << "\t#{gt}"
            end
            hybrid_indices.each do |x|
              gt = line_ary[x].split(":")[0]
              string << "\t#{gt}"
            end
            parent2_indices.each do |x|
              gt = line_ary[x].split(":")[0]
              string << "\t#{gt}"
            end 
            string << "\n"
            output_table_file.write(string)
          end
        end
      end
    end
  end
end
