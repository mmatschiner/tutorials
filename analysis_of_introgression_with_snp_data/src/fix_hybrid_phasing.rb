# m_matschiner Thu May 3 18:11:54 CEST 2018

# Get the command-line arguments.
invcf_header_file_name = ARGV[0]
invcf_body_file_name = ARGV[1]
outvcf_body_file_name = ARGV[2]
samples_table_file_name = ARGV[3]
trios = ARGV[4..-1]

# Read the table with sample information.
samples_table_file = File.open(samples_table_file_name)
samples_table_lines = samples_table_file.readlines
samples_table_species_ids = []
samples_table_specimen_ids = []
samples_table_lines.each do |l|
  line_ary = l.split
  samples_table_specimen_ids << line_ary[0]
  samples_table_species_ids << line_ary[1]
end

# Read the header line to get the indices of samples.
invcf_header_file = File.open(invcf_header_file_name)
header_line_ary = invcf_header_file.readlines.last.split

# Get the indices of samples that are part of trios.
parent1_indices = []
parent2_indices = []
hybrid_indices = []
trios.each do |t|
  trio_ary = t.split(",")
  parent1_species_id = trio_ary[0]
  parent2_species_id = trio_ary[2]
  hybrid_species_id = trio_ary[1]
  parent1_specimen_ids_this_trio = []
  parent2_specimen_ids_this_trio = []
  hybrid_specimen_ids_this_trio = []
  parent1_indices_this_trio = []
  parent2_indices_this_trio = []
  hybrid_indices_this_trio = []
  samples_table_species_ids.size.times do |x|
    if samples_table_species_ids[x] == parent1_species_id
      parent1_specimen_id = samples_table_specimen_ids[x]
      parent1_specimen_ids_this_trio << parent1_specimen_id
      if header_line_ary.include?(parent1_specimen_id)
        parent1_indices_this_trio << header_line_ary.index(parent1_specimen_id)
      else
        puts "ERROR: Specimen id #{parent1_specimen_id} could not be found in the vcf header!"
        exit
      end
    elsif samples_table_species_ids[x] == parent2_species_id
      parent2_specimen_id = samples_table_specimen_ids[x]
      parent2_specimen_ids_this_trio << parent2_specimen_id
      if header_line_ary.include?(parent2_specimen_id)
        parent2_indices_this_trio << header_line_ary.index(parent2_specimen_id)
      else
        puts "ERROR: Specimen id #{parent2_specimen_id} could not be found in the vcf header!"
        exit
      end
    elsif samples_table_species_ids[x] == hybrid_species_id
      hybrid_specimen_id = samples_table_specimen_ids[x]
      hybrid_specimen_ids_this_trio << hybrid_specimen_id
      if header_line_ary.include?(hybrid_specimen_id)
        hybrid_indices_this_trio << header_line_ary.index(hybrid_specimen_id)
      else
        puts "ERROR: Specimen id #{hybrid_specimen_id} could not be found in the vcf header!"
        exit
      end
    end
  end
  parent1_indices << parent1_indices_this_trio
  parent2_indices << parent2_indices_this_trio
  hybrid_indices << hybrid_indices_this_trio
end

# Open the input and output vcf files.
invcf_body_file = File.open(invcf_body_file_name,"r")
outvcf_body_file = File.open(outvcf_body_file_name,"w")

# Read each line of the input vcf file, modify it, and write it to the output vcf file.
n_modified_lines = 0
invcf_body_file.each do |l|
  line_ary = l.split
  # First check if any of the hybrids is heterozygous, if not, just write the line to the output as it is.
  heterozygous_hybrid_found = false
  trios.size.times do |z|
    hybrid_indices_this_trio = hybrid_indices[z]
    hybrid_indices_this_trio.each do |i|
      hybrid_gt = line_ary[i]
      # Make sure the gt format is as expected (e.g. '1|1').
      if hybrid_gt.include?(":")
        puts "ERROR: Expected genotype to contain no colons, but found #{hybrid_gt.include}!"
        exit
      end
      if hybrid_gt.include?("|") == false
        puts "ERROR: Expected genotype to contain the pipe symbol as a separator, but found #{hybrid_gt.include}!"
        exit
      end
      hybrid_gt_ary = hybrid_gt.split("|")
      hybrid_gt1 = hybrid_gt_ary[0]
      hybrid_gt2 = hybrid_gt_ary[1]
      if hybrid_gt1 != hybrid_gt2
        heterozygous_hybrid_found = true
        break
      end
    end
  end
  if heterozygous_hybrid_found    
    trios.size.times do |z|
      parent1_indices_this_trio = parent1_indices[z]
      parent2_indices_this_trio = parent2_indices[z]
      hybrid_indices_this_trio = hybrid_indices[z]
      # Get the most frequent allele in parent1.
      parent1_gts = []
      parent1_indices_this_trio.each do |i|
        parent1_gts << line_ary[i]
      end
      parent1_alleles = []
      parent1_gts.each do |g|
        parent1_gt_ary = g.split("|")
        parent1_alleles << parent1_gt_ary[0] unless parent1_gt_ary[0] == "."
      end
      next if parent1_alleles == []
      parent1_majority_allele = parent1_alleles[0]
      parent1_alleles.each do |a|
        if parent1_alleles.count(a) > parent1_alleles.count(parent1_majority_allele)
          parent1_majority_allele = a
        end
      end
      # Get the most frequent allele in parent2.
      parent2_gts = []
      parent2_indices_this_trio.each do |i|
        parent2_gts << line_ary[i]
      end
      parent2_alleles =[]
      parent2_gts.each do |g|
        parent2_gt_ary = g.split("|")
        parent2_alleles<< parent2_gt_ary[0] unless parent2_gt_ary[0] == "."
      end
      next if parent2_alleles == []
      parent2_majority_allele = parent2_alleles[0]
      parent2_alleles.each do |a|
        if parent2_alleles.count(a) > parent2_alleles.count(parent2_majority_allele)
          parent2_majority_allele = a
        end
      end
      # Don't change hybrid gts for this trio if the parent's most frequent alleles are identical.
      next if parent1_majority_allele == parent2_majority_allele
      # If parent's most frequent alleles differ, check the hybrid gts.
      hybrid_indices_this_trio.each do |i|
        hybrid_gt = line_ary[i]
        hybrid_gt_ary = hybrid_gt.split("|")
        hybrid_allele1 = hybrid_gt_ary[0]
        hybrid_allele2 = hybrid_gt_ary[1]
        next if hybrid_allele1 == hybrid_allele2
        next if [parent1_majority_allele,parent2_majority_allele].include?(hybrid_allele1) == false
        next if [parent1_majority_allele,parent2_majority_allele].include?(hybrid_allele2) == false
        # If none of the above criteria were met, this must mean that the hybrid is heterozygous for the two
        # parental majority alleles. If so, just make the first allele the allele of parent1 and the second
        # allele the allele of parent2.
        n_modified_lines += 1
        line_ary[i] = "#{parent1_majority_allele}|#{parent2_majority_allele}"
      end
    end
    # Compose a new string from the modified array.
    new_line = ""
    line_ary.each {|i| new_line << "#{i}\t"}
    new_line.strip!
    new_line << "\n"
    # Write the modified string to the output vcf file.
    outvcf_body_file.write(new_line)
  else
    # Write the unmodified string to the output vcf file.
    outvcf_body_file.write(l)
  end
end

# Feedback.
puts "Modified #{n_modified_lines} lines of the vcf file."
