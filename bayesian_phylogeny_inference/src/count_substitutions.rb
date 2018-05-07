# m_matschiner Mon May 7 16:21:54 CEST 2018

# Get the command-line argument.
input_file_name = ARGV[0]

# Read the input file.
file = File.open(input_file_name)
lines = file.readlines
unless lines[0].strip.downcase == "#nexus"
    puts "ERROR: File is not in NEXUS format: #{a}!"
    exit 1
end

# Analyze the input file.
seqs = []
in_matrix = false
lines.each do |l|
    l.strip!
    if l.downcase == "matrix"
        in_matrix = true
    elsif l == ";"
        in_matrix = false
    elsif l != "" and in_matrix
        seqs << l.split[1]
    end
end
[["A","C"],["A","G"],["A","T"],["C","G"],["C","T"],["G","T"]].each do |alleles|
    allele1 = alleles[0]
    allele2 = alleles[1]
    n_substitution_sites = 0
    seqs[0].size.times do |x|
        site_has_allele1 = false
        site_has_allele2 = false
        seqs.each do |s|
            if s[x].upcase == allele1
                site_has_allele1 = true
            elsif s[x].upcase == allele2
                site_has_allele2 = true
            end
        end
        n_substitution_sites += 1 if site_has_allele1 and site_has_allele2
    end

    # Feedback.
    puts "Alleles #{allele1} and #{allele2} were found at #{n_substitution_sites} sites."
end