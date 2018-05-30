#!/usr/bin/ruby

## Michael Matschiner 2017-07-04.
## This script produces XML files for BEAST2.

require 'fileutils'

# If the arguments include "-h", display the help text.
help_text = ""
help_text << "\n"
help_text << "beauti.rb\n"
help_text << "\n"
help_text << "Produces XML input files for BEAST2 from a directory containing\n"
help_text << "alignments in nexus format."
help_text << "\n"
help_text << "Usage: beauti.rb -id analysis_id -n input_directory_name\n"
help_text << "       [-o output_directory] [-l run_length] [-ps [run_length_for_PS]]\n"
help_text << "       [-psd step_directory_for_PS] [-c constraints_file_name]\n"
help_text << "       [-t starting_tree_file_name] [-m substitution_model] [-g]\n"
help_text << "       [-bd] [-e] [-u] [-*] [-*tl tree_linkage_file_name_for_*BEAST]\n"
help_text << "       [-*sl species_linkage_file_name_for_*BEAST] [-q] [-s]\n"
help_text << "\n"
help_text << "\n"
help_text << "Available options:\n"
help_text << "\n"
help_text << "    -id      Specify an analysis id (required).\n"
help_text << "             This id will be used as a name for output files from the BEAST2 analysis\n"
help_text << "             (e.g. the .trees and .log files).\n"
help_text << "\n"
help_text << "    -n       Specify the name of an input file directory (required).\n"
help_text << "             All nexus files with ending \".nex\" will be read from this directory.\n"
help_text << "             Note that each nexus file will be treated as a separate partition\n"
help_text << "             with unlinked substitution and clock models, and if the -* is specified\n"
help_text << "             for *BEAST analyses, also with unlinked gene trees unless two or more\n"
help_text << "             partitions are linked with command -*tl (see below).\n"
help_text << "\n"
help_text << "    -o       Specify the name of an output directory (default: .).\n"
help_text << "             The output XML file will be named after the specified analysis id and.\n"
help_text << "             written to the directory specified here.\n"
help_text << "\n"
help_text << "    -l       Specify a run length for the BEAST2 analysis in number of MCMC\n"
help_text << "             generations (default: 50000000).\n"
help_text << "\n"
help_text << "    -ps      Use Path Sampling to determine the marginal likelihood (defaut: off).\n"
help_text << "\n"
help_text << "    -psd     Specify a directory in which to store the step files for Path Sampling\n"
help_text << "             (default: .). When using Path Sampling, a large number of step files will\n"
help_text << "             be written during the BEAST2 analysis. To avoid clogging you XML file\n"
help_text << "             directory, you can use this option and save these files elsewhere.\n"
help_text << "\n"
help_text << "    -c       Specify the name of a file with a list of fossil constraints\n"
help_text << "             (default: off).\n"
help_text << "             This file should include nothing except the fossil constraints in BEAST2\n"
help_text << "             XML format. The string read from this file will be inserted into the XML\n"
help_text << "             without modification.\n"
help_text << "\n"
help_text << "    -t       Specify the name of a file with a starting tree (default: off).\n"
help_text << "             If this option is not specified, a random starting tree will be used.\n"
help_text << "\n"
help_text << "    -i       Specify a number by which the age of the starting tree will be increased\n"
help_text << "             (default: off). This may help if the tree conflicts with fossil constraints.\n"
help_text << "\n"
help_text << "    -m       Specify a substitution model other than bModelTest (default: off).\n"
help_text << "             If this option is not used, the bModelTest model will be applied.\n"
help_text << "             The following alternative substitution models can be chosen: HKY, GTR, RB.\n"
help_text << "\n"
help_text << "    -g       Use a gamma distribution of rate heterogeneity, with 4 rate categories.\n"
help_text << "             (default: off).\n"
help_text << "\n"
help_text << "    -bd      Use a birth death tree prior instead of the Yule tree prior\n"
help_text << "             (default: off).\n"
help_text << "\n"
help_text << "    -e       Estimate base frequencies (default: off)\n"
help_text << "\n"
help_text << "    -u       Use the UCLN clock instead of the strict molecular clock (default: off)\n"
help_text << "\n"
help_text << "    -*       Use *BEAST (default off).\n"
help_text << "             In most cases when this option is used, options -*tl and -*sl should also\n"
help_text << "             be specified (see below).\n"
help_text << "\n"
help_text << "    -usd     Specify a mean of the exponential prior on the standard deviation of the.\n"
help_text << "             UCLN clock (default: 0.1)\n"
help_text << "\n"
help_text << "    -*tl     Specifies a file name with a command how to link gene trees for *BEAST\n"
help_text << "             analyses (default: off).\n"
help_text << "             By default, individual gene trees are allowed for all partitions. If however\n"
help_text << "             partitions are not defined per molecular marker, but per codon position,\n"
help_text << "             or separately for intron and exon sequences of the same marker, these\n"
help_text << "             partitions should be assumed to have the same gene tree (albeit with\n"
help_text << "             different substitution and clock model parameters). In this case, a command\n"
help_text << "             string in the file specified here can be used to link gene trees for\n"
help_text << "             particular partitions. The format for this command is as follows:\n"
help_text << "               mtdna:mtdna_cp12,mtdna_cp3\n"
help_text << "               nuclear:nuclear_cp12,nuclear_cp3\n"
help_text << "             (the id of the linked gene tree should be specified before a colon, and\n"
help_text << "             comma-separated alignment ids should follow after the colon).\n"
help_text << "             This command also invokes -*.\n"
help_text << "\n"
help_text << "    -*sl     Specifies a file name with a command how to group taxa into species for\n"
help_text << "             *BEAST analyses (default: off).\n"
help_text << "             By default, each taxon is considered as a separate species. In *BEAST\n"
help_text << "             analyses, however, multiple individuals may have been sequenced of the same\n"
help_text << "             species. In this case, the command string in the file specified here can be\n"
help_text << "             used to group taxa into species. The command format is similar to the one of\n"
help_text << "             option -*tl:\n"
help_text << "               human:human_ind1,human_ind2,human_ind3\n"
help_text << "               chimp:chimp_ind1,chimp_ind2,chimp_ind3\n"
help_text << "             (the id of the species should be specified before a colon, and comma-\n"
help_text << "             separated taxon ids should follow after the colon).\n"
help_text << "             This command also invokes -*.\n"
help_text << "\n"
help_text << "    -q       Write a qsub file for job submission to a Sun Grid Engine cluster\n"
help_text << "             (default: off).\n"
help_text << "             If specified, the qsub file will be named start.sh and written to the same\n"
help_text << "             directory as the BEAST2 XML file.\n"
help_text << "\n"
help_text << "    -s       Write a slurm script file for job submission to a Linux cluster\n"
help_text << "             (default: off).\n"
help_text << "             If specified, the slurm file will be named start.sh and written to the same\n"
help_text << "             directory as the BEAST2 XML file.\n"
help_text << "\n"
help_text << "\n"
if ARGV.include?("-h") or ARGV == []
    puts help_text
    exit
end

# Determine the analysis id.
if ARGV.include?("-id")
    analysis_id = ARGV[ARGV.index("-id")+1]
else
    raise "Please specify an analysis id with option -id!"
end

# Determine the input file directory.
if ARGV.include?("-n")
    nexus_dir = ARGV[ARGV.index("-n")+1]
else
    raise "Please specify a directory with alignments in nexus format with the option -n!"
end
nexus_dir_entries = Dir.entries(nexus_dir)
nexus_file_names = []
nexus_dir_entries.each {|e| nexus_file_names << e if e.include?(".nex")}

# Make sure some nexus files were found.
raise "ERROR: No nexus files were found!" unless nexus_file_names.size > 0

# Determine the output file name.
output_file_name = "#{analysis_id}.xml"
output_file_directory = ""
if ARGV.include?("-o")
    output_file_directory = ARGV[ARGV.index("-o")+1]
    output_file_directory = output_file_directory.chomp("/")
    output_file_directory << "/"
end

# Determine whether or not the marginal likelihood should be calculated by Path Sampling.
path_sampling = false
first_tab = "\t"
steps_directory = ""
path_sampling = true if ARGV.include?("-ps")
if path_sampling == true
    first_tab = "\t\t"
    if ARGV.index("-psd") == ARGV.size - 1
        steps_directory = "tmp"
    elsif ARGV[ARGV.index("-psd")+1][0] == "-"
        steps_directory = "tmp"
    else
        steps_directory = ARGV[ARGV.index("-psd")+1]
    end
end

# Determine the run length.
run_length = 50000000
run_length = ARGV[ARGV.index("-l")+1].to_i if ARGV.include?("-l")
store_length = run_length/2000

# Determine the starting tree if there is any.
starting_tree_lines = ""
if ARGV.include?("-t")
    starting_tree_file_name = ARGV[ARGV.index("-t")+1]
    starting_tree_file = File.open(starting_tree_file_name)
    starting_tree_lines = starting_tree_file.readlines
end

# Determine the increase number for the age of the starting tree if one has been specified.
age_increase = 0
if ARGV.include?("-i")
    age_increase = ARGV[ARGV.index("-i")+1].to_i
end

# Determine the constraints if there are any.
constraints = ""
if ARGV.include?("-c")
    constraints_file_name = ARGV[ARGV.index("-c")+1]
    constraints_file = File.open(constraints_file_name)
    constraints = constraints_file.read
end

# Determine the substitution model.
substitution_model = "bmodeltest"
substitution_model = ARGV[ARGV.index("-m")+1].downcase if ARGV.include?("-m")
valid_substitution_models = ["bmodeltest","hky","gtr","rb"]
raise "Please specify a substitution model with the option -m (Available options: BMODELTEST, HKY, GTR, and RB)!" unless valid_substitution_models.include?(substitution_model)

# Determine site rate heterogeneity.
plus_gamma = false
plus_gamma = true if ARGV.include?("-g")

# Determine the tree prior.
tree_model = "yule"
tree_model = "birthdeath" if ARGV.include?("-bd")

# Determine whether base frequencies should be estimated.
estimate_freqs = false
estimate_freqs = true if ARGV.include?("-e")

# Determine the clock model.
ucln_clock = false
ucln_clock = true if ARGV.include?("-u")

# Determine the mean of the exponential distribution on the standard deviation of the UCLN clock.
ucln_clock_stdev = 0.1
ucln_clock_stdev = ARGV[ARGV.index("-usd")+1].to_f if ARGV.include?("-usd")

# Determine whether *BEAST should be used, and if so, determine whether a tree linkage and/or species linkage command has been given.
starbeast = false
starbeast_tree_linkage = false
starbeast = true if ARGV.include?("-*")
if ARGV.include?("-*tl")
    starbeast = true
    starbeast_tree_linkage = true
end
if ARGV.include?("-*sl")
    starbeast = true
    starbeast_species_linkage = true
end
if starbeast_tree_linkage
    tree_linkage_command_file_name = ARGV[ARGV.index("-*tl")+1]
    tree_linkage_command_file = File.open(tree_linkage_command_file_name)
    tree_linkage_command = tree_linkage_command_file.read
    tree_linkage_command_parts = tree_linkage_command.split("\n")
    tree_linkage_command_gene_tree_ids = []
    tree_linkage_command_alignment_ids = []
    tree_linkage_command_parts.each do |p|
        if p.count(":") == 1
            tree_linkage_command_parts_ary = p.split(":")
            tree_linkage_command_gene_tree_ids << tree_linkage_command_parts_ary[0]
            tree_linkage_command_alignments_ids_string = tree_linkage_command_parts_ary[1]
            if tree_linkage_command_alignments_ids_string.include?(",")
                tree_linkage_command_alignment_ids << tree_linkage_command_alignments_ids_string.split(",")
            else
                tree_linkage_command_alignment_ids << [tree_linkage_command_alignments_ids_string]
            end
        else
            raise "Could not understand the tree linkage command!"
        end
    end
end
if starbeast_species_linkage
    species_linkage_command_file_name = ARGV[ARGV.index("-*sl")+1]
    species_linkage_command_file = File.open(species_linkage_command_file_name)
    species_linkage_command = species_linkage_command_file.read
    species_linkage_command_parts = species_linkage_command.split("\n")
    species_linkage_command_species_ids = []
    species_linkage_command_taxon_ids = []
    species_linkage_command_parts.each do |p|
        if p.count(":") == 1
            species_linkage_command_parts_ary = p.split(":")
            species_linkage_command_species_ids << species_linkage_command_parts_ary[0]
            species_linkage_command_alignments_ids_string = species_linkage_command_parts_ary[1]
            if species_linkage_command_alignments_ids_string.include?(",")
                species_linkage_command_taxon_ids << species_linkage_command_alignments_ids_string.split(",")
            else
                species_linkage_command_taxon_ids << [species_linkage_command_alignments_ids_string]
            end
        end
    end
end

# Determine whether a qsub file should be written.
write_qsub = false
write_qsub = true if ARGV.include?("-q")

# Determine whether a slurm script should be written.
write_slurm = false
write_slurm = true if ARGV.include?("-s")

# Inititate arrays for alignment ids and taxon ids.
alignment_ids = []
taxon_ids = []

# Parse the starting tree if one has been specified.
starting_tree = ""
if starting_tree_lines != ""

    # Determine whether the starting tree is given in nexus format (if not it is
    # assumed to be plain newick format).
    if starting_tree_lines[0].downcase.strip == "#nexus"
        # Get the first tree string (there could be more than one trees in this file).
        tree_line = ""
        translate = false
        starting_tree_lines.each do |l|
            if l.downcase.strip[0..3] == "tree"
                tree_line = l
                break
            elsif l.downcase.strip == "translate"
                translate = true
            end
        end
        raise "ERROR: No tree found in file #{starting_tree_file_name}!" if tree_line == ""
        starting_tree = tree_line.gsub(/\[.*?\]/,"").split("=")[1].strip.chomp(";")

        # If there was a translate block, reread it.
        if translate
            translate_lines = []
            in_translate_block = false
            starting_tree_lines.each do |l|
                if l.downcase.strip == "translate"
                    in_translate_block = true
                elsif l.downcase.strip == ";"
                    in_translate_block = false
                elsif in_translate_block
                    translate_lines << l
                end
            end
            raise "ERROR: Translate block could not be read!" if translate_lines == []
            translate_numbers = []
            translate_names = []
            translate_lines.each do |l|
                translate_line_ary = l.split
                translate_numbers << translate_line_ary[0]
                translate_names << translate_line_ary[1].chomp(",")
            end
            translate_numbers.size.times do |x|
                starting_tree.sub!("(#{translate_numbers[x]}:","(#{translate_names[x]}:")
                starting_tree.sub!(",#{translate_numbers[x]}:",",#{translate_names[x]}:")
            end
        end
    else
        starting_tree = starting_tree_lines[0].strip.chomp(";").gsub(/\[.*?\]/,"")
    end
end

# If the age of the starting tree should be increased, do so.
if age_increase > 0
    taxon_names_raw = starting_tree.scan(/\(.+?:/) + starting_tree.scan(/\,.+?:/)
    taxon_names = []
    taxon_names_raw.each do |n|
        taxon_names << n.gsub(",","").gsub(":","").gsub("(","")
    end
    taxon_names.uniq!
    taxon_names.each do |n|
        if starting_tree.match(/\(#{n}:(\d+\.\d+),/)
            starting_tree.sub!(/\(#{n}:\d+\.\d+,/,"(#{n}:#{$1.to_f+age_increase},")
        elsif starting_tree.match(/,#{n}:(\d+\.\d+)\)/)
            starting_tree.sub!(/,#{n}:\d+\.\d+\)/,",#{n}:#{$1.to_f+age_increase})")
        end
    end
end

# If a starting tree has been specified, starbeast is used, starbeast_species_linkage is false, and starting tree taxa do not contain "_spc", add it.
if starbeast and starting_tree != ""
    if starbeast_species_linkage == false
        taxon_names_raw = starting_tree.scan(/\(.+?:/) + starting_tree.scan(/\,.+?:/)
        taxon_names = []
        taxon_names_raw.each do |n|
            taxon_names << n.gsub(",","").gsub(":","").gsub("(","")
        end
        taxon_names.uniq!
        taxon_names.each do |n|
            unless n.include?("_spc")
                if starting_tree.match(/\(#{n}:/)
                    starting_tree.sub!(/\(#{n}:/,"(#{n}_spc:")
                elsif starting_tree.match(/,#{n}:/)
                    starting_tree.sub!(/,#{n}:/,",#{n}_spc:")
                end
            end
        end
    end
end

# Prepare the xml string and write the header.
xml_string = ""
xml_string << "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?>\n"
if starbeast
    if substitution_model == "bmodeltest"
        xml_string << "<beast beautitemplate=\"StarBeast\" beautistatus=\"\" namespace=\"beast.core:beast.evolution.alignment:beast.evolution.tree.coalescent:beast.core.util:beast.evolution.nuc:beast.evolution.operators:beast.evolution.sitemodel:beast.evolution.substitutionmodel:beast.evolution.likelihood\" required=\"bModelTest v1.0.3\" version=\"2.4\">\n"
    else
        xml_string << "<beast beautitemplate=\"StarBeast\" beautistatus=\"\" namespace=\"beast.core:beast.evolution.alignment:beast.evolution.tree.coalescent:beast.core.util:beast.evolution.nuc:beast.evolution.operators:beast.evolution.sitemodel:beast.evolution.substitutionmodel:beast.evolution.likelihood\" version=\"2.4\">\n"
    end
else
    if substitution_model == "bmodeltest"
        xml_string << "<beast beautitemplate=\"Standard\" beautistatus=\"\" namespace=\"beast.core:beast.evolution.alignment:beast.evolution.tree.coalescent:beast.core.util:beast.evolution.nuc:beast.evolution.operators:beast.evolution.sitemodel:beast.evolution.substitutionmodel:beast.evolution.likelihood\" required=\"bModelTest v1.0.3\" version=\"2.4\">\n"
    else
        xml_string << "<beast beautitemplate=\"Standard\" beautistatus=\"\" namespace=\"beast.core:beast.evolution.alignment:beast.evolution.tree.coalescent:beast.core.util:beast.evolution.nuc:beast.evolution.operators:beast.evolution.sitemodel:beast.evolution.substitutionmodel:beast.evolution.likelihood\" version=\"2.4\">\n"
    end
end
xml_string << "\n"

# Inform about the model.
xml_string << "<!--\n"
xml_string << "XML produced with beauti.rb:\n"
xml_string << "beauti.rb "
ARGV.each do |i|
    xml_string << "#{i} "
end
xml_string.chomp!(" ")
xml_string << "\n"
xml_string << "-->\n"
xml_string << "\n"
xml_string << "<!--\n"
xml_string << "Models implemented in this XML:\n"
if starbeast
    xml_string << "Inference tool:          *BEAST\n"
else
    xml_string << "Inference tool:          BEAST\n"
end
if tree_model == "yule"
    xml_string << "Tree model:              Yule\n"
elsif tree_model == "birthdeath"
    xml_string << "Tree model:              birth-death\n"
end
if ucln_clock
    xml_string << "Molecular clock model:   uncorrelated lognormal\n"
else
    xml_string << "Molecular clock model:   strict clock\n"
end
if substitution_model == "bmodeltest"
    xml_string << "Substitution model:      bModelTest\n"
elsif substitution_model == "hky"
    xml_string << "Substitution model:      HKY\n"
elsif substitution_model == "gtr"
    xml_string << "Substitution model:      GTR\n"
elsif substitution_model == "rb"
    xml_string << "Substitution model:      RB\n"
end
unless substitution_model == "bmodeltest"
    if plus_gamma
        xml_string << "Site rate heterogeneity: gamma distribution\n"
    else
        xml_string << "Site rate heterogeneity: none\n"
    end
    if estimate_freqs
        xml_string << "Base frequencies:        estimated\n"
    else
        xml_string << "Base frequencies:        empirical\n"
    end
end
xml_string << "-->\n"
xml_string << "\n"

# Do a first reading of all nexus_files to determine all taxon ids.
nexus_file_names.each do |f|
    nexus_file_name = f
    nexus_file = File.open("#{nexus_dir}/#{nexus_file_name}")
    nexus_lines = nexus_file.readlines
    in_matrix = false
    nexus_lines.each do |l|
        if l.strip.downcase == "matrix"
            in_matrix = true
        elsif l.strip.downcase == "end;" or l.strip == ";"
            in_matrix = false
        elsif in_matrix and l.strip != ""
            line_ary = l.strip.split(" ")
            id = line_ary[0]
            taxon_ids << id unless taxon_ids.include?(id)
        end
    end
end

# Read and write the sequences.
xml_string << "\t<!-- Sequence alignments -->\n"
nexus_file_names.each do |f|
    nexus_file_name = f
    alignment_id = nexus_file_name.chomp(".nex")
    alignment_ids << alignment_id
    nexus_file = File.open("#{nexus_dir}/#{nexus_file_name}")
    nexus_lines = nexus_file.readlines
    ids = []
    seqs = []
    in_matrix = false
    nexus_lines.each do |l|
        if l.strip.downcase == "matrix"
            in_matrix = true
        elsif l.strip.downcase == "end;" or l.strip == ";"
            in_matrix = false
        elsif in_matrix and l.strip != ""
            line_ary = l.strip.split(" ")
            id = line_ary[0]
            ids << id
            seqs << line_ary[1].downcase
        end
    end
    alignment_length = seqs[0].size
    seqs.each do |s|
        raise "Some sequences of file #{f} have different lengths!" if s.size != alignment_length
    end
    checked_ids = []
    checked_seqs = []
    taxon_ids.each do |t|
        if ids.include?(t)
            checked_ids << t
            checked_seqs << seqs[ids.index(t)]
        # else
        #     checked_ids << t
        #     empty_seq = ""
        #     alignment_length.times {|x| empty_seq << "-"}
        #     checked_seqs << empty_seq
        end
    end
    alignment_length = seqs[0].size
    seqs.each do |s|
        raise "Some sequences have different lengths in file #{f}!" unless s.size == alignment_length
    end
    xml_string << "\t<data id=\"#{alignment_id}\" name=\"alignment\">\n"
    checked_ids.size.times do |x|
        xml_string << "\t\t<sequence id=\"#{alignment_id}_#{checked_ids[x]}\" taxon=\"#{checked_ids[x]}\" totalcount=\"4\" value=\"#{checked_seqs[x]}\"/>\n"
    end
    xml_string << "\t</data>\n"
end
xml_string << "\n"

# Now that all alignment ids are read, see for which alignment ids the the trees should be linked for *BEAST.
if starbeast
    gene_tree_ids = []
    gene_tree_alignment_ids = []
    if starbeast_tree_linkage
        alignments_to_be_distributed = []
        alignment_ids.each do |a|
            alignments_to_be_distributed << a
        end
        while alignments_to_be_distributed.size > 0
            alignment_id = alignments_to_be_distributed.shift
            unlinked = true
            tree_linkage_command_alignment_ids.each do |ca|
                unlinked = false if ca.include?(alignment_id)
            end
            if unlinked
                gene_tree_ids << alignment_id
                gene_tree_alignment_ids << [alignment_id]
            else
                gene_tree_id_for_this_alignment = nil
                tree_linkage_command_gene_tree_ids.size.times do |x|
                    if tree_linkage_command_alignment_ids[x].include?(alignment_id)
                        gene_tree_id_for_this_alignment = tree_linkage_command_gene_tree_ids[x]
                        break
                    end
                end
                raise "Gene tree id was not found!" if gene_tree_id_for_this_alignment == nil
                if gene_tree_ids.include?(gene_tree_id_for_this_alignment)
                    gene_tree_alignment_ids[gene_tree_ids.index(gene_tree_id_for_this_alignment)] << alignment_id
                else
                    gene_tree_ids << gene_tree_id_for_this_alignment
                    gene_tree_alignment_ids << [alignment_id]
                end
            end
        end

    else
        alignment_ids.each do |a|
            gene_tree_ids << a
            gene_tree_alignment_ids << [a]
        end
    end
end

# Write the maps part.
xml_string << "\t<map name=\"Beta\">beast.math.distributions.Beta</map>\n"
xml_string << "\t<map name=\"Exponential\">beast.math.distributions.Exponential</map>\n"
xml_string << "\t<map name=\"InverseGamma\">beast.math.distributions.InverseGamma</map>\n"
xml_string << "\t<map name=\"LogNormal\">beast.math.distributions.LogNormalDistributionModel</map>\n"
xml_string << "\t<map name=\"Gamma\">beast.math.distributions.Gamma</map>\n"
xml_string << "\t<map name=\"Uniform\">beast.math.distributions.Uniform</map>\n"
xml_string << "\t<map name=\"prior\">beast.math.distributions.Prior</map>\n"
xml_string << "\t<map name=\"LaplaceDistribution\">beast.math.distributions.LaplaceDistribution</map>\n"
xml_string << "\t<map name=\"OneOnX\">beast.math.distributions.OneOnX</map>\n"
xml_string << "\t<map name=\"Normal\">beast.math.distributions.Normal</map>\n"
xml_string << "\n"

if starbeast and starting_tree != ""
    xml_string << "#{first_tab}<init estimate=\"false\" birthRate=\"@birthRate.t:Species\" popMean=\"@popMean\" spec=\"beast.evolution.speciation.StarBeastStartState\" speciesTree=\"@tree.t:Species\">\n"
    gene_tree_ids.size.times do |x|
        xml_string << "#{first_tab}\t<tree idref=\"tree.t:#{gene_tree_ids[x]}\" name=\"gene\"/>\n"
    end
    xml_string << "#{first_tab}\t<speciesTreePrior bottomPopSize=\"@popSize\" gammaParameter=\"@popMean\" spec=\"beast.evolution.speciation.SpeciesTreePrior\" taxonset=\"@allTaxa\" tree=\"@tree.t:Species\"/>\n"
    xml_string << "#{first_tab}</init>\n"
    xml_string << "\n"
end

# Open the run part if using path sampling.
if path_sampling
    nrofsteps = 64
    xml_string << "\t<run spec=\"beast.inference.PathSampler\"\n"
    xml_string << "\t\tchainLength=\"400000\"\n"
    xml_string << "\t\talpha=\"0.3\"\n"
    xml_string << "\t\trootdir=\"#{steps_directory}\"\n"
    xml_string << "\t\tburnInPercentage=\"50\"\n"
    xml_string << "\t\tpreBurnin=\"50000\"\n"
    xml_string << "\t\tdeleteOldLogs=\"true\"\n"
    xml_string << "\t\tnrofsteps=\"#{nrofsteps}\">\n"
    nrofsteps.times do |x|
        xml_string << "\t\tcp -r BEASTii tmp/step#{x}\n"
    end
    if substitution_model == "bmodeltest"
        nrofsteps.times do |x|
            xml_string << "\t\tcp -r bModelTest tmp/step#{x}\n"
            xml_string << "\t\tcp -r BEASTLabs tmp/step#{x}\n"
        end        
    elsif substitution_model == "rb"
        nrofsteps.times do |x|
            xml_string << "\t\tcp -r RBS tmp/step#{x}\n"
        end
    end
    xml_string << "\t\tcd $(dir)\n"
    # In the next line it makes sense to use -threads 1 even if mulitple cores are available, as multithreading
    # is more efficient when the number of cpus is specified when starting the XML produced here.
    xml_string << "\t\tjava -cp $(java.class.path) -Xmx4098m -jar ../../beast.jar -threads 1 $(resume/overwrite) -java -seed $(seed) beast.xml\n"
    xml_string << "\n"
    xml_string << "#{first_tab}<mcmc chainLength=\"#{run_length}\" id=\"mcmc\" spec=\"MCMC\" storeEvery=\"#{store_length}\">\n"
else
    xml_string << "\t<run chainLength=\"#{run_length}\" id=\"mcmc\" spec=\"MCMC\" storeEvery=\"#{store_length}\">\n"
end
xml_string << "\n"

# Open the state part.
xml_string << "#{first_tab}\t<!-- Set up all parameters as part of the state -->\n"
xml_string << "#{first_tab}\t<state id=\"state\" storeEvery=\"#{store_length}\">\n"
xml_string << "\n"

# A parameter for the species tree, and a taxon set in which all taxon ids are defined.
xml_string << "#{first_tab}\t\t<!-- The species tree as a parameter, and defining a taxon set of species -->\n"
xml_string << "#{first_tab}\t\t<tree id=\"tree.t:Species\" name=\"stateNode\">\n"
xml_string << "#{first_tab}\t\t\t<taxonset id=\"allTaxa\" spec=\"TaxonSet\">\n"
if starbeast
    if starbeast_species_linkage
        # Make sure that all taxon_ids are found in all species_linkage_command_taxon_ids.
        unless taxon_ids.sort == species_linkage_command_taxon_ids.flatten.sort
            raise "Taxon ids in alignments and taxon ids in the spceis linkage command could not be mapped 1:1."
        end
        species_linkage_command_species_ids.size.times do |x|
            xml_string << "#{first_tab}\t\t\t\t<taxon id=\"#{species_linkage_command_species_ids[x]}\" spec=\"TaxonSet\">\n"
            species_linkage_command_taxon_ids[x].each do |t|
                xml_string << "#{first_tab}\t\t\t\t\t<taxon id=\"#{t}\" spec=\"Taxon\"/>\n"
            end
            xml_string << "#{first_tab}\t\t\t\t</taxon>\n"
        end
    else
        taxon_ids.each do |t|
            xml_string << "#{first_tab}\t\t\t\t<taxon id=\"#{t}_spc\" spec=\"TaxonSet\">\n"
            xml_string << "#{first_tab}\t\t\t\t\t<taxon id=\"#{t}\" spec=\"Taxon\"/>\n"
            xml_string << "#{first_tab}\t\t\t\t</taxon>\n"
        end
    end
else
    taxon_ids.each do |t|
        xml_string << "#{first_tab}\t\t\t\t<taxon id=\"#{t}\" spec=\"Taxon\"/>\n"
    end
end
xml_string << "#{first_tab}\t\t\t</taxonset>\n"
xml_string << "#{first_tab}\t\t</tree>\n"
xml_string << "\n"

# Parameters of individual gene trees for *BEAST.
if starbeast
    if gene_tree_ids.size == 1
        xml_string << "#{first_tab}\t\t<!-- The gene tree as a parameter, with taxon sets taken from its alignment -->\n"
    else
        xml_string << "#{first_tab}\t\t<!-- The #{gene_tree_ids.size} gene trees as parameters, with taxon sets taken from their alignments -->\n"
    end
    gene_tree_ids.size.times do |x|
        xml_string << "#{first_tab}\t\t<tree id=\"tree.t:#{gene_tree_ids[x]}\" name=\"stateNode\">\n"
        xml_string << "#{first_tab}\t\t\t<taxonset id=\"taxonSet.#{gene_tree_ids[x]}\" spec=\"TaxonSet\">\n"
        xml_string << "#{first_tab}\t\t\t\t<data idref=\"#{gene_tree_alignment_ids[x][0]}\" name=\"alignment\"/>\n"
        xml_string << "#{first_tab}\t\t\t</taxonset>\n"
        xml_string << "#{first_tab}\t\t</tree>\n"
    end
    xml_string << "\n"
end

# Parameters of the tree.
if tree_model == "yule"
    xml_string << "#{first_tab}\t\t<!-- A parameters to be used as birth rate. -->\n"
    xml_string << "#{first_tab}\t\t<parameter id=\"birthRate.t:Species\" value=\"0.01\" name=\"stateNode\"/>\n"
elsif tree_model == "birthdeath"
    xml_string << "#{first_tab}\t\t<!-- Two parameters to be used as birth and relative death rates. -->\n"
    xml_string << "#{first_tab}\t\t<parameter id=\"birthRate.t:Species\" value=\"0.01\" name=\"stateNode\"/>\n"
    xml_string << "#{first_tab}\t\t<parameter id=\"relativeDeathRate.t:Species\" value=\"0.5\" name=\"stateNode\"/>\n"
end
xml_string << "\n"

# Population parameters for *BEAST.
if starbeast
    xml_string << "#{first_tab}\t\t<!-- A population size parameter -->\n"
    xml_string << "#{first_tab}\t\t<parameter id=\"popSize\" value=\"1\" name=\"stateNode\"/>\n"
    xml_string << "\n"
    xml_string << "#{first_tab}\t\t<!-- A population mean parameter (this is the mean of the gamma distribution for population sizes) -->\n"
    xml_string << "#{first_tab}\t\t<parameter id=\"popMean\" value=\"1\" name=\"stateNode\"/>\n"
    xml_string << "\n"
end

# Parameters of the molecular clock.
if ucln_clock
    xml_string << "#{first_tab}\t\t<!-- Parameters of the UCLN clock for each gene -->\n"
    if starbeast
        gene_tree_ids.each do |g|
            xml_string << "#{first_tab}\t\t<parameter id=\"ucldStdev.c:#{g}\" value=\"0.005\" name=\"stateNode\"/>\n"
            xml_string << "#{first_tab}\t\t<parameter dimension=\"10\" id=\"rateCategories.c:#{g}\" spec=\"parameter.IntegerParameter\" value=\"1\" name=\"stateNode\"/>\n"
        end
    else
        xml_string << "#{first_tab}\t\t<parameter id=\"ucldStdev.c:Species\" value=\"0.005\" name=\"stateNode\"/>\n"
        xml_string << "#{first_tab}\t\t<parameter dimension=\"10\" id=\"rateCategories.c:Species\" spec=\"parameter.IntegerParameter\" value=\"1\" name=\"stateNode\"/>\n"
    end
end
xml_string << "\n"

# Parameters of the substitution rates.
xml_string << "#{first_tab}\t\t<!-- Parameters of the substitution rates -->\n"
alignment_ids.each do |a|
    xml_string << "#{first_tab}\t\t<parameter id=\"mutationRate.s:#{a}\" value=\"1\" name=\"stateNode\"/>\n"
end
xml_string << "\n"

# Parameters of the bModelTest substitution model.
if substitution_model == "bmodeltest"
    xml_string << "#{first_tab}\t\t<!-- Parameters of the bModelTest substitution model for each alignment -->\n"
    alignment_ids.each do |a|
        xml_string << "#{first_tab}\t\t<stateNode id=\"BMT_ModelIndicator.s:#{a}\" spec=\"parameter.IntegerParameter\" lower=\"0\" upper=\"5\">5</stateNode>\n"
        xml_string << "#{first_tab}\t\t<parameter id=\"BMT_Rates.s:#{a}\" dimension=\"6\" lower=\"0.01\" name=\"stateNode\" upper=\"100.0\">1.0</parameter>\n"
        xml_string << "#{first_tab}\t\t<parameter id=\"BMT_gammaShape.s:#{a}\" name=\"stateNode\">1.0</parameter>\n"
        xml_string << "#{first_tab}\t\t<parameter id=\"BMT_ProportionInvariable.s:#{a}\" lower=\"0.0\" name=\"stateNode\" upper=\"1.0\">0.1</parameter>\n"
        xml_string << "#{first_tab}\t\t<stateNode id=\"hasInvariableSites.s:#{a}\" spec=\"parameter.IntegerParameter\">1</stateNode>\n"
        xml_string << "#{first_tab}\t\t<stateNode id=\"hasGammaRates.s:#{a}\" spec=\"parameter.IntegerParameter\">1</stateNode>\n"
        xml_string << "#{first_tab}\t\t<stateNode id=\"hasEqualFreqs.s:#{a}\" spec=\"parameter.BooleanParameter\">false</stateNode>\n"
        xml_string << "#{first_tab}\t\t<parameter id=\"BMT_frequencies.s:#{a}\" dimension=\"4\" name=\"stateNode\">0.25 0.25 0.25 0.25</parameter>\n"
    end
else
    # Parameters of other substitution models.
    if substitution_model == "hky"
        xml_string << "#{first_tab}\t\t<!-- Parameters of the HKY substitution model for each alignment -->\n"
        alignment_ids.each do |a|
            xml_string << "#{first_tab}\t\t<parameter id=\"kappa.s:#{a}\" value=\"2.0\" name=\"stateNode\"/>\n"
        end
    elsif substitution_model == "gtr"
        xml_string << "#{first_tab}\t\t<!-- Parameters of the GTR substitution model for each alignment -->\n"
        alignment_ids.each do |a|
            xml_string << "#{first_tab}\t\t<parameter id=\"rateAC.s:#{a}\" value=\"1.0\" name=\"stateNode\"/>\n"
            xml_string << "#{first_tab}\t\t<parameter id=\"rateAG.s:#{a}\" value=\"1.0\" name=\"stateNode\"/>\n"
            xml_string << "#{first_tab}\t\t<parameter id=\"rateAT.s:#{a}\" value=\"1.0\" name=\"stateNode\"/>\n"
            xml_string << "#{first_tab}\t\t<parameter id=\"rateCG.s:#{a}\" value=\"1.0\" name=\"stateNode\"/>\n"
            xml_string << "#{first_tab}\t\t<parameter id=\"rateGT.s:#{a}\" value=\"1.0\" name=\"stateNode\"/>\n"
        end
    elsif substitution_model == "rb"
        xml_string << "#{first_tab}\t\t<!-- Parameters of the RB substitution model for each alignment -->\n"
        alignment_ids.each do |a|
            xml_string << "#{first_tab}\t\t<parameter id=\"RBcount.s:#{a}\" spec=\"parameter.IntegerParameter\" lower=\"0\" upper=\"5\" value=\"5\" name=\"stateNode\"/>\n"
            xml_string << "#{first_tab}\t\t<parameter dimension=\"5\" id=\"RBrates.s:#{a}\" lower=\"0.01\" upper=\"100.0\" value=\"1\" name=\"stateNode\"/>\n"
        end
    end
    xml_string << "\n"

    # Parameters of the gamma distribution of site rate heterogeneity
    if plus_gamma
        xml_string << "#{first_tab}\t\t<!-- Parameters of the gamma distribution of site rate heterogeneity for each alignment -->\n"
        alignment_ids.each do |a|
            xml_string << "#{first_tab}\t\t<parameter id=\"gammaShape.s:#{a}\" value=\"1.0\" name=\"stateNode\"/>\n"
        end
        xml_string << "\n"
    end

    # Parameters of the base frequencies.
    if estimate_freqs
        xml_string << "#{first_tab}\t\t<!-- Parameters of the base frequencies for each alignment -->\n"
        alignment_ids.each do |a|
            xml_string << "#{first_tab}\t\t<parameter id=\"freqParameter.s:#{a}\" dimension=\"4\" lower=\"0.0\" upper=\"1.0\" value=\"0.25\" name=\"stateNode\"/>\n"
        end
        xml_string << "\n"
    end
end

# Close the state part.
xml_string << "#{first_tab}\t</state>\n"
xml_string << "\n"

# Open the initiation part.
xml_string << "#{first_tab}\t<!-- Intitators -->\n"

# Inititate the starting trees.
if starbeast
    if path_sampling
        # XXX TODO.
        raise "Use of path sampling is not completely implemented. Check how to initiate trees with path sampling."
    else
        if starting_tree == ""
            # *BEAST and random tree.
            xml_string << "#{first_tab}\t<init estimate=\"false\" birthRate=\"@birthRate.t:Species\" popMean=\"@popMean\" spec=\"beast.evolution.speciation.StarBeastStartState\" speciesTree=\"@tree.t:Species\">\n"
            gene_tree_ids.size.times do |x|
                xml_string << "#{first_tab}\t\t<tree idref=\"tree.t:#{gene_tree_ids[x]}\" name=\"gene\"/>\n"
            end
            xml_string << "#{first_tab}\t\t<speciesTreePrior bottomPopSize=\"@popSize\" gammaParameter=\"@popMean\" spec=\"beast.evolution.speciation.SpeciesTreePrior\" taxonset=\"@allTaxa\" tree=\"@tree.t:Species\"/>\n"
            xml_string << "#{first_tab}\t</init>\n"
        else
            # *BEAST and starting_tree.
            xml_string << "#{first_tab}\t<init spec=\"beast.util.TreeParser\" id=\"NewickTree.t:Species\"\n"
            xml_string << "#{first_tab}\t\tinitial=\"@tree.t:Species\" taxonset=\"@allTaxa\"\n"
            xml_string << "#{first_tab}\t\tIsLabelledNewick=\"true\"\n"
            xml_string << "#{first_tab}\t\tnewick=\"#{starting_tree}\"/>\n"
            xml_string << "\n"
            # The first gene tree intitiator.
            xml_string << "#{first_tab}\t<init spec=\"beast.evolution.speciation.RandomGeneTree\"\n"
            xml_string << "#{first_tab}\t\tid=\"randomGeneTree1\" initial=\"@tree.t:#{gene_tree_ids[0]}\"\n"
            xml_string << "#{first_tab}\t\tspeciesTree=\"@tree.t:Species\" taxa=\"@#{alignment_ids[0]}\">\n"
            xml_string << "#{first_tab}\t\t<populationModel id=\"popmodel\" spec=\"ConstantPopulation\" popSize=\"1.0\"/>\n"
            xml_string << "#{first_tab}\t</init>\n"
            # All other gene tree initiators.
            1.upto(gene_tree_ids.size-1) do |x|
                xml_string << "#{first_tab}\t<init spec=\"beast.evolution.speciation.RandomGeneTree\"\n" 
                xml_string << "#{first_tab}\t\tid=\"randomGeneTree#{x+1}\" initial=\"@tree.t:#{gene_tree_ids[x]}\"\n"
                xml_string << "#{first_tab}\t\tspeciesTree=\"@tree.t:Species\" taxa=\"@#{alignment_ids[x]}\"\n"
                xml_string << "#{first_tab}\t\tpopulationModel=\"@popmodel\"/>\n"
            end
        end
    end
else
    if path_sampling
        # XXX TODO.
        raise "Use of path sampling is not completely implemented. Check how to initiate trees with path sampling."
    else
        if starting_tree == ""
            # standard beast and random tree.
            xml_string << "#{first_tab}\t<init estimate=\"false\" id=\"randomTree.t:Species\" spec=\"beast.evolution.tree.RandomTree\" taxa=\"@#{alignment_ids[0]}\" initial=\"@tree.t:Species\">\n"
            xml_string << "#{first_tab}\t\t<populationModel id=\"constantPopulation.t:Species\" spec=\"ConstantPopulation\">\n"
            xml_string << "#{first_tab}\t\t\t<parameter id=\"randomPopSize.t:Species\" name=\"popSize\">1.0</parameter>\n"
            xml_string << "#{first_tab}\t\t</populationModel>\n"
            xml_string << "#{first_tab}\t</init>\n"
        else
            # standard beast and starting tree.
            xml_string << "#{first_tab}\t<init spec=\"beast.util.TreeParser\" id=\"NewickTree.t:#{alignment_ids[0]}\"\n" 
            xml_string << "#{first_tab}\t\tinitial=\"@tree.t:Species\" taxa=\"@#{alignment_ids[0]}\"\n"
            xml_string << "#{first_tab}\t\tIsLabelledNewick=\"true\"\n"
            xml_string << "#{first_tab}\tnewick=\"#{starting_tree}\"/>\n"
            xml_string << "\n"
        end
    end
end
xml_string << "\n"

# Open the posterior part.
xml_string << "#{first_tab}\t<!-- The posterior -->\n"
xml_string << "#{first_tab}\t<distribution id=\"posterior\" spec=\"util.CompoundDistribution\">\n"
xml_string << "\n"

# The species coalescent for *BEAST.
if starbeast
    xml_string << "#{first_tab}\t\t<!-- Part 1 of the posterior: the species coalescent prior -->\n"
    if path_sampling
        xml_string << "#{first_tab}\t\t<distribution id=\"prior\" spec=\"util.CompoundDistribution\">\n"
    else
        xml_string << "#{first_tab}\t\t<distribution id=\"speciescoalescent\" spec=\"util.CompoundDistribution\">\n"
    end
    xml_string << "\n"
    xml_string << "#{first_tab}\t\t\t<!-- Coalescent population size parameter -->\n"
    xml_string << "#{first_tab}\t\t\t<distribution bottomPopSize=\"@popSize\" gammaParameter=\"@popMean\" id=\"speciesTreePopSize.Species\" spec=\"beast.evolution.speciation.SpeciesTreePrior\" taxonset=\"@allTaxa\" tree=\"@tree.t:Species\">\n"
    xml_string << "#{first_tab}\t\t\t\t<parameter id=\"popSizeTop\" name=\"topPopSize\" value=\"1\"/>\n"
    xml_string << "#{first_tab}\t\t\t</distribution>\n"
    gene_tree_ids.each do |i|
        if i.downcase.include?("mtdna") or i.downcase.include?("mitoc")
            xml_string << "#{first_tab}\t\t\t<distribution id=\"treePrior.t:#{i}\" ploidy=\"0.5\" spec=\"beast.evolution.speciation.GeneTreeForSpeciesTreeDistribution\" speciesTree=\"@tree.t:Species\" speciesTreePrior=\"@speciesTreePopSize.Species\" tree=\"@tree.t:#{i}\"/>\n"
        else
            xml_string << "#{first_tab}\t\t\t<distribution id=\"treePrior.t:#{i}\" spec=\"beast.evolution.speciation.GeneTreeForSpeciesTreeDistribution\" speciesTree=\"@tree.t:Species\" speciesTreePrior=\"@speciesTreePopSize.Species\" tree=\"@tree.t:#{i}\"/>\n"
        end
    end
    xml_string << "\n"
    xml_string << "#{first_tab}\t\t<!-- End of: Part 1 of the posterior: the species coalescent prior -->\n"
    unless path_sampling
        xml_string << "#{first_tab}\t\t</distribution>\n"
    end
    xml_string << "\n"
end

# Open the prior part.
if starbeast
    xml_string << "#{first_tab}\t\t<!-- Part 2 of the posterior: the prior -->\n"
else
    xml_string << "#{first_tab}\t\t<!-- Part 1 of the posterior: the prior -->\n"
end
unless path_sampling
    xml_string << "#{first_tab}\t\t<distribution id=\"prior\" spec=\"util.CompoundDistribution\">\n"
end
xml_string << "\n"

# Divergence age priors.
xml_string << "#{first_tab}\t\t\t<!-- Divergence age priors -->\n"
if constraints == ""
    xml_string << "#{first_tab}\t\t\t<!-- Add divergence age priors here -->\n"
else
    xml_string << constraints
end
xml_string << "\n"

# The tree prior.
if tree_model == "yule"
    xml_string << "#{first_tab}\t\t\t<!-- The Yule tree prior -->\n"
    xml_string << "#{first_tab}\t\t\t<distribution id=\"yuleModel.t:Species\" birthDiffRate=\"@birthRate.t:Species\" spec=\"beast.evolution.speciation.YuleModel\" tree=\"@tree.t:Species\"/>\n"
    xml_string << "\n"
    xml_string << "#{first_tab}\t\t\t<!-- Prior on the birth rate -->\n"
    xml_string << "#{first_tab}\t\t\t<prior id=\"birthRatePrior\" name=\"distribution\" x=\"@birthRate.t:Species\">\n"
    xml_string << "#{first_tab}\t\t\t\t<Exponential name=\"distr\" mean=\"0.1\"/>\n"
    xml_string << "#{first_tab}\t\t\t</prior>\n"
elsif tree_model == "birthdeath"
    xml_string << "#{first_tab}\t\t\t<!-- The Birth Death tree prior -->\n"
    xml_string << "#{first_tab}\t\t\t<distribution id=\"birthDeath.t:Species\" birthDiffRate=\"@birthRate.t:Species\" relativeDeathRate=\"@relativeDeathRate.t:Species\" spec=\"beast.evolution.speciation.BirthDeathGernhard08Model\" tree=\"@tree.t:Species\"/>\n"
    xml_string << "\n"
    xml_string << "#{first_tab}\t\t\t<!-- Priors on the birth and death rates -->\n"
    xml_string << "#{first_tab}\t\t\t<prior id=\"birthRatePrior\" name=\"distribution\" x=\"@birthRate.t:Species\">\n"
    xml_string << "#{first_tab}\t\t\t\t<Exponential name=\"distr\" mean=\"0.1\"/>\n"
    xml_string << "#{first_tab}\t\t\t</prior>\n"
    xml_string << "#{first_tab}\t\t\t<prior id=\"relativeDeathRatePrior\" name=\"distribution\" x=\"@relativeDeathRate.t:Species\">\n"
    xml_string << "#{first_tab}\t\t\t\t<Uniform name=\"distr\" lower=\"0.0001\" upper=\"0.99\"/>\n"
    xml_string << "#{first_tab}\t\t\t</prior>\n"
end
xml_string << "\n"

# A prior on population sizes for *BEAST.
if starbeast
    xml_string << "#{first_tab}\t\t\t<!-- A one-over-x hyperprior on population sizes-->\n"
    xml_string << "#{first_tab}\t\t\t<prior id=\"popMean.prior\" name=\"distribution\" x=\"@popMean\">\n"
    xml_string << "#{first_tab}\t\t\t\t<OneOnX name=\"distr\"/>\n"
    xml_string << "#{first_tab}\t\t\t</prior>\n"
    xml_string << "\n"
end

# Priors on the molecular clock.
if ucln_clock
    xml_string << "#{first_tab}\t\t\t<!-- An exponential prior on the standard deviation of the rate of the each uncorrelated lognormal clock -->\n"
    if starbeast
        gene_tree_ids.each do |g|
            xml_string << "#{first_tab}\t\t\t<prior id=\"ucldStdevPrior.c:#{g}\" name=\"distribution\" x=\"@ucldStdev.c:#{g}\">\n"
            xml_string << "#{first_tab}\t\t\t\t<Exponential name=\"distr\" mean=\"#{ucln_clock_stdev}\"/>\n"
            xml_string << "#{first_tab}\t\t\t</prior>\n"
        end
    else
        xml_string << "#{first_tab}\t\t\t<prior id=\"ucldStdevPrior.c:Species\" name=\"distribution\" x=\"@ucldStdev.c:Species\">\n"
        xml_string << "#{first_tab}\t\t\t\t<Exponential name=\"distr\" mean=\"#{ucln_clock_stdev}\"/>\n"
        xml_string << "#{first_tab}\t\t\t</prior>\n"
    end
end
xml_string << "\n"

# Priors on the substitution rate.
xml_string << "#{first_tab}\t\t\t<!-- One-over-x priors on the substitution rates for each alignment -->\n"
alignment_ids.each do |a|
    xml_string << "#{first_tab}\t\t\t<prior id=\"MutationRatePrior.s:#{a}\" name=\"distribution\" x=\"@mutationRate.s:#{a}\">\n"
    xml_string << "#{first_tab}\t\t\t\t<OneOnX name=\"distr\"/>\n"
    xml_string << "#{first_tab}\t\t\t</prior>\n"
end
xml_string << "\n"

# Priors of the bmodeltest substitution model.
if substitution_model == "bmodeltest"
    xml_string << "#{first_tab}\t\t\t<!-- Priors for the bModelTest substitution models for each alignment -->\n"
    alignment_ids.each do |a|

        xml_string << "#{first_tab}\t\t\t<!-- A distribution for the model itself, for alignment #{a} -->\n"
        xml_string << "#{first_tab}\t\t\t<distribution id=\"BMT_RatesPrior.s:#{a}\" spec=\"beast.math.distributions.NucleotideRevJumpSubstModelRatePrior\" modelIndicator=\"@BMT_ModelIndicator.s:#{a}\" x=\"@BMT_Rates.s:#{a}\">\n"
        xml_string << "#{first_tab}\t\t\t\t<substModel id=\"RevJump.s:#{a}\" spec=\"NucleotideRevJumpSubstModel\" modelIndicator=\"@BMT_ModelIndicator.s:#{a}\" rates=\"@BMT_Rates.s:#{a}\">\n"
        xml_string << "#{first_tab}\t\t\t\t\t<frequencies id=\"BMTfreqs.s:#{a}\" spec=\"ModelFrequencies\" empirical=\"false\" frequencies=\"@BMT_frequencies.s:#{a}\" hasEqualFreqs=\"@hasEqualFreqs.s:#{a}\"/>\n"
        xml_string << "#{first_tab}\t\t\t\t</substModel>\n"
        xml_string << "#{first_tab}\t\t\t\t<Exponential id=\"BMT_RatesPrior.s:#{a}x\" name=\"distr\"/>\n"
        xml_string << "#{first_tab}\t\t\t</distribution>\n"

        xml_string << "#{first_tab}\t\t\t<!-- An exponential prior for the alpha parameter of the distribution of site rate heterogeneity, for alignment #{a} -->\n"
        xml_string << "#{first_tab}\t\t\t<distribution id=\"BMT_GammaShapePrior.s:#{a}\" spec=\"beast.math.distributions.BMTPrior\" count=\"@hasGammaRates.s:#{a}\" x=\"@BMT_gammaShape.s:#{a}\">\n"
        xml_string << "#{first_tab}\t\t\t\t<Exponential id=\"BMT_GammaShapePrior.s:#{a}x\" name=\"distr\">\n"
        xml_string << "#{first_tab}\t\t\t\t\t<parameter name=\"mean\">0.1</parameter>\n"
        xml_string << "#{first_tab}\t\t\t\t</Exponential>\n"
        xml_string << "#{first_tab}\t\t\t</distribution>\n"

        xml_string << "#{first_tab}\t\t\t<!-- A beta prior for invariable sites, for alignment #{a} -->\n"
        xml_string << "#{first_tab}\t\t\t<distribution id=\"BMT_PropInvariablePrior.s:#{a}\" spec=\"beast.math.distributions.BMTPrior\" count=\"@hasInvariableSites.s:#{a}\" x=\"@BMT_ProportionInvariable.s:#{a}\">\n"
        xml_string << "#{first_tab}\t\t\t\t<Beta id=\"BMT_PropInvariablePrior.s:#{a}b\" name=\"distr\">\n"
        xml_string << "#{first_tab}\t\t\t\t\t<parameter name=\"alpha\">1.0</parameter>\n"
        xml_string << "#{first_tab}\t\t\t\t\t<parameter name=\"beta\">4.0</parameter>\n"
        xml_string << "#{first_tab}\t\t\t\t</Beta>\n"
        xml_string << "#{first_tab}\t\t\t</distribution>\n"

        xml_string << "#{first_tab}\t\t\t<!-- A Dirichlet prior for frequency distributions, for alignment #{a} -->\n"
        xml_string << "#{first_tab}\t\t\t<prior id=\"BMT_freqsPrior.s:#{a}\" name=\"distribution\" x=\"@BMT_frequencies.s:#{a}\">\n"
        xml_string << "#{first_tab}\t\t\t\t<distr id=\"BMT_freqsPrior.s:#{a}d\" spec=\"beast.math.distributions.Dirichlet\">\n"
        xml_string << "#{first_tab}\t\t\t\t\t<parameter dimension=\"4\" lower=\"0.0\" name=\"alpha\" upper=\"0.0\">1.0 1.0 1.0 1.0</parameter>\n"
        xml_string << "#{first_tab}\t\t\t\t</distr>\n"
        xml_string << "#{first_tab}\t\t\t</prior>\n"
        xml_string << "\n"
    end
# Priors of other substitution models.
else
    if substitution_model == "hky"
        xml_string << "#{first_tab}\t\t\t<!-- Priors for the HKY substitution models for each alignment -->\n"
        alignment_ids.each do |a|
            xml_string << "#{first_tab}\t\t\t<prior id=\"KappaPrior.s:#{a}\" name=\"distribution\" x=\"@kappa.s:#{a}\">\n"
            xml_string << "#{first_tab}\t\t\t\t<LogNormal name=\"distr\" M=\"1.0\" S=\"1.25\"/>\n"
            xml_string << "#{first_tab}\t\t\t</prior>\n"
        end
    elsif substitution_model == "gtr"
        xml_string << "#{first_tab}\t\t\t<!-- Priors for the GTR substitution models for each alignment -->\n"
        alignment_ids.each do |a|
            xml_string << "#{first_tab}\t\t\t<prior id=\"rateACPrior.s:#{a}\" name=\"distribution\" x=\"@rateAC.s:#{a}\">\n"
            xml_string << "#{first_tab}\t\t\t\t<Gamma name=\"distr\" alpha=\"0.05\" beta=\"10.0\"/>\n"
            xml_string << "#{first_tab}\t\t\t</prior>\n"
            xml_string << "#{first_tab}\t\t\t<prior id=\"rateAGPrior.s:#{a}\" name=\"distribution\" x=\"@rateAG.s:#{a}\">\n"
            xml_string << "#{first_tab}\t\t\t\t<Gamma name=\"distr\" alpha=\"0.05\" beta=\"20.0\"/>\n"
            xml_string << "#{first_tab}\t\t\t</prior>\n"
            xml_string << "#{first_tab}\t\t\t<prior id=\"rateATPrior.s:#{a}\" name=\"distribution\" x=\"@rateAT.s:#{a}\">\n"
            xml_string << "#{first_tab}\t\t\t\t<Gamma name=\"distr\" alpha=\"0.05\" beta=\"10.0\"/>\n"
            xml_string << "#{first_tab}\t\t\t</prior>\n"
            xml_string << "#{first_tab}\t\t\t<prior id=\"rateCGPrior.s:#{a}\" name=\"distribution\" x=\"@rateCG.s:#{a}\">\n"
            xml_string << "#{first_tab}\t\t\t\t<Gamma name=\"distr\" alpha=\"0.05\" beta=\"10.0\"/>\n"
            xml_string << "#{first_tab}\t\t\t</prior>\n"
            xml_string << "#{first_tab}\t\t\t<prior id=\"rateGTPrior.s:#{a}\" name=\"distribution\" x=\"@rateGT.s:#{a}\">\n"
            xml_string << "#{first_tab}\t\t\t\t<Gamma name=\"distr\" alpha=\"0.05\" beta=\"10.0\"/>\n"
            xml_string << "#{first_tab}\t\t\t</prior>\n"
        end
    elsif substitution_model == "rb"
        xml_string << "#{first_tab}\t\t\t<!-- Priors for the RB substitution models for each alignment -->\n"
        alignment_ids.each do |a|
            xml_string << "#{first_tab}\t\t\t<distribution id=\"RBprior.s:#{a}\" count=\"@RBcount.s:#{a}\" spec=\"beast.math.distributions.RBPrior\" x=\"@RBrates.s:#{a}\">\n"
            xml_string << "#{first_tab}\t\t\t\t<Gamma name=\"distr\" alpha=\"0.2\" beta=\"5.0\"/>\n"
            xml_string << "#{first_tab}\t\t\t</distribution>\n"
        end
    end
    xml_string << "\n"

    if plus_gamma
        xml_string << "#{first_tab}\t\t\t<!-- Priors for the alpha parameter of the gamma distribution of site rate heterogeneity for each alignment -->\n"
        alignment_ids.each do |a|
            xml_string << "#{first_tab}\t\t\t<prior id=\"GammaShapePrior.s:#{a}\" name=\"distribution\" x=\"@gammaShape.s:#{a}\">\n"
            xml_string << "#{first_tab}\t\t\t\t<Exponential name=\"distr\" mean=\"1\"/>\n"
            xml_string << "#{first_tab}\t\t\t</prior>\n"
        end
        xml_string << "\n"
    end
end

# Close the prior part.
if starbeast
    xml_string << "#{first_tab}\t\t<!-- End of part 2 of the posterior: the prior -->\n"
else
    xml_string << "#{first_tab}\t\t<!-- End of part 1 of the posterior: the prior -->\n"
end
xml_string << "#{first_tab}\t\t</distribution>\n"
xml_string << "\n"

# Open the likelihood part.
if starbeast
    xml_string << "#{first_tab}\t\t<!-- Part 3 of the posterior: the Likelihood -->\n"
else
    xml_string << "#{first_tab}\t\t<!-- Part 2 of the posterior: the Likelihood -->\n"
end
xml_string << "#{first_tab}\t\t<distribution id=\"likelihood\" spec=\"util.CompoundDistribution\" useThreads=\"true\">\n"
xml_string << "\n"

# The likelihoods.
xml_string << "#{first_tab}\t\t\t<!-- Likelihoods for each alignment -->\n"
alignment_ids.each do |a|
    if starbeast
        tree_id = nil
        gene_tree_ids.size.times do |x|
            if gene_tree_alignment_ids[x].include?(a)
                tree_id = gene_tree_ids[x]
                break
            end
        end
        raise "Gene tree id could not be found for aligment id #{a}!" if tree_id == nil
    else
        tree_id = "Species"
    end
    xml_string << "#{first_tab}\t\t\t<distribution id=\"treeLikelihood:#{a}\" data=\"@#{a}\" spec=\"TreeLikelihood\" tree=\"@tree.t:#{tree_id}\">\n"
    if substitution_model == "bmodeltest"
        xml_string << "#{first_tab}\t\t\t\t<siteModel id=\"BEASTModelTest.s:#{a}\" spec=\"BEASTModelTestSiteModel\" gammaCategoryCount=\"4\" hasGammaRates=\"@hasGammaRates.s:#{a}\" hasInvariantSites=\"@hasInvariableSites.s:#{a}\" proportionInvariant=\"@BMT_ProportionInvariable.s:#{a}\" shape=\"@BMT_gammaShape.s:#{a}\" substModel=\"@RevJump.s:#{a}\">\n"
        xml_string << "#{first_tab}\t\t\t\t\t<parameter idref=\"mutationRate.s:#{a}\" name=\"mutationRate\">1.0</parameter>\n"
        xml_string << "#{first_tab}\t\t\t\t</siteModel>\n"
    else
        if plus_gamma
            xml_string << "#{first_tab}\t\t\t\t<siteModel id=\"SiteModel.s:#{a}\" gammaCategoryCount=\"4\" shape=\"@gammaShape.s:#{a}\" proportionInvariant=\"0.0\" mutationRate=\"@mutationRate.s:#{a}\" spec=\"SiteModel\">\n"
        else
            xml_string << "#{first_tab}\t\t\t\t<siteModel id=\"SiteModel.s:#{a}\" proportionInvariant=\"0.0\" mutationRate=\"@mutationRate.s:#{a}\" spec=\"SiteModel\">\n"
        end
        if substitution_model == "hky"
            xml_string << "#{first_tab}\t\t\t\t\t<substModel id=\"hky.s:#{a}\" kappa=\"@kappa.s:#{a}\" spec=\"HKY\">\n"
        elsif substitution_model == "gtr"
            xml_string << "#{first_tab}\t\t\t\t\t<substModel id=\"gtr.s:#{a}\" rateAC=\"@rateAC.s:#{a}\" rateAG=\"@rateAG.s:#{a}\" rateAT=\"@rateAT.s:#{a}\" rateCG=\"@rateCG.s:#{a}\" rateGT=\"@rateGT.s:#{a}\" spec=\"GTR\">\n"
        elsif substitution_model == "rb"
            xml_string << "#{first_tab}\t\t\t\t\t<substModel id=\"RB.s:#{a}\" count=\"@RBcount.s:#{a}\" rates=\"@RBrates.s:#{a}\" spec=\"RB\">\n"
        end
        if estimate_freqs
            xml_string << "#{first_tab}\t\t\t\t\t\t<frequencies frequencies=\"@freqParameter.s:#{a}\" id=\"estimatedFreqs.s:#{a}\" spec=\"Frequencies\"/>\n"
        else
            xml_string << "#{first_tab}\t\t\t\t\t\t<frequencies data=\"@#{a}\" id=\"freqs.s:#{a}\" spec=\"Frequencies\"/>\n"
        end
        xml_string << "#{first_tab}\t\t\t\t\t</substModel>\n"
        xml_string << "#{first_tab}\t\t\t\t</siteModel>\n"
    end
    if ucln_clock
        xml_string << "#{first_tab}\t\t\t\t<branchRateModel id=\"RelaxedClock.c:#{a}\" clock.rate=\"1.0\" rateCategories=\"@rateCategories.c:#{tree_id}\" spec=\"beast.evolution.branchratemodel.UCRelaxedClockModel\" tree=\"@tree.t:#{tree_id}\">\n"
        xml_string << "#{first_tab}\t\t\t\t\t<LogNormal id=\"LogNormalDistributionModel.c:#{a}\" name=\"distr\" M=\"1.0\" meanInRealSpace=\"true\" S=\"@ucldStdev.c:#{tree_id}\"/>\n"
        xml_string << "#{first_tab}\t\t\t\t</branchRateModel>\n"
    else
        xml_string << "#{first_tab}\t\t\t\t<branchRateModel id=\"strictClock.c:#{a}\" clock.rate=\"1.0\" spec=\"beast.evolution.branchratemodel.StrictClockModel\"/>\n"
    end
    xml_string << "#{first_tab}\t\t\t</distribution>\n"
end
xml_string << "\n"

# Close the likelihoods part.
if starbeast
    xml_string << "#{first_tab}\t\t<!-- End of part 3 of the posterior: the Likelihood -->\n"
else
    xml_string << "#{first_tab}\t\t<!-- End of part 2 of the posterior: the Likelihood -->\n"
end
xml_string << "#{first_tab}\t\t</distribution>\n"
xml_string << "\n"

# Close the posterior part.
xml_string << "#{first_tab}\t<!-- End of the posterior -->\n"
xml_string << "#{first_tab}\t</distribution>\n"
xml_string << "\n"

# Open the operators part.
xml_string << "#{first_tab}\t<!-- Operators -->\n"
xml_string << "\n"

# Operators on the species tree height for *BEAST.
if starbeast
    xml_string << "#{first_tab}\t<!-- A scale operator on the node heights of all gene trees -->\n"
    if taxon_ids.size < 150
        xml_string << "#{first_tab}\t<operator id=\"speciesTreeReheight\" spec=\"NodeReheight\" taxonset=\"@allTaxa\" tree=\"@tree.t:Species\" weight=\"100.0\">\n"
    else
        xml_string << "#{first_tab}\t<operator id=\"speciesTreeReheight\" spec=\"NodeReheight\" taxonset=\"@allTaxa\" tree=\"@tree.t:Species\" weight=\"400.0\">\n"
    end
    gene_tree_ids.each do |i|
        xml_string << "#{first_tab}\t\t<tree idref=\"tree.t:#{i}\" name=\"genetree\"/>\n"
    end
    xml_string << "#{first_tab}\t</operator>\n"
    xml_string << "\n"
end

# Operators on the tree height for *BEAST.
if starbeast
    gene_tree_ids.each do |i|
        xml_string << "#{first_tab}\t<!-- Operators on the height of the #{i} gene tree -->\n"
        xml_string << "#{first_tab}\t<operator id=\"treeScale:#{i}\" scaleFactor=\"0.9\" spec=\"ScaleOperator\" tree=\"@tree.t:#{i}\" weight=\"10.0\"/>\n"
        xml_string << "#{first_tab}\t<operator id=\"treeRootOnlyScale:#{i}\" rootOnly=\"true\" scaleFactor=\"0.9\" spec=\"ScaleOperator\" tree=\"@tree.t:#{i}\" weight=\"10.0\"/>\n"
        xml_string << "#{first_tab}\t<operator id=\"treeUniform:#{i}\" spec=\"Uniform\" tree=\"@tree.t:#{i}\" weight=\"30.0\"/>\n"
        xml_string << "\n"
    end
end

# Operators on the tree height.
unless starbeast
    xml_string << "#{first_tab}\t<!-- Operators on the height of the tree -->\n"
    xml_string << "#{first_tab}\t<operator id=\"treeScale:Species\" scaleFactor=\"0.9\" spec=\"ScaleOperator\" tree=\"@tree.t:Species\" weight=\"3.0\"/>\n"
    xml_string << "#{first_tab}\t<operator id=\"treeRootOnlyScale:Species\" rootOnly=\"true\" scaleFactor=\"0.9\" spec=\"ScaleOperator\" tree=\"@tree.t:Species\" weight=\"3.0\"/>\n"
    xml_string << "#{first_tab}\t<operator id=\"treeUniform:Species\" spec=\"Uniform\" tree=\"@tree.t:Species\" weight=\"30.0\"/>\n"
    xml_string << "\n"
end

# Operators on all gene tree topologies for *BEAST.
if starbeast
    gene_tree_ids.each do |i|
        xml_string << "#{first_tab}\t<!-- Operators on the topology of the #{i} gene tree -->\n"
        xml_string << "#{first_tab}\t<operator id=\"treeSubtreeSlide:#{i}\" spec=\"SubtreeSlide\" tree=\"@tree.t:#{i}\" weight=\"15.0\"/>\n"
        xml_string << "#{first_tab}\t<operator id=\"treeExchange:#{i}\" spec=\"Exchange\" tree=\"@tree.t:#{i}\" weight=\"15.0\"/>\n"
        xml_string << "#{first_tab}\t<operator id=\"treeNarrowExchange:#{i}\" isNarrow=\"false\" spec=\"Exchange\" tree=\"@tree.t:#{i}\" weight=\"3.0\"/>\n"
        xml_string << "#{first_tab}\t<operator id=\"treeWilsonBalding:#{i}\" spec=\"WilsonBalding\" tree=\"@tree.t:#{i}\" weight=\"3.0\"/>\n"
        xml_string << "\n"
    end
end

# Operators on the tree topology.
unless starbeast
    xml_string << "#{first_tab}\t<!-- Operators on the topology of the tree -->\n"
    xml_string << "#{first_tab}\t<operator id=\"treeSubtreeSlide:Species\" spec=\"SubtreeSlide\" tree=\"@tree.t:Species\" weight=\"5.0\"/>\n"
    xml_string << "#{first_tab}\t<operator id=\"treeExchange:Species\" spec=\"Exchange\" tree=\"@tree.t:Species\" weight=\"5.0\"/>\n"
    xml_string << "#{first_tab}\t<operator id=\"treeNarrowExchange:Species\" isNarrow=\"false\" spec=\"Exchange\" tree=\"@tree.t:Species\" weight=\"5.0\"/>\n"
    xml_string << "#{first_tab}\t<operator id=\"treeWilsonBalding:Species\" spec=\"WilsonBalding\" tree=\"@tree.t:Species\" weight=\"3.0\"/>\n"
    xml_string << "\n"
end

# Operators on population parameters for *BEAST.
if starbeast
    xml_string << "#{first_tab}\t<!-- Scale operators on the population size and the population mean hyperprior -->\n"
    xml_string << "#{first_tab}\t<operator id=\"popSizeScale\" degreesOfFreedom=\"1\" parameter=\"@popSize\" scaleFactor=\"0.5\" spec=\"ScaleOperator\" weight=\"5.0\"/>\n"
    xml_string << "#{first_tab}\t<operator id=\"popMeanScale\" parameter=\"@popMean\" scaleFactor=\"0.75\" spec=\"ScaleOperator\" weight=\"3.0\"/>\n"
    xml_string << "\n"
end

# Operators on the tree model.
if tree_model == "yule"
    xml_string << "#{first_tab}\t<!-- Scale operator on the birth rate of Yule model -->\n"
    xml_string << "#{first_tab}\t<operator id=\"birthRateScale\" parameter=\"@birthRate.t:Species\" scaleFactor=\"0.75\" spec=\"ScaleOperator\" weight=\"3.0\"/>\n"
elsif tree_model == "birthdeath"
    xml_string << "#{first_tab}\t<!-- Scale operators on the birth and death rates -->\n"
    xml_string << "#{first_tab}\t<operator id=\"birthRateScale\" parameter=\"@birthRate.t:Species\" scaleFactor=\"0.75\" spec=\"ScaleOperator\" weight=\"3.0\"/>\n"
    xml_string << "#{first_tab}\t<operator id=\"relativeDeathRateScale\" parameter=\"@relativeDeathRate.t:Species\" scaleFactor=\"0.8\" spec=\"ScaleOperator\" weight=\"3.0\"/>\n"
end
xml_string << "\n"

# Operators on the molecular clock.
if ucln_clock
    xml_string << "#{first_tab}\t<!-- Operators on the molecular clock for each alignment -->\n"
    if starbeast
        gene_tree_ids.each do |g|
            xml_string << "#{first_tab}\t<operator id=\"relaxedStdevScaler:#{g}\" parameter=\"@ucldStdev.c:#{g}\" scaleFactor=\"0.8\" spec=\"ScaleOperator\" weight=\"3.0\"/>\n"
            xml_string << "#{first_tab}\t<operator id=\"relaxedCategoriesRandomWalk:#{g}\" parameter=\"@rateCategories.c:#{g}\" spec=\"IntRandomWalkOperator\" weight=\"10.0\" windowSize=\"1\"/>\n"
            xml_string << "#{first_tab}\t<operator id=\"relaxedCategoriesSwap:#{g}\" intparameter=\"@rateCategories.c:#{g}\" spec=\"SwapOperator\" weight=\"10.0\"/>\n"
            xml_string << "#{first_tab}\t<operator id=\"relaxedCategoriesUniform:#{g}\" parameter=\"@rateCategories.c:#{g}\" spec=\"UniformOperator\" weight=\"10.0\"/>\n"
        end
    else
        xml_string << "#{first_tab}\t<operator id=\"relaxedStdevScaler:Species\" parameter=\"@ucldStdev.c:Species\" scaleFactor=\"0.8\" spec=\"ScaleOperator\" weight=\"3.0\"/>\n"
        xml_string << "#{first_tab}\t<operator id=\"relaxedCategoriesRandomWalk:Species\" parameter=\"@rateCategories.c:Species\" spec=\"IntRandomWalkOperator\" weight=\"10.0\" windowSize=\"1\"/>\n"
        xml_string << "#{first_tab}\t<operator id=\"relaxedCategoriesSwap:Species\" intparameter=\"@rateCategories.c:Species\" spec=\"SwapOperator\" weight=\"10.0\"/>\n"
        xml_string << "#{first_tab}\t<operator id=\"relaxedCategoriesUniform:Species\" parameter=\"@rateCategories.c:Species\" spec=\"UniformOperator\" weight=\"10.0\"/>\n"
    end
    xml_string << "\n"
end

# Operators on the substitution rates.
xml_string << "#{first_tab}\t<!-- Operators on the substitution rates for each alignment -->\n"
alignment_ids.each do |a|
    xml_string << "#{first_tab}\t<operator id=\"mutationRateScaler.s:#{a}\" parameter=\"@mutationRate.s:#{a}\" spec=\"ScaleOperator\" scaleFactor=\"0.9\" weight=\"5\"/>\n"
end
xml_string << "\n"

# Up-down operators on individual substitution rates and either the species tree or the corresponding gene trees.
if starbeast
    xml_string << "#{first_tab}\t<!-- Up-down operators on the substitution rates for each alignment and the corresponding gene tree -->\n"
else
    xml_string << "#{first_tab}\t<!-- Up-down operators on the substitution rates for each alignment and the species tree -->\n"
end
if ucln_clock
    if starbeast
        gene_tree_ids.size.times do |x|
            if taxon_ids.size < 150
                xml_string << "#{first_tab}\t<operator id=\"relaxedUpDownRateAndTree.c:#{gene_tree_ids[x]}\" scaleFactor=\"0.9\" spec=\"UpDownOperator\" weight=\"3.0\">\n"
            else
                xml_string << "#{first_tab}\t<operator id=\"relaxedUpDownRateAndTree.c:#{gene_tree_ids[x]}\" scaleFactor=\"0.9\" spec=\"UpDownOperator\" weight=\"20.0\">\n"
            end
            gene_tree_alignment_ids[x].each do |a|
                xml_string << "#{first_tab}\t\t<parameter idref=\"mutationRate.s:#{a}\" name=\"up\"/>\n"
            end
            xml_string << "#{first_tab}\t\t<tree idref=\"tree.t:#{gene_tree_ids[x]}\" name=\"down\"/>\n"
            xml_string << "#{first_tab}\t</operator>\n"
        end
    else
        alignment_ids.each do |a|
            xml_string << "#{first_tab}\t<operator id=\"relaxedUpDownRateAndTree.c:#{a}\" scaleFactor=\"0.9\" spec=\"UpDownOperator\" weight=\"3.0\">\n"
            xml_string << "#{first_tab}\t\t<parameter idref=\"mutationRate.s:#{a}\" name=\"up\"/>\n"
            xml_string << "#{first_tab}\t\t<tree idref=\"tree.t:Species\" name=\"down\"/>\n"
            xml_string << "#{first_tab}\t</operator>\n"
        end
    end
else
    if starbeast
        gene_tree_ids.size.times do |x|
            if taxon_ids.size < 150
                xml_string << "#{first_tab}\t<operator id=\"strictClockUpDownOperator.c:#{gene_tree_ids[x]}\" scaleFactor=\"0.75\" spec=\"UpDownOperator\" weight=\"3.0\">\n"
            else
                xml_string << "#{first_tab}\t<operator id=\"strictClockUpDownOperator.c:#{gene_tree_ids[x]}\" scaleFactor=\"0.75\" spec=\"UpDownOperator\" weight=\"20.0\">\n"
            end
            gene_tree_alignment_ids[x].each do |a|
                xml_string << "#{first_tab}\t\t<parameter idref=\"mutationRate.s:#{a}\" name=\"up\"/>\n"
            end
            xml_string << "#{first_tab}\t\t<tree idref=\"tree.t:#{gene_tree_ids[x]}\" name=\"down\"/>\n"
            xml_string << "#{first_tab}\t</operator>\n"
        end
    else
        alignment_ids.each do |a|
            xml_string << "#{first_tab}\t<operator id=\"strictClockUpDownOperator.c:#{a}\" scaleFactor=\"0.75\" spec=\"UpDownOperator\" weight=\"3.0\">\n"
            xml_string << "#{first_tab}\t\t<parameter idref=\"mutationRate.s:#{a}\" name=\"up\"/>\n"
            xml_string << "#{first_tab}\t\t<tree idref=\"tree.t:Species\" name=\"down\"/>\n"
            xml_string << "#{first_tab}\t</operator>\n"
        end
    end
end
xml_string << "\n"

# A joint up-down operator on all clocks, trees, and other parameters.
if starbeast
    if ucln_clock
        xml_string << "#{first_tab}\t<!-- Joint up-down operator on the means of all substitution rates and all gene trees -->\n"
    else
        xml_string << "#{first_tab}\t<!-- Joint up-down operator on all substitution rates and all gene trees -->\n"
    end
else
    if ucln_clock
        xml_string << "#{first_tab}\t<!-- Joint up-down operator on the means of all substitution rates and the tree -->\n"
    else
        xml_string << "#{first_tab}\t<!-- Joint up-down operator on all substitution rates and the tree -->\n"
    end
end
if taxon_ids.size < 150
    xml_string << "#{first_tab}\t<operator id=\"jointUpDown:Species\" scaleFactor=\"0.9\" spec=\"UpDownOperator\" weight=\"10.0\">\n"
else
    xml_string << "#{first_tab}\t<operator id=\"jointUpDown:Species\" scaleFactor=\"0.9\" spec=\"UpDownOperator\" weight=\"100.0\">\n"
end
xml_string << "#{first_tab}\t\t<parameter idref=\"birthRate.t:Species\" name=\"up\"/>\n"
if ucln_clock
    alignment_ids.each do |a|
        xml_string << "#{first_tab}\t\t<parameter idref=\"mutationRate.s:#{a}\" name=\"up\"/>\n"
    end
else
    alignment_ids.each do |a|
        xml_string << "#{first_tab}\t\t<parameter idref=\"mutationRate.s:#{a}\" name=\"up\"/>\n"
    end
end
if starbeast
    xml_string << "#{first_tab}\t\t<parameter idref=\"popMean\" name=\"down\"/>\n"
    xml_string << "#{first_tab}\t\t<parameter idref=\"popSize\" name=\"down\"/>\n"
end
xml_string << "#{first_tab}\t\t<tree idref=\"tree.t:Species\" name=\"down\"/>\n"
if starbeast
    gene_tree_ids.each do |i|
        xml_string << "#{first_tab}\t\t<tree idref=\"tree.t:#{i}\" name=\"down\"/>\n"
    end
end
xml_string << "#{first_tab}\t</operator>\n"
xml_string << "\n"

# Operators on the bModelTest subsitution model.
if substitution_model == "bmodeltest"
    xml_string << "#{first_tab}\t<!-- Operators on the bModelTest substitution model for each alignment -->\n"
    alignment_ids.each do |a|
        xml_string << "#{first_tab}\t<operator id=\"BMT_ModelTestOperator.s:#{a}\" spec=\"BMTMergeSplitOperator\" modelIndicator=\"@BMT_ModelIndicator.s:#{a}\" rates=\"@BMT_Rates.s:#{a}\" substModel=\"@RevJump.s:#{a}\" weight=\"1.0\"/>\n"
        xml_string << "#{first_tab}\t<operator id=\"BMT_Ratescaler.s:#{a}\" spec=\"BMTExchangeOperator\" delta=\"0.15\" modelIndicator=\"@BMT_ModelIndicator.s:#{a}\" rates=\"@BMT_Rates.s:#{a}\" substModel=\"@RevJump.s:#{a}\" weight=\"1.0\"/>\n"
        xml_string << "#{first_tab}\t<operator id=\"BMT_gammaShapeScaler.s:#{a}\" spec=\"BMTScaleOperator\" count=\"@hasGammaRates.s:#{a}\" parameter=\"@BMT_gammaShape.s:#{a}\" scaleFactor=\"0.5\" weight=\"0.5\"/>\n"
        xml_string << "#{first_tab}\t<operator id=\"BMT_ProportionInvariableScaler.s:#{a}\" spec=\"BMTScaleOperator\" count=\"@hasInvariableSites.s:#{a}\" parameter=\"@BMT_ProportionInvariable.s:#{a}\" scaleFactor=\"0.5\" weight=\"0.5\"/>\n"
        xml_string << "#{first_tab}\t<operator id=\"BMT_hasGammaRatesFlipper.s:#{a}\" spec=\"BMTBirthDeathOperator\" count=\"@hasGammaRates.s:#{a}\" rates=\"@BMT_gammaShape.s:#{a}\" weight=\"0.1\"/>\n"
        xml_string << "#{first_tab}\t<operator id=\"BMT_hasInvariableSitesFlipper.s:#{a}\" spec=\"BMTBirthDeathOperator\" count=\"@hasInvariableSites.s:#{a}\" rates=\"@BMT_ProportionInvariable.s:#{a}\" weight=\"0.1\"/>\n"
        xml_string << "#{first_tab}\t<operator id=\"BMT_FreqsFlipOperator.s:#{a}\" spec=\"BitFlipOperator\" parameter=\"@hasEqualFreqs.s:#{a}\" weight=\"0.1\"/>\n"
        xml_string << "#{first_tab}\t<operator id=\"BMT_FrequenciesExchanger.s:#{a}\" spec=\"DeltaExchangeOperator\" delta=\"0.01\" weight=\"0.1\">\n"
        xml_string << "#{first_tab}\t\t<parameter idref=\"BMT_frequencies.s:#{a}\"/>\n"
        xml_string << "#{first_tab}\t</operator>\n"
    end
    xml_string << "\n"
else
    # Operators on other subsitution models.
    if substitution_model == "hky"
        xml_string << "#{first_tab}\t<!-- Operators on the HKY substitution model for each alignment -->\n"
        alignment_ids.each do |a|
            xml_string << "#{first_tab}\t<operator id=\"kappaScale.s:#{a}\" parameter=\"@kappa.s:#{a}\" scaleFactor=\"0.5\" spec=\"ScaleOperator\" weight=\"1\"/>\n"
        end
    elsif substitution_model == "gtr"
        xml_string << "#{first_tab}\t<!-- Operators on the GTR substitution model for each alignment -->\n"
        alignment_ids.each do |a|
            xml_string << "#{first_tab}\t<operator id=\"rateACScaler.s:#{a}\" parameter=\"@rateAC.s:#{a}\" scaleFactor=\"0.5\" spec=\"ScaleOperator\" weight=\"0.1\"/>\n"
            xml_string << "#{first_tab}\t<operator id=\"rateAGScaler.s:#{a}\" parameter=\"@rateAG.s:#{a}\" scaleFactor=\"0.5\" spec=\"ScaleOperator\" weight=\"0.1\"/>\n"
            xml_string << "#{first_tab}\t<operator id=\"rateATScaler.s:#{a}\" parameter=\"@rateAT.s:#{a}\" scaleFactor=\"0.5\" spec=\"ScaleOperator\" weight=\"0.1\"/>\n"
            xml_string << "#{first_tab}\t<operator id=\"rateCGScaler.s:#{a}\" parameter=\"@rateCG.s:#{a}\" scaleFactor=\"0.5\" spec=\"ScaleOperator\" weight=\"0.1\"/>\n"
            xml_string << "#{first_tab}\t<operator id=\"rateGTScaler.s:#{a}\" parameter=\"@rateGT.s:#{a}\" scaleFactor=\"0.5\" spec=\"ScaleOperator\" weight=\"0.1\"/>\n"
        end
    elsif substitution_model == "rb"
        xml_string << "#{first_tab}\t<!-- Operators on the RB substitution model for each alignment -->\n"
        alignment_ids.each do |a|
            xml_string << "#{first_tab}\t<operator id=\"RBOperator.s:#{a}\" count=\"@RBcount.s:#{a}\" rates=\"@RBrates.s:#{a}\" spec=\"RBOperator\" weight=\"1.0\"/>\n"
            xml_string << "#{first_tab}\t<operator id=\"RBratescaler.s:#{a}\" count=\"@RBcount.s:#{a}\" parameter=\"@RBrates.s:#{a}\" scaleFactor=\"0.5\" spec=\"RBScaleOperator\" weight=\"1.0\"/>\n"
        end
    end
    xml_string << "\n"

    # Operators on the gamma distribution of site rate heterogeneity.
    if plus_gamma
        xml_string << "#{first_tab}\t<!-- Operators on the gamma distribution of site rate heterogeneity for each alignment -->\n"
        alignment_ids.each do |a|
            xml_string << "#{first_tab}\t<operator id=\"gammaScale.s:#{a}\" parameter=\"@gammaShape.s:#{a}\" scaleFactor=\"0.5\" spec=\"ScaleOperator\" weight=\"0.1\"/>\n"
        end
        xml_string << "\n"
    end

    # Operators on the base frequencies.
    if estimate_freqs
        xml_string << "#{first_tab}\t<!-- Operators on the base frequencies for each alignment -->\n"
        alignment_ids.each do |a|
            xml_string << "#{first_tab}\t<operator id=\"FrequenciesExchanger.s:#{a}\" parameter=\"@freqParameter.s:#{a}\" delta=\"0.01\" spec=\"DeltaExchangeOperator\" weight=\"0.1\"/>\n"
        end
        xml_string << "\n"
    end
end

# Open the log part.
xml_string << "#{first_tab}\t<!-- Loggers -->\n"
xml_string << "\n"

# The parameter logger.
xml_string << "#{first_tab}\t<!-- The parameter logger -->\n"
xml_string << "#{first_tab}\t<logger fileName=\"#{analysis_id}.log\" logEvery=\"#{store_length}\" model=\"@posterior\" sort=\"smart\">\n"
xml_string << "#{first_tab}\t\t<log idref=\"posterior\"/>\n"
xml_string << "#{first_tab}\t\t<log idref=\"likelihood\"/>\n"
xml_string << "#{first_tab}\t\t<log idref=\"prior\"/>\n"
if starbeast
    unless path_sampling
        xml_string << "#{first_tab}\t\t<log idref=\"speciescoalescent\"/>\n"
    end
end
xml_string << "#{first_tab}\t\t<log id=\"treeHeight:Species\" spec=\"beast.evolution.tree.TreeHeightLogger\" tree=\"@tree.t:Species\"/>\n"
if starbeast
    gene_tree_ids.each do |i|
        xml_string << "#{first_tab}\t\t<log idref=\"treePrior.t:#{i}\"/>\n"
    end
end
alignment_ids.each do |a|
    xml_string << "#{first_tab}\t\t<log idref=\"treeLikelihood:#{a}\"/>\n"
end
if starbeast
    gene_tree_ids.each do |i|
        xml_string << "#{first_tab}\t\t<log id=\"treeHeight.t:#{i}\" spec=\"beast.evolution.tree.TreeHeightLogger\" tree=\"@tree.t:#{i}\"/>\n"
    end
end
if tree_model == "yule"
    xml_string << "#{first_tab}\t\t<log idref=\"yuleModel.t:Species\"/>\n"
    xml_string << "#{first_tab}\t\t<log idref=\"birthRate.t:Species\"/>\n"
elsif tree_model == "birthdeath"
    xml_string << "#{first_tab}\t\t<log idref=\"birthDeath.t:Species\"/>\n"
    xml_string << "#{first_tab}\t\t<log idref=\"birthRate.t:Species\"/>\n"
    xml_string << "#{first_tab}\t\t<log idref=\"relativeDeathRate.t:Species\"/>\n"
end
if ucln_clock
    if starbeast
        gene_tree_ids.each do |g|
            xml_string << "#{first_tab}\t\t<log idref=\"ucldStdev.c:#{g}\"/>\n"
        end
    else
        xml_string << "#{first_tab}\t\t<log idref=\"ucldStdev.c:Species\"/>\n"
    end
end
alignment_ids.each do |a|
    xml_string << "#{first_tab}\t\t<log idref=\"mutationRate.s:#{a}\"/>\n"
end

if substitution_model == "bmodeltest"
    alignment_ids.each do |a|
        xml_string << "#{first_tab}\t\t<log idref=\"RevJump.s:#{a}\"/>\n"
        xml_string << "#{first_tab}\t\t<log idref=\"BMT_ModelIndicator.s:#{a}\"/>\n"
        xml_string << "#{first_tab}\t\t<log idref=\"BMT_Rates.s:#{a}\"/>\n"
        xml_string << "#{first_tab}\t\t<log idref=\"BMT_gammaShape.s:#{a}\"/>\n"
        xml_string << "#{first_tab}\t\t<log idref=\"BMT_ProportionInvariable.s:#{a}\"/>\n"
        xml_string << "#{first_tab}\t\t<log idref=\"hasGammaRates.s:#{a}\"/>\n"
        xml_string << "#{first_tab}\t\t<log idref=\"hasInvariableSites.s:#{a}\"/>\n"
        xml_string << "#{first_tab}\t\t<log id=\"ActivePropInvariable.s:#{a}\" spec=\"beast.util.Script\" argnames=\"BMT_ProportionInvariable hasInvariableSites\" expression=\"BMT_ProportionInvariable * hasInvariableSites\">\n"
        xml_string << "#{first_tab}\t\t\t<x idref=\"BMT_ProportionInvariable.s:#{a}\"/>\n"
        xml_string << "#{first_tab}\t\t\t<x idref=\"hasInvariableSites.s:#{a}\"/>\n"
        xml_string << "#{first_tab}\t\t</log>\n"
        xml_string << "#{first_tab}\t\t<log id=\"ActiveGammaShape.s:#{a}\" spec=\"beast.util.Script\" argnames=\"BMT_gammaShape hasGammaRates\" expression=\"BMT_gammaShape * hasGammaRates\">\n"
        xml_string << "#{first_tab}\t\t\t<x idref=\"BMT_gammaShape.s:#{a}\"/>\n"
        xml_string << "#{first_tab}\t\t\t<x idref=\"hasGammaRates.s:#{a}\"/>\n"
        xml_string << "#{first_tab}\t\t</log>\n"
        xml_string << "#{first_tab}\t\t<log idref=\"hasEqualFreqs.s:#{a}\"/>\n"
        xml_string << "#{first_tab}\t\t<log idref=\"BMT_frequencies.s:#{a}\"/>\n"
    end
else
    if substitution_model == "hky"
        alignment_ids.each do |a|
            xml_string << "#{first_tab}\t\t<log idref=\"kappa.s:#{a}\"/>\n"
        end
    elsif substitution_model == "gtr"
        alignment_ids.each do |a|
            xml_string << "#{first_tab}\t\t<log idref=\"rateAC.s:#{a}\"/>\n"
            xml_string << "#{first_tab}\t\t<log idref=\"rateAG.s:#{a}\"/>\n"
            xml_string << "#{first_tab}\t\t<log idref=\"rateAT.s:#{a}\"/>\n"
            xml_string << "#{first_tab}\t\t<log idref=\"rateCG.s:#{a}\"/>\n"
            xml_string << "#{first_tab}\t\t<log idref=\"rateGT.s:#{a}\"/>\n"
        end
    elsif substitution_model == "rb"
        alignment_ids.each do |a|
            xml_string << "#{first_tab}\t\t<log idref=\"RBrates.s:#{a}\"/>\n"
            xml_string << "#{first_tab}\t\t<log idref=\"RBcount.s:#{a}\"/>\n"
        end
    end
    if plus_gamma
        alignment_ids.each do |a|
            xml_string << "#{first_tab}\t\t<log idref=\"gammaShape.s:#{a}\"/>\n"
        end
    end
    if estimate_freqs
        alignment_ids.each do |a|
            xml_string << "#{first_tab}\t\t<log idref=\"freqParameter.s:#{a}\"/>\n"
        end
    end
end
xml_string << "#{first_tab}\t</logger>\n"
xml_string << "\n"

# The species tree logger.
if starbeast
    xml_string << "#{first_tab}\t<!-- The species tree log file output -->\n"
    xml_string << "#{first_tab}\t<logger fileName=\"#{analysis_id}_species.trees\" id=\"treelog.t:Species\" logEvery=\"#{store_length}\" mode=\"tree\">\n"
    xml_string << "#{first_tab}\t\t<log id=\"speciesTreeLogger\" popSize=\"@popSize\" popSizeTop=\"@popSizeTop\" spec=\"beast.evolution.speciation.SpeciesTreeLogger\" speciesTreePrior=\"@speciesTreePopSize.Species\" tree=\"@tree.t:Species\">\n"
    xml_string << "#{first_tab}\t\t\t<treetop id=\"treeTopFinder\" spec=\"beast.evolution.speciation.TreeTopFinder\">\n"
    gene_tree_ids.each do |i|
        xml_string << "#{first_tab}\t\t\t\t<tree idref=\"tree.t:#{i}\"/>\n"
    end
    xml_string << "#{first_tab}\t\t\t</treetop>\n"
    xml_string << "#{first_tab}\t\t</log>\n"
else
    xml_string << "#{first_tab}\t<!-- The species tree log file output -->\n"
    xml_string << "#{first_tab}\t<logger fileName=\"#{analysis_id}.trees\" id=\"treelog.t:Species\" logEvery=\"#{store_length}\" mode=\"tree\">\n"
    if ucln_clock
        xml_string << "#{first_tab}\t\t<log id=\"treeWithMetaDataLogger.t:Species\" branchratemodel=\"@RelaxedClock.c:#{alignment_ids[0]}\" spec=\"beast.evolution.tree.TreeWithMetaDataLogger\" tree=\"@tree.t:Species\"/>\n"
    else
        xml_string << "#{first_tab}\t\t<log id=\"treeWithMetaDataLogger.t:Species\" spec=\"beast.evolution.tree.TreeWithMetaDataLogger\" tree=\"@tree.t:Species\"/>\n"
    end
end
xml_string << "#{first_tab}\t</logger>\n"
xml_string << "\n"

# Loggers for individual gene trees of *BEAST.
if starbeast
    gene_tree_ids.each do |i|
        xml_string << "#{first_tab}\t<!-- The #{i} gene tree log file output -->\n"
        xml_string << "#{first_tab}\t<logger fileName=\"#{analysis_id}_#{i}.trees\" id=\"treelog.t:#{i}\" logEvery=\"#{store_length}\" mode=\"tree\">\n"
        xml_string << "#{first_tab}\t\t<log spec=\"beast.evolution.tree.TreeWithMetaDataLogger\" tree=\"@tree.t:#{i}\"/>\n"
        xml_string << "#{first_tab}\t</logger>\n"
        xml_string << "\n"
    end
end

# The logger to screen.
xml_string << "#{first_tab}\t<!-- The screen output -->\n"
xml_string << "#{first_tab}\t<logger logEvery=\"#{store_length/10}\" model=\"@posterior\">\n"
xml_string << "#{first_tab}\t\t<log idref=\"posterior\"/>\n"
xml_string << "#{first_tab}\t\t<log arg=\"@posterior\" spec=\"util.ESS\"/>\n"
xml_string << "#{first_tab}\t\t<log idref=\"likelihood\"/>\n"
xml_string << "#{first_tab}\t\t<log idref=\"prior\"/>\n"
xml_string << "#{first_tab}\t\t<log idref=\"treeHeight:Species\"/>\n"
xml_string << "#{first_tab}\t</logger>\n"
xml_string << "\n"

if path_sampling
    xml_string << "\t\t</mcmc>\n"
    xml_string << "\n"
end

# Close the run part.
xml_string << "\t</run>\n"
xml_string << "\n"

# Close the XML.
xml_string << "</beast>\n"

# Make the directory if it doesn't yet exist.
if output_file_directory != ""
    FileUtils.mkdir_p(output_file_directory) unless Dir.exists?(output_file_directory)
end

# Write the XML file.
xml_file = File.new("#{output_file_directory}#{output_file_name}","w")
xml_file.write(xml_string)
xml_file.close

# Feedback.
puts "Wrote file #{output_file_name}."

# Write a qsub file in the same directory as the xml output file.
if write_qsub
    qsub_string_part1 = ""
    qsub_string_part1 << "#!/bin/bash\n"
    qsub_string_part1 << "#$ -N B#{analysis_id[0..7]}\n"
    qsub_string_part1 << "#$ -cwd\n"
    qsub_string_part1 << "#$ -S /bin/bash\n"
    qsub_string_part1 << "\n"
    qsub_string_part1 << "SUBMITDIR=`pwd`\n"
    qsub_string_part1 << "\n"
    qsub_string_part1 << "cp -a $SUBMITDIR/* $TMPDIR\n"
    qsub_string_part1 << "cd $TMPDIR\n"
    qsub_string_start = "java -jar -Dbeast.user.package.dir=\".\" -Xmx4g beast.jar -seed #{rand(100000)} #{analysis_id}.xml\n"
    qsub_string_resume = "java -jar -Dbeast.user.package.dir=\".\" -Xmx4g beast.jar -seed #{rand(100000)} -resume #{analysis_id}.xml\n"
    qsub_string_part2 = "\n"
    qsub_string_part2 << "cp -a $TMPDIR/* $SUBMITDIR/\n"
    qsub_start_file_name = "#{output_file_directory}start.sh"
    qsub_start_file = File.new(qsub_start_file_name,"w")
    qsub_start_file.write(qsub_string_part1 + qsub_string_start + qsub_string_part2)
    qsub_start_file.close
    qsub_resume_file_name = "#{output_file_directory}resume.sh"
    qsub_resume_file = File.new(qsub_resume_file_name,"w")
    qsub_resume_file.write(qsub_string_part1 + qsub_string_resume + qsub_string_part2)
    qsub_resume_file.close
end

# Write slurm scripts to start and to resume the same directory as the XML output file.
if write_slurm
    slurm_string1 = ""
    slurm_string1 << "#!/bin/bash\n"
    slurm_string1 << "\n"
    slurm_string1 << "# Job name:\n"
    slurm_string1 << "#SBATCH --job-name=#{analysis_id[0..7]}\n"
    slurm_string1 << "#\n"
    slurm_string1 << "# Project:\n"
    slurm_string1 << "#SBATCH --account=nn9244k\n"
    slurm_string1 << "#\n"
    slurm_string1 << "# Wall clock limit:\n"
    slurm_string1 << "#\n"
    slurm_string1 << "#SBATCH --time=168:00:00\n"
    slurm_string1 << "#\n"
    slurm_string1 << "# Processor and memory usage:\n"
    slurm_string1 << "#SBATCH --ntasks-per-node=4\n"
    slurm_string1 << "#SBATCH --nodes=1\n"
    slurm_string1 << "#SBATCH --mem-per-cpu=5G\n"
    slurm_string1 << "#\n"
    slurm_string1 << "# Outfile:\n"
    slurm_string1 << "#SBATCH --output=#{analysis_id}.out\n"
    slurm_string1 << "\n"
    slurm_string1 << "## Set up job environment:\n"
    slurm_string1 << "source /cluster/bin/jobsetup\n"
    slurm_string1 << "module load beagle\n"
    slurm_string1 << "\n"
    slurm_string1 << "## Copy input files to the work directory:\n"
    slurm_string1 << "cp #{analysis_id}.xml* $SCRATCH\n"
    slurm_string1 << "if [ -f #{analysis_id}.log ]; then\n"
    slurm_string1 << "  cp #{analysis_id}.log $SCRATCH\n"
    slurm_string1 << "fi\n"
    slurm_string1 << "if [ -f #{analysis_id}.trees ]; then\n"
    slurm_string1 << "  cp #{analysis_id}.trees $SCRATCH\n"
    slurm_string1 << "fi\n"
    slurm_string1 << "for i in #{analysis_id}_*; do\n"
    slurm_string1 << "  if [ -f $i ]; then\n"
    slurm_string1 << "    cp $i $SCRATCH\n"
    slurm_string1 << "  fi\n"
    slurm_string1 << "done\n"
    slurm_string1 << "cp beast.jar $SCRATCH\n"
    slurm_string1 << "if [ -d RBS ]; then\n"
    slurm_string1 << "  cp -r RBS $SCRATCH\n"
    slurm_string1 << "fi\n"
    slurm_string1 << "if [ -d CA ]; then\n"
    slurm_string1 << "  cp -r CA $SCRATCH\n"
    slurm_string1 << "fi\n"
    slurm_string1 << "if [ -d BEASTLabs ]; then\n"
    slurm_string1 << "  cp -r BEASTLabs $SCRATCH\n"
    slurm_string1 << "fi\n"
    slurm_string1 << "if [ -d bModelTest ]; then\n"
    slurm_string1 << "  cp -r bModelTest $SCRATCH\n"
    slurm_string1 << "fi\n"
    slurm_string1 << "\n"
    slurm_string1 << "## Run BEAST:\n"
    slurm_string1 << "cd $SCRATCH\n"
    random = rand(100000)
    slurm_string_start = "java -jar -Dbeast.user.package.dir=\".\" -Xmx4g beast.jar -threads 3 -seed #{random} -beagle #{analysis_id}.xml\n"
    slurm_string_resume = "java -jar -Dbeast.user.package.dir=\".\" -Xmx4g beast.jar -threads 3 -seed #{random} -beagle -resume #{analysis_id}.xml\n"
    slurm_string2 = ""
    slurm_string2 << "\n"
    slurm_string2 << "# Copy files back to the submission directory\n"
    slurm_string2 << "cp #{analysis_id}* $SUBMITDIR\n"
    slurm_string2 << "\n"
    slurm_start_file_name = "#{output_file_directory}start.slurm"
    slurm_start_file = File.new(slurm_start_file_name,"w")
    slurm_start_file.write(slurm_string1 + slurm_string_start + slurm_string2)
    slurm_start_file.close
    slurm_resume_file_name = "#{output_file_directory}resume.slurm"
    slurm_resume_file = File.new(slurm_resume_file_name,"w")
    slurm_resume_file.write(slurm_string1 + slurm_string_resume + slurm_string2)
    slurm_resume_file.close
end
