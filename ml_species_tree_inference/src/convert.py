#!/usr/local/bin/python3

# Michael Matschiner, 2016-06-19
# michaelmatschiner@mac.com

# Import libraries and make sure we're on python 3.
import sys
if sys.version_info[0] < 3:
    print('Python 3 is needed to run this script!')
    sys.exit(0)
import argparse, textwrap, random, os, re
from subprocess import call

# Parse the command line arguments.
parser = argparse.ArgumentParser(
    formatter_class=argparse.RawDescriptionHelpFormatter,
    description=textwrap.dedent('''\
      %(prog)s
    -----------------------------------------
      Converts sequence data into multiple formats.
      Accepted input formats: phylip, fasta, nexus, ped
      Accepted output formats: phylip, fasta, nexus, genepop,
       stampp, treemix, treemix_bi (as treemix, but only
       biallelic SNPs).
    '''))
parser.add_argument(
    '-v', '--version',
    action='version',
    version='%(prog)s 0.95'
    )
parser.add_argument(
    '-p', '--populations',
    nargs='*',
    type=str,
    help="One or more population identifiers. If specified, only sequences\
    matching a population identifier will be included.")
parser.add_argument(
    '-x', '--exclude',
    nargs=1,
    type=str,
    default=[-1],
    dest='exclude',
    help="Name of a single individual to be excluded from the data set."
    )
parser.add_argument(
    '-f', '--format',
    nargs=1,
    type=str,
    default=[-1],
    dest='output_format',
    help="Ouput format."
    )
parser.add_argument(
    '--phased',
    action='store_true',
    dest='phased',
    help="Assume that each pair of sequences are from the same individual \
    (with phylip, fasta, nexus input formats), or that each consecutive \
    pair of bases is sorted according to haplotype (ped input format). \
    Only has an effect with genepop, stampp, or treemix output formats."
    )
parser.add_argument(
    '-m', '--mind',
    nargs=1,
    type=float,
    default=[-1],
    dest='missingness_per_individual',
    help="Allowed missingness per individual."
    )
parser.add_argument(
    '-g', '--geno',
    nargs=1,
    type=float,
    default=[-1],
    dest='missingness_per_marker',
    help="Allowed missingness per marker."
    )
parser.add_argument(
    'infile',
    nargs='?',
    type=argparse.FileType('r'),
    default='-',
    help='The input file name.'
    )
parser.add_argument(
    'outfile',
    nargs='?',
    type=argparse.FileType('w'),
    default=sys.stdout,
    help='The output file name.'
    )

# Get the command line arguments.
args = parser.parse_args()
pops = args.populations
if pops == None:
    pops = []
exclude = args.exclude[0]
output_format = args.output_format[0]
if output_format == -1:
    print("ERROR: No output format specified!")
    sys.exit(1)
phased = args.phased
missingness_per_individual = args.missingness_per_individual[0]
missingness_per_marker = args.missingness_per_marker[0]
infile = args.infile
outfile = args.outfile

# Read the input.
if infile.isatty():
    print('No input file specified, and no input piped through stdin!')
    sys.exit(0)
instring = infile.read()
inlines = instring.split('\n')

# Determine the input format.
input_format = None
if inlines[0].strip().lower() == "#nexus":
    input_format = "nexus"
elif inlines[0][0] == ">":
    input_format = "fasta"
elif len(inlines[0].split()) == 2 and int(inlines[0].split()[0]) > 0 and int(inlines[0].split()[1]) > 0:
    input_format = "phylip"
elif "[[Samples]]" in instring:
    input_format = "arlequin"
else:
    # Check wether it could be ped format.
    is_ped = True
    first_line_num_cols = len(inlines[0].strip().split())
    if first_line_num_cols < 10:
        is_ped = False
    for inline in inlines:
        if inline.strip() != "":
            if len(inline.strip().split()) != first_line_num_cols:
                is_ped = False
    if is_ped:
        input_format = "ped"
    else:
        print("ERROR: Input format could not be recognized!")
        sys.exit(1)

# Parse the input lines.
record_ids = []
record_seqs = []
if input_format == "nexus":
    in_matrix = False
    for line in inlines:
        if line.strip().lower() == 'matrix':
            in_matrix = True
        elif line.strip() == ';':
            in_matrix = False
            in_tree = False
        elif in_matrix and line.strip() is not '':
            record_ary = line.split()
            record_ids.append(record_ary[0])
            record_seqs.append(record_ary[1].upper())

elif input_format == "fasta":
    for line in inlines:
        if line.strip() != "":
            if line[0] == ">":
                record_ids.append(line[1:].strip())
                record_seqs.append("")
            elif line.strip() != "":
                record_seqs[-1] = record_seqs[-1].upper() + line.strip()

elif input_format == "phylip":
    ntax = int(inlines[0].split()[0])
    nchar = int(inlines[0].split()[1])
    for line in inlines[1:]:
        if line.strip() != "":
            record_ary = line.split()
            record_ids.append(record_ary[0])
            record_seqs.append(record_ary[1].upper())
    if len(record_ids) < ntax:
        print("ERROR: Found less sequences than specified in the phylip header!")
        sys.exit(1)
    elif len(record_ids) > ntax:
        print("ERROR: Found more sequences than specified in the phylip header!")
        sys.exit(1)
    if len(record_seqs[0]) != nchar:
        print("ERROR: The sequence length does not match the length specified in the phylip header!")

elif input_format == "ped":
    for line in inlines:
        line = line.strip()
        if line != "" and line[0] != "#":
            record_ary = line.split()
            record_seq_ary = record_ary[6:]
            if len(record_seq_ary)%2 != 0:
                print("ERROR: With ped input format, the number of alleles per individual should be an even number (is " + str(len(record_seq_ary)) + ")!")
                sys.exit(1)
            for x in range(round(len(record_seq_ary)/2)):
                allele_a_missing = True
                allele_b_missing = True
                if record_seq_ary[2*x] == "A":
                    allele_a_missing = False
                elif record_seq_ary[2*x] == "C":
                    allele_a_missing = False
                elif record_seq_ary[2*x] == "G":
                    allele_a_missing = False
                elif record_seq_ary[2*x] == "T":
                    allele_a_missing = False
                if record_seq_ary[2*x+1] == "A":
                    allele_b_missing = False
                elif record_seq_ary[2*x+1] == "C":
                    allele_b_missing = False
                elif record_seq_ary[2*x+1] == "G":
                    allele_b_missing = False
                elif record_seq_ary[2*x+1] == "T":
                    allele_b_missing = False
                if allele_a_missing == True or allele_b_missing == True:
                    record_seq_ary[2*x] = "N"
                    record_seq_ary[2*x+1] = "N"
            if phased == True:
                record_ids.append(record_ary[1] + "_a")
                record_ids.append(record_ary[1] + "_b")
                record_seq_a = ""
                record_seq_b = ""
                for x in range(round(len(record_seq_ary)/2)):
                    record_seq_a += record_seq_ary[2*x]
                    record_seq_b += record_seq_ary[2*x+1]
                record_seqs.append(record_seq_a)
                record_seqs.append(record_seq_b)
            else:
                record_ids.append(record_ary[1])
                record_seq = ""
                for x in range(round(len(record_seq_ary)/2)):
                    allele_a = record_seq_ary[2*x]
                    allele_b = record_seq_ary[2*x+1]
                    if allele_a == "A" and allele_b == "A":
                        record_seq += "A"
                    elif allele_a == "A" and allele_b == "C":
                        record_seq += "M"
                    elif allele_a == "A" and allele_b == "G":
                        record_seq += "R"
                    elif allele_a == "A" and allele_b == "T":
                        record_seq += "W"
                    elif allele_a == "C" and allele_b == "A":
                        record_seq += "M"
                    elif allele_a == "C" and allele_b == "C":
                        record_seq += "C"
                    elif allele_a == "C" and allele_b == "G":
                        record_seq += "S"
                    elif allele_a == "C" and allele_b == "T":
                        record_seq += "Y"
                    elif allele_a == "G" and allele_b == "A":
                        record_seq += "R"
                    elif allele_a == "G" and allele_b == "C":
                        record_seq += "S"
                    elif allele_a == "G" and allele_b == "G":
                        record_seq += "G"
                    elif allele_a == "G" and allele_b == "T":
                        record_seq += "K"
                    elif allele_a == "T" and allele_b == "A":
                        record_seq += "W"
                    elif allele_a == "T" and allele_b == "C":
                        record_seq += "Y"
                    elif allele_a == "T" and allele_b == "G":
                        record_seq += "K"
                    elif allele_a == "T" and allele_b == "T":
                        record_seq += "T"
                    else:
                        record_seq += "N"
                record_seqs.append(record_seq)

elif input_format == 'arlequin':
    in_alignment = False
    for inline in inlines:
        if len(inline) > 0:
            if inline.strip() == "SampleData= {":
                in_alignment = True
            elif inline.strip() == "}":
                in_alignment = False
            elif in_alignment == True:
                record_id = inline.split()[0]
                record_seq = inline.split()[2].replace("0","A").replace("1","C")
                record_ids.append(record_id)
                record_seqs.append(record_seq)
else:
    print("ERROR: Unsupported input format: " + str(input_format) + "!")
    sys.exit(1)

# Make sure the same number of records and sequences are found.
if len(record_ids) != len(record_seqs):
    print("ERROR: The number of record ids and sequences differ!")
    sys.exit(1)

# Remove records specified for exclusion.
if exclude  != -1:
    exclude_taxon_found = False
    reduced_record_ids = []
    reduced_record_seqs = []
    for x in range(len(record_ids)):
        if exclude in record_ids[x]:
            exclude_taxon_found = True
            print("INFO: Excluding " + record_ids[x])
        else:
            reduced_record_ids.append(record_ids[x])
            reduced_record_seqs.append(record_seqs[x])
    record_ids = reduced_record_ids
    record_seqs = reduced_record_seqs
    if exclude_taxon_found == False:
        print("WARNING: The taxon specified to be exluded (" + exclude + ") could not be found.")

# Make sure that all sequences have the same length.
lengths_differ = False
for x in range(1,len(record_seqs)):
    if len(record_seqs[0]) != len(record_seqs[x]):
        lengths_differ = True
if lengths_differ == True:
    print("WARNING: Not all sequences have the same length!")
    longest_length = 0
    for record_seq in record_seqs:
        if len(record_seq) > longest_length:
            longest_length = len(record_seq)
    for x in range(len(record_seqs)):
        record_length = len(record_seqs[x])
        if record_length != longest_length:
            for z in range(longest_length-record_length):
                record_seqs[x] += "-"

# Remove records that do not fit population identifiers, if population identifiers
# have been specified.
if pops != []:
    filtered_record_ids = []
    filtered_record_seqs = []
    for x in range(len(record_ids)):
        record_should_be_included = False
        for pop in pops:
            if pop in record_ids[x]:
                record_should_be_included = True
        if record_should_be_included == True:
            filtered_record_ids.append(record_ids[x])
            filtered_record_seqs.append(record_seqs[x])
    record_ids = filtered_record_ids
    record_seqs = filtered_record_seqs

# Remove individuals based on missing data, if the "mind" option is used.
if missingness_per_individual != -1:
    filtered_record_ids = []
    filtered_record_seqs = []
    for x in range(len(record_ids)):
        missing = (record_seqs[x].count("N") + record_seqs[x].count("?") + record_seqs[x].count("-"))/len(record_seqs[x])
        if missing <= missingness_per_individual:
            filtered_record_ids.append(record_ids[x])
            filtered_record_seqs.append(record_seqs[x])
    record_ids = filtered_record_ids
    record_seqs = filtered_record_seqs
    if len(record_ids) == 0:
        print("ERROR: No individuals left after filtering!")
        sys.exit(1)

# Remove sites based on missing data, if the "geno" option is used.
if missingness_per_marker != -1:
    filtered_record_seqs = []
    for x in range(len(record_ids)):
        filtered_record_seqs.append("")
    for pos in range(len(record_seqs[0])):
        missing = 0
        for x in range(len(record_ids)):
            if record_seqs[x][pos] == "N":
                missing += 1
            elif record_seqs[x][pos] == "?":
                missing += 1
            elif record_seqs[x][pos] == "-":
                missing += 1
        if missing/len(record_ids) <= missingness_per_marker:
            for x in range(len(record_ids)):
                filtered_record_seqs[x] += record_seqs[x][pos]
    record_seqs = filtered_record_seqs
    if len(record_seqs[0]) == 0:
        print("ERROR: No markers left after filtering!")
        sys.exit(0)

# Prepare the output string.
output_string = ""
if output_format == "phylip":
    output_string += str(len(record_ids)) + " " + str(len(record_seqs[0])) + "\n"
    max_id_length = 0
    for record_id in record_ids:
        if len(record_id) > max_id_length:
            max_id_length = len(record_id)
    for x in range(len(record_ids)):
        output_string += record_ids[x].ljust(max_id_length+2) + record_seqs[x] + "\n"

elif output_format == "fasta":
    for x in range(len(record_ids)):
        output_string += ">" + record_ids[x] + "\n"
        output_string += record_seqs[x] + "\n"

elif output_format == "nexus":
    output_string += "#nexus\n"
    output_string += "begin data;\n"
    output_string += "  dimensions ntax=" + str(len(record_ids)) + " nchar=" + str(len(record_seqs[0])) + ";\n"
    output_string += "  format datatype=dna missing=? gap=-;\n"
    output_string += "  matrix\n"
    max_id_length = 0
    for record_id in record_ids:
        if len(record_id) > max_id_length:
            max_id_length = len(record_id)
    for x in range(len(record_ids)):
        output_string += "  " + record_ids[x].ljust(max_id_length+2) + record_seqs[x] + "\n"
    output_string += "  ;\n"
    output_string += "end;\n"

elif output_format == "genepop":
    variable_positions = []
    for pos in range(len(record_seqs[0])):
        bases_at_this_site = []
        for record_seq in record_seqs:
            if record_seq[pos] in ["A","C","G","T"]:
                bases_at_this_site.append(record_seq[pos])
        bases_at_this_site = list(set(bases_at_this_site))
        if len(bases_at_this_site) > 1:
            variable_positions.append(pos)

    genepop_record_ids = []
    genepop_record_seqs = []
    if phased == False:
        for x in range(len(record_ids)):
            genepop_record_id = record_ids[x]
            genepop_record_seq = ""
            for variable_position in variable_positions:
                if record_seqs[x][variable_position] == "A":
                    genepop_record_seq += "\t01"
                elif record_seqs[x][variable_position] == "C":
                    genepop_record_seq += "\t02"
                elif record_seqs[x][variable_position] == "G":
                    genepop_record_seq += "\t03"
                elif record_seqs[x][variable_position] == "T":
                    genepop_record_seq += "\t04"
                else:
                    genepop_record_seq += "\t00"
            genepop_record_ids.append(genepop_record_id)
            genepop_record_seqs.append(genepop_record_seq)
    else:
        for x in range(round(len(record_ids)/2)):
            genepop_record_id = record_ids[2*x]
            genepop_record_seq = ""
            for variable_position in variable_positions:
                if record_seqs[2*x][variable_position] == "A":
                    genepop_record_seq += "\t01"
                elif record_seqs[2*x][variable_position] == "C":
                    genepop_record_seq += "\t02"
                elif record_seqs[2*x][variable_position] == "G":
                    genepop_record_seq += "\t03"
                elif record_seqs[2*x][variable_position] == "T":
                    genepop_record_seq += "\t04"
                else:
                    genepop_record_seq += "00"
                if record_seqs[(2*x)+1][variable_position] == "A":
                    genepop_record_seq += "01"
                elif record_seqs[(2*x)+1][variable_position] == "C":
                    genepop_record_seq += "02"
                elif record_seqs[(2*x)+1][variable_position] == "G":
                    genepop_record_seq += "03"
                elif record_seqs[(2*x)+1][variable_position] == "T":
                    genepop_record_seq += "04"
                else:
                    genepop_record_seq += "00"
            genepop_record_ids.append(genepop_record_id)
            genepop_record_seqs.append(genepop_record_seq)

    output_string += "converted from " + input_format + " format by converter.py" + "\n"
    for variable_position in variable_positions: 
        output_string += "SNP" + str(variable_position).rjust(len(str(variable_positions[-1]))).replace(" ","0") + "\n"
    max_id_length = 0
    for genepop_record_id in genepop_record_ids:
        if len(genepop_record_id) > max_id_length:
            max_id_length = len(genepop_record_id)
    if pops == []:
        output_string += "pop\n"
        for x in range(len(genepop_record_ids)):
            output_string += genepop_record_ids[x].ljust(max_id_length) + "\t," + genepop_record_seqs[x] + "\n"
    else:
        for pop in pops:
            output_string += "pop\n"
            for x in range(len(genepop_record_ids)):
                if pop in genepop_record_ids[x]:
                    output_string += genepop_record_ids[x].ljust(max_id_length) + "\t," + genepop_record_seqs[x] + "\n"

elif output_format == "stampp":

    if phased == False:
        print("ERROR: When output format stampp is chosen, phased sequences are required. Specify with \'--phased\'.")
        sys.exit(1)

    variable_positions = []
    alleles_at_variable_positions = []
    for pos in range(len(record_seqs[0])):
        bases_at_this_site = []
        for record_seq in record_seqs:
            if record_seq[pos] in ["A","C","G","T"]:
                bases_at_this_site.append(record_seq[pos])
        bases_at_this_site = list(set(bases_at_this_site))
        if len(bases_at_this_site) == 2:
            variable_positions.append(pos)
            alleles_at_variable_positions.append(sorted(bases_at_this_site))

    stampp_record_ids = []
    stampp_record_seqs = []
    for x in range(round(len(record_ids)/2)):
        stampp_record_id = record_ids[2*x]
        stampp_record_seq = ""
        for z in range(len(variable_positions)):
            if record_seqs[2*x][variable_positions[z]] == alleles_at_variable_positions[z][0]:
                if record_seqs[2*x+1][variable_positions[z]] == alleles_at_variable_positions[z][0]:
                    stampp_record_seq += ",AA"
                elif record_seqs[2*x+1][variable_positions[z]] == alleles_at_variable_positions[z][1]:
                    stampp_record_seq += ",AB"
                else:
                    stampp_record_seq += ",A-"
            elif record_seqs[2*x][variable_positions[z]] == alleles_at_variable_positions[z][1]:
                if record_seqs[2*x+1][variable_positions[z]] == alleles_at_variable_positions[z][0]:
                    stampp_record_seq += ",AB"
                elif record_seqs[2*x+1][variable_positions[z]] == alleles_at_variable_positions[z][1]:
                    stampp_record_seq += ",BB"
                else:
                    stampp_record_seq += ",B-"
            else:
                if record_seqs[2*x+1][variable_positions[z]] == alleles_at_variable_positions[z][0]:
                    stampp_record_seq += ",A-"
                elif record_seqs[2*x+1][variable_positions[z]] == alleles_at_variable_positions[z][1]:
                    stampp_record_seq += ",B-"
                else:
                    stampp_record_seq += ",--"
        stampp_record_ids.append(stampp_record_id)
        stampp_record_seqs.append(stampp_record_seq)

    output_string += "Inds,Pop,Ploidy,Format"
    for variable_position in variable_positions: 
        output_string += ",SNP" + str(variable_position).rjust(len(str(variable_positions[-1]))).replace(" ","0")
    output_string += "\n"

    if pops == []:
        for x in range(len(stampp_record_ids)):
            output_string += stampp_record_ids[x] + ","
            output_string += "all,2,BiA"
            output_string += stampp_record_seqs[x] + "\n"
    else:
        for pop in pops:
            for x in range(len(stampp_record_ids)):
                if pop in stampp_record_ids[x]:
                    output_string += stampp_record_ids[x] + ","
                    output_string += pop + ",2,BiA"
                    output_string += stampp_record_seqs[x] + "\n"

elif output_format == "treemix" or output_format == "treemix_bi":
    if pops == []:
        print("ERROR: Population assignments must be specified with output format 'treemix' (option -p)!")
        sys.exit(0)
    first_pop = True
    for pop in pops:
        if first_pop == True:
            output_string += pop
            first_pop = False
        else:
            output_string += " " + pop
    output_string += "\n"
    for pos in range(len(record_seqs[0])):
        # First find out what the two alleles are at this pos, for which individuals from all populations must be read.
        alleles_at_this_pos = []
        for x in range(len(record_ids)):
            for pop in pops:
                if pop in record_ids[x]:
                    if record_seqs[x][pos] == "A":
                        alleles_at_this_pos.append("A")
                        alleles_at_this_pos.append("A")
                    elif record_seqs[x][pos] == "C":
                        alleles_at_this_pos.append("C")
                        alleles_at_this_pos.append("C")
                    elif record_seqs[x][pos] == "G":
                        alleles_at_this_pos.append("G")
                        alleles_at_this_pos.append("G")
                    elif record_seqs[x][pos] == "T":
                        alleles_at_this_pos.append("T")
                        alleles_at_this_pos.append("T")
                    elif record_seqs[x][pos] == "M":
                        alleles_at_this_pos.append("A")
                        alleles_at_this_pos.append("C")
                    elif record_seqs[x][pos] == "R":
                        alleles_at_this_pos.append("A")
                        alleles_at_this_pos.append("G")
                    elif record_seqs[x][pos] == "W":
                        alleles_at_this_pos.append("A")
                        alleles_at_this_pos.append("T")
                    elif record_seqs[x][pos] == "S":
                        alleles_at_this_pos.append("C")
                        alleles_at_this_pos.append("G")
                    elif record_seqs[x][pos] == "Y":
                        alleles_at_this_pos.append("C")
                        alleles_at_this_pos.append("T")
                    elif record_seqs[x][pos] == "K":
                        alleles_at_this_pos.append("G")
                        alleles_at_this_pos.append("T")
        unique_alleles_at_this_pos = sorted(list(set(alleles_at_this_pos)))
        # Only use biallelic SNPs when treemix_bi is chosen, or both mono-
        # and biallelic SNPs when treemix is chosen.
        use_this_marker = False
        if len(unique_alleles_at_this_pos) == 2:
            use_this_marker = True
        elif len(unique_alleles_at_this_pos) == 1 and output_format == "treemix":
            use_this_marker = True
        if use_this_marker == True:
            first_pop = True
            for pop in pops:
                alleles_of_this_pop_at_this_pos = []
                if phased == True:
                    for x in range(len(record_ids)):
                        if pop in record_ids[x]:
                            alleles_of_this_pop_at_this_pos.append(record_seqs[x][pos])
                else:
                    for x in range(len(record_ids)):
                        if pop in record_ids[x]:
                            if record_seqs[x][pos] == "A":
                                alleles_of_this_pop_at_this_pos.append("A")
                                alleles_of_this_pop_at_this_pos.append("A")
                            elif record_seqs[x][pos] == "C":
                                alleles_of_this_pop_at_this_pos.append("C")
                                alleles_of_this_pop_at_this_pos.append("C")
                            elif record_seqs[x][pos] == "G":
                                alleles_of_this_pop_at_this_pos.append("G")
                                alleles_of_this_pop_at_this_pos.append("G")
                            elif record_seqs[x][pos] == "T":
                                alleles_of_this_pop_at_this_pos.append("T")
                                alleles_of_this_pop_at_this_pos.append("T")
                            elif record_seqs[x][pos] == "M":
                                alleles_of_this_pop_at_this_pos.append("A")
                                alleles_of_this_pop_at_this_pos.append("C")
                            elif record_seqs[x][pos] == "R":
                                alleles_of_this_pop_at_this_pos.append("A")
                                alleles_of_this_pop_at_this_pos.append("G")
                            elif record_seqs[x][pos] == "W":
                                alleles_of_this_pop_at_this_pos.append("A")
                                alleles_of_this_pop_at_this_pos.append("T")
                            elif record_seqs[x][pos] == "S":
                                alleles_of_this_pop_at_this_pos.append("C")
                                alleles_of_this_pop_at_this_pos.append("G")
                            elif record_seqs[x][pos] == "Y":
                                alleles_of_this_pop_at_this_pos.append("C")
                                alleles_of_this_pop_at_this_pos.append("T")
                            elif record_seqs[x][pos] == "K":
                                alleles_of_this_pop_at_this_pos.append("G")
                                alleles_of_this_pop_at_this_pos.append("T")
                num_occurrences_allele_a = 0
                num_occurrences_allele_b = 0
                for allele in alleles_of_this_pop_at_this_pos:
                    if allele == unique_alleles_at_this_pos[0]:
                        num_occurrences_allele_a += 1
                    else:
                        num_occurrences_allele_b += 1

                allele_frequency_string = str(num_occurrences_allele_a) + "," + str(num_occurrences_allele_b)
                if first_pop == True:
                    output_string += allele_frequency_string
                    first_pop = False
                else:
                    output_string += " " + allele_frequency_string
            output_string += "\n"

else:
    print("ERROR: Unsupported output format: " + str(output_format) + "!")
    sys.exit(1)

# Write the output string to file or STDOUT.
outfile.write(output_string)
