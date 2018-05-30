# m_matschiner Wed May 30 00:56:28 CEST 2018

# Import libraries and make sure we're on python 3.
import sys
if sys.version_info[0] < 3:
    print('Python 3 is needed to run this script!')
    sys.exit(0)
import argparse, textwrap, random, os, re
import dendropy
import itertools
import copy

# Print a first empty line.
print("")

# Parse the command line arguments.
parser = argparse.ArgumentParser(
    formatter_class=argparse.RawDescriptionHelpFormatter,
    description=textwrap.dedent('''\
      %(prog)s
    -----------------------------------------
      Counts the frequencies of the three possible
      rooted topologies in one or more files in Nexus
      or Newick format, each containing one or more
      phylogenetic trees.
    '''))
parser.add_argument(
    '-v', '--version',
    action='version',
    version='%(prog)s 0.1'
    )
parser.add_argument(
    '-s', '--species',
    nargs='*',
    type=str,
    help="The IDs of multiple species included in the trees or the name of a file assigning sample IDs to species IDs.")
parser.add_argument(
    '-t', '--tree_files',
    nargs='*',
    type=str,
    help="One or more input tree files in Nexus format.")
args = parser.parse_args()
tree_file_names = args.tree_files
taxon_strings = []
for species_id in args.species:
    taxon_strings.append(species_id)

# Get the species and sample IDs from the taxon strings. If a single taxon string is given and it is the name of a file,
# read species and sample IDs from that file.
species_ids = []
sample_ids = []
if len(taxon_strings) == 1 and os.path.isfile(taxon_strings[0]):
    species_file = open(taxon_strings[0])
    species_lines = species_file.readlines()
    for line in species_lines:
        line_ary = line.split()
        sample_ids.append(line_ary[0])
        species_ids.append(line_ary[1])
    if len(set(sample_ids)) < len(set(species_ids)):
        print("WARNING: Expected sample IDs on the left and species IDs on the right in the")
        print("    table in file ", taxon_strings[0], " but apparently the columns are swapped.", sep="")
        print("    Sample IDs are thus now read from the right column and species IDs are read")
        print("    from the left column.")
        print("")
        species_ids, sample_ids = sample_ids, species_ids
else:
    for taxon_string in taxon_strings:
        species_ids.append(taxon_string)
        sample_ids.append(taxon_string)

# Make a list of unique species ids.
uniq_species_ids = sorted(list(set(species_ids)))
if len(uniq_species_ids) < 3 or len(uniq_species_ids) > 10:
    print("ERROR: Expected between 3 and 10 unique species IDs but found ", len(uniq_species_ids), "!", sep="")
    print("")
    sys.exit(0)

# Sort the sample ids according to unique species id.
sample_ids_per_unique_species = []
for uniq_species_id in uniq_species_ids:
    sample_ids_for_this_unique_species = []
    for x in range(len(sample_ids)):
        if species_ids[x] == uniq_species_id:
            sample_ids_for_this_unique_species.append(sample_ids[x])
    sample_ids_per_unique_species.append(sample_ids_for_this_unique_species)

# Get all possible combinations of one sample id per unique species id.
if len(sample_ids_per_unique_species) == 3:
    sample_id_combinations = list(itertools.product(sample_ids_per_unique_species[0],sample_ids_per_unique_species[1],sample_ids_per_unique_species[2]))
elif len(sample_ids_per_unique_species) == 4:
    sample_id_combinations = list(itertools.product(sample_ids_per_unique_species[0],sample_ids_per_unique_species[1],sample_ids_per_unique_species[2],sample_ids_per_unique_species[3]))
elif len(sample_ids_per_unique_species) == 5:
    sample_id_combinations = list(itertools.product(sample_ids_per_unique_species[0],sample_ids_per_unique_species[1],sample_ids_per_unique_species[2],sample_ids_per_unique_species[3],sample_ids_per_unique_species[4]))
elif len(sample_ids_per_unique_species) == 6:
    sample_id_combinations = list(itertools.product(sample_ids_per_unique_species[0],sample_ids_per_unique_species[1],sample_ids_per_unique_species[2],sample_ids_per_unique_species[3],sample_ids_per_unique_species[4],sample_ids_per_unique_species[5]))
elif len(sample_ids_per_unique_species) == 7:
    sample_id_combinations = list(itertools.product(sample_ids_per_unique_species[0],sample_ids_per_unique_species[1],sample_ids_per_unique_species[2],sample_ids_per_unique_species[3],sample_ids_per_unique_species[4],sample_ids_per_unique_species[5],sample_ids_per_unique_species[6]))
elif len(sample_ids_per_unique_species) == 8:
    sample_id_combinations = list(itertools.product(sample_ids_per_unique_species[0],sample_ids_per_unique_species[1],sample_ids_per_unique_species[2],sample_ids_per_unique_species[3],sample_ids_per_unique_species[4],sample_ids_per_unique_species[5],sample_ids_per_unique_species[6],sample_ids_per_unique_species[7]))
elif len(sample_ids_per_unique_species) == 9:
    sample_id_combinations = list(itertools.product(sample_ids_per_unique_species[0],sample_ids_per_unique_species[1],sample_ids_per_unique_species[2],sample_ids_per_unique_species[3],sample_ids_per_unique_species[4],sample_ids_per_unique_species[5],sample_ids_per_unique_species[6],sample_ids_per_unique_species[7],sample_ids_per_unique_species[8]))
elif len(sample_ids_per_unique_species) == 10:
    sample_id_combinations = list(itertools.product(sample_ids_per_unique_species[0],sample_ids_per_unique_species[1],sample_ids_per_unique_species[2],sample_ids_per_unique_species[3],sample_ids_per_unique_species[4],sample_ids_per_unique_species[5],sample_ids_per_unique_species[6],sample_ids_per_unique_species[7],sample_ids_per_unique_species[8],sample_ids_per_unique_species[9]))

# Read the input tree files.
n_trees_per_file = []
n_trees_total = 0
print("Reading tree files...", end='')
sys.stdout.flush()
tree_lists = []
for tree_file_name in tree_file_names:
    tree_list = dendropy.TreeList.get(path=tree_file_name, schema="nexus", preserve_underscores=True)
    tree_lists.append(tree_list)
    n_trees_per_file.append(len(tree_list))
    n_trees_total += len(tree_list)
print(" done.")
sys.stdout.flush()

# Give a warning if not all input files have the same number of trees
average_n_trees_per_file = n_trees_total/len(n_trees_per_file)
unique_n_trees_per_file = sorted(list(set(n_trees_per_file)))
if len(n_trees_per_file) > 1 and len(unique_n_trees_per_file) > 1:
    print("WARNING: Expected all input tree files to contain the same number of trees,")
    print("    but found different numbers. An average of ", n_trees_total/len(n_trees_per_file), " will be ", sep="")
    print("    assumed when weighing topologies.")
    print("")

# Get all topology strings.
print("Analyzing topology frequencies...", end='')
sys.stdout.flush()
topology_strings = []
for sample_id_combination in sample_id_combinations:
    sorted_sample_id_combination = sorted(list(sample_id_combination))
    for tree_list in tree_lists:
        for tree in tree_list:
            # Make sure that all sample ids are included in the tree:
            tree_taxon_ids = []
            for taxon_id in tree.taxon_namespace:
                tree_taxon_ids.append(taxon_id.label)
            for sample_id in sorted_sample_id_combination:
                if sample_id not in tree_taxon_ids:
                    print("ERROR: Could not find sample ", sample_id, " in  tree ", tree.as_string(schema="newick").rstrip(), "!", sep = "")
                    sys.exit(1)
            pruned_tree = copy.copy(tree)
            pruned_tree.retain_taxa_with_labels(sorted_sample_id_combination)
            for edge in pruned_tree.postorder_edge_iter():
                edge.length = None
            pruned_tree.ladderize(ascending=True)
            topology_string = pruned_tree.as_string(schema="newick").rstrip().rstrip(";").replace("'","")
            # Replace sample ids with species ids in the topology string.
            for sample_id in sample_id_combination:
                # Get the species id for this sample id.
                for x in range(len(species_ids)):
                    if sample_id == sample_ids[x]:
                        species_id = species_ids[x]
                        break
                topology_string = topology_string.replace(sample_id,species_id)
            # Sort cherries in the topology string alphabetically.
            topology_string_is_sorted = False
            while not topology_string_is_sorted:
                topology_string_is_sorted = True
                matches = re.findall('\([\w0-9_]+?,[\w0-9_]+?\)', topology_string)
                for match in matches:
                    match_ary = match.split(",")
                    species_id1 = match_ary[0][1:]
                    species_id2 = match_ary[1][0:-1]
                    if species_id1 > species_id2:
                        replace_string = "(" + species_id2 + "," + species_id1 + ")"
                        topology_string = topology_string.replace(match,replace_string)
                        topology_string_is_sorted = False
            topology_strings.append(topology_string)
print(" done.")
print("")
sys.stdout.flush()

# Count the frequencies of unique topologies and report them.
print("Found these topology frequencies:")
uniq_topology_strings = sorted(list(set(topology_strings)))
for uniq_topology_string in uniq_topology_strings:
    print("    ", uniq_topology_string, ": ", "{0:.4f}".format(topology_strings.count(uniq_topology_string)/(average_n_trees_per_file*len(sample_id_combinations))), sep="")
print("")
