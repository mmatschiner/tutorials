#!/usr/local/bin/python3

# Michael Matschiner, 2017-04-04
# michaelmatschiner@mac.com

# Import libraries and make sure we're on python 3.
import sys
if sys.version_info[0] < 3:
    print('Python 3 is needed to run this script!')
    sys.exit(0)
import argparse, textwrap, random, os, re
from subprocess import call
from Bio import AlignIO

# Parse the command line arguments.
parser = argparse.ArgumentParser(
    formatter_class=argparse.RawDescriptionHelpFormatter,
    description=textwrap.dedent('''\
      %(prog)s
    -----------------------------------------
      Combines .log and .trees files produced by BEAST.
      The input file can be a list of names of .log or
      .trees files.
    '''))
parser.add_argument(
    '-v', '--version',
    action='version',
    version='%(prog)s 0.94'
    )
parser.add_argument(
    '-b', '--burnin',
    nargs=1,
    type=float,
    default=[20],
    dest='burnin',
    help="Percentage of states to be removed as burnin per file (default = 20)."
    )
parser.add_argument(
    '-n', '--number-of-states',
    nargs=1,
    type=int,
    default=[2000],
    dest='number_of_states',
    help="Number of states to sample post-burnin per file (default = 2000 or all,\
    whatever is smaller). Specify \'0\' to use all states."
    )
parser.add_argument(
    '-s', '--sampling-scheme',
    nargs=1,
    type=str,
    default=["even"],
    dest='sampling_scheme',
    help="Sampling scheme to sample post-burnin states. Available options:\
     1.) even - sample at regular intervals between the first and the last state;\
     2.) head - sample from the first states;\
     3.) tail - sample from the last states\
     (default = even)"
    )
parser.add_argument(
    '--remove-comments',
    action='store_true',
    dest='remove_comments',
    help='When file type is \'trees\': Remove all comments in square brackets from \
    the tree string.'
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
args = parser.parse_args()
burnin = args.burnin[0]
number_of_states = args.number_of_states[0]
sampling_scheme = args.sampling_scheme[0]
remove_comments = args.remove_comments
infile = args.infile
outfile = args.outfile

if infile.isatty():
    print('No input file specified, and no input piped through stdin!')
    sys.exit(0)
instring = infile.read()
inlines = instring.split('\n')

# Determine the file type (list, log, or trees).
filetype = "list"
if inlines[0].lower().strip() == "#nexus":
    filetype = "trees"
else:
    for inline in inlines:
        if inline[0] != "#":
            if inline[0:6].lower() == "sample":
                filetype = "log"
            break

# Compile a regexp to remove comments from tree string if needed.
if filetype == "trees" and remove_comments == True:
    pattern = re.compile("\[.*?\]")

# If the file is a list, read the first file, and store the names of all
# files that need to be read later.
names_of_files_to_be_read = []
if filetype == "list":
    for inline in inlines[1:]:
        if inline.strip() != "":
            names_of_files_to_be_read.append(inline.strip())
    infile = open(inlines[0],"r")
    instring = infile.read()
    infile.close()
    inlines = instring.split('\n')
    if inlines[0].lower().strip() == "#nexus":
        filetype = "trees"
    else:
        filetype = "log"

# Get the comment part of the first log or trees file.
# For log files, use all lines before and including a lines that begins with "Sample",
# for trees files, use all lines before and including a lines that begins with "tree",
comments_string = ""
for inline in inlines:
    if filetype == "log":
        comments_string += inline + "\n"
        if inline.strip()[:6].lower() == "sample":
            break
    elif filetype == "trees":
        if inline.strip()[:4].lower() == "tree":
            break
        if inline.strip()[:1] != "[": # this is a quick and dirty way to exclude comments, but works only if the comment is only on a single line and starts with "["
            comments_string += inline + "\n"
    else:
        print("ERROR: Unexpected file type: " + filetype + "!")

# Analyse all files one by one.
body_string = ""
state_count = 0
for x in range(len(names_of_files_to_be_read)+1):

    # Read the file if this is not the first file anymore (only applies if a list has been specified).
    if x > 0:
        infile = open(names_of_files_to_be_read[x-1],"r")
        instring = infile.read()
        inlines = instring.split('\n')

    # Get the body part of this log or trees file.
    body_lines = []
    in_body = False
    if filetype == "log":
        for inline in inlines:
            if len(inline) > 0:
                if inline.strip()[:6].lower() == "sample":
                    in_body = True
                elif in_body == True:
                    body_lines.append(inline)
    elif filetype == "trees":
        # First read the comments line of this file to confirm that the comment part is identical.
        # If it isn't the files may not be combinable.
        comments_string_this_file = ""
        for inline in inlines:
            if inline.strip()[:4].lower() == "tree":
                break
            if inline.strip()[:1] != "[": # this is a quick and dirty way to exclude comments, but works only if the comment is only on a single line and starts with "["
                comments_string_this_file += inline + "\n"
        if comments_string_this_file != comments_string:
            # print(comments_string)
            print("WARNING: File " + names_of_files_to_be_read[x-1] + " has a different header than the first file in the list. These may not be combinable.")
            # sys.exit(1)
        for inline in inlines:
            if len(inline) > 0 and inline.strip() != ";" and inline.lower().strip() != "end;":
                if inline.strip()[:4].lower() == "tree":
                    in_body = True
                if in_body == True:
                    if remove_comments:
                        inline = pattern.sub("",inline)
                        body_lines.append(inline)
                    else:
                        body_lines.append(inline)
    else:
        print("ERROR: Unexpected file type: " + filetype + "!")

    # Remove the burnin from the body lines.
    if burnin > 0:
        number_of_states_to_remove = int(len(body_lines)*(burnin/100))
        body_lines = body_lines[number_of_states_to_remove:]
    
    # Downsample the states of each file.
    if number_of_states == 0 or len(body_lines) < number_of_states:
        sampled_body_lines = body_lines
    else:
        states_to_sample = []
        if sampling_scheme == "even":
            for n in range(number_of_states):
                decimal_state = n * ((len(body_lines)-1)/(number_of_states-1))
                rounded_state = round(decimal_state)
                states_to_sample.append(rounded_state)
        elif sampling_scheme == "head":
            states_to_sample = range(0,number_of_states)
        elif sampling_scheme == "tail":
            states_to_sample = range((len(body_lines)-number_of_states),len(body_lines))
        else:
            print("ERROR: Unknonw sampling scheme: " + sampling_scheme)
            sys.exit(1)
        sampled_body_lines = []
        for state in states_to_sample:
            sampled_body_lines.append(body_lines[state])

    # Add the sampled body lines to the body string.
    if filetype == "log":
        for body_line in sampled_body_lines:
            body_line_ary = body_line.split("\t")
            body_string += str(state_count)
            for item in body_line_ary[1:]:
                body_string += "\t" + item
            body_string += "\n"
            state_count += 1
    elif filetype == "trees":
        for body_line in sampled_body_lines:
            body_line_ary = body_line.split("=")
            body_string += "tree " + str(state_count) + " "
            for item in body_line_ary[1:]:
                body_string += "=" + item
            body_string += "\n"
            state_count += 1
    else:
        print("ERROR: Unexpected file type: " + filetype + "!")

# Prepare the output string.
out_string = comments_string + body_string
if filetype == "trees":
    out_string += "End;\n"

# Write the output string to file or STDOUT.
outfile.write(out_string)
