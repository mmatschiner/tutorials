#!/usr/local/bin/python3

# Michael Matschiner, 2015-01-24
# michaelmatschiner@mac.com

# Import libraries and make sure we're on python 3.
import sys
if sys.version_info[0] < 3:
    print('Python 3 is needed to run this script!')
    sys.exit(0)
import argparse, textwrap, re

# Parse the command line arguments.
parser = argparse.ArgumentParser(
    formatter_class=argparse.RawDescriptionHelpFormatter,
    description=textwrap.dedent('''\
      %(prog)s
    -----------------------------------------
      Get the mean bootstrap support of a phylogenetic tree.
    '''))
parser.add_argument('-v', '--version', action='version', version='%(prog)s 1.0')
parser.add_argument(
    'infile',
    nargs='?',
    type=argparse.FileType('r'),
    default='-',
    help='the input file name.'
    )
parser.add_argument(
    'outfile',
    nargs='?',
    type=argparse.FileType('w'),
    default=sys.stdout,
    help='the output file name.'
    )
args = parser.parse_args()
infile = args.infile
outfile = args.outfile
if infile.isatty():
    print('No input file specified, and no input piped through stdin!')
    sys.exit(0)
instring = infile.read()

support_values = []
if "posterior=" in instring:
    pattern=re.compile('posterior=\d+\.*\d*')
    matches = pattern.findall(instring)
    for match in matches:
        support_values.append(float(match[10:]))
else:
    pattern = re.compile('\)\d+\.*\d*:')
    matches = pattern.findall(instring)
    for match in matches:
        support_values.append(float(match[1:-1]))
if len(support_values) > 1:
    outfile.write(str(sum(support_values)/len(support_values)) + "\n")
else:
    outfile.write("None\n")
