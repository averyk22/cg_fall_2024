# Used code from class notebook: https://colab.research.google.com/github/BenLangmead/comp-genomics-class/blob/master/notebooks/CG_BWT_SimpleBuild.ipynb#scrollTo=kzOHVO7ujlrX

import sys

# Function taken from notebook
def suffixArray(s):
    satups = sorted([(s[i:], i) for i in range(len(s))])
    return map(lambda x: x[1], satups)

if len(sys.argv) != 3:
    print("Usage: python3 hw3q1.py <input_file> <output_file>")
    sys.exit(1)

# Assign the input and output filenames from command-line arguments
input_file = sys.argv[1] # First argument: input filename
output_file = sys.argv[2] # Second argument: output filename

string = ''

# Read in the input
with open(input_file, 'r') as in_file:
    string = in_file.read().strip()

suffix_array = list(suffixArray(string))

# Output suffix array
with open(output_file, 'w') as out_file:
    out_file.write(' '.join(map(str, suffix_array)))
