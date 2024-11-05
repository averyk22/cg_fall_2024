# Used code from class notebook: https://colab.research.google.com/github/BenLangmead/comp-genomics-class/blob/master/notebooks/CG_BWT_SimpleBuild.ipynb#scrollTo=kzOHVO7ujlrX
import sys

if len(sys.argv) != 3:
    print("Usage: python3 hw3q3.py <input_file> <output_file>")
    sys.exit(1)

# Assign the input and output filenames from command-line arguments
input_file = sys.argv[1] # First argument: input filename
output_file = sys.argv[2] # Second argument: output filename

# Read in the input
with open(input_file, 'r') as in_file:
    string = in_file.read().strip()

# sort all the suffixes
suffixes = sorted([string[i:] for i in range(len(string))])

# get the index for each of the suffixes from the string
indices = []
for suffix in suffixes:
    indices.append(len(string)-len(suffix))

# Create permutations
permutations = []
for i in range(len(indices)):
    if indices[i] == 0:
        permutations.append(len(suffixes)-1)
    else:
        permutations.append(indices[i] - 1)

# Write to output file
with open(output_file, 'w') as out_file:
    out_file.write(' '.join(map(str, permutations)) + '\n')
