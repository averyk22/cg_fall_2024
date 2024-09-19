# Build a k-mer index of a genome for a given value of k from a FASTA file
# Input: FASTA file and k.txt which contains a single positive integer value
# Output: two values separated by space to the output file: total # of distinct keys and the number of keys that appear once exactly 

import sys
import collections
from collections import defaultdict

if len(sys.argv) != 4:
    print("Usage: python3 hw2q2.py <input_file> <input_file> <output_file>")
    sys.exit(1)

# Assign the input and output filenames from command-line arguments
fasta_file = sys.argv[1] # First argument: input filename
k_file = sys.argv[2] # Second argument: output filename
output_file = sys.argv[3]

# Handle the fasta file and collect the valid bases in the file
with open(fasta_file, 'r') as fasta_file:
    fasta_input = ''
    for line in fasta_file:
        line = line.strip()
        if not line.startswith('>'):
            for c in line:
                if c in ('ACGT'):
                    fasta_input += c

# Read the value of k that is needed to build the k-mer index
with open(k_file, 'r') as k_file:
    k_index = k_file.readline().strip()
    k_index = int(k_index)

# Code modified from k-mer implementation from notebook (https://nbviewer.org/github/BenLangmead/comp-genomics-class/blob/master/notebooks/CG_KmerIndexHash.ipynb)
index = collections.defaultdict()
for i in range(len(fasta_input) - k_index + 1):  # for each k-mer
    kmer = fasta_input[i:i+k_index]
    if kmer not in index:
        index[kmer] = [i]
    else:
        index[kmer].append(i)

# Counting how many keys appear exactly once
exactly_one = 0
for value in index.values():
    if len(value) == 1:
        exactly_one += 1

# Counting total number of keys
total_keys = len(index.keys())
with open(output_file, 'w') as out_file:
    out_file.write(f"{total_keys} {exactly_one}")