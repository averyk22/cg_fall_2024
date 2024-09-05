# Translation
# Takes an RNA string and outputs the counts of the amino acids in the amino acid sequence by translating the RNA codons
# Output should be a series of 20 numbers separated by commas
# Ignore any non-uppercase, non-A/C/G/U characters in the input and write the output to a single line
# Ex: CGCCAGAUCGCCCCCGGCCAG (input) 1,0,0,0,0,1,0,1,0,0,0,0,1,2,1,0,0,0,0,0 (output)

import sys
import collections

# Codon table
table = {'UUU': 'F', 'UUC': 'F', 'UUA': 'L', 'UUG': 'L',
    'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S',
    'UAU': 'Y', 'UAC': 'Y', 'UAA': '*', 'UAG': '*',
    'UGU': 'C', 'UGC': 'C', 'UGA': '*', 'UGG': 'W',
    'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L',
    'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'CAU': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AUU': 'I', 'AUC': 'I', 'AUA': 'I', 'AUG': 'M',
    'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'AAU': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGU': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'V',
    'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'GAU': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'}

# Amino acid labels
amino_acids = 'ACDEFGHIKLMNPQRSTVWY'

# Convert RNA 
def convert_rna_to_aa_sequence(rna):
    aa_sequence = ''
    for i in range(0, len(rna)- 2, 3):
        codon = rna[i: i + 3]
        aa = table.get(codon, '*')
        if aa != '*':
            aa_sequence += aa
    return aa_sequence

if len(sys.argv) != 3:
    print("Usage: python3 hw1q5b.py <input_file> <output_file>")
    sys.exit(1)

# Assign the input and output filenames from command-line arguments
input_file = sys.argv[1] # First argument: input filename
output_file = sys.argv[2] # Second argument: output filename
rna_sequence = ''
with open(input_file, 'r') as in_file, open(output_file, 'w') as out_file:
    for line in in_file:
        for c in line:
            if c in 'ACGU':
                rna_sequence += c
    protein_sequence = convert_rna_to_aa_sequence(rna_sequence)

    # Create a dictionary to be able to count aa for sequence
    aa_dict = collections.Counter(protein_sequence)
    result = [aa_dict.get(aa, 0) for aa in amino_acids]
    out_file.write(','.join(map(str, result)) + '\n')
