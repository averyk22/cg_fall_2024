# mRNA vaccine design
# Input: Given an amino-acid sequence
# Ouput: Nucleotide Sequence that Maximizes the combined number of Cs and Gs, ties resolved with lexico
# make another codon table with max C and G, but choose the one lexicographically greater if the are equal
# run through codon table once for each key value 
# pull all the codons and build the sequence

table = {
    'A': ['GCU', 'GCC', 'GCA', 'GCG'],
    'C': ['UGU','UGC'],
    'D': ['GAU', 'GAC'],
    'E': ['GAA', 'GAG'],  
    'F': ['UUU', 'UUC'],   
    'G': ['GGU', 'GGC', 'GGA', 'GGG'],
    'H': ['CAU', 'CAC'],  
    'I': ['AUU', 'AUC', 'AUA'],
    'K': ['AAA', 'AAG'],
    'L': ['UUA', 'UUG', 'CUU', 'CUC', 'CUA', 'CUG'],
    'M': ['AUG'], 
    'N': ['AAU', 'AAC'],  
    'P': ['CCU', 'CCC', 'CCA', 'CCG'],
    'Q': ['CAA', 'CAG'],  
    'R': ['CGU', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
    'S': ['UCU', 'UCC', 'UCA', 'UCG', 'AGU', 'AGC'],
    'T': ['ACU', 'ACC', 'ACA', 'ACG'],
    'V': ['GUU', 'GUC', 'GUA', 'GUG'],
    'W': ['UGG'],                       
    'Y': ['UAU', 'UAC'],     
    '*': ['UAA', 'UAG', 'UGA'] 
}
def count_cg(codon):
    return codon.count('C') + codon.count('G')

def modify_table(table):
    modified_table = {}
    for aa, codons in table.items():
        sorted_codons = sorted(codons, key=lambda x: (-count_cg(x), x))
        modified_table[aa] = sorted_codons[0]
    return modified_table

import sys
if len(sys.argv) != 3:
    print("Usage: python3 script.py <input_file> <output_file>")
    sys.exit(1)

input_file = sys.argv[1] # First argument: input filename
output_file = sys.argv[2] # Second argument: output filename

codon_table = modify_table(table)
aa_sequence = ''
with open(input_file, 'r') as in_file:
    for line in in_file:
        # Filter valid amino acid characters and concatenate
        aa_sequence += ''.join(char for char in line if char in 'ACDEFGHIKLMNPQRSTVWY*')


result = []
for aa in aa_sequence:
    codon = codon_table[aa]
    if codon == '*':
        break
    result.append(codon)

with open(output_file, 'w') as out_file:
    out_file.write(''.join(result))