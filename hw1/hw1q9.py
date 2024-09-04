# mRNA vaccine design
# Input: Given an amino-acid sequence
# Ouput: Nucleotide Sequence that Maximizes the combined number of Cs and Gs, ties resolved with lexico

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

def generate_sequences(aa_sequence):
    def generate_helper(curr_ind, curr_sequence):
        if curr_ind == len(aa_sequence):
            sequences.append(curr_sequence)
            return
        aa = aa_sequence[curr_ind]
        if aa not in table:
            return
        for codon in table[aa]:
            generate_helper(curr_ind+1, curr_sequence + codon)

    sequences = []
    generate_helper(0,"")
    return sequences

import sys
if len(sys.argv) != 3:
    print("Usage: python3 script.py <input_file> <output_file>")
    sys.exit(1)

input_file = sys.argv[1] # First argument: input filename
output_file = sys.argv[2] # Second argument: output filename

with open(input_file, 'r') as in_file:
    aa_sequence = in_file.readline().strip()
    aa_sequence = ''.join(char for char in aa_sequence if char in 'ACDEFGHIKLMNPQRSTVWY')

sequences = generate_sequences(aa_sequence)
sequences.sort(key=lambda seq: (-seq.count('C'), -seq.count('G'), seq))

with open(output_file, 'w') as out_file:
        out_file.write(sequences[0])
