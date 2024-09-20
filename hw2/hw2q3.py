# Build 6mer index and perform exact matching
# Report the key that got the most index hits, number of times any reads occur in reference genome, smallest 0-based index in the reference genome
# Notebook links references:
# - https://nbviewer.org/github/BenLangmead/comp-genomics-class/blob/master/notebooks/CG_Naive.ipynb
# - https://nbviewer.org/github/BenLangmead/comp-genomics-class/blob/master/notebooks/CG_KmerIndexHash.ipynb
# - https://nbviewer.org/github/BenLangmead/comp-genomics-class/blob/master/notebooks/FASTQ.ipynb

import sys
import collections

# Parses the fasta, code taken from previous questions/notebook
def parse_fasta(fasta_file):
    with open(fasta_file, 'r') as file:
        fasta_input = ''
        for line in file:
            line = line.strip()
            if not line.startswith('>'):
                fasta_input += ''.join(c for c in line if c in ('A', 'C', 'G', 'T'))
    return fasta_input

# Create the index of kmers given fasta input, code taken from previous questions
def create_index(fasta_input):
    index = collections.defaultdict()
    for i in range(len(fasta_input) - 6 + 1):  # For each 6-mer
        kmer = fasta_input[i:i+6]
        if kmer not in index:
            index[kmer] = [i]
        else:
            index[kmer].append(i)
    return index

# Parse fastq file, code taken from previous questions
def parse_fastq(fastq_file):
    with open(fastq_file, 'r') as file:
        # Don't look at duplicate reads
        reads = set()
        while True:
            first_line = file.readline()
            if not first_line:
                break  # End of file
            seq = file.readline().rstrip()
            file.readline()
            file.readline()
            reads.add(seq)
        return list(reads)


# Find exact matches and the 3 summary info, looked at jupyter notebook for reference
def exact_matching(reads, index, fasta_input):
    positions = collections.defaultdict(int)
    total_hits_dict = collections.defaultdict(list)
    max_indices_set = set()
    for read in reads:
        read_len = len(read)
        # Checking to make sure there is a p to compare
        if read_len < 6:
            continue
        # First 6 characters in read
        first_6 = read[:6]
        # Leftover characters in read
        offset_length = len(read) - 6
        offset_chars_first = read[6:]
        if first_6 in index:
            for idx in index[first_6]:
                # Nothing to even check if part of fasta will not be in it
                if idx + read_len > len(fasta_input):
                    continue
                # There is a match at idx, so look at the rest of the characters
                if fasta_input[idx + 6:idx+read_len] == offset_chars_first:
                    total_hits_dict[read].append(idx)
            
    # Find the k-mer(s) with the most index hits
    max_hits = max(len(index[kmer]) for kmer in index)
    most_hits_key = [kmer for kmer in index if len(index[kmer]) == max_hits]
    most_hits_key.sort()
    
    # Find the total number of matches
    total_hits = sum(len(positions) for positions in total_hits_dict.values())
    
    # Find the smallest index with the maximum number of read alignments
    max_positions_count = max(len(indices) for indices in total_hits_dict.values())
    max_keys = [key for key, positions in total_hits_dict.items() if len(positions) == max_positions_count]
    for key in max_keys:
        max_indices_set.update(total_hits_dict[key])
    smallest_index = min(max_indices_set)
    return ','.join(most_hits_key), total_hits, smallest_index

# Main Program
if len(sys.argv) != 4:
    print("Usage: python3 hw2q3.py <fasta_file> <fastq_file> <output_file>")
    sys.exit(1)

fasta_file = sys.argv[1]
fastq_file = sys.argv[2]
output_file = sys.argv[3]

fasta_input = parse_fasta(fasta_file)
index = create_index(fasta_input)
reads = parse_fastq(fastq_file)

most_hits_key, total_hits, min_index_max_read = exact_matching(reads, index, fasta_input)

with open(output_file, 'w') as out_file:
    out_file.write(f"{most_hits_key} {total_hits} {min_index_max_read}\n")