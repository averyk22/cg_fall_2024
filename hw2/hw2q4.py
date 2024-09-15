# Builds a k-mer index of a genome from a FASTA file then finds approximate matches with ≤ 4 mismatches for each of the reads in a FASTQ file
# Each partition is k = |P |/(x + 1) long, we want x = 4 in this method
# Output: Write one line of output per read, in the same order as the input reads. 
# Let x be the number of mismatches allowed and k be the length of a partition. 
# Each line must have 2(x+1) items each separated by a single space. The only spaces in the output should be the ones separating items on a line.

# AI Usage: I used ChatGPT to help with the output string formatting and fed it the mismatch part of the question and asked it to explain it in different 
# terms because I was a little confused on what was being asked. 

import sys
import collections
def build_kmer_index(partition_length, fasta_file):
    with open(fasta_file, 'r') as file:
        fasta_input = ''
        for line in file:
            line = line.strip()
            if not line.startswith('>'):
                fasta_input += ''.join(c for c in line if c in ('A', 'C', 'G', 'T'))

    index = collections.defaultdict(list)
    for i in range(len(fasta_input) - partition_length + 1):  # For each 6-mer
        kmer = fasta_input[i:i+partition_length]
        index[kmer].append(i)
    return index, fasta_input

# Parse fastq file
def parse_fastq(fastq_file):
    with open(fastq_file, 'r') as file:
        reads = []
        while True:
            first_line = file.readline()
            if not first_line:
                break  # End of file
            seq = file.readline().strip()
            file.readline()  
            file.readline() 
            reads.append(seq)
        return reads

# Function for approximate matching and generating the two kinds of outputs
def approximate_matching(read, index, fasta_input, num_mismatches, kmer_len):
    partitions = []
    index_hits_counts = []
    num_partitions = num_mismatches + 1
    total_hits_dict = collections.defaultdict(list)
    result = ''

    # Creating the partitions
    for i in range(0, len(read) - kmer_len + 1, kmer_len):
        partition = read[i:i + kmer_len]
        partitions.append(partition)
    
    print(partitions)
    for p, string in enumerate(partitions):
        # Checking for condition 1 about num index hits
        index_hits = index[string]
        index_hits_counts.append(len(index_hits))

        # Checking for condition 2 mismatches
        for index_hit in index_hits:
            start_pos = index_hit - p * kmer_len
            # Check the strings that are not in bound
            if start_pos < 0 or start_pos + len(read) > len(fasta_input):
                continue 
            mismatches = 0
            for i in range(len(read)):
                if fasta_input[start_pos + i] != read[i]:
                    mismatches += 1
                    if mismatches > num_mismatches:
                        break
            print(mismatches)
            if mismatches <= num_mismatches:
                total_hits_dict[mismatches].append(start_pos)
    
    # Remove duplicates and sort offsets for each mismatch count
    for mismatches in total_hits_dict:
        total_hits_dict[mismatches] = sorted(set(total_hits_dict[mismatches]))
    # Formatting of result
    # First 5 elements that are being printed
    result = ' '.join(map(str, index_hits_counts))
    for mismatches in range(num_mismatches + 1):
        offsets = total_hits_dict.get(mismatches, [])
        if offsets:
            offsets_str = ','.join(map(str, offsets))
            result += f' {mismatches}:{offsets_str}'
        else:
            result += f' {mismatches}:'
    return result

if len(sys.argv) != 4:
    print("Usage: python3 hw2q3.py <fasta_file> <fastq_file> <output_file>")
    sys.exit(1)
    
num_mismatches = 4
kmer_len = int(30/(4 + 1))
fasta_file = sys.argv[1]
fastq_file = sys.argv[2]
output_file = sys.argv[3]
index, fasta_input = build_kmer_index(int(kmer_len), fasta_file)
print(index)
reads = parse_fastq(fastq_file)
with open(output_file, 'w') as out_file:
    for read in reads:
        result = approximate_matching(read, index, fasta_input, num_mismatches, kmer_len)
        out_file.write(result + '\n')
