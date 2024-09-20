
# build a k-mer index from a FASTA file, 
# locate approximate matches for each read in a FASTQ, 
# scan the reference FASTA to categorize specific loci
# use the quality scores of the reference alternatives to identify the reads which consistently differ from the reference
# For each read that had an occurrence overlapping this reference offset, extract the specific read base and quality value overlapping the offset.
# If no reads overlapped the offset, or if all read bases overlapping the offset match the reference, continue without printing anything
# If any of the read bases overlapping the offset do not match the reference- calculate a “total weight” for each base, equal to the sum of all the associated quality values
# If any non-reference base has total weight is greater than 20:
# Output: The offset, 0-based. The reference base. The non-reference base with the greatest total weight. That greatest total weight.
# The non-reference base with the second-greatest total weight (‘‘-’’ if there were no other bases or the second-greatest weight not greater than 20).
# That second-greatest total weight (0 if there were no other bases or the second- greatest weight not greater than 20).

# AI Usage: I was not passing the non-public case on Gradescope, but my q4 approximate matching was passing. I fed in my q4 approximate matching function and the code I had written to ask it if there were any differences.
# GPT pointed out I forgot an indent for the index hits line, and this fix allowed me to pass the test. 

import sys
import collections

# Function taken from Notebook in the HW directions
def phred33_to_q(qual):
  """ Turn Phred+33 ASCII-encoded quality into Phred-scaled integer """
  return ord(qual)-33

# Building a kmer index based off of partition length fed in and fasta file
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
        if kmer not in index:
            index[kmer] = [i]
        else:
            index[kmer].append(i)
    return index, fasta_input

# Parse fastq file but want to save name, sequence, and quality -- taken from jupyter notebook
def parse_fastq(fastq_file):
    with open(fastq_file, 'r') as file:
        reads = []
        while True:
            name = file.readline().rstrip()
            if len(name) == 0:
                break  # end of file
            seq = file.readline().strip()
            file.readline()  
            qual = file.readline().strip()
            reads.append((name, seq, qual))
        return reads

# Function for approximate matching 
def approximate_matching(name, read, qual, index, fasta_input, num_mismatches, kmer_len, pos_data):
    partitions = []
    num_partitions = num_mismatches + 1
    unique_start_positions = set()
    
    # Creating the partitions
    for i in range(0, len(read) - kmer_len + 1, kmer_len):
        partition = read[i:i + kmer_len]
        partitions.append((partition, i))
    
    for p, (string, p_offset) in enumerate(partitions):
        index_hits = index[string]
        for index_hit in index_hits:
            start_pos = index_hit - p * kmer_len
            if start_pos < 0 or start_pos + len(read) > len(fasta_input):
                continue 
            mismatches = 0
            for i in range(len(read)):
                if fasta_input[start_pos + i] != read[i]:
                    mismatches += 1
                    if mismatches > num_mismatches:
                        break
            # Collecting the data that has base in the read, base in the reference, and the quality score at the base
            if mismatches <= num_mismatches:  
                unique_start_positions.add(start_pos)
                for i in range(len(read)):
                    ref_pos = start_pos + i 
                    base = read[i]
                    qual_score = phred33_to_q(qual[i])
                    ref_base = fasta_input[ref_pos]
                    if name not in pos_data[ref_pos]:
                        pos_data[ref_pos][name].append((base, qual_score, ref_base))

if len(sys.argv) != 4:
    print("Usage: python3 hw2q3.py <fasta_file> <fastq_file> <output_file>")
    sys.exit(1)
    
num_mismatches = 4
kmer_len = int(30/(4 + 1))
fasta_file = sys.argv[1]
fastq_file = sys.argv[2]
output_file = sys.argv[3]
index, fasta_input = build_kmer_index(int(kmer_len), fasta_file)
reads = parse_fastq(fastq_file)
pos_data = collections.defaultdict(lambda:collections.defaultdict(list))
for name, read, qual in reads:
    approximate_matching(name, read, qual, index, fasta_input, num_mismatches, kmer_len, pos_data)
with open(output_file, 'w') as out_file:
    for pos in range(len(fasta_input)):
        if pos in pos_data:
            ref_base = fasta_input[pos]
            base_weights = collections.defaultdict(int)
            for name in pos_data[pos]:
                for base, qual_score, ref_base in pos_data[pos][name]:
                    if base != ref_base:
                        base_weights[base] += qual_score
            if base_weights:
                sig_bases = {base: weight for base, weight in base_weights.items() if weight > 20}
                if sig_bases:
                    sig_bases = sorted(sig_bases.items(), key=lambda x: (-x[1], x[0]))
                    f_base, f_weight = sig_bases[0]
                    if len(sig_bases) > 1:
                        s_base, s_weight = sig_bases[1]
                    else:
                        s_base = '-'
                        s_weight = 0    
                    out_file.write(f"{pos} {ref_base} {f_base} {f_weight} {s_base} {s_weight}\n")
