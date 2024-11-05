import sys
import collections

# Code from Notebook
def phred33_to_q(qual):
  """ Turn Phred+33 ASCII-encoded quality into Phred-scaled integer """
  return ord(qual)-33

# Parse the fasta file, code taken from previous questions/notebook
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
        reads = []
        while True:
            first_line = file.readline()
            if not first_line:
                break  # End of file
            seq = file.readline().rstrip()
            file.readline()  
            qual = file.readline().rstrip() 
            reads.append((seq, qual))
        return reads

# Function for approximate matching
def approximate_matching(read, qual, index, fasta_input, num_mismatches, kmer_len, pos_data):
    partitions = []
    num_partitions = num_mismatches + 1
    unique_start_positions = set()

    for i in range(0, len(read) - kmer_len + 1, kmer_len):
        partition = read[i:i + kmer_len]
        partitions.append((partition, i))
    
    for partition, i in partitions:
        if partition not in index:
            continue
        index_hits = index[partition]
        for index_hit in index_hits:
            start_pos = index_hit-i*kmer_len
            if start_pos < 0 or start_pos + len(read) > len(fasta_input):
                continue 
            nmm = 0
            for i in range(len(read)):
                if fasta_input[start_pos + i] != read[i]:
                    nmm += 1
                    if nmm > num_mismatches:
                        break
            if nmm <= num_mismatches:
                unique_start_positions.add(start_pos)
                for j in range(len(read)):
                    ref_pos = start_pos + j
                    if ref_pos >= len(fasta_input):
                        continue
                    ref_base = fasta_input[ref_pos]
                    read_base = read[j]
                    qual_score = phred33_to_q(qual[j])
                    if read_base != ref_base:
                        pos_data[ref_pos][read_base] += qual_score

# Main Program
if len(sys.argv) != 4:
    print("Usage: python3 hw2q4.py <fasta_file> <fastq_file> <output_file>")
    sys.exit(1)
num_mismatches = 4
kmer_len = int(30/(4 + 1))
fasta_file = sys.argv[1]
fastq_file = sys.argv[2]
output_file = sys.argv[3]
reads = parse_fastq(fastq_file)
fasta_input = parse_fasta(fasta_file)
index = create_index(fasta_input)

pos_data = collections.defaultdict(lambda:collections.defaultdict(int))
for read,qual in reads:
    approximate_matching(read, qual, index, fasta_input, num_mismatches, kmer_len, pos_data)
with open(output_file, 'w') as out_file:
    for ref_pos, base_weights in sorted(pos_data.items()):
        ref_base = fasta_input[ref_pos]
        sig_bases = {base: weight for base, weight in base_weights.items() if weight > 20}
        if sig_bases:
            sig_bases = sorted(sig_bases.items(), key=lambda x: (-x[1], x[0]))
            f_base, f_weight = sig_bases[0]
            if len(sig_bases) > 1:
                s_base, s_weight = sig_bases[1]
            else:
                s_base = '-'
                s_weight = 0    
            out_file.write(f"{ref_pos} {ref_base} {f_base} {f_weight} {s_base} {s_weight}\n")


#             base_weights = collections.defaultdict(int)
#             for name in pos_data[pos]:
#                 for base, qual_score, ref_base in pos_data[pos][name]:
#                     if base != ref_base:
#                         base_weights[base] += qual_score
#             if base_weights:
#                 sig_bases = {base: weight for base, weight in base_weights.items() if weight > 20}
#                 if sig_bases:
#                     sig_bases = sorted(sig_bases.items(), key=lambda x: (-x[1], x[0]))
#                     f_base, f_weight = sig_bases[0]
#                     if len(sig_bases) > 1:
#                         s_base, s_weight = sig_bases[1]
#                     else:
#                         s_base = '-'
#                         s_weight = 0    
#                     out_file.write(f"{pos} {ref_base} {f_base} {f_weight} {s_base} {s_weight}\n")
