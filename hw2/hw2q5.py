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

# Notebook code used: 
# - https://nbviewer.org/github/BenLangmead/comp-genomics-class/blob/master/notebooks/CG_KmerIndexHash.ipynb
# - https://nbviewer.org/github/BenLangmead/comp-genomics-class/blob/master/notebooks/FASTQ.ipynb
# - http://nbviewer.ipython.org/github/BenLangmead/comp-genomics-class/blob/master/notebooks/CG_NaiveApprox.ipynb

# Function taken from Notebook in the HW directions
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

# Function for approximate matching, code taken from some of the previous questions
def approximate_matching(name, read, qual, index, fasta_input, num_mismatches, kmer_len, pos_data):
    partitions = []
    num_partitions = num_mismatches + 1
    # Making sure there is no overriding
    unique_indexes = set()
    
    # Creating the partitions
    for i in range(0, len(read) - kmer_len + 1, kmer_len):
        partition = read[i:i + kmer_len]
        partitions.append(partition)
    # print(partitions)
    
    for p, string in enumerate(partitions):
        if string in index:
            for offset in index[string]:
                start_pos = offset - p * kmer_len
                if start_pos < 0 or start_pos + len(read) > len(fasta_input):
                    continue 
                nmm = 0
                for i in range(len(read)):
                    if fasta_input[start_pos + i] != read[i]:
                        nmm += 1
                        if nmm > num_mismatches:
                            break
                # Adding quality score if the base doesn't match ref base
                if nmm <= num_mismatches:  
                    if start_pos not in unique_indexes:
                        unique_indexes.add(start_pos)
                        for j in range(len(read)):
                            ref_pos = start_pos + j
                            # Checking bounds if reference position goes out 
                            if ref_pos < len(fasta_input):
                                ref_base = fasta_input[ref_pos]
                                base = read[j]
                                score = phred33_to_q(qual[j])
                                # Need to do this check before adding score, can just handle base weights
                                if base != ref_base:
                                    pos_data[ref_pos][base] += score

if len(sys.argv) != 4:
    print("Usage: python3 hw2q3.py <fasta_file> <fastq_file> <output_file>")
    sys.exit(1)
    
num_mismatches = 4
kmer_len = int(30/(4 + 1))
fasta_file = sys.argv[1]
fastq_file = sys.argv[2]
output_file = sys.argv[3]
fasta_input = parse_fasta(fasta_file)
index = create_index(fasta_input)
reads = parse_fastq(fastq_file)
# This holds the reference position, the non-reference base, and its quality score
pos_data = collections.defaultdict(lambda:collections.defaultdict(int))
for name, read, qual in reads:
    approximate_matching(name, read, qual, index, fasta_input, num_mismatches, kmer_len, pos_data)
#print (pos_data.items())
with open(output_file, 'w') as out_file:
    # Looking at all the the position and the weights -- sorting and outputting
    for ref_pos, weights in sorted(pos_data.items()):
        ref_base = fasta_input[ref_pos]
        if weights:
            sig_bases = {base: weight for base, weight in weights.items() if weight > 20}
            if sig_bases:
                sig_bases = sorted(sig_bases.items(), key=lambda x: (-x[1], x[0]))
                f_base, f_weight = sig_bases[0]
                if len(sig_bases) > 1:
                    s_base, s_weight = sig_bases[1]
                else:
                    s_base, s_weight  = '-', 0
                out_file.write(f"{ref_pos} {ref_base} {f_base} {f_weight} {s_base} {s_weight}\n")
