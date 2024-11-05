# Used code from class notebook
import sys

def suffix_prefix_match(str1, str2, min_overlap):
    if len(str2) < min_overlap:
        return 0
    str2_prefix = str2[:min_overlap]
    str1_pos = -1
    while True:
        str1_pos = str1.find(str2_prefix, str1_pos + 1)
        if str1_pos == -1:
            return 0
        str1_suffix = str1[str1_pos:]
        if str2.startswith(str1_suffix):
            return len(str1_suffix)



def make_kmer_table(seqs, k):
    """ Given read dictionary and integer k, return a dictionary that
    maps each k-mer to the set of names of reads containing the k-mer. """
    table = {}
    for name, seq in seqs.items():
        for i in range(0, len(seq) - k + 1):
            kmer = seq[i:i+k]
            if kmer not in table:
                table[kmer] = set()
            table[kmer].add(name)
    return table

# Parse fastq file
def parse_fastq(fastq_file):
    with open(fastq_file, 'r') as file:
        reads = {}
        while True:
            first_line = file.readline()
            if not first_line:
                break  # End of file
            name = first_line[1:].strip()
            seq = file.readline().strip()
            file.readline()  
            file.readline() 
            reads[name] = seq
        return reads

def find_best_match(reads, k_value):
    kmer_table = make_kmer_table(reads, k_value)
    best_matches = {}

    for id, seq in reads.items():
        best_match = None
        best_overlap_len = 0
        suffix = seq[-k_value:]
        possible_reads = set()
        if suffix in kmer_table:
            possible_reads.update(kmer_table[suffix])
        possible_reads.discard(id)

        for b_read in possible_reads:
            curr_overlap_len = suffix_prefix_match(seq, reads[b_read], k_value)
            if curr_overlap_len > best_overlap_len:
                best_overlap_len = curr_overlap_len
                best_match = b_read
            elif curr_overlap_len == best_overlap_len:
                best_match = None
        
        if best_match and best_overlap_len >= k_value:
            best_matches[id] = (best_overlap_len, best_match)
            
    return best_matches

if len(sys.argv) != 4:
    print("Usage: python3 hw4q2.py <fastq_file> <k> <output_file>")
    sys.exit(1)

fastq_file = sys.argv[1]
k_value = int(sys.argv[2])
output_file = sys.argv[3]

reads = parse_fastq(fastq_file)
best_matches = find_best_match(reads, k_value)
sorted_best_matches = dict(sorted(best_matches.items()))
with open(output_file, 'w') as out_file:
    for id, (overlap, best_match) in sorted_best_matches.items():
        out_file.write(f"{id} {overlap} {best_match}\n")