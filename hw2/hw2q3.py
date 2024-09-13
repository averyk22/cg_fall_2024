import sys
import collections

def build_6mer(fasta_file):
    """Build a 6-mer index from the FASTA file."""
    with open(fasta_file, 'r') as file:
        fasta_input = ''
        for line in file:
            line = line.strip()
            if not line.startswith('>'):
                fasta_input += ''.join(c for c in line if c in ('A', 'C', 'G', 'T'))

    index = collections.defaultdict(list)
    for i in range(len(fasta_input) - 6 + 1):  # For each 6-mer
        kmer = fasta_input[i:i+6]
        index[kmer].append(i)
    return index, fasta_input

def parse_fastq(fastq_file):
    """Parse reads from a FASTQ filehandle."""
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

# Find exact matches and the 3 summary info
def exact_matching(reads, index, fasta_input):
    positions = collections.defaultdict(int)
    total_hits_dict = set()
    total_hits = 0
    for read in reads:
        read_len = len(read)
        if read_len < 6:
            continue
        for i in range(len(fasta_input) - read_len + 1):
            if fasta_input[i:i + read_len] == read:
                total_hits_dict.add(i)
        for i in range(read_len - 6 + 1):
            kmer = read[i:i+6]
            if kmer in index:
                for position in index[kmer]:
                    positions[position] += 1
    
    # Find the k-mer(s) with the most index hits
    max_hits = max(len(index[kmer]) for kmer in index)
    most_hits_key = [kmer for kmer in index if len(index[kmer]) == max_hits]
    most_hits_key.sort()
    
    # Find the total number of matches
    total_hits = len(total_hits_dict)

    # Find the smallest index with the maximum number of read alignments
    if positions:
        max_read_align = max(positions.values())
        min_index_max_read = min(index for index, count in positions.items() if count == max_read_align)
    else:
        min_index_max_read = -1
    
    return ','.join(most_hits_key), total_hits, min_index_max_read

def main():
    if len(sys.argv) != 4:
        print("Usage: python3 hw2q3.py <fasta_file> <fastq_file> <output_file>")
        sys.exit(1)

    fasta_file = sys.argv[1]
    fastq_file = sys.argv[2]
    output_file = sys.argv[3]

    index, fasta_input = build_6mer(fasta_file)
    reads = parse_fastq(fastq_file)

    most_hits_key, total_hits, min_index_max_read = exact_matching(reads, index, fasta_input)

    with open(output_file, 'w') as out_file:
        out_file.write(f"{most_hits_key} {total_hits} {min_index_max_read}\n")

if __name__ == "__main__":
    main()
