# canonical_kmers.py

import sys

def read_fasta(filepath):
    """
    Reads a FASTA file and returns the concatenated DNA sequence.
    Ignores any characters that are not A, C, G, or T.
    """
    sequence = []
    with open(filepath, 'r') as f:
        for line in f:
            if not line.startswith('>'):  # Skip header lines
                # Filter out non-ACGT characters and convert to uppercase
                filtered_line = ''.join([c.upper() for c in line.strip() if c.upper() in {'A', 'C', 'G', 'T'}])
                sequence.append(filtered_line)
    return ''.join(sequence)

def reverse_complement(seq):
    """
    Returns the reverse complement of a DNA sequence.
    """
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    try:
        rc = ''.join([complement[base] for base in reversed(seq)])
    except KeyError as e:
        print(f"Invalid nucleotide found: {e}")
        sys.exit(1)
    return rc

def get_canonical_kmers(sequence, k):
    """
    Generates all canonical k-mers from the given sequence.
    A k-mer is canonical if it is lexicographically less than or equal to its reverse complement.
    """
    kmers = set()
    for i in range(len(sequence) - k + 1):
        kmer = sequence[i:i+k]
        rc_kmer = reverse_complement(kmer)
        canonical_kmer = min(kmer, rc_kmer)
        kmers.add(canonical_kmer)
    return kmers

def main():
    if len(sys.argv) != 2:
        print("Usage: python3 canonical_kmers.py <input_fasta>")
        sys.exit(1)

    input_fasta = sys.argv[1]

    # Read and preprocess the sequence
    sequence = read_fasta(input_fasta)

    # Define the values of k for parts (a) and (b)
    ks = [4, 5]

    for k in ks:
        canonical_kmers = get_canonical_kmers(sequence, k)
        print (canonical_kmers)
        print(f"Number of distinct canonical {k}-mers: {len(canonical_kmers)}")

if __name__ == "__main__":
    main()
