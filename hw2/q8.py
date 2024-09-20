from collections import defaultdict

# Function to read FASTA files and extract the genome sequence
def read_fasta(filepath):
    sequence = []
    with open(filepath, 'r') as f:
        for line in f:
            if not line.startswith(">"):  # Ignore header lines
                sequence.append(line.strip())  # Collect sequence lines
    return ''.join(sequence)

# Function to generate all k-mers from a genome
def generate_kmers(genome, k):
    kmers = set()
    for i in range(len(genome) - k + 1):
        kmers.add(genome[i:i + k])
    return kmers

# Part (a): Find unique 8-mers in omicron
def unique_8mers_in_omicron(wuhan_kmers, alpha_kmers, delta_kmers, omicron_kmers):
    unique_8mers = omicron_kmers - (wuhan_kmers | alpha_kmers | delta_kmers)
    return len(unique_8mers)

# Part (b): Find the minimum k for a k-mer that is unique to omicron
def find_min_unique_k(wuhan_genome, alpha_genome, delta_genome, omicron_genome):
    for k in range(1, len(omicron_genome)):
        wuhan_kmers = generate_kmers(wuhan_genome, k)
        alpha_kmers = generate_kmers(alpha_genome, k)
        delta_kmers = generate_kmers(delta_genome, k)
        omicron_kmers = generate_kmers(omicron_genome, k)
        
        unique_kmers = omicron_kmers - (wuhan_kmers | alpha_kmers | delta_kmers)
        if unique_kmers:
            return k

# Part (c): Find shared 8-mers between omicron, delta, alpha but absent in wuhan
def shared_8mers_omicron_delta_alpha(wuhan_kmers, alpha_kmers, delta_kmers, omicron_kmers):
    shared_kmers = omicron_kmers & delta_kmers & alpha_kmers
    unique_shared_kmers = shared_kmers - wuhan_kmers
    return len(unique_shared_kmers)

# Load the genome sequences
wuhan_genome = read_fasta('sarscov2_wuhan.fa')
alpha_genome = read_fasta('sarscov2_alpha.fa')
delta_genome = read_fasta('sarscov2_delta.fa')
omicron_genome = read_fasta('sarscov2_omicron.fa')

# Generate 8-mers for each genome (for parts a and c)
wuhan_8mers = generate_kmers(wuhan_genome, 8)
alpha_8mers = generate_kmers(alpha_genome, 8)
delta_8mers = generate_kmers(delta_genome, 8)
omicron_8mers = generate_kmers(omicron_genome, 8)

# Part (a): Unique 8-mers in omicron
result_a = unique_8mers_in_omicron(wuhan_8mers, alpha_8mers, delta_8mers, omicron_8mers)
print(f"Total unique 8-mers in omicron: {result_a}")

# Part (b): Minimum k for unique k-mer in omicron
result_b = find_min_unique_k(wuhan_genome, alpha_genome, delta_genome, omicron_genome)
print(f"Minimum value of k with unique k-mer in omicron: {result_b}")

# Part (c): 8-mers shared by omicron, delta, and alpha but absent in wuhan
result_c = shared_8mers_omicron_delta_alpha(wuhan_8mers, alpha_8mers, delta_8mers, omicron_8mers)
print(f"Total 8-mers in omicron, delta, and alpha but absent in wuhan: {result_c}")
