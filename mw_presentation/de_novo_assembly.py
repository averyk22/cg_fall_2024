import sys

# This function checks to see if the suffix of str1 overlaps the prefix of str2 for min values
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


# This kmer table will store all the names of reads associated with a specific kmer of size k
# seqs: names and reads
def make_kmer_table(seqs, k):
    table = collections.defaultdict()
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

import sys
from collections import defaultdict

def find_best_right_matches(old_output_file):
    right_matches = defaultdict(list)
    with open(old_output_file, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) == 3:
                read_a, overlap_len, read_b = parts
                overlap_len = int(overlap_len)
                right_matches[read_a] = (read_b, overlap_len)
    return right_matches

def find_best_left_matches(right_matches):
    best_left = defaultdict(list)
    for read_a, (read_b, overlap_length) in right_matches.items():
        best_left[read_b].append((read_a, overlap_length))
    checked_best_l = {}
    for read_b, candidates in best_left.items():
        if len(candidates) == 1:
            checked_best_l[read_b] = candidates[0]
        else:
            candidates.sort(key=lambda x: x[1], reverse=True)
            if candidates[0][1] > candidates[1][1]:
                checked_best_l[read_b] = candidates[0]               
    return checked_best_l

def find_and_write_unitigs(best_right, best_left):
    visited = set()
    unitigs = []
    for read in best_right:
        unitig = []
        curr_read, curr_path = read, set()
        if read in visited:
            continue
        while curr_read in best_left:
            prev_read, o_len = best_left[curr_read]
            if prev_read not in visited and prev_read in best_right and best_right[prev_read][0] == curr_read:
                if prev_read in curr_path:
                    break
                curr_path.add(prev_read)
                curr_read = prev_read
        curr_path.clear()

        while curr_read in best_right and curr_read not in visited:
            unitig.append(curr_read)
            visited.add(curr_read)
            curr_path.add(curr_read)
            next_read, o_len = best_right[curr_read]
            if best_left.get(next_read, (None,))[0] == curr_read: 
                if next_read in curr_path:
                    break
                unitig.append((o_len, next_read))
                curr_read = next_read
        if len(unitig) > 1:
            unitigs.append(unitig)
    unitigs.sort(key=lambda x: x[0]) 
    return unitigs

# Main portion
if len(sys.argv) != 3:
    print("Usage: python3 hw4q3.py <old_output_file> <output_file>")
    sys.exit(1)

old_output_file = sys.argv[1]
output_file = sys.argv[2]

best_right = find_best_right_matches(old_output_file)
best_left = find_best_left_matches(best_right)

unitigs = find_and_write_unitigs(best_right, best_left)

with open(output_file, 'w') as f:
    for unitig in unitigs:
        f.write(unitig[0] + "\n")
        for item in unitig[1:]:
            if isinstance(item, tuple) and len(item) == 2:
                overlap_length, read = item
                f.write(f"{overlap_length} {read}\n")

import sys
from collections import defaultdict
import itertools

def read_unitigs(unitig_file):
    unitigs = []
    with open(unitig_file, 'r') as f:
        lines = f.read().splitlines()
    i = 0
    while i < len(lines):
        read_ids = [lines[i]]
        i += 1
        while i < len(lines) and lines[i]:
            parts = lines[i].split()
            if len(parts) == 2:
                overlap_len, read_id = parts
                read_ids.append(read_id)
            i += 1
        unitigs.append(read_ids)
    return unitigs

def load_reads(fastq_file):
    reads = {}
    with open(fastq_file, 'r') as f:
        while True:
            name = f.readline().strip()
            seq = f.readline().strip()
            plus = f.readline()
            qual = f.readline()
            if not qual:
                break
            read_id = name[1:].split()[0]
            reads[read_id] = seq
    return reads

def build_unitig_sequences(unitigs, reads, k):
    unitig_sequences = {}
    for idx, unitig in enumerate(unitigs):
        seq = reads[unitig[0]]
        for read_id in unitig[1:]:
            overlap_len = k - 1  # Since overlap length is k-1
            seq += reads[read_id][overlap_len:]
        unitig_name = f"Unitig_{idx}"
        unitig_sequences[unitig_name] = seq
    return unitig_sequences

def compute_overlaps(unitig_sequences, min_overlap):
    overlaps = defaultdict(dict)
    unitig_names = list(unitig_sequences.keys())
    for u1, u2 in itertools.permutations(unitig_names, 2):
        seq1 = unitig_sequences[u1]
        seq2 = unitig_sequences[u2]
        max_overlap_len = min(len(seq1), len(seq2))
        for olen in range(max_overlap_len, min_overlap - 1, -1):
            if seq1[-olen:] == seq2[:olen]:
                overlaps[u1][u2] = olen
                break
    return overlaps

def find_best_path(unitig_sequences, overlaps):
    max_total_overlap = -1
    best_path = None
    unitig_names = list(unitig_sequences.keys())

    def dfs(path, used, total_overlap):
        nonlocal max_total_overlap, best_path
        last_unitig = path[-1]
        if len(used) == len(unitig_names):
            if total_overlap > max_total_overlap:
                max_total_overlap = total_overlap
                best_path = path[:]
            return
        for next_unitig in overlaps[last_unitig]:
            overlap_len = overlaps[last_unitig][next_unitig]
            path.append(next_unitig)
            used.add(next_unitig)
            dfs(path, used, total_overlap + overlap_len)
            path.pop()
            used.discard(next_unitig)

    for start_unitig in unitig_names:
        dfs([start_unitig], {start_unitig}, 0)

    return best_path

def assemble_sequence(best_path, unitig_sequences, overlaps):
    assembled_seq = unitig_sequences[best_path[0]]
    for i in range(len(best_path) - 1):
        u1 = best_path[i]
        u2 = best_path[i+1]
        olen = overlaps[u1][u2]
        assembled_seq += unitig_sequences[u2][olen:]
    return assembled_seq

# Main portion
if len(sys.argv) != 5:
    print("Usage: python3 assembler.py <unitig_file> <reads_fastq> <k-mer_size> <output_file>")
    sys.exit(1)

unitig_file = sys.argv[1]
reads_fastq = sys.argv[2]
kmer_size = int(sys.argv[3])
output_file = sys.argv[4]

unitigs = read_unitigs(unitig_file)
reads = load_reads(reads_fastq)
unitig_sequences = build_unitig_sequences(unitigs, reads, kmer_size)
min_overlap = kmer_size - 1
overlaps = compute_overlaps(unitig_sequences, min_overlap)
best_path = find_best_path(unitig_sequences, overlaps)
assembled_sequence = assemble_sequence(best_path, unitig_sequences, overlaps)

with open(output_file, 'w') as f:
    f.write(f">{''.join(best_path)}\n")
    f.write(assembled_sequence + "\n")
