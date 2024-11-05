import sys
from collections import defaultdict

def load_unitigs(unitig_file):
    unitigs = []
    with open(unitig_file, 'r') as f:
        current_unitig = []
        for line in f:
            line = line.strip()
            if line.isdigit() or " " not in line:  # New unitig start
                if current_unitig:
                    unitigs.append(current_unitig)
                current_unitig = [line]
            else:
                overlap_len, read = line.split()
                current_unitig.append((int(overlap_len), read))
        if current_unitig:
            unitigs.append(current_unitig)
    return unitigs

def compose_sequences(unitigs, read_sequences):
    sequences = {}
    for unitig in unitigs:
        sequence = read_sequences[unitig[0]]
        for overlap_len, read in unitig[1:]:
            sequence += read_sequences[read][overlap_len:]
        sequences[unitig[0]] = sequence
    return sequences

def find_longest_overlap(seq_a, seq_b, min_overlap=1):
    max_overlap_len = min(len(seq_a), len(seq_b))
    for overlap_len in range(max_overlap_len, min_overlap - 1, -1):
        if seq_a[-overlap_len:] == seq_b[:overlap_len]:
            return overlap_len
    return 0

def assemble_final_path(sequences):
    path = []
    used = set()
    
    # Start with the longest sequence as the initial path
    sorted_unitigs = sorted(sequences.keys(), key=lambda x: -len(sequences[x]))
    path.append(sequences[sorted_unitigs[0]])
    used.add(sorted_unitigs[0])

    while len(used) < len(sequences):
        best_overlap_len, best_seq = 0, None
        for unitig in sequences:
            if unitig in used:
                continue
            for i in range(len(path)):
                overlap_len = find_longest_overlap(path[i], sequences[unitig])
                if overlap_len > best_overlap_len:
                    best_overlap_len, best_seq = overlap_len, (i, unitig)
        
        if best_seq:
            i, unitig = best_seq
            path[i] += sequences[unitig][best_overlap_len:]
            used.add(unitig)
        else:
            # No overlap found, start a new path entry
            for unitig in sequences:
                if unitig not in used:
                    path.append(sequences[unitig])
                    used.add(unitig)
                    break

    return ''.join(path)

# Main code block
if len(sys.argv) != 4:
    print("Usage: python3 assembler.py <unitigs_file> <reads_file> <final_output>")
    sys.exit(1)

unitigs_file = sys.argv[1]
reads_file = sys.argv[2]
output_file = sys.argv[3]

# Load the reads into a dictionary for quick lookup
read_sequences = {}
with open(reads_file, 'r') as f:
    for line in f:
        parts = line.strip().split()
        if len(parts) == 2:
            read_id, sequence = parts
            read_sequences[read_id] = sequence

# Load the unitigs and compose sequences
unitigs = load_unitigs(unitigs_file)
composed_sequences = compose_sequences(unitigs, read_sequences)

# Find the best path and write it to the output file
final_path = assemble_final_path(composed_sequences)
with open(output_file, 'w') as f:
    f.write(final_path)
