from io import StringIO
import sys

# Function to convert Phred+33 ASCII-encoded quality into Phred-scaled integer
def phred33_to_q(qual):
    return ord(qual) - 33

# Function to parse FASTQ file and find summary statistics
def parse_fastq_and_find_summary(fh):
    """ Parse reads from a FASTQ filehandle. For each read, we return a name, nucleotide-string, quality-string triple. """
    
    read_index, reads = 1, []
    tot_less_than_10, tot_ge_30 = 0, 0
    num_not_chars = 0
    lowest_score, highest_score = float('inf'), float('-inf')
    lowest_index, highest_index = -1, -1

    # Dictionary to track low-quality base count per cycle
    cycle_low_quality_counts = {}

    while True:
        first_line = fh.readline()
        if len(first_line) == 0:
            break  # end of file
        name = first_line[1:].rstrip()
        seq = fh.readline().rstrip()
        
        # Count characters that are not A, C, G, T
        for char in seq:
            if char not in ('A', 'C', 'G', 'T'):
                num_not_chars += 1
                
        fh.readline()  # ignore line starting with +
        qual = fh.readline().rstrip()
        quality_scores = [phred33_to_q(char) for char in qual]

        # Track low-quality base count per cycle
        for cycle_index, score in enumerate(quality_scores):
            if score < 10:
                if cycle_index not in cycle_low_quality_counts:
                    cycle_low_quality_counts[cycle_index] = 0
                cycle_low_quality_counts[cycle_index] += 1

        total_quality_score = sum(quality_scores)
        reads.append((read_index, total_quality_score))
        
        # Check for lowest and highest quality score reads
        for read_idx, score in reads:
            if score < lowest_score:
                lowest_score = score
                lowest_index = read_idx
            if score > highest_score:
                highest_score = score
                highest_index = read_idx
        
        # Calculate total number of low and high-quality bases
        for score in quality_scores:
            if score < 10:
                tot_less_than_10 += 1
            elif score >= 30:
                tot_ge_30 += 1
        
        read_index += 1
    
    # Compare each cycle with its neighbors to find the most anomalous cycle
    most_anomalous_cycle = None
    max_deviation = float('-inf')
    cycle_keys = sorted(cycle_low_quality_counts.keys())

    for i in range(1, len(cycle_keys) - 1):
        prev_cycle = cycle_low_quality_counts[cycle_keys[i - 1]]
        curr_cycle = cycle_low_quality_counts[cycle_keys[i]]
        next_cycle = cycle_low_quality_counts[cycle_keys[i + 1]]
        # Calculate average of previous and next cycles
        neighbor_avg = (prev_cycle + next_cycle) / 2

        # Calculate the deviation from neighbor average
        deviation = curr_cycle - neighbor_avg
        if deviation > max_deviation:
            max_deviation = deviation
            most_anomalous_cycle = cycle_keys[i]

    # If there's no clear anomaly, use the cycle with the highest low-quality bases
    if most_anomalous_cycle is None:
        most_anomalous_cycle = max(cycle_low_quality_counts, key=cycle_low_quality_counts.get)

    most_anomalous_cycle_count = cycle_low_quality_counts[most_anomalous_cycle]
    
    # Return the summary data along with the most anomalous cycle information
    return (lowest_index, highest_index, tot_less_than_10, tot_ge_30, 
            num_not_chars, most_anomalous_cycle + 1, most_anomalous_cycle_count)

# Main Program
if len(sys.argv) != 3:
    print("Usage: python3 hw2q1.py <input_file> <output_file>")
    sys.exit(1)

# Assign the input and output filenames from command-line arguments
fastq_file = sys.argv[1] # First argument: input filename
summary_file = sys.argv[2] # Second argument: output filename

# Process the FASTQ file and gather the summary statistics
with open(fastq_file, 'r') as in_file:
    (lowest_index, highest_index, tot_less_than_10, tot_ge_30, 
     num_not_chars, most_anomalous_cycle, most_anomalous_cycle_count) = parse_fastq_and_find_summary(in_file)

# Write the summary data to the output file
with open(summary_file, 'w') as out_file:
    out_file.write(f"{lowest_index} {highest_index} {tot_less_than_10} {tot_ge_30} {num_not_chars} {most_anomalous_cycle} {most_anomalous_cycle_count}")
