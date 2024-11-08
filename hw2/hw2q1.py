# Read a FASTQ file and summarize info
# Two command-line arguments: FASTQ file and output summary file
# Input: FASTQ file
# Output: Summary 5 integers separated by spaces, L to R (
# 1. read w lowest total quality score, 
# 2. read w highest total quality score
# more than one read is tied choose the one the occurs first
# 3. total number base quality values < 10
# 4. total number base quality values >= 30
# 5. total number of chars not A, C, G, T

# Notebook code used: 
# - https://nbviewer.org/github/BenLangmead/comp-genomics-class/blob/master/notebooks/FASTQ.ipynb

from io import StringIO
import sys

# Function from notebook
def phred33_to_q(qual):
  """ Turn Phred+33 ASCII-encoded quality into Phred-scaled integer """
  return ord(qual)-33

# Function taken from notebook, but made modifications to it for the specific conditions
def parse_fastq_and_find_summary(fh):
    """ Parse reads from a FASTQ filehandle.  For each read, we
        return a name, nucleotide-string, quality-string triple. """
    
    read_index, reads = 1, []
    tot_less_than_10, tot_ge_30 = 0, 0
    num_not_chars = 0
    lowest_score, highest_score = float('inf'), float('-inf')
    lowest_index, highest_index = -1, -1

    while True:
        first_line = fh.readline()
        if len(first_line) == 0:
            break  # end of file
        name = first_line[1:].rstrip()
        seq = fh.readline().rstrip()
        # Condition 5: how many characters that are not in ACGT
        for char in seq:
            if char not in ('ACGT'):
                num_not_chars += 1
        fh.readline()  # ignore line starting with +
        qual = fh.readline().rstrip()
        quality_scores = [phred33_to_q(char) for char in qual]
        total_quality_score = sum(quality_scores)/len(quality_scores)
        reads.append((read_index, total_quality_score))
        for read_index, score in reads:
            # Condition 1: Finding the index of the lowest quality score
            if score < lowest_score:
                lowest_score = score
                lowest_index = read_index
            # Condition 2: Finding the index of the highest quality score
            if score > highest_score:
                highest_score = score
                highest_index = read_index
        for score in quality_scores:
            # Condition 3: number of scores less than 10
            if score < 10:
                tot_less_than_10 += 1
            # Condition 4: number of scores greater than or equal to 30
            elif score >= 30:
                tot_ge_30 += 1
        
        # Increment the read index, because next read will be looked at now
        read_index += 1
    return lowest_index, highest_index, tot_less_than_10, tot_ge_30, num_not_chars

# Main Program
if len(sys.argv) != 3:
    print("Usage: python3 hw2q1.py <input_file> <output_file>")
    sys.exit(1)

# Assign the input and output filenames from command-line arguments
fastq_file = sys.argv[1] # First argument: input filename
summary_file = sys.argv[2] # Second argument: output filename

with open(fastq_file, 'r') as in_file:
    lowest_index, highest_index, tot_less_than_10, tot_ge_30, num_not_chars = parse_fastq_and_find_summary(in_file)
with open(summary_file, 'w') as out_file:
    out_file.write(f"{lowest_index} {highest_index} {tot_less_than_10} {tot_ge_30} {num_not_chars}")