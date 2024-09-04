# Repeat detection
# Input: DNA string
# Output: Prints the indices of the 6 letter long substring that occurs most frequently in the DNA string
# Print the base-0 index values separated by commas. To resolve ties select the string that comes first in lexico 

import sys
from collections import defaultdict

if len(sys.argv) != 3:
    print("Usage: python3 script.py <input_file> <output_file>")
    sys.exit(1)

input_file = sys.argv[1] # First argument: input filename
output_file = sys.argv[2] # Second argument: output filename

with open(input_file, 'r') as in_file:
    dna_string = in_file.readline().strip()

# Dictionary to store counts and all the indices that have substrings of that count
count_dict = defaultdict(lambda: [0, []])
for i in range(len(dna_string) - 5):
    substring = dna_string[i:i+6]
    count_dict[substring][0] += 1
    count_dict[substring][1].append(i)

most_freq_ss = None
most_freq_ind = None

# Now need to locate most frequent
for substring, (count, index) in count_dict.items():
    if most_freq_ss is None or count > count_dict[most_freq_ss][0] or (count == count_dict[most_freq_ss][0] and substring < most_freq_ss):
        most_freq_ss = substring
        most_freq_ind = index

with open(output_file, 'w') as out_file:
    if most_freq_ind is not None: 
        out_file.write(','.join(map(str, most_freq_ind)))
    else:
        out_file.write('There were no 6-letter substrings found in the DNA sequence.')

