# Transcription
# Input: takes a DNA string. Output is a RNA version where T replaced by U
# Dash between every three characters
# INPUT in input.txt file: AcTGAAC
# Expected OUTPUT in output.txt file: AUG-AAC

import sys

if len(sys.argv) != 3:
    print("Usage: python3 hw1q3.py <input_file> <output_file>")
    sys.exit(1)

# Assign the input and output filenames from command-line arguments
input_file = sys.argv[1] # First argument: input filename
output_file = sys.argv[2] # Second argument: output filename
result = []
# Open both input and output files - input file for reading and output file for writing
with open(input_file, 'r') as in_file, open(output_file, 'w') as out_file:
    # Iterate through each line in the input file
    for line in in_file:
        # Iterate through each character in the line
        for c in line:
            # Check if the character is one of ’A’, ’C’, ’G’, or ’T’ and append the complement
            if c in "ACG":
                result.append(c)
            if c == "T":
                result.append("U")
    # Reverse the result string

    result_str = ''.join(result)
    split = [result_str[i:i+3] for i in range(0, len(result_str), 3)]
    split = '-'.join(split)
    out_file.write(split)