# Reverse Complement
# Input: should read from input file, Output: written on a single line to output file.
# Any input that is Not A, C, G, T should be ignored. 
import sys

if len(sys.argv) != 3:
    print("Usage: python3 hw1q2.py <input_file> <output_file>")
    sys.exit(1)
# Assign the input and output filenames from command-line arguments
input_file = sys.argv[1] # First argument: input filename
output_file = sys.argv[2] # Second argument: output filename
complement_dict = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'} 
result = []
# Open both input and output files - input file for reading and output file for writing
with open(input_file, 'r') as in_file, open(output_file, 'w') as out_file:
    # Iterate through each line in the input file
    for line in in_file:
        # Iterate through each character in the line
        for c in line:
            # Check if the character is one of ’A’, ’C’, ’G’, or ’T’ and append the complement
            if c in complement_dict:
                result.append(complement_dict[c])
    # Reverse the result string
    result.reverse()
    result_str = ''.join(result)
    out_file.write(result_str)