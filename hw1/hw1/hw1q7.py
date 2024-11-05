# Hamming Distance
# Input: two DNA strings of the same length
# Output: reports how many positions differ between the two

import sys

if len(sys.argv) != 3:
    print("Usage: python3 script.py <input_file> <output_file>")
    sys.exit(1)

input_file = sys.argv[1] # First argument: input filename
output_file = sys.argv[2] # Second argument: output filename

with open(input_file, 'r') as in_file:
    # Read the first line, remove whitespace from both ends
    a = in_file.readline().strip()
    # Read the second line, remove whitespace from both ends
    b = in_file.readline().strip()

count = 0
for i in range(0, len(a)):
    if a[i] != b[i]:
        count += 1

with open(output_file, 'w') as out_file:
    out_file.write(str(count))