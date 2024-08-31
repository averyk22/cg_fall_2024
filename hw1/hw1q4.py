# Palindromes 
# takes an integer n and outputs all distinct DNA strings that are equal to their own reverse complements with length ≤ n. 

import sys
if len(sys.argv) != 3:
    print("Usage: python3 hw1q4.py <input_file> <output_file>")
    sys.exit(1)

input_file = sys.argv[1] # First argument: input filename
output_file = sys.argv[2] # Second argument: output filename
n = 0
with open(input_file, 'r') as in_file, open(output_file, 'w') as out_file:
    for line in input_file:
        for c in line:
            n = c

   
# Checks if a string is equal with their reverse complement
def isPalindrome(str):
    complement_dict = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'} 
    for char in str:
        if c in complement_dict:
                result.append(complement_dict[c])
    result.reverse()
    result_str = ''.join(result)
    if result_str == str:
        return True
    return False
