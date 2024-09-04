# Palindromes 
# takes an integer n and outputs all distinct DNA strings that are equal to their own reverse complements with length ≤ n. 
# Checks if a string is equal with their reverse complement
def is_palindrome(str):
    complement_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    reverse_comp = []
    # Generate the reverse complement
    for c in str: 
        if c in complement_dict:
            reverse_comp.append(complement_dict[c])
    # Reverse the result string
    reverse_comp.reverse()
    reverse_comp_str = ''.join(reverse_comp)
    
    # Check if the original string is equal to its reverse complement
    return str == reverse_comp_str

# Generate all the possible DNA strings
def generate_valid_dna_strings(n):
    dna_bases = ['A', 'T', 'G', 'C']
    result = set()
    def generate(curr_str, len):
        if len > n:
            return
        if len > 0 and is_palindrome(curr_str):
            result.add(curr_str)
        for base in dna_bases:
            generate(curr_str + base, len + 1)
    result = set()
    generate('', 0)
    return result

import sys
if len(sys.argv) != 3:
    print("Usage: python3 hw1q4.py <input_file> <output_file>")
    sys.exit(1)

input_file = sys.argv[1] # First argument: input filename
output_file = sys.argv[2] # Second argument: output filename

with open(input_file, 'r') as in_file, open(output_file, 'w') as out_file:
    n = int(in_file.read())
    palindromes = sorted(generate_valid_dna_strings(n))
    for palindrome in palindromes:
        out_file.write(palindrome + '\n')
   
