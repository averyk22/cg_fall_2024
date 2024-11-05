# Used code from class notebook: https://colab.research.google.com/github/BenLangmead/comp-genomics-class/blob/master/notebooks/CG_BWT_SimpleBuild.ipynb#scrollTo=kzOHVO7ujlrX
import sys

# Function taken from notebook
def suffixArray(s):
    satups = sorted([(s[i:], i) for i in range(len(s))])
    return map(lambda x: x[1], satups)

# Function taken from notebook
def bwtViaSa(t):
    # Given T, returns BWT(T) by way of the suffix array
    bw = []
    for si in suffixArray(t):
        if si == 0:
            bw.append('$')
        else:
            bw.append(t[si-1])
    return ''.join(bw) # return string-ized version of list bw

# Count number of same-character runs
def btw_run_length(s):
    counter = 1  
    for i in range(1, len(s)):
        if s[i] != s[i - 1]:
            counter += 1
    return counter

# Find the length of the longest run in BWT
def longest_run_len(s):
    max_len = 1
    curr_len = 1
    for i in range(1, len(s)):
        if s[i] == s[i - 1]:
            curr_len += 1
        else:
            max_len = max(max_len, curr_len)
            curr_len = 1
    max_len = max(max_len, curr_len) 
    return max_len


if len(sys.argv) != 3:
    print("Usage: python3 hw3q4.py <input_file> <output_file>")
    sys.exit(1)

# Assign the input and output filenames from command-line arguments
input_file = sys.argv[1] # First argument: input filename
output_file = sys.argv[2] # Second argument: output filename

string = ''

# Read in the input
with open(input_file, 'r') as in_file:
    string = in_file.read().strip()

# create suffix array and bwt
suffix_array = suffixArray(string)
bwt = bwtViaSa(string)

# values we are looking to output
compression_n = len(string)
compression_r = btw_run_length(bwt)
longest_run_len = longest_run_len(bwt)

# Write to output file
with open(output_file, 'w') as out_file:
    out_file.write(f"{bwt}\n")
    out_file.write(f"{compression_n}:{compression_r}\n")
    out_file.write(f"{longest_run_len}\n")