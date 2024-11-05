import sys

if len(sys.argv) != 3:
    print("Usage: python3 hw4q1.py <input_file> <output_file>")
    sys.exit(1)

# Assign the input and output filenames from command-line arguments
input_file = sys.argv[1] # First argument: input filename
output_file = sys.argv[2] # Second argument: output filename

with open(input_file, 'r') as in_file:
    t_string = in_file.readline().strip()
    p_string = in_file.readline().strip()

t_len = len(t_string)
p_len = len(p_string)

prev_row = [0] * (p_len + 1)
curr_row = [0] * (p_len + 1)

for i in range(1, t_len + 1):
    for j in range(1, p_len + 1):
        if t_string[i - 1] == p_string[j - 1]:
            curr_row[j] = prev_row[j - 1] + 1
        else:
            curr_row[j] = max(prev_row[j], curr_row[j - 1])
    prev_row, curr_row = curr_row, prev_row

result = prev_row[p_len]
result_string = str(result)

with open(output_file, 'w') as out_file:
    out_file.write(result_string)

