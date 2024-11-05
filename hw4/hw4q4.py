import sys
 
def parse_fastq(fn):
    """ Parse FASTA into a dictionary """
    with open(fn) as myFile:
        lines = myFile.readlines()
    nreads = int(len(lines) / 4)
    fa = {}
    for i in range(nreads):
        nm = lines[i*4][1:].rstrip()
        seq = lines[(i*4)+1].rstrip().upper()
        fa[nm] = seq
    return fa
 
def write_solution(unitigID, unitigSequence, n, out_fh, per_line=60):
    offset = 0
    out_fh.write(f">{unitigID} {n}\n")
    while offset < len(unitigSequence):
        line = unitigSequence[offset:offset + per_line]
        offset += per_line
        out_fh.write(line + "\n")
 
sequences = parse_fastq(sys.argv[1])
 
unitigs = {}
cur_unitig = None
from collections import defaultdict
unitigs_num = defaultdict(int)
for line in open(sys.argv[2]):
    fields = line.split()
    if len(fields) == 1:
        cur_unitig = fields[0]
        unitigs[cur_unitig] = sequences[fields[0]]
    else:
        unitigs[cur_unitig] += sequences[fields[1]][int(fields[0]):]
    unitigs_num[cur_unitig] += 1
 
out_fh = open(sys.argv[3], 'w')
for u, v in unitigs.items():
    write_solution(u, v, unitigs_num[u], out_fh)
out_fh.close()