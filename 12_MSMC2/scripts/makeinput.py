import sys
import numpy as np
size_chrom_file, bedfile, tped = sys.argv[1:]

with open(size_chrom_file, 'r') as fh:
    size_chrom = int(fh.read().rstrip())

keep = np.zeros(int(size_chrom), dtype=np.int32)

with open(bedfile, 'r') as fh:
    for line in fh:
        _, start, stop = line.rstrip().split()
        keep[int(start):int(stop)] = 1

with open(tped, 'r') as fh:
    start=0
    for line in fh:
        fields = line.rstrip().split()
        chrom, _, _, pos = fields[:4]
        pos = int(pos)
        nhom = np.sum(keep[start:pos])
        data = "".join(fields[4:])
        data = data.replace("0", "?")
        print(chrom, pos, nhom, data)
        start = pos
