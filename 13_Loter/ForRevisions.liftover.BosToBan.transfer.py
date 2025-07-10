#from liftover import get_lifter
#converter = get_lifter('hg19', 'hg38')

from liftover import ChainFile
converter = ChainFile('./bosTau9ToGCF_032452875.1.over.chain.gz', '', '', one_based=True)
#converter[chrom][pos]


# chrom = '1'
# pos = 103786442
# o = converter[chrom][pos]

with open("./N935") as f:
    bed = f.readlines()

# out = open("GCST90018947_buildGRCh38_signif_chrompos.bed", "w")
out = open("./BosToBanCoordinate.N935.pos.clean", "w")

# import sys
# sys.stdout = out

for line in bed:
    t = line.strip().split("\t")
    chrom = t[0]
    pos = int(t[2])
    o = converter[chrom][pos]
    if len(o) > 0:
        s = t[0] + '\t' + t[2] + '\t'
        out.write(s + '\t'.join(map(str,o[0])) + '\n')
