import sys

## ONLY FOR DIALLELIC SNPS

CONV = {
    "0|0" : [0, 0],
    "0|1" : [0, 1],
    "1|0" : [1, 0],
    "1|1" : [1, 1],
    ".|0" : [2, 0],
    ".|1" : [2, 1],
    "0|." : [0, 2],
    "1|." : [1, 2],
    ".|." : [2, 2]

}

CONV_UNPHASED = {
    "0/0" : [0, 0],
    "0/1" : [0, 1],
    "1/0" : [1, 0],
    "1/1" : [1, 1],
    "./0" : [2, 0],
    "./1" : [2, 1],
    "0/." : [0, 2],
    "1/." : [1, 2],
    "./." : [2, 2]

}

def base(gt, ref_alt, phase):

    try:
        if phase:
            alleles = CONV[gt]
        else:
            alleles = CONV_UNPHASED[gt]
        return " ".join(ref_alt[a] for a in alleles)
    except KeyError:
        print(f"problems at gt {gt}. ref: {' '.join(ref_alt[:2])}. put to missing", file=sys.stderr)
        return "0 0"

phase = int(sys.argv[1])
for line in sys.stdin:
    if line.startswith("#"):
        continue
    fields = line.rstrip().split()
    (chrom, pos, snpid, ref,alt)  = fields[:5]
    ref_alt = (ref, alt, "0")
    gts = fields[9:]
    haplo = " ".join(base(gt, ref_alt, phase) for gt in gts)

    print(chrom, snpid, 0, pos, haplo)
