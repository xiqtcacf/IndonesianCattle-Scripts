import os
import itertools as it
from pathlib import Path
import scripts.helpers
from scripts.helpers import add_ext

CURR_DIR = os.getcwd()
MAKEIN = os.path.join(CURR_DIR, "scripts", "makeinput.py")
FAKEPHASE = os.path.join(CURR_DIR, "scripts", "fake_tped_phase.py")
PLOT_COALEST = os.path.join(CURR_DIR, "scripts", "coolbeans.R")
BCFTOOLS = config["BCFTOOLS"]
PLINK = config["PLINK"]
PYTHON = config["PYTHON"]
MSMC2 = config["MSMC2"]
CHROMS = config["CHROMS"]
POPS = config["POPS"]

BED = config["GOODBED"]
FAI = config["FAI"]
MUT = config["MUT"]
GEN = config["GEN"]

if not config["PHASED"]:
    print("SOFTWARE requires phased VCF")
    sys.exit(1)
PHASED_BOOL = config["PHASED"]

OD = Path(config["OUTMAIN"])
OD_plink_chr = OD / "plink" / "{p1}.{p2}.{chrom}"
OD_ref_chr = OD / "ref" / "{chrom}"
OD_msmcin_chr = OD / "msmc_input" / "{p1}.{p2}.{chrom}"
OD_msmcout = OD / "msmc_output" / "{p1}.{p2}"
OD_msmcoutplot = OD / "msmc_output_plots" / "res"

## TAKES TWO DIPLOID SAMPLES FROM EACH POP
N = 2

if ("PLINKBED" in config.keys()) or \
   ("VCF" not in config.keys()) :
    print("\n\t -> The software can ONLY take a PHASED 'VCF' as input :)\n")
    sys.exit(1)

(p1s, p2s) = list(zip(*(it.combinations(POPS.keys(), 2))))

wildcard_constraints:
    chrom = "|".join(map(str, CHROMS)),
    p1 = "|".join(map(str, POPS.keys())),
    p2 = "|".join(map(str, POPS.keys())),
    first_second = "1|2"

rule all:
    input:
        add_ext(OD_msmcoutplot, "coalest", "png"),
        expand(add_ext(OD_msmcout, "coalest.txt", t=str), zip, p1=p1s, p2=p2s),
        expand(add_ext(OD_msmcout, "im_input.txt", t=str), zip, p1=p1s, p2=p2s)

rule samples_chrom_to_plink2:
    output:
        temp(add_ext(OD_plink_chr, "tped")),
        add_ext(OD_plink_chr, "samples.kept")
    log:
        add_ext(OD_plink_chr, "log"),
    threads: 3
    run:

        baseout = output[0][:-5]
        samplelist_file1 = POPS[wildcards.p1]
        samplelist_file2 = POPS[wildcards.p2]
        with open(samplelist_file1, 'r') as fh:
            samplelist1 = [x.rstrip() for x in fh]
        with open(samplelist_file2, 'r') as fh:
            samplelist2 = [x.rstrip() for x in fh]
        (tped, samplekeep) = output

        samplelist = []
        with open(samplelist_file1, 'r') as fh:
            for idx, line in enumerate(fh):
                if idx==N:
                    break
                samplelist.append(line.rstrip())
        with open(samplelist_file2, 'r') as fh:
            for idx, line in enumerate(fh):
                if idx==N:
                    break
                samplelist.append(line.rstrip())

        with open(samplekeep, 'w') as fh:
            for x in samplelist:
                print(x, file=fh)

        DATA = config["VCF"]
        shell("{BCFTOOLS} view --threads {threads} -r {wildcards.chrom} -S {samplekeep} {DATA} | "
              "{BCFTOOLS} view -H -i 'INFO/AC>0 && INFO/AC<INFO/AN' | "
              "python3 {FAKEPHASE} {PHASED_BOOL} > {tped} 2> {log}")

rule chrom_size:
    input:
        FAI
    output:
        add_ext(OD_ref_chr, "size")
    params:
        awk = "{print $2}"
    shell:
        "awk '$1==\"{wildcards.chrom}\" {params.awk}' {input} > {output}" #### have to change based on chromosome name

rule bed_chrom:
    input:
        BED
    output:
        add_ext(OD_ref_chr, "bed")
    shell:
        "awk '$1==\"{wildcards.chrom}\"' {input} > {output}" #### have to change based on chromosome name
 
rule msmc_input:
    input:
        add_ext(OD_ref_chr, "size"),
        add_ext(OD_ref_chr, "bed"),
        add_ext(OD_plink_chr, "tped"),
    output:
        add_ext(OD_msmcin_chr, "input")
    shell:
        "{PYTHON} {MAKEIN} {input} > {output}"

rule msmc_run_between:
    input:
        expand(add_ext(OD_msmcin_chr, "input", t=str), chrom=CHROMS, allow_missing=True)
    output:
        multiext(add_ext(OD_msmcout, "merge", t=str), ".final.txt", ".loop.txt")
    log:
        add_ext(OD_msmcout, "merge", "log")
    threads: 10
    run:
        baseout = output[0][:-10]
        idxs = []
        for x in range(N*2):
            for y in range(N*2, N*2*2):
                idxs.append(f"{x}-{y}")
        idxs = ",".join(idxs)

        shell("{MSMC2} -t {threads} -o {baseout} -I {idxs} {input}")

## run within as well
rule msmc_run_within:
    input:
        expand(add_ext(OD_msmcin_chr, "input", t=str), chrom=CHROMS, allow_missing=True)
    output:
        multiext(add_ext(OD_msmcout, "{first_second}", t=str), ".final.txt", ".loop.txt")
    log:
        add_ext(OD_msmcout, "{first_second}", "log")
    threads: 10
    run:
        baseout = output[0][:-10]
        if int(wildcards.first_second) == 1:
            idxs = list(map(str, range(N*2)))
        elif int(wildcards.first_second) == 2:
            idxs = list(map(str, range(N*2, N*2*2)))
        idxs = ",".join(idxs)
        shell("{MSMC2} -t {threads} -o {baseout} -I {idxs} {input}")

## https://github.com/stschiff/msmc-tools/blob/master/msmc-tutorial/guide.md
## coalest and coalest2 now returns the same results
rule coalest:
    input:
        expand(add_ext(OD_msmcout, "{first_second}", "final.txt", t=str), first_second = [1,2], allow_missing=True),
        add_ext(OD_msmcout, "merge", "final.txt", t=str),
    output:
        add_ext(OD_msmcout, "coalest.txt")
    run:
        with open(output[0], 'w') as fh:
            (x,y) = helpers.crossCoalPlotCombined(*input, mu=MUT, gen=GEN)
            for n in range(len(x)):
                print(wildcards.p1, wildcards.p2, round(x[n], 3), round(y[n],3), file=fh)

# rule coalest2:
#     input:
#         expand(add_ext(OD_msmcout, "{first_second}", "final.txt", t=str), first_second = [1,2], allow_missing=True),
#         add_ext(OD_msmcout, "merge", "final.txt", t=str),
#     output:
#         add_ext(OD_msmcout, "coalest2.txt"),
#         add_ext(OD_msmcout, "coalest2.years_rel.txt"),
#     run:
#         with open(output[0], 'w') as fh:
#             helpers.combinecoal(fh, *input)

#         with open(output[0], 'r') as fh, open(output[1], 'w') as fhout:
#             header = next(fh)
#             for line in fh:
#                 fields = line.rstrip().split()
#                 lftime = float(fields[1])
#                 lambda00 = float(fields[3])
#                 lambda01 = float(fields[4])
#                 lambda11 = float(fields[5])

#                 years = lftime / MUT * GEN
#                 relcoal = 2 * lambda01 / (lambda00 + lambda11)

#                 print(wildcards.p1, wildcards.p2, years, relcoal, file=fhout)

rule plot_coalest:
    input:
        expand(add_ext(OD_msmcout, "coalest.txt", t=str), zip, p1=p1s, p2=p2s),
    output:
        add_ext(OD_msmcoutplot, "coalest", "files"),
        add_ext(OD_msmcoutplot, "coalest", "png"),
    run:
        with open(output[0], 'w') as fh:
            print("\n".join(input), file=fh)

        shell("Rscript {PLOT_COALEST} {output[0]} {output[1]}")

rule im_input:
    input:
        expand(add_ext(OD_msmcout, "{first_second}", "final.txt", t=str), first_second = [1,2], allow_missing=True),
        add_ext(OD_msmcout, "merge", "final.txt", t=str),
    output:
        add_ext(OD_msmcout, "im_input.txt"),
    run:
        with open(output[0], 'w') as fh:
            helpers.combinecoal(fh, *input)


# rule plot_im:
#     input:
#         add_ext(OD_msmcout, "im_input.txt"),
#     output:
#         add_ext(OD_msmcoutplot, "b1_1e-08.b2_1e-06.MSMC_IM.fittingdetails.xlog.pdf")
#     shell:
#         pass
