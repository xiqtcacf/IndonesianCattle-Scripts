import itertools as it

BCFTOOLS="/isdata/hellergrp/nzg134/programs/bcftools-1.10.2/bcftools"
PLINK="/isdata/hellergrp/nzg134/programs/bin/plink"
PLINK2="/isdata/hellergrp/nzg134/programs/bin/plink2"
SAMTOOLS="/isdata/hellergrp/nzg134/programs/samtools-1.12/samtools"
FASTQC="/isdata/hellergrp/nzg134/programs/FastQC/fastqc"
MULTIQC="/isdata/hellergrp/nzg134/programs/anaconda3/envs/snakenv/bin/multiqc"
SATC="/isdata/hellergrp/nzg134/programs/SATC/satc.R"
RSCRIPT="/isdata/hellergrp/nzg134/programs/anaconda3/bin/Rscript"
ANGSDDIR="/isdata/hellergrp/rheller/software_scratch/angsdv0.929"
ANGSD=os.path.join(ANGSDDIR, "angsd")
EST_ERROR=os.path.join(ANGSDDIR,"R/estError.R")
REALSFS=os.path.join(ANGSDDIR, "misc", "realSFS")
SCRIPT_DIR=config["script_dir"]
OUTMAIN=config["outmain"]


MINQ=config["minQ"]
MINMAPQ=config["minmapQ"]

BEDFILE_SPLIT=20

whatrefs = ["close", "distant"]

wildcard_constraints:
    sample = "|".join(config["samples"]),
    s1 = "|".join(config["samples"]),
    s2 = "|".join(config["samples"]),
    whatref = "|".join(whatrefs),

rule all:
    shell:
        "echo targets: RUN_json RUN_satc RUN_multiqc_pre RUN_multiqc_post RUN_errorrate RUN_depths RUN_hetero RUN_scaf RUN_related RUN_pca RUN_ALL"


rule RUN_json:
    input:
        ## jsons
        os.path.join(OUTMAIN, "plots", "matches_pct.passed.png"),
        os.path.join(OUTMAIN, "plots", "insert_sizes.passed.png"),
        os.path.join(OUTMAIN, "plots", "mapping_quality.passed.png"),
        os.path.join(OUTMAIN, "plots", "query_lengths.passed.png"),
        os.path.join(OUTMAIN, "plots", "stats_passed_failed.png"),
        os.path.join(OUTMAIN, "plots", "stats_exclusion_freq.pdf"),
        os.path.join(OUTMAIN, "plots", "stats_exclusion_count.pdf"),

rule RUN_satc:
    input:
        os.path.join(OUTMAIN, "SATC", config["species"] + "_sampleSex.tsv"),

rule RUN_multiqc_pre:
    input:
        pre=os.path.join(OUTMAIN, "multiqc", "premap.html"),
        pre2=os.path.join(OUTMAIN, "multiqc", "premap_perlane.html"),

rule RUN_multiqc_post:
    input:
        postclose=os.path.join(OUTMAIN, "multiqc", "postmap_closeref.html"),
        postdist=os.path.join(OUTMAIN, "multiqc", "postmap_distantref.html"),
        postclose_minmapq=os.path.join(OUTMAIN, "multiqc", "postmap_closeref_minmapQ.html"),
        postdist_minmapq=os.path.join(OUTMAIN, "multiqc", "postmap_distantref_minmapQ.html"),
        stats=os.path.join(OUTMAIN, "multiqc", "stats.html"),

rule RUN_errorrate:
    input:
        os.path.join(OUTMAIN, "errorrate", "angsdErrorEst.txt"),

rule RUN_depths:
    input:
        os.path.join(OUTMAIN, "plots", "depths_zoom.png"),
        os.path.join(OUTMAIN, "plots", "depths_global.png"),



rule RUN_hetero:
    input:
        os.path.join(OUTMAIN, "hetero", "collect_het.png")

rule RUN_scaf:
    input:
        os.path.join(OUTMAIN, "plots", "cumsum_scafs.png"),

rule RUN_related:
    input:
        ## expand(os.path.join(OUTMAIN, "sfs_2d", "{whatref}", "collected.txt"), whatref=whatrefs)
        expand(os.path.join(OUTMAIN, "sfs_2d", "{whatref}_king.png"), whatref="close")

rule RUN_pca:
    input:
        expand(os.path.join(OUTMAIN, "pca", "{whatref}_pop.png"), whatref=whatrefs)

rule RUN_ALL:
    input:
        rules.RUN_json.input,
        rules.RUN_satc.input,
        rules.RUN_multiqc_pre.input,
        rules.RUN_multiqc_post.input,
        rules.RUN_errorrate.input,
        rules.RUN_depths.input,
        rules.RUN_hetero.input,
        rules.RUN_scaf.input,
        rules.RUN_related.input,
        rules.RUN_pca.input,

rule fastqc_premap:
    ## https://stackoverflow.com/a/50882104
    input:
        fq1 = lambda wildcards: config["fastq"][1][wildcards.sample],
        fq2 = lambda wildcards: config["fastq"][2][wildcards.sample]
    output:
        directory(os.path.join(OUTMAIN, "fastqc", "premap", "{sample}"))
    threads: 1
    log:
        os.path.join(OUTMAIN, "fastqc", "premap", "{sample}.log")
    run:
        shell("mkdir -p {output}")
        shell("zcat {input.fq1} | {FASTQC} -f fastq -o {output} -t 1 -q stdin:{wildcards.sample}_1 &> {log}")
        shell("zcat {input.fq2} | {FASTQC} -f fastq -o {output} -t 1 -q stdin:{wildcards.sample}_2 &> {log}")


rule fastqc_premap_perlane:
    ## https://stackoverflow.com/a/50882104
    input:
        fq1 = lambda wildcards: config["fastq"][1][wildcards.sample],
        fq2 = lambda wildcards: config["fastq"][2][wildcards.sample]
    output:
        directory(os.path.join(OUTMAIN, "fastqc", "premap_perlane", "{sample}"))
    threads: 1
    log:
        os.path.join(OUTMAIN, "fastqc", "premap_perlane", "{sample}.log")
    run:
        shell("mkdir -p {output}")
        for f in input.fq1:
            bf = os.path.basename(f).replace(".fq.gz", "")
            shell("zcat {f} | {FASTQC} -f fastq -o {output} -t 1 -q stdin:{wildcards.sample}_{bf} &> {log}")
        for f in input.fq2:
            bf = os.path.basename(f).replace(".fq.gz", "")
            shell("zcat {f} | {FASTQC} -f fastq -o {output} -t 1 -q stdin:{wildcards.sample}_{bf} &> {log}")

# rule fastqc_postmap_close:
#     ## https://stackoverflow.com/a/50882104
#     input:
#         bam = lambda wildcards: config["bams"]["close"][wildcards.sample]
#     output:
#         directory(os.path.join(OUTMAIN, "fastqc", "postmap_closeref", "{sample}"))
#     threads: 1
#     log:
#         os.path.join(OUTMAIN, "fastqc", "postmap_closeref", "{sample}.log")
#     run:
#         shell("mkdir -p {output}")
#         shell("{FASTQC} -f bam -o {output} -t 1 {input.bam} &> {log}")

# rule fastqc_postmap_distant:
#     ## https://stackoverflow.com/a/50882104
#     input:
#         bam = lambda wildcards: config["bams"]["distant"][wildcards.sample]
#     output:
#         directory(os.path.join(OUTMAIN, "fastqc", "postmap_distantref", "{sample}"))
#     threads: 1
#     log:
#         os.path.join(OUTMAIN, "fastqc", "postmap_distantref", "{sample}.log")
#     run:
#         shell("mkdir -p {output}")
#         shell("{FASTQC} -f bam -o {output} -t 1 {input.bam} &> {log}")

rule fastqc_postmap:
    ## https://stackoverflow.com/a/50882104
    input:
        bam = lambda wildcards: config["bams"][wildcards.whatref][wildcards.sample]
    output:
        directory(os.path.join(OUTMAIN, "fastqc", "postmap_{whatref}ref", "{sample}"))
    threads: 2
    log:
        os.path.join(OUTMAIN, "fastqc", "postmap_{whatref}ref", "{sample}.log")
    run:
        shell("mkdir -p {output}")
        shell("{FASTQC} -f bam -o {output} -t {threads} {input.bam} &> {log}")

rule fastqc_postmap_minmapQ:
    ## https://stackoverflow.com/a/50882104
    input:
        bam = lambda wildcards: config["bams"][wildcards.whatref][wildcards.sample]
    output:
        directory(os.path.join(OUTMAIN, "fastqc", "postmap_{whatref}ref_minmapQ", "{sample}"))
    threads: 2
    log:
        os.path.join(OUTMAIN, "fastqc", "postmap_{whatref}ref_minmapQ", "{sample}.log")
    run:
        shell("mkdir -p {output}")
        f = input.bam
        bf = os.path.basename(f)[:-4]
        shell("{SAMTOOLS} view -h -q {MINMAPQ} {f} | {SAMTOOLS} fastq /dev/stdin | {FASTQC} -t {threads}  -o {output} -q stdin:{bf}")

rule multiqc_pre:
    input:
        pre=expand(os.path.join(OUTMAIN, "fastqc", "premap", "{sample}"), sample=config["samples"]),
        prelane=expand(os.path.join(OUTMAIN, "fastqc", "premap_perlane", "{sample}"), sample=config["samples"]),
    output:
        pre=os.path.join(OUTMAIN, "multiqc", "premap.html"),
        prelane=os.path.join(OUTMAIN, "multiqc", "premap_perlane.html"),
    run:
        predir, prename = os.path.split(output.prelane)
        prename = prename[:-5]
        shell("{MULTIQC} --interactive -f -o {predir} -n {prename} {input.prelane}")

        predir, prename = os.path.split(output.pre)
        prename = prename[:-5]
        shell("{MULTIQC} --interactive -f -o {predir} -n {prename} {input.pre}")


rule multiqc_post:
    input:
        postclose=expand(os.path.join(OUTMAIN, "fastqc", "postmap_closeref", "{sample}"), sample=config["samples"]),
        postdist=expand(os.path.join(OUTMAIN, "fastqc", "postmap_distantref", "{sample}"), sample=config["samples"]),
        stats=config["stats"]["close"].values(),
        stats2=config["stats"]["distant"].values(),
    output:
        postclose=os.path.join(OUTMAIN, "multiqc", "postmap_closeref.html"),
        postdist=os.path.join(OUTMAIN, "multiqc", "postmap_distantref.html"),
        stats=os.path.join(OUTMAIN, "multiqc", "stats.html"),
    run:
        postdir, postname = os.path.split(output.postclose)
        postname = postname[:-5]
        shell("{MULTIQC} -f -o {postdir} -n {postname} {input.postclose}")

        postdir, postname = os.path.split(output.postdist)
        postname = postname[:-5]
        shell("{MULTIQC} -f -o {postdir} -n {postname} {input.postdist}")

        statsdir, statsname = os.path.split(output.stats)
        statsname = statsname[:-5]
        shell("{MULTIQC} -f -o {statsdir} -n {statsname} {input.stats} {input.stats2}")

rule multiqc_post_minmapQ:
    input:
        postclose=expand(os.path.join(OUTMAIN, "fastqc", "postmap_closeref_minmapQ", "{sample}"), sample=config["samples"]),
        postdist=expand(os.path.join(OUTMAIN, "fastqc", "postmap_distantref_minmapQ", "{sample}"), sample=config["samples"]),
    output:
        postclose=os.path.join(OUTMAIN, "multiqc", "postmap_closeref_minmapQ.html"),
        postdist=os.path.join(OUTMAIN, "multiqc", "postmap_distantref_minmapQ.html"),
    run:
        postdir, postname = os.path.split(output.postclose)
        postname = postname[:-5]
        shell("{MULTIQC} -f -o {postdir} -n {postname} {input.postclose}")

        postdir, postname = os.path.split(output.postdist)
        postname = postname[:-5]
        shell("{MULTIQC} -f -o {postdir} -n {postname} {input.postdist}")

rule collect_jsons:
    input:
        jsons = config["jsons"]
    output:
        directory(os.path.join(OUTMAIN, "jsons_statistics"))
    log:
        os.path.join(OUTMAIN, "jsons_statistics.log")
    shell:
        "{SCRIPT_DIR}/collect_jsons.py --jsons {input} --output_dir {output} > {log}"


rule SATC:
    input:
        idxs = config["idxstats"]["close"].values()
    output:
        out = os.path.join(OUTMAIN, "SATC", config["species"] + "_sampleSex.tsv"),
        t = os.path.join(OUTMAIN, "SATC", config["species"] + ".idx")
    log:
        os.path.join(OUTMAIN, "SATC", config["species"] + ".log"),
    run:
        with open(output.t, 'w') as fh:
            for x in input.idxs:
                print(x, file=fh)

        species = config["species"]
        outdir = os.path.dirname(output.out)
        shell("{RSCRIPT} {SATC} {species} {output.t} {outdir}  > {log} 2>&1")


rule exclude_small_scaffolds:
    input:
        c = lambda wildcards: config["refs"][wildcards.whatref] + ".fai",
    output:
        c = os.path.join(OUTMAIN, "scaffold_size", "{whatref}.txt"),
    log:
        c = os.path.join(OUTMAIN, "scaffold_size", "{whatref}.log"),
    params:
        minsize = 100000
    shell: """
    awk '{{OFS="\t"}} $2>={params.minsize}{{print $1,1,$2}}' {input.c} > {output.c}
    {ANGSD} sites index {output.c} 2> {log.c}
    """

rule plot_cumsum_scafolds:
    input:
        config["refs"]["close"] + ".fai",
        config["refs"]["distant"] + ".fai",
    output:
         os.path.join(OUTMAIN, "plots", "cumsum_scafs.png"),
    params:
        config["refs"]["close_name"],
        config["refs"]["distant_name"],
    shell:
        "{RSCRIPT} {SCRIPT_DIR}/plotting/plot_cumsum_scaf.R {input} {params} {output}"


rule subsample_genome:
    input:
        r = os.path.join(OUTMAIN, "scaffold_size", "{whatref}.txt"),
    output:
        c = os.path.join(OUTMAIN, "scaffold_size", "subsample", "{whatref}.txt"),
        s = multiext(os.path.join(OUTMAIN, "scaffold_size", "subsample", "{whatref}"), ".sites", ".sites.bin", ".sites.idx"),
        b = os.path.join(OUTMAIN, "scaffold_size", "subsample", "{whatref}.bed"),
        chroms = os.path.join(OUTMAIN, "scaffold_size", "subsample", "{whatref}.chrom"),
    log:
        l=os.path.join(OUTMAIN, "scaffold_size", "subsample", "{whatref}.log"),
    params:
        blocksize = 100000,
        blocks = 1000
    run:
        import numpy as np
        np.random.seed(666)
        with open(input.r, 'r') as fh:
            data = [line.rstrip().split() for line in fh]
        ## data = sorted(data, key=lambda x: int(x[2]), reverse=True)
        abc = []
        for x in data:
            if ":" in x[0] or "-" in x[0] or x[0] == "X":
                continue
            for start in range(1,int(x[2]), params.blocksize):
                end = start+params.blocksize-1
                if end>int(x[2]):
                    continue
                abc.append((x[0], start, end))
        a = np.random.choice(len(abc),size=params.blocks, replace=False)
        a.sort()
        keep = []
        with open(output.c, 'w') as fh, open(output.b, 'w') as fh2, open(output.s[0], 'w') as fh_sites, open(output.chroms, 'w') as fh_chrom:
            for x in a:
                if abc[x][0] not in keep:
                    print(abc[x][0], file=fh_chrom)
                    keep.append(abc[x][0])
                print(f"{abc[x][0]}\t{abc[x][1]}\t{abc[x][2]}", file=fh_sites)
                print(f"{abc[x][0]}:{abc[x][1]}-{abc[x][2]}", file=fh)
                print(f"{abc[x][0]}\t{abc[x][1]-1}\t{abc[x][2]}", file=fh2)

        shell("{ANGSD} sites index {output.s[0]} 2> {log.l}")

rule make_bamlist:
    input:
        bamsc = config["bams"]["close"].values(),
        bamsd = config["bams"]["distant"].values(),
    output:
        c1 = os.path.join(OUTMAIN, "bamlists", "close.bamlist"),
        c2 = os.path.join(OUTMAIN, "bamlists", "close.names"),
        d1 = os.path.join(OUTMAIN, "bamlists", "distant.bamlist"),
        d2 = os.path.join(OUTMAIN, "bamlists", "distant.names"),
    run:
        with open(output.c1, 'w') as fh1, open(output.c2, 'w') as fh2:
            for x in input.bamsc:
                print(x, file=fh1)
                name, ref, _ = os.path.basename(x).split(".")
                print(name, file=fh2)

        with open(output.d1, 'w') as fh1, open(output.d2, 'w') as fh2:
            for x in input.bamsd:
                print(x, file=fh1)
                name, ref, _ = os.path.basename(x).split(".")
                print(name, file=fh2)


## Call SNPS and GT
rule split_bed:
    input:
        regions = os.path.join(OUTMAIN, "scaffold_size", "subsample", "{whatref}.bed")
    output:
        expand(os.path.join(OUTMAIN, "gts", "bedsplit", "{{whatref}}_{idx}.bed"), idx=range(BEDFILE_SPLIT))
    run:
        with open(input.regions, 'r') as fh:
            beds = [line.rstrip() for line in fh]
        fs = [open(f, 'w') for f in output]
        MAX_PER_FILE = len(beds) // BEDFILE_SPLIT
        for idx, bed in enumerate(beds):
            f_idx = min(BEDFILE_SPLIT, idx//MAX_PER_FILE)
            print(bed, file=fs[f_idx])
        for f in fs:
            f.close()

rule get_GT:
    input:
        bamlist=os.path.join(OUTMAIN, "bamlists", "{whatref}.bamlist"),
        regions=os.path.join(OUTMAIN, "gts", "bedsplit", "{whatref}_{idx}.bed"),
        ref=lambda wildcards: config["refs"][wildcards.whatref],
    output:
        vcf=temp(os.path.join(OUTMAIN, "gts", "{whatref}_{idx}.bcf.gz")),
        idx=temp(os.path.join(OUTMAIN, "gts", "{whatref}_{idx}.bcf.gz.csi")),
    threads: 1
    shell: """
    {BCFTOOLS} mpileup --bam-list {input.bamlist} --annotate FORMAT/DP,FORMAT/SP,FORMAT/AD  --fasta-ref {input.ref} --min-BQ {MINQ} --min-MQ {MINMAPQ} --output-type u --per-sample-mF --regions-file {input.regions}  --threads {threads} | {BCFTOOLS} call  --multiallelic-caller --variants-only --output-type u --threads {threads} | {BCFTOOLS} view --threads {threads} -i 'STRLEN(REF)==1 & STRLEN(ALT)==1' -O b -o {output.vcf}
    {BCFTOOLS} index --threads {threads} {output.vcf}
    """

rule concat_gt:
    input:
        expand(os.path.join(OUTMAIN, "gts", "{{whatref}}_{idx}.bcf.gz"), idx=range(BEDFILE_SPLIT))
    output:
        vcf=os.path.join(OUTMAIN, "gts", "{whatref}.bcf.gz"),
        idx=os.path.join(OUTMAIN, "gts", "{whatref}.bcf.gz.csi"),
    threads: 10
    shell: """
    {BCFTOOLS} concat --threads {threads} -O b -o {output.vcf} {input}
    {BCFTOOLS} index --threads {threads} {output.vcf}
    """

rule vcf_plink:
    input:
        os.path.join(OUTMAIN, "gts", "{whatref}.bcf.gz")
    output:
        multiext(os.path.join(OUTMAIN, "gts", "plink", "{whatref}"),".bim", ".fam", ".bed")
    log:
        os.path.join(OUTMAIN, "gts", "plink", "{whatref}.log")
    params:
        b = lambda wildcards, output: output[0][:-4]
    shell: """
    {BCFTOOLS} view {input} | {PLINK} --vcf /dev/stdin  --out {params.b}  --double-id --allow-extra-chr
    """

rule plink_pca:
    input:
        multiext(os.path.join(OUTMAIN, "gts", "plink", "{whatref}"),".bim", ".fam", ".bed")
    output:
        multiext(os.path.join(OUTMAIN, "pca", "{whatref}"),".eigenvec", ".eigenval")
    log:
        os.path.join(OUTMAIN, "pca", "{whatref}.log")
    params:
        bo = lambda wildcards, output: output[0][:-9],
        bi = lambda wildcards, input: input[0][:-4]
    shell:
        "plink2 --bfile {params.bi} --out {params.bo}  --maf 0.05 --pca --allow-extra-chr"

rule plot_plink_pca:
    input:
        multiext(os.path.join(OUTMAIN, "pca", "{whatref}"),".eigenvec", ".eigenval")
    output:
        multiext(os.path.join(OUTMAIN, "pca", "{whatref}"),"_pop.png", "_nation.png", "_species.png")
    shell:
        "{RSCRIPT} {SCRIPT_DIR}/plotting/plot_pca.R {input} {output}"


## FROM GENIS IN QCSeq (the reads mapped to distant refs)
rule do_perfect_fasta:
    input:
        bam = config["perfect_bam"],
        ref=config["refs"]["distant"],
        regions=os.path.join(OUTMAIN, "scaffold_size", "subsample", "distant.txt"),
    output:
        fasta = os.path.join(OUTMAIN, "errorrate", "perfect.fa.gz")
    params:
        outprefix=lambda wildcards, output: output.fasta[:-6],
    log:
        fasta = os.path.join(OUTMAIN, "errorrate", "perfect.log")
    threads: 2
    shell: """
    {ANGSD} -rf {input.regions} -P {threads} -i {input.bam} -ref {input.ref} -doCounts 1 -doFasta 2 -out {params.outprefix} -minMapQ {MINMAPQ} -minQ {MINQ} -howoften 10000000 2> {log}
"""

rule do_error_rates:
    input:
        perfect = rules.do_perfect_fasta.output.fasta,
        anc=config["refs"]["distant"],
        bams = rules.make_bamlist.output.d1,
        regions=os.path.join(OUTMAIN, "scaffold_size", "subsample", "distant.txt"),
    output:
        ancerror = os.path.join(OUTMAIN, "errorrate", "angsderrates.ancError"),
        ancerrorchr = os.path.join(OUTMAIN, "errorrate", "angsderrates.ancErrorChr")
    params:
        outprefix=lambda wildcards, output: output.ancerror[:-9]
    threads: 2
    log:
        ancerror = os.path.join(OUTMAIN, "errorrate", "angsderrates.log")
    run:
        shell("{ANGSD} -P {threads} -rf {input.regions} -doAncError 1 -bam {input.bams} -out {params.outprefix} -minMapQ {MINMAPQ} -minQ {MINQ} -anc {input.anc} -ref {input.perfect}  -howoften 10000000  2> {log} ")

rule est_error_rates:
    input:
        ancerror = rules.do_error_rates.output.ancerror,
        indlist = rules.make_bamlist.output.d2
    output:
        e=os.path.join(OUTMAIN, "errorrate", "angsdErrorEst.txt"),
    params:
        outprefix=lambda wildcards, output: output.e[:-4]
    shell: """
    {RSCRIPT} {EST_ERROR} file={input.ancerror} out={params.outprefix} indNames={input.indlist}
    """

rule do_depth:
    input:
        bams = os.path.join(OUTMAIN, "bamlists", "{whatref}.bamlist"),  ## whatref = close, distant
        regions = os.path.join(OUTMAIN, "scaffold_size", "subsample", "{whatref}.txt"),
    output:
        os.path.join(OUTMAIN, "depths", "{whatref}.arg"),
        os.path.join(OUTMAIN, "depths", "{whatref}.pos.gz"),
        os.path.join(OUTMAIN, "depths", "{whatref}.counts.gz"),
        os.path.join(OUTMAIN, "depths", "{whatref}.depthSample"),
        os.path.join(OUTMAIN, "depths", "{whatref}.depthGlobal"),
    params:
        outbase = lambda wildcards, output: output[0][:-4]
    log:
        os.path.join(OUTMAIN, "depths", "{whatref}.log")
    threads: 2
    shell: """
    {ANGSD} -P {threads} -rf {input.regions} -howoften 1000000 -minMapQ {MINMAPQ} -minQ {MINQ}  -doCounts 1 -doDepth 1 -dumpCounts 2 -maxdepth 3000 -b {input.bams} -out {params.outbase} 2> {log}
    """

rule plot_whatever:
    input:
        os.path.join(OUTMAIN, "jsons_statistics")
    output:
        os.path.join(OUTMAIN, "plots", "matches_pct.passed.png"),
        os.path.join(OUTMAIN, "plots", "insert_sizes.passed.png"),
        os.path.join(OUTMAIN, "plots", "mapping_quality.passed.png"),
        os.path.join(OUTMAIN, "plots", "query_lengths.passed.png"),
    run:
        f = os.path.join(input[0], "matches_pct.passed.txt")
        shell("{RSCRIPT} {SCRIPT_DIR}/plotting/whatever.R {f} {output[0]} pct_matches bar")
        f = os.path.join(input[0], "insert_sizes.passed.txt")
        shell("{RSCRIPT} {SCRIPT_DIR}/plotting/whatever.R {f} {output[1]} insert_size line")
        f = os.path.join(input[0], "mapping_quality.passed.txt")
        shell("{RSCRIPT} {SCRIPT_DIR}/plotting/whatever.R {f} {output[2]} mapping_quality bar")
        f = os.path.join(input[0], "query_lengths.passed.txt")
        shell("{RSCRIPT} {SCRIPT_DIR}/plotting/whatever.R {f} {output[3]} query_length bar")

rule plot_whatever_stats_excl:
    input:
        os.path.join(OUTMAIN, "jsons_statistics")
    output:
        os.path.join(OUTMAIN, "plots", "stats_exclusion_freq.pdf"),
        os.path.join(OUTMAIN, "plots", "stats_exclusion_count.pdf"),
    run:
        f = os.path.join(input[0], "stats.txt")
        shell("{RSCRIPT} {SCRIPT_DIR}/plotting/whatever_stats_exclusioncrit.R {f} {output}")

rule plot_whatever_stats_pass:
    input:
        os.path.join(OUTMAIN, "jsons_statistics")
    output:
        os.path.join(OUTMAIN, "plots", "stats_passed_failed.png"),
    run:
        f = os.path.join(input[0], "stats.txt")
        shell("{RSCRIPT} {SCRIPT_DIR}/plotting/whatever_stats_passed.R {f} {output}")


rule plot_whatever_depths:
    input:
        os.path.join(OUTMAIN, "depths", "close.depthSample"),
        os.path.join(OUTMAIN, "depths", "distant.depthSample"),
        os.path.join(OUTMAIN, "bamlists", "close.names"),
        os.path.join(OUTMAIN, "bamlists", "distant.names"),
    params:
        config["refs"]["close_name"],
        config["refs"]["distant_name"],
    output:
        os.path.join(OUTMAIN, "plots", "depths.png"),
        os.path.join(OUTMAIN, "plots", "depths_zoom.png"),
    shell:
        "{RSCRIPT} {SCRIPT_DIR}/plotting/whatever_depths.R {input} {params} {output}"

rule plot_whatever_depths_global:
    input:
        os.path.join(OUTMAIN, "depths", "close.depthGlobal"),
        os.path.join(OUTMAIN, "depths", "distant.depthGlobal"),
    params:
        config["refs"]["close_name"],
        config["refs"]["distant_name"],
    output:
        os.path.join(OUTMAIN, "plots", "depths_global.png"),
    shell:
        "{RSCRIPT} {SCRIPT_DIR}/plotting/whatever_depths_global.R {input} {params} {output}"


##########
# HETERO #
##########

rule angsd_hetero:
    input:
        bam = lambda wildcards: config["bams"][wildcards.whatref][wildcards.sample],
        regions = os.path.join(OUTMAIN, "scaffold_size", "subsample", "{whatref}.chrom"),
        sites = os.path.join(OUTMAIN, "scaffold_size", "subsample", "{whatref}.sites"),
        ref = lambda wildcards: config["refs"][wildcards.whatref],
    output:
        os.path.join(OUTMAIN, "hetero", "{whatref}", "{sample}.saf.idx"),
        os.path.join(OUTMAIN, "hetero", "{whatref}", "{sample}.saf.pos.gz"),
        os.path.join(OUTMAIN, "hetero", "{whatref}", "{sample}.saf.gz"),
    params:
        outbase = lambda wildcards, output: output[0][:-8]
    log:
        os.path.join(OUTMAIN, "hetero", "{whatref}", "{sample}.log")
    shell:
        "{ANGSD} -sites {input.sites} -anc {input.ref} -rf {input.regions} -howoften 1000000 -minMapQ {MINMAPQ} -minQ {MINQ}  -i {input.bam} -out {params.outbase} -gl 2 -dosaf 1 2> {log}"

rule realsfs_hetero:
    input:
        os.path.join(OUTMAIN, "hetero", "{whatref}", "{sample}.saf.idx")
    output:
        os.path.join(OUTMAIN, "hetero", "{whatref}", "{sample}.sfs")
    log:
        os.path.join(OUTMAIN, "hetero", "{whatref}", "{sample}.realsfs.log")
    threads: 10
    shell:
        "{REALSFS} -P {threads} -fold 1 {input} > {output} 2> {log}"

rule collect_hetero:
    input:
        expand(os.path.join(OUTMAIN, "hetero", "{whatref}", "{sample}.sfs"), whatref=whatrefs, sample=config["samples"])
    output:
        os.path.join(OUTMAIN, "hetero", "collect_het.txt")
    run:
        import os
        with open(output[0], 'w') as fh:
            for x in input:
                fields = x.split("/")
                sample = fields[-1][:-4]
                whatref = fields[-2]
                with open(x, 'r') as fhin:
                    data = fhin.read().rstrip().split()
                    data = [float(x) for x in data]
                    het = data[1] / (data[0]+data[1])
                print(f"{sample} {whatref} {het}", file=fh)

rule plot_hetero:
    input:
        os.path.join(OUTMAIN, "hetero", "collect_het.txt")
    output:
        os.path.join(OUTMAIN, "hetero", "collect_het.png")
    shell:
        "{RSCRIPT} {SCRIPT_DIR}/plotting/hetero_plot.R {input} {output}"




######################
# 2D-SFS RELATEDNESS #
######################


rule sfs_2d:
    input:
        saf_idx1 = os.path.join(OUTMAIN, "hetero", "{whatref}", "{s1}.saf.idx"),
        saf_idx2 = os.path.join(OUTMAIN, "hetero", "{whatref}", "{s2}.saf.idx"),
    output:
        os.path.join(OUTMAIN, "sfs_2d", "{whatref}", "{s1}-{s2}.sfs")
    threads: 3
    log:
        os.path.join(OUTMAIN, "sfs_2d", "{whatref}", "{s1}-{s2}.log")
    shell:
        "{REALSFS} -P {threads} {input.saf_idx1} {input.saf_idx2} > {output} 2> {log}"


rule collect_2dsfs:
    input:
        expand(os.path.join(OUTMAIN, "sfs_2d", "{{whatref}}", "{s[0]}-{s[1]}.sfs"), s=it.combinations(config["samples"], 2))
    output:
        f=os.path.join(OUTMAIN, "sfs_2d", "{whatref}", "collected.txt")
    run:
        import os
        import pandas as pd
        data = []
        names = []
        for x in input:
            name = os.path.basename(x).replace(".sfs", "")
            names.append(name)
            with open(x, 'r') as fh:
                t = fh.read()
                data.append([float(x) for x in t.rstrip().split()])
        a = pd.DataFrame(data, index=names, columns=["aaAA","aaAD","aaDD","adAA","adAD","adDD","ddAA","ddAD","ddDD"])

        a["r0"] = (a["aaDD"]+a["ddAA"])/a["adAD"]
        a["r1"] = a["adAD"] / (a.iloc[:,[1,2,3,5,6,7]].sum(1))
        a["king"] = (a["adAD"] - 2*a[["aaDD", "ddAA"]].sum(1)) / (a.iloc[:,[1,3,5,7]].sum(1) + 2*a["adAD"])
        a.to_csv(output.f, index=True, header=True, index_label="id", sep=" ")


rule plot_high_king:
    input:
        os.path.join(OUTMAIN, "sfs_2d", "{whatref}", "collected.txt")
    output:
        os.path.join(OUTMAIN, "sfs_2d", "{whatref}_king.png")
    shell:
        "{RSCRIPT} {SCRIPT_DIR}/plotting/related_king.R {input} {output}"
