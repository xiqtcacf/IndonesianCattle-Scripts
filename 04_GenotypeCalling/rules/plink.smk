rule vcf_to_plink:
    """ Convert VCF to plink """
    input:
        vcf=VCF_DIR / "{prefix}.bcf.gz",
        index=VCF_DIR / "{prefix}.bcf.gz.csi",
    output:
        plink_files=multiext(
            str(PLINK_DIR / "{prefix}"),
            ".bed",
            ".bim",
            ".fam",
            ".frq",
            ".imiss",
            ".lmiss",
            ".log",
            ".nosex",
        ),
    params:
        stem=lambda wildcards, output: output.plink_files[0][:-len(".bed")],
    threads: 1
    shell:
        """
        ( \
            n_chr=$({BCFTOOLS} index --stats {input.index} | wc -l); \
            {BCFTOOLS} view \
                --output-type v \
                {input.vcf} | \
            {PLINK} \
                --vcf /dev/stdin \
                --allow-extra-chr \
                --chr-set ${{n_chr}} \
                --out {params.stem} \
                --double-id \
                --make-bed \
                --freq \
                --missing \
        ) > /dev/null
        """


rule plink_filter_maf:
    """ Convert VCF to plink """
    input:
        plink_files=multiext(str(PLINK_DIR / "{prefix}"), ".bed", ".bim", ".fam"),
        index=VCF_DIR / "{prefix}.bcf.gz.csi",
    output:
        plink_files=multiext(
            str(PLINK_DIR / "{prefix}_maf{maf}pct"),
            ".bed",
            ".bim",
            ".fam",
            ".frq",
            ".imiss",
            ".lmiss",
            ".log",
            ".nosex",
        ),
    params:
        input_stem=lambda wildcards, input: input.plink_files[0][:-len(".bed")],
        output_stem=lambda wildcards, output: output.plink_files[0][:-len(".bed")],
        maf=lambda wildcards: int(wildcards.maf) / 100,
    threads: 1
    shell:
        """
        ( \
            n_chr=$({BCFTOOLS} index --stats {input.index} | wc -l); \
            {PLINK} \
                --bfile {params.input_stem} \
                --allow-extra-chr \
                --chr-set ${{n_chr}} \
                --maf {params.maf} \
                --out {params.output_stem} \
                --make-bed \
                --freq \
                --missing \
        ) > /dev/null
        """
