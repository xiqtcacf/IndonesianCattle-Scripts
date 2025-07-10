#!/bin/bash
##bash vcf_splitter.sh vcfinput outputdir prefix_output threads
echo "VCF input is "$1 # holds the current script
echo "Output dir is " $2 
echo "Prefix for output is "$3
echo "Threads set to "$4

input=$1
outdir=$2
prefix=$3
t=$4

bcftools index -s $input | cut -f 1 | while read C; do bcftools view --threads $t -O z -o ${outdir}/${prefix}_${C}.vcf.gz $input "${C}" ; done
