#!/bin/bash

pop=$1 ### give a list of individuals for specific population
#window=$2 ### fix with window Size 10kbp but mutation rate size fix with 1Mbp
#outgroup=$3 ###fix with SouthAsianZebu Reference,

export DIR=/isdata/hellergrp/wlk579/Indonesia_Bos_project/5.Analyses_King015/hmmix/$pop
cd $DIR

###modified jason file,get average Banteng ancestry for each admixed population, e.g.
grep 'Madura' Casia.bantengPro.txt | cut -f2 | awk  'BEGIN { total = 0; count = 0 } { total += $1; count += 1; } END { avg = total / count; print avg}'
###get per chr bcf
#bcftools view 314inds.imputated.BosTau9.bcf --regions NC_037328.1 -O b -o Chr1.314inds.imputated.BosTau9.bcf
#vcf=../314inds.imputated.BosTau9.bcf
#vcf=../Chr1.314inds.imputated.BosTau9.bcf

vcf_input=/isdata/hellergrp/wlk579/Indonesia_Bos_project/5.Analyses_King015/hmmix/bcf_perChr
BosTau9_fasta=../BosTau9.fasta
BosTau9_fai=../BosTau9.fasta.fai
Sites_filtering=/isdata/hellergrp/wlk579/Indonesia_Bos_project/Good_sites/autosome_rep_het_dep_map.bed

Zebu="Nelore_BINE1,Nelore_BINE2,Bhagnari_23,Gir_BIGI3,Gir_BIGI4,Hariana_Har03,Kangayam_SAMN10131257,Lohani_18,Red_Sindhi_303,Sahiwal_Sha3b,Tharparkar_Thar1,N_22,N_23,N_3,N_45,N_6,N_1,N_11,N_14,N_16,N_17,N_24,N_38_copy,N_8,N_9"

for i in {1..29} ###29chromosomes
do
     while read sample ; do
#######use weights, using modified jason
      hmmix create_outgroup -ind=$Zebu -vcf=$vcf_input/Chr${i}.314inds.imputated.BosTau9.bcf -out=Chr${i}.Weights.outgroup.txt -weights=$Sites_filtering
      hmmix mutation_rate -outgroup=Chr${i}.Weights.outgroup.txt -window_size=1000000 -out Chr${i}.Weights.mutationrate.bed -weights=$Sites_filtering
      hmmix create_ingroup -ind=$sample -vcf=$vcf_input/Chr${i}.314inds.imputated.BosTau9.bcf -outgroup=Chr${i}.Weights.outgroup.txt -out=Chr${i}.obs.Weights -weights=$Sites_filtering

      hmmix train -obs=Chr${i}.obs.Weights.$sample.txt -mutrates=Chr${i}.Weights.mutationrate.bed -param=10000.$pop.Cattle.Initialguesses.json -window_size=10000 -out=Chr${i}.Weights.trained.$sample.json -weights=$Sites_filtering
      hmmix decode -obs=Chr${i}.obs.Weights.$sample.txt -mutrates=Chr${i}.Weights.mutationrate.bed -param=Chr${i}.Weights.trained.$sample.json -window_size=10000 -out=Chr${i}.Weights.$sample.decoded -weights=$Sites_filtering
     done < ../$pop.inds
done

