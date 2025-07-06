#!/bin/bash

batch=$1 ####pop name
export DIR=/isdata/hellergrp/wlk579/Indonesia_Bos_project/5.Analyses_King015/AncestryHMM/$batch
cd $DIR

#input_file=/isdata/hellergrp/wlk579/Bos_project/5.Analyses/AdmixtureTime/Ancestry_HMM/imputed.bcf
admixture_proportion=/isdata/hellergrp/wlk579/Indonesia_Bos_project/5.Analyses_King015/AncestryHMM/format.Casia.bantengPro.txt

#### step1 get the correct way vcf: all of chromosomes
bcftools view $input_file -o imputed.vcf -O v
bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\t%INFO\tGT\t[%GT\t]\n' imputed.vcf > format.imputed.vcf
grep '#' imputed.vcf > HEADER.txt
cat  HEADER.txt format.imputed.vcf > new.format.imputed.vcf

#### step2 get the specific pop input file
grep "$batch" $admixture_proportion | cut -f1,3 > $batch.samples
sed "s/$batch/admixed/g" "$batch.samples" > $batch.inds
sed "s/$batch//g" "$batch.samples" > $batch.inds.forBCFTOOLS
sed  "s/\t//g" "$batch.inds.forBCFTOOLS" -i
sed "s/$batch/2/g" "$batch.samples" > $batch.inds.sample_file
cat ../2ref.inds $batch.inds > 2ref.$batch.inds
cat ../2ref.inds.forBCFTOOLS $batch.inds.forBCFTOOLS > 2ref.$batch.inds.forBCFTOOLS

bcftools view -S 2ref.$batch.inds.forBCFTOOLS ../new.format.imputed.vcf -o 2ref.$batch.new.format.imputed.vcf -O v
bcftools view -g ^miss -q 0.01 -Q 0.99 -o qc.2ref.$batch.new.format.imputed.vcf -Ov 2ref.$batch.new.format.imputed.vcf
python3 ../vcf2ahmm.py -g 1 -v qc.2ref.$batch.new.format.imputed.vcf -s 2ref.$batch.inds > 2ref.$batch.new.format.imputed.vcf.ahmm.input ###vcf2ahmm.py is from AncestryHMM 

####running programme notice here samples should only have admixed samples
grep 'admixed' 2ref.$batch.inds > 2ref.$batch.sample_file
sed 's/admixed/2/g' "2ref.$batch.sample_file" -i
grep "$batch" $admixture_proportion > $batch.loter
Banteng_Proportion=$(awk '{ total += $2; count++ } END { print total/count }' $batch.loter) #### admixture proportion per pop
Zebu_Proportion=$(bc <<< "1 - $Banteng_Proportion")
ancestry_hmm -i 2ref.$batch.new.format.imputed.vcf.ahmm.input -s 2ref.$batch.sample_file -a 2 $Banteng_Proportion $Zebu_Proportion -p 0 -500 -0.4 -p 1 100000 -0.6 -g -b 100 5000 > $batch.allChr.one-pulse.out
