#!/bin/bash
pop=$1 ### give a list of individuals for specific population e.g. Madura Pesisir Here should be All of Pops:allPop

export DIR=/home/wlk579/2.0Indonesia_Bos_project/Indonesia_Bos_project/5.Analyses_King015/emu_PCA
cd $DIR

module load bedops
module load R/4.0.3
bedtools=/home/wlk579/2.0Indonesia_Bos_project/Indonesia_Bos_project/5.Analyses_King015/hmmix0.6.9_AllResultsPlot/heatmap/bedtools

source /maps/projects/seqafrica/apps/modules/activate.sh
module load plink
bcftools=/projects/popgen/apps/Modules/software/bcftools/1.16/bin/bcftools

angsd=/home/wlk579/Server_bos/apps/angsd/angsd

####First get all pure Banteng reference bcf file and with maf5%: BaliBali, BaliAustralia, CaptiveBanteng
$bcftools view -i 'AF>0.05 & AF<0.95' -s a1,a12,a47,a49,a50,a8,a10,a45,a48,a5,a6,a7,N_2hmr,N_5hmr_1,N_6hmr,N_6t_11,N_Gt_25,N_gt15,N_gt16,N_gt6,N_gt8,N_kb_3,N_mk01,N_mk1,N_tk_1,N_tk5,N_tk7,N_wn_2,N_wn_4,N_y_13,N_6t_7,LIB112407_Banteng_85B_Texas,LIB112408_Banteng_82B_Texas,LIB112409_Banteng_86B1_Texas,B._javanicus_W.Zoo_DB,B._javanicus_SD.Zoo_OR206 314inds.imputated.BosTau9.bcf -Ob -o Maf5.PureBantengRef.BaliAusiCaptive.bcf.gz
$bcftools index Maf5.PureBantengRef.bcf.gz

###Second get bed file for each individual for all of 29 chromosomes for Archaic regions
while read sample ; do
# get correct end position which need end + 10k and different quantile for ratio of BantengAncestryAdmixSNPs
   Rscript scripts/get_BantengAdmixRation.R input_admixedPop/Archaic.Pro0.9.Java.Gaur.$sample.$pop
## get the quantile file proportion
   cut -f14 q20_data.NewFormat.txt > BantengAdmixedSNPs_PCAemu/quantile_proportion/q20_data.$sample
#only output bed file format for each quantile
   cut -f1,2,12 q20_data.NewFormat.txt > Pro0.9.Archaic.q20_data.$sample.bed
   sed -i '1d' Pro0.9.Archaic.q20_data.$sample.bed
   rm q20_data.NewFormat.txt
## get bcf file for this specific sample and BantengRation regions, maf0.05
   $bcftools view -i 'AF>0.05 & AF<0.95' -R Pro0.9.Archaic.q20_data.$sample.bed -s $sample 314inds.imputated.BosTau9.bcf -Ob -o BantengAdmixedSNPs_PCAemu/results/Pro0.9.Archaic.q20_data.$sample.bcf.gz
   $bcftools index BantengAdmixedSNPs_PCAemu/results/Pro0.9.Archaic.q20_data.$sample.bcf.gz
   rm  Pro0.9.Archaic.q20_data.$sample.bed
done < input/$pop.pure.inds

###merge all admixed samples and  PureBantengBcf
$bcftools merge BantengAdmixedSNPs_PCAemu/results/*.bcf.gz -Ob -o Hmmix.q20BantengRegion.All.merged.bcf.gz
$bcftools index Hmmix.q20BantengRegion.All.merged.bcf.gz
###do plink version
plink --bcf Hmmix.q20BantengRegion.All.merged.bcf.gz \
      --chr-set 29 \
      --out Hmmix.q20BantengRegion.All.merged \
      --double-id \
      --allow-extra-chr

###Banteng=8, Zebu=7
emu --bfile Hmmix.q20BantengRegion.All.merged --n_eig 8 --n_out 8 --iter 1000 --threads 64 --out Hmmix.q20BantengRegion.All.merged.emu -f 0.05

