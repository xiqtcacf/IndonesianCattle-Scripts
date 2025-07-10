#!/bin/bash
#name=$1 ###Zebu file name
export DIR=/home/wlk579/2.0Indonesia_Bos_project/Revision_NC/EmuZebu_SouthernIndian/emuPCA/
cd $DIR

module load perl
module unload gsl
module load gsl/2.5
module load bcftools
module load gcc
module load R/4.0.3
module load bedops

bedtools=/home/wlk579/2.0Indonesia_Bos_project/Indonesia_Bos_project/5.Analyses_King015/hmmix0.6.9_AllResultsPlot/heatmap/bedtools
plink=/home/wlk579/Server_seqafrica/apps/modules/software/plink/1.90b6.21/bin/plink

angsd=/home/wlk579/Server_bos/apps/angsd/angsd

###maf filter for new India samples
bcftools view -i 'AF>0.05 & AF<0.95' autosome29.all.imputed.vcf.gz -o Maf5.autosome29.all.imputed.bcf.gz -Ob
bcftools index Maf5.autosome29.all.imputed.bcf.gz

bcftools view -i 'AF>0.05 & AF<0.95' -S sample.list Maf5.autosome29.all.imputed.bcf.gz -o SouthIndia.bcf.gz -Ob
bcftools index SouthIndia.bcf.gz
bcftools merge NoK1.loter.SNP.ZebuAnc.bcfAll.bcf.gz SouthIndia.bcf.gz -o NoK1.loter.SNP.ZebuAnc.bcfAll.SouthIndia.bcf.gz -Ob
bcftools index NoK1.loter.SNP.ZebuAnc.bcfAll.SouthIndia.bcf.gz

###do emu K-1, so Banteng=8, Zebu=7
emu --bfile NoK1.loter.SNP.ZebuAnc.bcfAll.SouthIndia --n_eig 7 --n_out 7 --iter 1000 --threads 64 --out NoK1.loter.SNP.ZebuAnc.bcfAll.SouthIndia.emu -f 0.05
