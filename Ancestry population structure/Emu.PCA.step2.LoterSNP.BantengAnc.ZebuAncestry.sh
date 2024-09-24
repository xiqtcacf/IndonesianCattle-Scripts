#!/bin/bash
pop=$1 ### give a list of individuals for specific population e.g. Madura Pesisir Here should be All of Pops:allPop

export DIR=/home/wlk579/2.0Indonesia_Bos_project/Indonesia_Bos_project/5.Analyses_King015/emu_PCA/0Using.Final_results_all/all_individuals
cd $DIR

module load bedops
module load R/4.0.3
bedtools=/home/wlk579/2.0Indonesia_Bos_project/Indonesia_Bos_project/5.Analyses_King015/hmmix0.6.9_AllResultsPlot/heatmap/bedtools

source /maps/projects/seqafrica/apps/modules/activate.sh
module load plink
bcftools=/projects/popgen/apps/Modules/software/bcftools/1.16/bin/bcftools

angsd=/home/wlk579/Server_bos/apps/angsd/angsd

###get Zebu/Banteng Ancestry loter specific bcf for each individual
while read sample ; do
### For loter, get Zebu/Banteng Ancestry bed file for each individual
  cut -f1,2 ../loter.SNP.bed_Ancestry_perInd/loter.ZebuAncestry.$sample | grep -v 'chr' > pos.loter.ZebuAncestry.$sample
  cut -f1,2 ../loter.SNP.bed_Ancestry_perInd/loter.BantengAncestry.$sample | grep -v 'chr' > pos.loter.BantengAncestry.$sample
### get Zebu/Banteng Ancestry loter specific bcf for each individual, maf0.05
  $bcftools view -i 'AF>0.05 & AF<0.95' -R pos.loter.ZebuAncestry.$sample -s $sample ../../314inds.imputated.BosTau9.bcf -Ob -o loter.SNP.ZebuAnc.bcfAll/$sample.bcf.gz
  $bcftools index loter.SNP.ZebuAnc.bcfAll/$sample.bcf.gz
  $bcftools view -i 'AF>0.05 & AF<0.95' -R pos.loter.BantengAncestry.$sample -s $sample ../../314inds.imputated.BosTau9.bcf -Ob -o loter.SNP.BantengAnc.bcfAll/$sample.bcf.gz
  $bcftools index loter.SNP.BantengAnc.bcfAll/$sample.bcf.gz
  rm pos.loter.ZebuAncestry.$sample
  rm pos.loter.BantengAncestry.$sample
done < ../$pop

##merge all admixed samples and PureZebuBCF
$bcftools merge results/*.bcf.gz -Ob -o loter.ZebuRegion.All.merged.bcf.gz
$bcftools index loter.ZebuRegion.All.merged.bcf.gz
###do plink version
plink --bcf loter.ZebuRegion.All.merged.bcf.gz \
      --chr-set 29 \
      --out loter.ZebuRegion.All.merged \
      --double-id \
      --allow-extra-chr

###merge all admixed samples and  PureBantengBcf
$bcftools merge results_BantengAncestry/*.bcf.gz -Ob -o loter.BantengRegion.All.merged.bcf.gz
$bcftools index loter.BantengRegion.All.merged.bcf.gz
###do plink version
plink --bcf loter.BantengRegion.All.merged.bcf.gz \
      --chr-set 29 \
      --out loter.BantengRegion.All.merged \
      --double-id \
      --allow-extra-chr

###Banteng=8, Zebu=7
emu --bfile loter.BantengRegion.All.merged --n_eig 8 --n_out 8 --iter 1000 --threads 64 --out loter.BantengRegion.All.merged.emu -f 0.05
emu --bfile loter.ZebuRegion.All.merged --n_eig 7 --n_out 7 --iter 1000 --threads 64 --out loter.ZebuRegion.All.merged.emu -f 0.05


