#!/bin/bash
batch=$1 ###pairwise files

export DIR=/home/wlk579/2.0Indonesia_Bos_project/Indonesia_Bos_project/5.Analyses_King015/hmmix0.6.9_AllResultsPlot/heatmap/Step2_results_IBS_matrix
cd $DIR

file1=/home/wlk579/2.0Indonesia_Bos_project/Indonesia_Bos_project/5.Analyses_King015/hmmix0.6.9_AllResultsPlot/heatmap/Step2_Pairwise_name/$batch.1
file2=/home/wlk579/2.0Indonesia_Bos_project/Indonesia_Bos_project/5.Analyses_King015/hmmix0.6.9_AllResultsPlot/heatmap/Step2_Pairwise_name/$batch.2

module load bedops
module load R/4.0.3
bedtools=/home/wlk579/2.0Indonesia_Bos_project/Indonesia_Bos_project/5.Analyses_King015/hmmix0.6.9_AllResultsPlot/heatmap/bedtools

source /maps/projects/seqafrica/apps/modules/activate.sh
module load plink
bcftools=/projects/popgen/apps/Modules/software/bcftools/1.16/bin/bcftools

angsd=/home/wlk579/Server_bos/apps/angsd/angsd

input=/home/wlk579/2.0Indonesia_Bos_project/Indonesia_Bos_project/5.Analyses_King015/hmmix0.6.9_AllResultsPlot/heatmap/Step1_results_ArchaicPro0.9.AllChr.10kWindow.PerSample

while IFS= read -r line1 && IFS= read -r line2 <&3; do

## get 10kWindow 100% overlapping archaic regions for each pairwise two individuals
    sort-bed $input/10kWindow.Archaic.Pro0.9.AllChr.$line1.allPop  > $input/Sorted.10kWindow.Archaic.Pro0.9.AllChr.$line1.allPop
    sort-bed $input/10kWindow.Archaic.Pro0.9.AllChr.$line2.allPop  > $input/Sorted.10kWindow.Archaic.Pro0.9.AllChr.$line2.allPop
    bedops --intersect $input/Sorted.10kWindow.Archaic.Pro0.9.AllChr.$line1.allPop $input/Sorted.10kWindow.Archaic.Pro0.9.AllChr.$line2.allPop > $input/ArchaicPro0.9.10kWindow.$line1.$line2.bed
    rm $input/Sorted.10kWindow.Archaic.Pro0.9.AllChr.$line1.allPop
    rm $input/Sorted.10kWindow.Archaic.Pro0.9.AllChr.$line2.allPop
## get two samples, maf0.05 and region bcf files by bcftools
    $bcftools view -i 'AF>0.05 & AF<0.95' -R $input/ArchaicPro0.9.10kWindow.$line1.$line2.bed -s $line1,$line2 ../314inds.imputated.BosTau9.bcf -Ob -o $line1.$line2.bcf.gz
    $bcftools index $line1.$line2.bcf.gz
#### do IBS tree for all samples, further be used as input of NJ tree
    plink --bcf $line1.$line2.bcf.gz \
          --chr-set 29 \
          --distance square ibs allele-ct \
          --out $line1.$line2.plink.ibs \
          --double-id \
          --allow-extra-chr
    rm $line1.$line2.bcf.gz
    rm $line1.$line2.bcf.gz.csi

done < $file1 3< $file2

