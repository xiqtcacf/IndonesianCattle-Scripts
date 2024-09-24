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

####First get pure Zebu bcf file and with maf5%
$bcftools view -i 'AF>0.05 & AF<0.95' -s Nelore_BINE1,Nelore_BINE2,Bhagnari_23,Gir_BIGI3,Gir_BIGI4,Hariana_Har03,Kangayam_SAMN10131257,Lohani_18,Red_Sindhi_303,Sahiwal_Sha3b,Tharparkar_Thar1,N_22,N_23,N_3,N_45,N_6,N_1,N_11,N_14,N_16,N_17,N_24,N_38_copy,N_8,N_9 314inds.imputated.BosTau9.bcf -Ob -o Maf5.PureZebuRef.bcf.gz
$bcftools index Maf5.PureZebuRef.bcf.gz

###Second get bed file for each individual for all of 29 chromosomes, no Archaic region, only cattle regions, no matter Probablity
while read sample ; do
   for i in {1..29}
   do
#only get Cattle regions
    grep -E 'state|Cattle' input/$pop/Chr${i}.Weights.$sample.decoded.diploid.txt > Cattle.Chr${i}.Weights.$sample.decoded.diploid.txt
# get correct end position which need end + 10k
    Rscript scripts/get_correctEndPosition.R Cattle.Chr${i}.Weights.$sample.decoded.diploid.txt
    rm Cattle.Chr${i}.Weights.$sample.decoded.diploid.txt
#only output bed file format
    cut -f1,2,8 NewFormat.txt > Cattle.Chr${i}.$sample.$pop
    sed -i '1d' Cattle.Chr${i}.$sample.$pop
    rm NewFormat.txt
    done
### get all chromosomes together
  cat Cattle.Chr1.$sample.$pop \
      Cattle.Chr2.$sample.$pop \
      Cattle.Chr3.$sample.$pop \
      Cattle.Chr4.$sample.$pop \
      Cattle.Chr5.$sample.$pop \
      Cattle.Chr6.$sample.$pop \
      Cattle.Chr7.$sample.$pop \
      Cattle.Chr8.$sample.$pop \
      Cattle.Chr9.$sample.$pop \
      Cattle.Chr10.$sample.$pop \
      Cattle.Chr11.$sample.$pop \
      Cattle.Chr12.$sample.$pop \
      Cattle.Chr13.$sample.$pop \
      Cattle.Chr14.$sample.$pop \
      Cattle.Chr15.$sample.$pop \
      Cattle.Chr16.$sample.$pop \
      Cattle.Chr17.$sample.$pop \
      Cattle.Chr18.$sample.$pop \
      Cattle.Chr19.$sample.$pop \
      Cattle.Chr20.$sample.$pop \
      Cattle.Chr21.$sample.$pop \
      Cattle.Chr22.$sample.$pop \
      Cattle.Chr23.$sample.$pop \
      Cattle.Chr24.$sample.$pop \
      Cattle.Chr25.$sample.$pop \
      Cattle.Chr26.$sample.$pop \
      Cattle.Chr27.$sample.$pop \
      Cattle.Chr28.$sample.$pop \
      Cattle.Chr29.$sample.$pop > Cattle.AllChr.$sample.$pop.bed
   rm Cattle.Chr*.$sample.$pop
## get bcf file for this specific sample and Cattles regions, maf0.05
   $bcftools view -i 'AF>0.05 & AF<0.95' -R Cattle.AllChr.$sample.$pop.bed -s $sample 314inds.imputated.BosTau9.bcf -Ob -o results/$sample.bcf.gz
   $bcftools index results/$sample.bcf.gz
   rm  Cattle.AllChr.$sample.$pop.bed
done < input/$pop.pure.inds

##merge all 90 admixed samples and PureZebuBCF
$bcftools merge results/*.bcf.gz -Ob -o CattlesRegion.All.merged.bcf.gz
$bcftools index CattlesRegion.All.merged.bcf.gz
###do plink version
plink --bcf CattlesRegion.All.merged.bcf.gz \
      --chr-set 29 \
      --out CattlesRegion.All.merged \
      --double-id \
      --allow-extra-chr
###do emu
emu --bfile CattlesRegion.All.merged --n_eig 7 --n_out 7 --iter 1000 --threads 64 --out CattlesRegion.All.merged.emu -f 0.05
