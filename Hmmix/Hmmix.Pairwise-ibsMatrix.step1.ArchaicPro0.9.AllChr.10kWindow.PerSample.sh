#!/bin/bash
pop=$1 ### give a list of individuals for specific population e.g. Madura Pesisir Here should be All of Pops:allPop
pair_inds=$2 ###each pair of individuals

export DIR=/home/wlk579/2.0Indonesia_Bos_project/Indonesia_Bos_project/5.Analyses_King015/hmmix0.6.9_AllResultsPlot/heatmap/
cd $DIR

module load bedops
module load R/4.0.3
bedtools=/home/wlk579/2.0Indonesia_Bos_project/Indonesia_Bos_project/5.Analyses_King015/hmmix0.6.9_AllResultsPlot/heatmap/bedtools

###How to use this script: bash Hmmix.Pairwise-ibsMatrix.step1.ArchaicPro0.9.AllChr.10kWindow.PerSample.sh allPop
###get bed file for each individual for all of 29 chromosomes
while read sample ; do
   for i in {1..29}
   do
#Probality > 0.9
    awk '$6 > 0.9' input/$pop/Chr${i}.Weights.$sample.decoded.diploid.txt > Pro0.9.Chr${i}.Weights.$sample.decoded.diploid.txt
#only get Archaic regions
    grep -E 'state|Archaic' Pro0.9.Chr${i}.Weights.$sample.decoded.diploid.txt > Archaic.Pro0.9.Chr${i}.Weights.$sample.decoded.diploid.txt
    rm Pro0.9.Chr${i}.Weights.$sample.decoded.diploid.txt
# get correct end position which need end + 10k
    Rscript scripts/get_correctEndPosition.R Archaic.Pro0.9.Chr${i}.Weights.$sample.decoded.diploid.txt
    rm Archaic.Pro0.9.Chr${i}.Weights.$sample.decoded.diploid.txt
#only output bed file format
    cut -f1,2,8 NewFormat.txt > Archaic.Pro0.9.Chr${i}.$sample.$pop
    sed -i '1d' Archaic.Pro0.9.Chr${i}.$sample.$pop
    rm NewFormat.txt

## merge Archric regions based on maximum distance of 10000bp
##  $bedtools merge -d 10000 -i a.Archaic.Pro0.9.Chr${i}.$sample.$pop > Archaic.Pro0.9.Chr${i}.$sample.$pop
##  rm a.Archaic.Pro0.9.Chr${i}.$sample.$pop

#get 10k window for each chr and each individual, this will be much more easier to get overlapping region
    $bedtools makewindows -b Archaic.Pro0.9.Chr${i}.$sample.$pop -w 10000 > 10kWindow.Archaic.Pro0.9.Chr${i}.$sample.$pop
    rm Archaic.Pro0.9.Chr${i}.$sample.$pop
    done
### get all chromosomes together
  cat 10kWindow.Archaic.Pro0.9.Chr1.$sample.$pop \
      10kWindow.Archaic.Pro0.9.Chr2.$sample.$pop \
      10kWindow.Archaic.Pro0.9.Chr3.$sample.$pop \
      10kWindow.Archaic.Pro0.9.Chr4.$sample.$pop \
      10kWindow.Archaic.Pro0.9.Chr5.$sample.$pop \
      10kWindow.Archaic.Pro0.9.Chr6.$sample.$pop \
      10kWindow.Archaic.Pro0.9.Chr7.$sample.$pop \
      10kWindow.Archaic.Pro0.9.Chr8.$sample.$pop \
      10kWindow.Archaic.Pro0.9.Chr9.$sample.$pop \
      10kWindow.Archaic.Pro0.9.Chr10.$sample.$pop \
      10kWindow.Archaic.Pro0.9.Chr11.$sample.$pop \
      10kWindow.Archaic.Pro0.9.Chr12.$sample.$pop \
      10kWindow.Archaic.Pro0.9.Chr13.$sample.$pop \
      10kWindow.Archaic.Pro0.9.Chr14.$sample.$pop \
      10kWindow.Archaic.Pro0.9.Chr15.$sample.$pop \
      10kWindow.Archaic.Pro0.9.Chr16.$sample.$pop \
      10kWindow.Archaic.Pro0.9.Chr17.$sample.$pop \
      10kWindow.Archaic.Pro0.9.Chr18.$sample.$pop \
      10kWindow.Archaic.Pro0.9.Chr19.$sample.$pop \
      10kWindow.Archaic.Pro0.9.Chr20.$sample.$pop \
      10kWindow.Archaic.Pro0.9.Chr21.$sample.$pop \
      10kWindow.Archaic.Pro0.9.Chr22.$sample.$pop \
      10kWindow.Archaic.Pro0.9.Chr23.$sample.$pop \
      10kWindow.Archaic.Pro0.9.Chr24.$sample.$pop \
      10kWindow.Archaic.Pro0.9.Chr25.$sample.$pop \
      10kWindow.Archaic.Pro0.9.Chr26.$sample.$pop \
      10kWindow.Archaic.Pro0.9.Chr27.$sample.$pop \
      10kWindow.Archaic.Pro0.9.Chr28.$sample.$pop \
      10kWindow.Archaic.Pro0.9.Chr29.$sample.$pop > 10kWindow.Archaic.Pro0.9.AllChr.$sample.$pop
   rm 10kWindow.Archaic.Pro0.9.Chr*.$sample.$pop

done < input/$pop.pure.inds
