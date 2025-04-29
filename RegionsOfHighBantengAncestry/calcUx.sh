#!/bin/sh

#Purpose: To calculate number of SNPs with certain filters in non-overlapping windows

#Make windows - already done ###/binf-isilon/hellergrp/wlk579/Bos_project/Josiah/pure_chinese_zebu/
#bedtools makewindows -g /home/klp437/bos/data/bosChrSizes.txt -w 50000 > /binf-isilon/hellergrp/wlk579/Bos_project/Josiah/pure_chinese_zebu/bosARS-UCD1.3.50kbwindows.bed

#Define input variables
outputdir="/maps/projects/bos/people/mdn487/ancestry_selection/UaLoterc"
inputtxt=$1 # phased233_maf05NEWnaming.abcfiltered.filltags.dAFs.txt #All sites describing AFs
countFilter="0" #Minimum AC
aFilter="0.05"  #Cattle
bFilter=$2   # top 5% loter threshold of each respective breeds
cFilter="0.9"   #Banteng
windowsfile=$3

#Output variables
filteredSites="$outputdir/`basename $inputtxt .txt`.${countFilter}countFilt_${aFilter}a_${bFilter}b_${cFilter}c.txt"
outputtxt="$outputdir/`basename $windowsfile .bed`.`basename $inputtxt .txt`.${countFilter}countFilt_${aFilter}a_${bFilter}b_${cFilter}c.abcSNPcounts.txt"
tmpdir=$outputdir/tmp
mkdir -p $tmpdir

#######
#All sites
#Filter for counts, with AC!=0
tail -n+2 $inputtxt | awk -v countFilter=$countFilter '$6>countFilter' | awk '{print $1,$2,($2+1)}' OFS="\t" > $tmpdir/allsites.$countFilter.tmp


#Count all sites per window
bedtools intersect -c -a $windowsfile -b $tmpdir/allsites.$countFilter.tmp > $tmpdir/allsites.$countFilter.windows.tmp
rm -f $tmpdir/allsites.$countFilter.tmp

#Filtered sites
#Filter for ≤0.05, ≥0.5/0.8, ≥0.9
tail -n+2 $inputtxt | awk -v countFilter=$countFilter '$6>countFilter'  | awk -v aFilter=$aFilter -v bFilter=$bFilter -v cFilter=$cFilter '$17<=aFilter && $18>=bFilter && $19>=cFilter' OFS="\t" > $filteredSites
