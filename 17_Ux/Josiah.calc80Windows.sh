#!/bin/bash

#Purpose: To calculate number of SNPs with certain filters in non-overlapping windows

#Make windows - already done ###/binf-isilon/hellergrp/wlk579/Bos_project/Josiah/pure_chinese_zebu/
#bedtools makewindows -g /home/klp437/bos/data/bosChrSizes.txt -w 50000 > /binf-isilon/hellergrp/wlk579/Bos_project/Josiah/pure_chinese_zebu/bosARS-UCD1.3.50kbwindows.bed

#Define input variables
outputdir="/binf-isilon/hellergrp/wlk579/Bos_project/Josiah/newAnalyses/madura"
inputtxt="/binf-isilon/hellergrp/wlk579/Bos_project/Josiah/newAnalyses/madura/phased233_maf05NEWnaming.abcfiltered.filltags.dAFs.txt" #All sites describing AFs
countFilter="0" #Minimum AC
aFilter="0.05"  #Cattle
bFilter="0.8"   #madura
cFilter="0.9"   #Banteng
windowsfile="/binf-isilon/hellergrp/wlk579/Bos_project/Josiah/madura/bosARS-UCD1.3.50kbwindows.bed"

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
awk '{print $1,$2,($2+1)}' OFS="\t" $filteredSites > $tmpdir/`basename $filteredSites`.bed.tmp

#Count filtered sites per window
bedtools intersect -c -a $windowsfile -b $tmpdir/`basename $filteredSites`.bed.tmp > $tmpdir/`basename $filteredSites`.tmp

#Combine output and remove temporary bed files
join -j1 -o1.2,1.3,1.4,1.5,2.5 <(<$tmpdir/`basename $filteredSites`.tmp awk '{print $1"-"$2"-"$3" "$0}' | sort -k1,1) <(< $tmpdir/allsites.$countFilter.windows.tmp  awk '{print $1"-"$2"-"$3" "$0}' | sort -k1,1) | awk '$1=$1' OFS="\t" | sort -k 1,1 -k2,2n > $outputtxt
rm -f $tmpdir/`basename $filteredSites`.tmp $tmpdir/allsites.$countFilter.windows.tmp $tmpdir/`basename $filteredSites`.bed.tmp
rmdir $tmpdir