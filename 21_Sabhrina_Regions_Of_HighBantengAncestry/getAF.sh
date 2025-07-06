#!/bin/sh

export BCFTOOLS=/isdata/hellergrp/wlk579/software/bcftools/bcftools
export BCFTOOLS_PLUGINS=/isdata/hellergrp/wlk579/software/bcftools/plugins

inputvcf="/binf-isilon/hellergrp/wlk579/Bos_project/Josiah/newAnalyses/phased233_maf05NEWnaming.abcfilte
red.bcf"
groupfile="/binf-isilon/hellergrp/wlk579/Bos_project/Josiah/newAnalyses/pasundan/newPasundan.txt"
tmpvcf="/binf-isilon/hellergrp/wlk579/Bos_project/Josiah/newAnalyses/pasundan/phased233_maf05NEWnaming.abcfiltered.tmp.bcf"
outputtmpvcf="/binf-isilon/hellergrp/wlk579/Bos_project/Josiah/newAnalyses/pasundan/phased233_maf05NEWnamin.abcfiltered.filltags.tmp.bcf"
outputvcf="/binf-isilon/hellergrp/wlk579/Bos_project/Josiah/newAnalyses/pasundan/phased233_maf05NEWnaming.abcfiltered.filltags.bcf"
outtxt="/binf-isilon/hellergrp/wlk579/Bos_project/Josiah/newAnalyses/pasundan/phased233_maf05NEWnaming.abcfiltered.filltags.txt"
outtxt2="/binf-isilon/hellergrp/wlk579/Bos_project/Josiah/newAnalyses/pasundan/phased233_maf05NEWnaming.abcfiltered.filltags.dAFs.txt"


#Rename a5, a6, a7 as bcftools won't take it -> a5a a6a a7a
$BCFTOOLS reheader --threads 20 -s /binf-isilon/hellergrp/wlk579/Bos_project/Josiah/newAnalyses/newsamplenames.txt $inputvcf > $tmpvcf

#Fill tags
$BCFTOOLS plugin fill-tags  --output-type b --threads 60 -o $outputtmpvcf $tmpvcf -- --samples-file $groupfile
$BCFTOOLS index --threads 60 $outputtmpvcf

#Reheader to match original sample names
$BCFTOOLS reheader --threads 60 -s /binf-isilon/hellergrp/wlk579/Bos_project/Josiah/newAnalyses/oldsamplenames.txt $outputtmpvcf > $outputvcf
$BCFTOOLS index --threads 60 $outputvcf

#Clean up
rm -f $outputtmpvcf $tmpvcf

#Get header
$BCFTOOLS query -H -f '%CHROM\t%POS\t%REF\t%ALT\t%AN\t%AC\t%AF\t%MAF\t%AC_Hom\t%AC_Het\t%AC_Cattle\t%AC_Pasundan\t%AC_Banteng\t%AN_Cattle\t%AN_Pasundan\t%AN_Banteng\t%AF_Cattle\t%AF_Pasundan\t%AF_Banteng\n' $outputvcf | head -n1  | sed 's/ # //' | sed 's/# //' | sed -e 's/\[[^][]*\]//g' > $outtxt2

#Generate dAF values
$BCFTOOLS query -H -f '%CHROM\t%POS\t%REF\t%ALT\t%AN\t%AC\t%AF\t%MAF\t%AC_Hom\t%AC_Het\t%AC_Cattle\t%AC_Pasundan\t%AC_Banteng\t%AN_Cattle\t%AN_Pasundan\t%AN_Banteng\t%AF_Cattle\t%AF_Pasundan\t%AF_Banteng\n' $outputvcf | grep -v '#' >> $outtxt2 ##Note non-divisible values are "." e.g. row 1586

$BCFTOOLS query -H -f '%CHROM\t%POS\t%REF\t%ALT\t%AN\t%AC\t%AF\t%MAF\t%AC_Hom\t%AC_Het\t%AC_Cattle\t%AC_Pasundan\t%AC_Banteng\t%AN_Cattle\t%AN_Pasundan\t%AN_Banteng\n' $outputvcf | grep -v '#' | awk '{print $0,$14!=0?$11/$14:"NA", $15!=0?$12/$15:"NA",$16!=0?$13/$16:"NA"}' OFS="\t" >> $outtxt2
