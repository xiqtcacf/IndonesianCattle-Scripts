#!/bin/bash

batch=$1 #### each individual as a batch
folder=$2 ###9 folder for each 26inds

export DIR=/home/wlk579/2.0Indonesia_Bos_project/Indonesia_Bos_project/5.Analyses_King015/2.Fig.Hete_Roh/Hete_NoROH/
cd $DIR

bedtools=/home/wlk579/Server_bos/apps/bedtools2/bin/bedtools
bcftools=/home/wlk579/Server_bos/apps/bcftools

input=/home/wlk579/2.0Indonesia_Bos_project/Indonesia_Bos_project/5.Analyses_King015/2.Fig.Hete_Roh/Hete_NoROH/

###get sites by removing ROH >0.5M from good_bed.sites per individual

while read ind ; do

grep "$ind" $input/Using.ForHeteNoRoh.roh_5M.txt | cut -f2,3,4 > $input/$ind.ROH.location
$bedtools subtract -a $input/Good_sites_Bostau9.bed -b $input/$ind.ROH.location > $input/Good_sites_NoRoh/$ind.Good_sites_Bostau9.bed
rm $input/$ind.ROH.location

###get depth for each individual as input of config file
grep -A 2 "$ind:" $input/config_all.yaml | sed 1,5d | sed 2d > $input/$ind.low.depth
grep -A 2 "$ind:" $input/config_all.yaml | sed 1,6d > $input/$ind.high.depth

### second step, prepare config file for each individual

echo "PSMC_DIR: "/home/wlk579/Server_bos/apps/psmc"" > $input/$folder/$ind.config
echo "BCFTOOLS: "/home/wlk579/Server_bos/apps/bcftools"" >> $input/$folder/$ind.config
echo "VCFUTILS: "/projects/popgen/apps/Modules/software/bcftools/1.16/bin/vcfutils.pl"" >> $input/$folder/$ind.config
echo "" >> $input/$folder/$ind.config
echo "" >> $input/$folder/$ind.config
echo "ALLELE_SUPPORT: 2 " >> $input/$folder/$ind.config
echo "MIN_MQ: 25" >> $input/$folder/$ind.config
echo "MIN_BQ: 30">> $input/$folder/$ind.config
echo "" >> $input/$folder/$ind.config
echo "" >> $input/$folder/$ind.config
echo "OUTDIR: "$ind"" >> $input/$folder/$ind.config
echo "REF: "/home/wlk579/Server_bos/data/reference/bosTau9.fasta"" >> $input/$folder/$ind.config
echo "GOODBED: "/home/wlk579/2.0Indonesia_Bos_project/Indonesia_Bos_project/5.Analyses_King015/2.Fig.Hete_Roh/Hete_NoROH/Good_sites_NoRoh/$ind.Good_sites_Bostau9.bed"" >> $input/$folder/$ind.config
echo "" >> $input/$folder/$ind.config
echo "" >> $input/$folder/$ind.config
echo "SAMPLES:" >> $input/$folder/$ind.config
echo "  $ind:  /home/wlk579/2.0Indonesia_Bos_project/Indonesia_Bos_project/data_King015/all_filtered_raw_bam/$ind.BosTau9.bam" >> $input/$folder/$ind.config
echo "" >> $input/$folder/$ind.config
echo "DEPTHS:" >> $input/$folder/$ind.config
echo "  $ind:" >> $input/$folder/$ind.config
echo "$(<$input/$ind.low.depth)" >> $input/$folder/$ind.config
echo "$(<$input/$ind.high.depth)" >> $input/$folder/$ind.config

rm $input/$ind.low.depth
rm $input/$ind.high.depth

cd $input/$folder
mkdir $input/$folder/$ind
mkdir $input/$folder/$ind/vcf
cp /home/wlk579/2.0Indonesia_Bos_project/Indonesia_Bos_project/5.Analyses_King015/2.Fig.Hete_Roh/Hete_Casia/vcf/$ind.di_mono_allelic.bcf.gz $input/$folder/$ind/vcf/$ind.bcf.gz
$bcftools index $input/$folder/$ind/vcf/$ind.bcf.gz

snakemake --snakefile Hete.psmc.snakefile --configfile $ind.config -c20
mv $ind/heterozygosity.txt $ind.heterozygosity.txt

done < $batch

