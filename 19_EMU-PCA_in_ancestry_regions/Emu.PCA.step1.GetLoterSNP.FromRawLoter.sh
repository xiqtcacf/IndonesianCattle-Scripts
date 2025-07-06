#!/bin/bash
#pop=$1 ### give a list of individuals for specific population e.g. Madura Pesisir Here should be All of Pops:allPop

export DIR=/home/wlk579/2.0Indonesia_Bos_project/Indonesia_Bos_project/5.Analyses_King015/emu_PCA/loter/SNP_position/bed_Ancestry/
cd $DIR

module load bedops
module load R/4.0.3
bedtools=/home/wlk579/2.0Indonesia_Bos_project/Indonesia_Bos_project/5.Analyses_King015/hmmix0.6.9_AllResultsPlot/heatmap/bedtools

source /maps/projects/seqafrica/apps/modules/activate.sh
module load plink
bcftools=/projects/popgen/apps/Modules/software/bcftools/1.16/bin/bcftools

angsd=/home/wlk579/Server_bos/apps/angsd/angsd

####Step1 get correct format results from raw loter
while read Chr ; do
zless by_Chr_input_loter_raw/internal_eastasian_ausie_${Chr}.vcf.gz | grep -v '#' | cut -f1,2 > internal_eastasian_ausie_${Chr}.vcf.gz.chr_pos
cat header_chr_pos internal_eastasian_ausie_${Chr}.vcf.gz.chr_pos > test1
rm  internal_eastasian_ausie_${Chr}.vcf.gz.chr_pos
paste 188D_names by_Chr_input_loter_raw/internal_eastasian_ausie_nophaseC_${Chr}_anc.txt > test2
Rscript scripts/Rawloter.R test2
paste test1 test.txt > format.internal_eastasian_ausie_nophaseC_${Chr}_anc.txt
rm test1
rm test2
rm test.txt
done < chr29_list
 sed -i '1d' format.internal_eastasian_ausie_nophaseC_${Chr}_anc.txt ###except the first chromosome
 cat format.internal_eastasian_ausie_nophaseC_${Chr}_anc.txt > format.internal_eastasian_ausie_nophaseC_all29Chr_anc.txt

####Step2 get bed file per individual per chromosome, Zebu Ancestry:00 Banteng Ancestry:11
while read Chr ; do
      Rscript scripts/get_column_1_0_position.R format.internal_eastasian_ausie_nophaseC_${Chr}_anc.txt
      cd bed_Ancestry
      for i in *.txt; do cat $i | sed 's/$/\t'$i'/' > $Chr.$i.new; done
      rm *.txt
      sed -i '1d' $Chr.$i.new
done < chr29_list

####Step3 get Ancestry SNP posotion per individual for 29 chrs
while read sample ; do
      cat NC_037328.1.ZebuAncestry.$sample.txt.new \
          NC_037329.1.ZebuAncestry.$sample.txt.new \
          NC_037330.1.ZebuAncestry.$sample.txt.new \
          NC_037331.1.ZebuAncestry.$sample.txt.new \
          NC_037332.1.ZebuAncestry.$sample.txt.new \
          NC_037333.1.ZebuAncestry.$sample.txt.new \
          NC_037334.1.ZebuAncestry.$sample.txt.new \
          NC_037335.1.ZebuAncestry.$sample.txt.new \
          NC_037336.1.ZebuAncestry.$sample.txt.new \
          NC_037337.1.ZebuAncestry.$sample.txt.new \
          NC_037338.1.ZebuAncestry.$sample.txt.new \
          NC_037339.1.ZebuAncestry.$sample.txt.new \
          NC_037340.1.ZebuAncestry.$sample.txt.new \
          NC_037341.1.ZebuAncestry.$sample.txt.new \
          NC_037342.1.ZebuAncestry.$sample.txt.new \
          NC_037343.1.ZebuAncestry.$sample.txt.new \
          NC_037344.1.ZebuAncestry.$sample.txt.new \
          NC_037345.1.ZebuAncestry.$sample.txt.new \
          NC_037346.1.ZebuAncestry.$sample.txt.new \
          NC_037347.1.ZebuAncestry.$sample.txt.new \
          NC_037348.1.ZebuAncestry.$sample.txt.new \
          NC_037349.1.ZebuAncestry.$sample.txt.new \
          NC_037350.1.ZebuAncestry.$sample.txt.new \
          NC_037351.1.ZebuAncestry.$sample.txt.new \
          NC_037352.1.ZebuAncestry.$sample.txt.new \
          NC_037353.1.ZebuAncestry.$sample.txt.new \
          NC_037354.1.ZebuAncestry.$sample.txt.new \
          NC_037355.1.ZebuAncestry.$sample.txt.new \
          NC_037356.1.ZebuAncestry.$sample.txt.new > ../bed_Ancestry_perInd/loter.ZebuAncestry.$sample

      cat NC_037328.1.BantengAncestry.$sample.txt.new \
          NC_037329.1.BantengAncestry.$sample.txt.new \
          NC_037330.1.BantengAncestry.$sample.txt.new \
          NC_037331.1.BantengAncestry.$sample.txt.new \
          NC_037332.1.BantengAncestry.$sample.txt.new \
          NC_037333.1.BantengAncestry.$sample.txt.new \
          NC_037334.1.BantengAncestry.$sample.txt.new \
          NC_037335.1.BantengAncestry.$sample.txt.new \
          NC_037336.1.BantengAncestry.$sample.txt.new \
          NC_037337.1.BantengAncestry.$sample.txt.new \
          NC_037338.1.BantengAncestry.$sample.txt.new \
          NC_037339.1.BantengAncestry.$sample.txt.new \
          NC_037340.1.BantengAncestry.$sample.txt.new \
          NC_037341.1.BantengAncestry.$sample.txt.new \
          NC_037342.1.BantengAncestry.$sample.txt.new \
          NC_037343.1.BantengAncestry.$sample.txt.new \
          NC_037344.1.BantengAncestry.$sample.txt.new \
          NC_037345.1.BantengAncestry.$sample.txt.new \
          NC_037346.1.BantengAncestry.$sample.txt.new \
          NC_037347.1.BantengAncestry.$sample.txt.new \
          NC_037348.1.BantengAncestry.$sample.txt.new \
          NC_037349.1.BantengAncestry.$sample.txt.new \
          NC_037350.1.BantengAncestry.$sample.txt.new \
          NC_037351.1.BantengAncestry.$sample.txt.new \
          NC_037352.1.BantengAncestry.$sample.txt.new \
          NC_037353.1.BantengAncestry.$sample.txt.new \
          NC_037354.1.BantengAncestry.$sample.txt.new \
          NC_037355.1.BantengAncestry.$sample.txt.new \
          NC_037356.1.BantengAncestry.$sample.txt.new > ../bed_Ancestry_perInd/loter.BantengAncestry.$sample
done < ../188_names

