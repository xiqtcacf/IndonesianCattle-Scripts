# run post QC

OUT=/isdata/hellergrp/wlk579/scripts/Data_PreProcess/bos_tmp/imputation/output
outdir=$OUT/imputed

echo "run post QC"

vcfin=$outdir/all.bcf
vcfphased=$outdir/all.phased.bcf
out=$outdir/all
stats=$out.imputed.sites.stats
bcftools +fill-tags $vcfin -- -t MAF |bcftools query -f "%ID\t%MAF\n" >$out.maf
# # merge r2 and maf together
echo -e "ID\tR2\tMAF" >$stats
paste $out.r2 <(awk '{print $2}' $out.maf) >>$stats
/isdata/hellergrp/wlk579/scripts/Data_PreProcess/bos_tmp/imputation/mafBinSmaller.plot-r2-vs-maf.R $stats $stats && echo plot done
# apply filters MAF>0.05 && R2 >0.95
# remove NaN sites in R2 file
#awk 'NR>1 && $3>0.05 && $2>0.95 && $2!="NaN"{print $1}' $stats >$out.maf0.05.r2.0.95.sites
#bcftools view -i ID=@$out.maf0.05.r2.0.95.sites -Ob -o $out.maf0.05.r2.0.95.bcf --threads 20 $vcfin && bcftools index -f $out.maf0.05.r2.0.95.bcf
# # for phased vcf
#bcftools view -i ID=@$out.maf0.05.r2.0.95.sites -Ob -o $out.phased.maf0.05.r2.0.95.bcf --threads 20 $vcfphased && bcftools index -f $out.phased.maf0.05.r2.0.95.bcf

#echo "done post QC"
