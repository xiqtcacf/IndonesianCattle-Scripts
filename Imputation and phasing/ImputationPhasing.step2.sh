#!/bin/bash

bcftools() {
	/isdata/hellergrp/wlk579/bos_tmp/imputation/bcftools-1.9/bcftools "$@"
}
export -f bcftools
export BCFTOOLS_PLUGINS="/isdata/hellergrp/wlk579/bos_tmp/imputation/bcftools-1.9/plugins"

beagle3=/isdata/hellergrp/wlk579/software/BioScripts/beagle3.jar
VCF=/isdata/hellergrp/wlk579/bos_tmp/imputation/input/bosTau9_bos_variable_sites_nomultiallelics_noindels.bcf.gz
OUT=/isdata/hellergrp/wlk579/bos_tmp/imputation/output
outdir=$OUT/imputed
indir=$OUT/vcf
mkdir -p $indir $outdir

#runBeagle3_by_region() {
#    # show verbose log
#    rg="$1"
#    chr=$(perl -e '@t=split(":", $ARGV[0]);print $t[0];' $rg)
#    echo "start running beagle3 for $rg"
#    # convert PL tag to GL in beagle format 
#    VCF=$2
#    OUT=$3
#    beagle3=$4
#
#    indir=$OUT/vcf
#    bname=$indir/`basename $VCF`.$rg
#    bcftools view -r $rg -O u $VCF | bcftools +tag2tag -Ov -o ${bname}.vcf -- -r --pl-to-gl && vcftools --vcf ${bname}.vcf --out ${bname} --BEAGLE-GL --chr $chr \
#      && rm -f ${bname}.vcf && echo 'convert to BEAGLE-GL done'
#
#    # run imputation
#    outdir=$OUT/imputed
#    out=$outdir/$rg
#    input=${bname}.BEAGLE.GL
#    /usr/bin/time -v java -Xss5m -Xmx20g -Djava.io.tmpdir=$outdir -jar $beagle3 seed=-99999 like=$input out=$out
#
#    echo "done running beagle3 for $rg"
#}
#export -f runBeagle3_by_region

#threads=20
regions=$(bcftools index -s $VCF |perl -slane '$size=10000000; $n=sprintf("%f", $F[1]/$size);for($i=0;$i<$n;$i++){$s=$i*$size+1;$e=($i+1)*$size;$e = ($e > $F[1]) ? $F[1] : $e;printf("%s:%d-%d\n", $F[0], $s, $e);}')
# echo $regions | tr " " "\n" | parallel -j $threads -k -I {} runBeagle3_by_region {} $VCF $OUT $beagle3
#wait
#echo "runBeagle3_by_region done"


echo "start concating files by chromos"
chroms=$(bcftools index -s $VCF | cut -f1)      # not array just string sep by space

for chr in $chroms;do
  {
    rgs=$(echo $regions | tr ' ' '\n' | grep -w $chr)
    out=$outdir/$chr
    #n=1
    #for i in $rgs;do
    #  gps=$outdir/`basename $VCF`.$i.BEAGLE.GL.gprobs.gz
    #  if [ $n == 1 ]; then
    #    zcat $gps ;
    #  else
    #    zcat $gps | awk "NR>1";
    #  fi
    #  n=$((n + 1))
    #done | gzip -c >$out.gprobs.gz
    #zcat $out.gprobs.gz | java -jar /isdata/hellergrp/wlk579/software/BioScripts/gprobs2beagle.jar 0.9 -1 | gzip -c >$out.bgl.gz && \
    #  zcat $out.gprobs.gz | awk 'NR>1{split($1,a,":");print $1,a[2],$2,$3}' >$out.bgl.sites && \
    #  java -jar /isdata/hellergrp/wlk579/software/BioScripts/beagle2vcf.jar $chr $out.bgl.sites $out.bgl.gz -1 |bgzip -c >$out.vcf.gz && \
    #  bcftools index -f $out.vcf.gz && echo "convert $chr to vcf done"
    #
    #for i in $rgs;do
    #  cat $outdir/`basename $VCF`.$i.BEAGLE.GL.r2
    #done >$out.r2

    n=1
    for i in $rgs;do
      phs=$outdir/`basename $VCF`.$i.BEAGLE.GL.phased.gz
      if [ -f $phs ];then
	      if [ $n == 1 ]; then
		zcat $phs ;
	      else
		zcat $phs | awk "NR>1";
	      fi
      fi
      n=$((n + 1))
    done | gzip -c >$out.phased.gz
    /isdata/hellergrp/wlk579/software/BioScripts/beagle_phased_to_hap_sample.py $chr $out.phased.gz $out.bgl.sites $out && bcftools convert --hapsample2vcf $out |bcftools annotate -I +'%CHROM:%POS' -Oz -o $out.phased.vcf.gz && bcftools index -f $out.phased.vcf.gz
  } &
done
wait 
echo "done concating files by chromos"

echo "start concating all"
out=$outdir/all
# bcftools concat --threads 20 -Ob -o $out.bcf `for i in $chroms; do echo $outdir/$i.vcf.gz;done` && bcftools index -f $out.bcf &
bcftools concat --threads 20 `for i in $chroms; do echo $outdir/$i.phased.vcf.gz;done` | bcftools annotate -I +'%CHROM:%POS' -Ob -o $out.phased.bcf --threads 20  && bcftools index -f $out.phased.bcf &
for i in $chroms;do
    cat $outdir/$i.r2
done >$out.r2
wait
echo "done concating all"

# run post QC

echo "run post QC"

vcfin=$outdir/all.bcf
vcfphased=$outdir/all.phased.bcf
out=$outdir/all
stats=$out.imputed.sites.stats
bcftools +fill-tags $vcfin -- -t MAF |bcftools query -f "%ID\t%MAF\n" >$out.maf
# # merge r2 and maf together
echo -e "ID\tR2\tMAF" >$stats
paste $out.r2 <(awk '{print $2}' $out.maf) >>$stats
/isdata/hellergrp/wlk579/software/BioScripts/plot-r2-vs-maf.R $stats $stats && echo plot done
# apply filters MAF>0.05 && R2 >0.95
# remove NaN sites in R2 file
awk 'NR>1 && $3>0.05 && $2>0.95 && $2!="NaN"{print $1}' $stats >$out.maf0.05.r2.0.95.sites
bcftools view -i ID=@$out.maf0.05.r2.0.95.sites -Ob -o $out.maf0.05.r2.0.95.bcf --threads 20 $vcfin && bcftools index -f $out.maf0.05.r2.0.95.bcf
# # for phased vcf
bcftools view -i ID=@$out.maf0.05.r2.0.95.sites -Ob -o $out.phased.maf0.05.r2.0.95.bcf --threads 20 $vcfphased && bcftools index -f $out.phased.maf0.05.r2.0.95.bcf

echo "done post QC"
