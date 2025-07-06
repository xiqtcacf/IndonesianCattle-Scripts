#!/bin/bash

beagle3=/isdata/hellergrp/wlk579/software/BioScripts/beagle3.jar
VCF=/isdata/hellergrp/wlk579/bos_tmp/imputation/input/bosTau9_bos_variable_sites_nomultiallelics_noindels.bcf.gz
OUT=/isdata/hellergrp/wlk579/bos_tmp/imputation/output
outdir=$OUT/imputed
indir=$OUT/vcf
mkdir -p $indir $outdir

runBeagle3_by_region() {
    # show verbose log
    rg="$1"
    chr=$(perl -e '@t=split(":", $ARGV[0]);print $t[0];' $rg)
    echo "start running beagle3 for $rg"
    # convert PL tag to GL in beagle format 
    VCF=$2
    OUT=$3
    beagle3=$4

    indir=$OUT/vcf
    mkdir -p $indir
    bname=$indir/`basename $VCF`.$rg
    bcftools view -r $rg -O u $VCF | bcftools +tag2tag -Ov -o ${bname}.vcf -- -r --pl-to-gl && vcftools --vcf ${bname}.vcf --out ${bname} --BEAGLE-GL --chr $chr \
      && rm -f ${bname}.vcf && echo 'convert to BEAGLE-GL done'

    # run imputation
    outdir=$OUT/imputed
    mkdir -p $outdir
    out=$outdir/$rg
    input=${bname}.BEAGLE.GL
    java -Xss5m -Xmx20g -Djava.io.tmpdir=$outdir -jar $beagle3 seed=-99999 omitprefix=true like=$input out=$out

    echo "done running beagle3 for $rg"
}
export -f runBeagle3_by_region

threads=20
regions=$(bcftools index -s $VCF |perl -slane '$size=10000000; $n=sprintf("%f", $F[1]/$size);for($i=0;$i<$n;$i++){$s=$i*$size+1;$e=($i+1)*$size;$e = ($e > $F[1]) ? $F[1] : $e;printf("%s:%d-%d\n", $F[0], $s, $e);}')
echo $regions | tr " " "\n" | parallel -j $threads -k -I {} runBeagle3_by_region {} $VCF $OUT $beagle3
wait
echo "runBeagle3_by_region done"
