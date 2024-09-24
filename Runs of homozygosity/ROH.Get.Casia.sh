num_procs=15
num_jobs="\j"

workdir="/isdata/hellergrp/asianBos/reRun_king015/roh/"
plinkFile="/isdata/hellergrp/asianBos/reRun_king015/data/phased231_nogaur_maf05_newNaming.vcf.gz"

while read sample; do

   while (( ${num_jobs@P} >= num_procs )); do
      #echo "${num_jobs@P} "
    wait -n
  done

#echo -e "$sample\t$sample" > ${workdir}/ROH/het5_min500kb_206samples/${sample}.list
nice -n10 plink --vcf ${plinkFile} \
--keep sample_data/${sample} \
--homozyg-window-het 5 --chr-set 29 --homozyg-kb 500 \
--geno 0 \
--out ${workdir}/het5_min500kb_231samples/${sample} \
--double-id  \
--allow-extra-chr &
done <  <(cat /isdata/hellergrp/asianBos/reRun_king015/data/231_noGaurNEWnaming.txt)

