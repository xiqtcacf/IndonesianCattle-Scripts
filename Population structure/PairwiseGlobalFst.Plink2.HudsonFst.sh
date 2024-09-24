#!/bin/bash

export DIR=/isdata/hellergrp/wlk579/Indonesia_Bos_project/5.Analyses_King015/Fst_globle
cd $DIR

bcf=/isdata/hellergrp/wlk579/Indonesia_Bos_project/5.Analyses_King015/Fst_globle/314inds.imputated.BosTau9.bcf
plink2=/isdata/hellergrp/wlk579/software/plink2/plink2

### must use the plink2 that I downloaded 2022.03 version
bcftools view -S purepop.inds $bcf -o purepop.imputated.BosTau9.bcf -O b
bcftools index purepop.imputated.BosTau9.bcf
$plink2 --bcf purepop.imputated.BosTau9.bcf --make-bed --out purepop.imputated.BosTau9 --allow-extra-chr
$plink2 --bfile breeds.imputated.BosTau9 --fst breeds.new.new.txt --family breeds.new.new.txt --out breeds.imputated.BosTau9.fst --allow-extra-chr

