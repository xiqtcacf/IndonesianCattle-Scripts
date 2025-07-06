##plink should be good population name! thin and maf! check also snp's name column
##make folder m_0 until m_10

### in bash

nohup plink2 --bfile /isdata/hellergrp/asianBos/imputed_analysis/treemix/treemix_rootWaterbuffalo/merge.261inds.waterbuffaloRefPlusBos --keep treemix_ind_2col_withGaur_mhss.txt --make-bed --out 155Gaur_plusWaterbuff_maf05 --allow-extra-chr --thin 0.01 --maf 0.05 &

###this steps might be different if you got diiferent plink (this one is from Xi)

awk -v OFS="\t" '{print $1,$2,$1}' 155Gaur_plusWaterbuff_maf05.fam > 155Gaur_plusWaterbuff_maf05.pop
awk -v OFS="\t" '{print $2,$1}' 155Gaur_plusWaterbuff_maf05.fam > treemixpop

###in bash
plink --bfile 155Gaur_plusWaterbuff_maf05 --freq --missing --within 155Gaur_plusWaterbuff_maf05.pop --allow-extra-chr --out 155Gaur_plusWaterbuff_maf05

gzip 155Gaur_plusWaterbuff_maf05.frq.strat

/isdata/hellergrp/asianBos/imputed_analysis/treemix/plink2treemix.py 155Gaur_plusWaterbuff_maf05.frq.strat.gz treemix.frq.gz

### then run treemix, careful to choose the roots!

