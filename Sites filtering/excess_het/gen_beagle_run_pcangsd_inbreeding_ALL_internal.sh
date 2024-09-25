BAM=/isdata/hellergrp/nzg134/asianbos/notes/internal_noDP_noFD_bamlist
OUT=/artemis-internal/casia16/asianbos/sites_QC/results
PCANGSD=/isdata/hellergrp/nzg134/programs/pcangsd/pcangsd.py
ANGSD=/isdata/hellergrp/rheller/software_scratch/angsdv0.929/angsd

$ANGSD -GL 2 -out ${OUT}/internal_noDP_noFD \
-nThreads 80 -doGlf 2 -doMajorMinor 1 \
-SNP_pval 1e-6 -doMaf 1 -minInd 75 -minQ 30 -minMapQ 30 -minMaf 0.05 \
-bam ${BAM}

python3 $PCANGSD -beagle ${OUT}/internal_noDP_noFD.beagle.gz \
-kinship -sites_save -inbreedSites -e 2 -o ${OUT}/internal_noDP_noFD_e2 \
-threads 80 > ${OUT}/internal_noDP_noFD_e2.log

