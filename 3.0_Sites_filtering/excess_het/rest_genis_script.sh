workdir="/artemis-internal/casia16/asianbos/sites_QC"
BEDTOOLS="/isdata/hellergrp/nzg134/programs/bedtools2/bin/bedtools"
${BEDTOOLS} slop -b 10000 -i ${workdir}/results/internal_noDP_noFD_e2.badsites.bed \
-g ${workdir}/note/genome_sizes.file > ${workdir}/results/internal_noDP_noFD_e2.temp_bed

${BEDTOOLS} merge -i ${workdir}/results/internal_noDP_noFD_e2.temp_bed \
> ${workdir}/results/internal_noDP_noFD_e2.badbed

Rscript summaryplot.R ${workdir}/results/internal_noDP_noFD_e2 ${workdir}/results/excessHet \
${workdir}/results/internal_noDP_noFD_e2.badbed ${workdir}/note/bostau_len_chr.txt -0.9 24

