BEDTOOLS=/isdata/hellergrp/nzg134/programs/bedtools2/bin/bedtools

auto=/artemis-internal/casia16/asianbos/sites_QC/note/bostau_len_chrY.txt
mapA=/isdata/hellergrp/nzg134/radseq/output/mappability/buffalo/GEM/mappability_K100_E2.bed
rep=/artemis-internal/casia16/asianbos/sites_QC/results/repeatMasked.bed
het=/artemis-internal/casia16/asianbos/sites_QC/results/internal_noDP_noFD_e2.badbed
dep=/isdata/hellergrp/asianBos/analysis_asianBos_Xi/3.site_filtering/depth/results/bed/all_keep.bed
mapY=/artemis-internal/casia16/asianbos/sites_QC/results/mappabilityYchr_m1_k100_e2.bed

#goodbed_het
${BEDTOOLS} subtract -a $auto -b $het > /artemis-internal/casia16/asianbos/sites_QC/results/het_goodregion.bed

cat $mapA $mapY > /artemis-internal/casia16/asianbos/sites_QC/results/mappability_all.bed
map=/artemis-internal/casia16/asianbos/sites_QC/results/mappability_all.bed


#mkdir beds
cd beds

$BEDTOOLS intersect -a $auto -b $rep > autosome_rep.bed
$BEDTOOLS subtract -a autosome_rep.bed -b $het > autosome_rep_het.bed
$BEDTOOLS intersect -a autosome_rep_het.bed -b $dep > autosome_rep_het_dep.bed
$BEDTOOLS intersect -a autosome_rep_het_dep.bed -b $map > autosome_rep_het_dep_map.bed

awk '{print $1"\t"$2+1"\t"$3}'  autosome_rep_het_dep_map.bed > autosome_rep_het_dep_map.regions
$ANGSD sites index autosome_rep_het_dep_map.regions


