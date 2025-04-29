This folder contains scripts to identify and analyse regions of extreme banteng ancestry.

### 1. Calculating Ux
We calculated Ux separately for each group of cattle (Aceh, Pesisir, Pasundan, Madura, and EastAsianZebu) as follows:
a. We prepare a sliding window to base our SNP counts using the following commands in bash:
```
awk -v OFS='\t' '{print $1,$2}' bosTau9.fasta.fai | head -n 29 > bosTau9_genome_autosome.txt
bedtools makewindows -g bosTau9_genome_autosome.txt -w 50000 > bosTau9_genome_autosome_50K_windows.txt
``` 
b. We generate allele frequency files using `getAF.sh`. This would generate per-population AF files. The script was hard-coded with the path to the .bcf file and sample names required to calculate allele frequency per cattle group.

c. After getting the output of `getAF.sh`, we feed it as the first input argument of `calcUx.sh` to get the Ux of each population. For example, for Aceh group it would be running the following command in bash:
```
bash calcUabc.sh phased233_maf05NEWnaming.abcfiltered.filltags.dAFs.Aceh.txt 0.265 bosTau9_genome_autosome_50K_windows.bed
```
where 0.265 is the top 5% allele frequency of banteng-specific allele across the genome obtained from LOTER for Aceh. This value changes with different cattle group, i.e. 0.438 for East Asian Zebu, 0.638 for Madura, 0.375 for Pasundan, and 0.333 for Pesisir.

### 2. Detecting regions of high banteng ancestry
To identify regions of banteng ancestry, we intersect banteng proportion from three methods, HMMIX, Loter, and Ux, using the following command (for Aceh):
```
bedtools intersect -a ../ancestry.Aceh.hmmix_Archaic_0.9_propPer50K.Uabc_Loter_50K.txt -b bosTau9_genome_autosome_50K_windows.phased233_maf05NEWnaming.abcfiltered.filltags.dAFs.Aceh.0countFilt_0.05a_0.265b_0.9c.abcSNPcounts.txt -wb | cut -f -21,25,26 > ../ancestry.Aceh.hmmix_Archaic_0.9_propPer50K.Uabc_Loter_50K.UaLoterc265.txt
```
Then we ran the output file on `R` to identify windows with the top 5% banteng ancestry.

### 3. Overlying regions of high banteng ancestry with gene and QTL traits
To do so while also getting the gene and the QTL traits, we first intersect the output file with gene annotation and QTL with bosTau9 coordinates:
```
# QTL reading and QCing
zcat Animal_QTLdb_release52_cattleARS_UCD1.gff.gz | grep -v "^#" | awk -v OFS='\t' -F '\t' '{print $1,$4,$5,$9}' > Animal_QTLdb_release52_cattleARS_UCD1.bed
sed -i 's/Chr.30/Chr.X/' Animal_QTLdb_release52_cattleARS_UCD1.bed 
grep "^Chr.[0-9]" Animal_QTLdb_release52_cattleARS_UCD1.bed | sed 's/Chr./chr/' > Animal_QTLdb_release52_cattleARS_UCD1.chrPREFIX.bed
zcat Animal_QTLdb_release52_cattleARS_UCD1.gff.gz | grep -v "^#" | awk -F'\t' '{print $1,$4,$5,$9}' | tail -n+3 | sed 's/;/\t/g' | grep "Map_Type=Genome" | grep "Significance=Significant" | sed 's/Chr.30/Chr.X/' | grep "^Chr.[0-9]" | sed 's/Chr./chr/' | sed 's/ /\t/g' | awk -F'\t' -v OFS='\t' '{print $1,$2,$3,$4}' > Animal_QTLdb_release52_cattleARS_UCD1.chrPREFIX.significant.mapTypeGenome.bed
zcat Animal_QTLdb_release52_cattleARS_UCD1.gff.gz | grep -v "^#" | awk -F'\t' -v OFS='\t' '{print $1,$4,$5,$9}' | tail -n+3 | sed 's/;/\t/g' | grep "Map_Type=Genome" | grep "Significance=Significant" | sed 's/Chr.30/Chr.X/' | grep "^Chr.[0-9]" | sed 's/Chr./chr/' > Animal_QTLdb_release52_cattleARS_UCD1.chrPREFIX.significant.mapTypeGenome.txt
cut -f 1-4,9 Animal_QTLdb_release52_cattleARS_UCD1.chrPREFIX.significant.mapTypeGenome.txt > Animal_QTLdb_release52_cattleARS_UCD1.chrPREFIX.significant.mapTypeGenome.bed 
grep -v "QTL_ID=282894" Animal_QTLdb_release52_cattleARS_UCD1.chrPREFIX.significant.mapTypeGenome.bed > Animal_QTLdb_release52_cattleARS_UCD1.chrPREFIX.significant.mapTypeGenome.no282894.bed

# Overlaying with ancestry file per cattle group, the following example is for Pesisir
join -1 1 -2 1 BosTau9_refSeq_ucsc_chrNames.txt ancestry.Pesisir.hmmix_Archaic_0.9_propPer50K.Uabc_Loter_50K.txt | sed 's/ /\t/g' > pesisir_test.txt
cut -f 2- pesisir_test.txt > ancestry.Pesisir.hmmix_Archaic_0.9_propPer50K.Uabc_Loter_50K.chrPrefix.bed
bedtools intersect -a ancestry.Pesisir.hmmix_Archaic_0.9_propPer50K.Uabc_Loter_50K.chrPrefix.bed -b Animal_QTLdb_release52_cattleARS_UCD1.chrPREFIX.significant.mapTypeGenome.no282894.bed -loj > ancestry.Pesisir.hmmix_Archaic_0.9_propPer50K.Uabc_Loter_50K.chrPrefix.QTL.bed
bedtools merge -i ancestry.Pesisir.hmmix_Archaic_0.9_propPer50K.Uabc_Loter_50K.chrPrefix.QTL.bed -c 25,26 -o collapse > ancestry.Pesisir.chrPrefix.QTL.collapse.bed
bedtools intersect -a ancestry.Pesisir.hmmix_Archaic_0.9_propPer50K.Uabc_Loter_50K.chrPrefix.bed -b ancestry.Pesisir.chrPrefix.QTL.collapse.bed -loj | cut -f -21,25,26 > ancestry.Peisisir.hmmix_Archaic_0.9_propPer50K.Uabc_Loter_50K.chrPrefix.QTL.collapse.bed 
```

### 4. GO-enrichment analyses
GO-enrichment analyses were conducted on genes that were in top 5% banteng ancestry proportion as identified by LOTER.
