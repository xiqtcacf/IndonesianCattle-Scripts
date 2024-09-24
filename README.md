# Scripts for
## Wang, X., Nursyifa, C., ... Heller, R. (2024). Multiple origins and extensive bovine introgression shaped the highly diverse Indonesian cattle

### Mapping, Sample filtering, and Sites filtering
A pipeline designed for mapping and post-mapping filtering of waterbuck.

A pipeline designed to remove problematic samples prior to downstream analyses, based on

mapping statistics;

statistics from MultipleQC;

exteremely high heterozygosity;

error rates;

relateness between pairwise samples: King calculated from global 2D-sfs;

relateness within each location by NgsRelate.

A pipeline designed to avoid biases from low-quality mapping, based on 

3.1.detection of problematic regions of the reference genomes:

3.1.1.mappability;

3.1.2.repeats;

3.1.3.sex-linked chromosomes or scaffolds;

3.2.sites that showed unusual depth; 

3.3.sites excess heterozygosity after mapping. 


### Genotype calling, imputation and phasing
Script to call genotypes;

Script to imputate and phase genotypes.

### Population structure
Script to infer principal component analysis (PCA) using HaploNet;

Script to estimate admixture proportions for each individual using HaploNet;

Script to infer Pairwise global Fst between each pair of populations;

### Genetic diversity and inbreeding
Script to calculate genome-wide heterozygosity;

Script to infer runs of homozygosity (ROH);

Script to calculate genome-wide heterozygosity without ROH;

### Admixture and origin of Indonesian cattle
Script to infer the population tree assuming different numbers of admixture events using TreeMix;

Script to infer the ancient admixture events, we calculated D statistics (ABBA-BABA) using the R package ADMIXTOOLS2;

Script to estimate genome-wide admixture proportions using local ancestry inference with LOTER;

Script to estimated genome-wide admixture proportions using local ancestry inference with F4-derived admixture ratio;

Script to identify regions in the genome that were introgressed from highly divergent source populations into the admixed breeds using Hmmix;

Script to calculate the genetic similarity of each introgressed region identified by Hmmix with two different potential source lineages, Javan banteng, and gaur;

Script to infer genetic distance shared between pairs of individuals within introgressed regions identified by Hmmix;

Script to inferred the introgression time using AncestryHMM.

Script to infer PCA for zebu-specific and banteng-specific ancestry regions identified by LOTER using EMU;

Script to infer PCA for zebu-specific and banteng-specific ancestry regions identified by Hmmix using EMU;

Script to calculate heterozygosity in each of the three ancestry states along the genome in admixed breeds: tracts homozygous for zebu ancestry, tracts heterozygous for zebu/banteng ancestry, and tracts homozygous for banteng ancestry;

### Genomic landscape of ancestry
Script to infer phylogenetic tree for mitochondrial DNA (mtDNA);

Script to infer phylogenetic tree for Y-chromosome (Ychr);

Script to estimate genome-wide admixture proportions using local ancestry inference on the X chromosome with LOTER;

Script to calculate the correlation between the mean banteng ancestry and three genomic features: recombination rate, conservation score, and coding region density in Madura cattle;

### Regions of high banteng ancestry in admixed cattle
Script to to identify putatively adaptively introgressed regions using a metric Ux;

Script to identify regions of high banteng ancestry (top 5%) within breeds using 50 kb non-overlapping windows across the genome with three complementary methods of local ancestry inference (LOTER, Hmmix, and Ux);

Script to performed GO-enrichment analyses for the outlier gene set from each breed;

### Shared outlier regions and convergent adaptive introgression
Script to tabulated the overlap between genes in the top 5% windows with banteng ancestry across the admixed breeds;

Script to investigate the haplotype structure in ASIP gene using Haplostrips;

### Divergence time
Script to infer divergence time between between zebu, banteng and Bali cattle using MSMC2;
