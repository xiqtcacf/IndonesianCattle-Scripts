# Scripts for
## Wang, X., Nursyifa, C., ... Heller, R. (2025). The genetic diversity of Indonesian cattle has been shaped by multiple introductions and adaptive introgression.

### 01_Mapping
A pipeline for mapping and post-mapping filtering using PALEOMIX BAM pipeline.

### 02_Sample_filtering
A pipeline designed to remove problematic samples prior to downstream analyses, based on

statistics from MultipleQC;

exteremely high heterozygosity;

relateness between pairwise samples: King calculated from global 2D-sfs;

### 03_Sites_filtering
A pipeline designed to avoid biases from low-quality mapping, based on

03_1_Reference_repeats: detection of problematic repeat regions of the reference genomes;

03_2_mappability: detection of problematic regions of the reference genomes based on mappability;

03_3_Depth: sites that showed unusual depth;

03_4_excess_het: sites excess heterozygosity after mapping.

### 04_GenotypeCalling
A pipeline designed to call genotypes, the origional pipeline was developed for the 'a1kg project'.

### 05_Imputation_Phasing
Scripts to imputate and phase genotypes.

### 06_Jonas_Haplonet
Script to infer principal component analysis (PCA) and to estimate admixture proportions for each individual using HaploNet.

### 07_GenomeWide_Heterozygosity
Script to calculate genome-wide heterozygosity per individual.

### 08_Runs_of_homozygosity
Script to infer runs of homozygosity (ROH) per individual;

Script to calculate genome-wide heterozygosity without ROH per individual.

### 09_Pairwise_globalFst
Script to infer Pairwise global Fst between each pair of populations using PLINK2.

### 10_Treemix
Script to infer the population tree assuming different numbers of admixture events using TreeMix.

### 11_Dstatistics
Script to infer the ancient admixture events, we calculated D statistics (ABBA-BABA) using the R package ADMIXTOOLS2.

### 12_MSMC2
Script to infer divergence time between between zebu, banteng and Bali cattle using MSMC2.

### 13_Loter
Script to estimate genome-wide admixture proportions using local ancestry inference with LOTER;

Script to estimate genome-wide admixture proportions using local ancestry inference on the X chromosome with LOTER;

Script to estimate genome-wide admixture proportions using local ancestry inference with LOTER, based on data mapped to Banteng reference during revisions.

### 14_Loter_Heterozygosity_in_ancestry_States
Script to calculate heterozygosity in each of the three ancestry states along the genome in admixed breeds: tracts homozygous for zebu ancestry, tracts heterozygous for zebu/banteng ancestry, and tracts homozygous for banteng ancestry.

### 15_F4_Ratio
Script to estimated genome-wide admixture proportions using local ancestry inference with F4 ratio.

### 16_Hmmix
Script to identify regions in the genome that were introgressed from highly divergent source populations into the admixed breeds using Hmmix;

Script to calculate the genetic similarity of each introgressed region identified by Hmmix with two different potential source lineages, Javan banteng, and gaur;

Script to infer genetic distance (ibs) shared between pairs of individuals within introgressed regions identified by Hmmix.

### 17_Ux
Script to to identify putatively adaptively introgressed regions using a metric Ux.

### 18_Thomas_Correlation_Ancestry_GenomicFeatures
Script to calculate the correlation between the mean banteng ancestry and three genomic features: recombination rate, conservation score, and coding region density in Madura cattle.

### 19_EMU-PCA_in_ancestry_regions
Script to infer PCA for zebu-specific and banteng-specific ancestry regions identified by LOTER using EMU;

Script to infer PCA for zebu-specific and banteng-specific ancestry regions identified by Hmmix using EMU.

### 20_AncestryHMM-AdmixtureTime
Script to inferred the introgression time using AncestryHMM.

### 21_Sabhrina_Regions_Of_HighBantengAncestry
Script to identify regions of high banteng ancestry (top 5%) within breeds using 50 kb non-overlapping windows across the genome with three complementary methods of local ancestry inference (LOTER, Hmmix, and Ux);

Script to performed GO-enrichment analyses for the outlier gene set from each breed.

### 22_HaplotypeStructure_AZIP
Script to investigate the haplotype structure in ASIP gene using Haplostrips.
