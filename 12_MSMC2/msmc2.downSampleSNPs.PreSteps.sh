####step1: generate samples targeted in : 2 for each pop, e.g. xx is sample name
bcftools view -s xx,xx,xx,xx filltag.phased.bosTau9_bos_variable_sites_nomultiallelics_noindels_10dp_3het.phased.bcf | bcftools view -i 'AC>0 & AC<AN' -o filltag.phased.bosTau9_bos_variable_sites_nomultiallelics_noindels_10dp_3het.phased.bcf -O b
bcftools index filltag.phased.bosTau9_bos_variable_sites_nomultiallelics_noindels_10dp_3het.phased.bcf
#### downsample snps to 10million
plink --bcf filltag.phased.bosTau9_bos_variable_sites_nomultiallelics_noindels_10dp_3het.phased.bcf --recode --out down_snp --double-id --allow-extra-chr
cut -f 2 down_snp.map > snps.map
shuf -n 10000000 snps.map > snps.subset10m.map
sed -i 's/:/\t/g' snps.subset10m.map

#### get 10million snps from targeted 24inds.bcf
bcftools view -R snps.subset10m.map filltag.phased.bosTau9_bos_variable_sites_nomultiallelics_noindels_10dp_3het.phased.bcf -o 10million.filltag.phased.bosTau9_bos_variable_sites_nomultiallelics_noindels_10dp_3het.phased.bcf -O b
### run snakemake file
