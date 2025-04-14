##### HaploNet analyses #####

### haplonet train (clustering) ###
for c in {28..56}
do
	haplonet train --vcf phased231_nogaur_maf05__NC_0373${c}.1.vcf.gz --threads 32 --batch 64 --out asian.bos.NC_0373${c}.1
	echo asian.bos.NC_0373${c}.1.loglike.npy >> asian.bos.filelist
done


### haplonet pca ###
haplonet pca --filelist asian.bos.filelist --threads 32 --out asian.bos.haplonet.pca


### haplonet admix (ancestry estimation) ###
for k in {3..12}
do
	for s in {1..50}
	do
		haplonet admix --filelist asian.bos.filelist --K $k --seed $s --threads 32 --out asian.bos.K${k}.s${s} 
	done
done
