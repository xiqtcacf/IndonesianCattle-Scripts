###Only Chr1
####only compare individual of 'LIB112407_Banteng_85B_Texas'
library(vcfppR)
######1X VS Ture (34.86X)
compare <- vcfcomp(test = '1X.vcf.gz', truth = 'real.bcf.gz', setid=T, samples='LIB112407_Banteng_85B_Texas', region="NC_037328.1", formats = c("GT","GT"), bin=c(0,0.01,0.05,0.1,0.2,0.3,0.4,0.5), af='MAF.real.bcf.NC_037328.1.txt')
######plotting
par(mar=c(5,5,2,2), cex.lab = 2)
vcfplot(compare, col = 2,cex = 2, lwd = 3, type = "b")
dev.off()

######5X VS Ture
compare <- vcfcomp(test = '5X.vcf.gz', truth = 'real.bcf.gz', setid=T, samples='LIB112407_Banteng_85B_Texas', region="NC_037328.1", formats = c("GT","GT"), bin=c(0,0.01,0.05,0.1,0.2,0.3,0.4,0.5), af='MAF.real.bcf.NC_037328.1.txt')
######plotting
par(mar=c(5,5,2,2), cex.lab = 2)
vcfplot(compare, col = 2,cex = 2, lwd = 3, type = "b")
dev.off()

######10X VS Ture
compare <- vcfcomp(test = '10X.vcf.gz', truth = 'real.bcf.gz', setid=T, samples='LIB112407_Banteng_85B_Texas', region="NC_037328.1", formats = c("GT","GT"), bin=c(0,0.01,0.05,0.1,0.2,0.3,0.4,0.5), af='MAF.real.bcf.NC_037328.1.txt')
######plotting
par(mar=c(5,5,2,2), cex.lab = 2)
vcfplot(compare, col = 2,cex = 2, lwd = 3, type = "b")
dev.off()
