args <- commandArgs(TRUE)

print(args[1])
setwd(args[1])
group<-args[1]
df=read.table("het_4cols.txt")
myhet=c(df[,2],df[,4],df[,6],df[,8])
mydat=data.frame(anc_reg=factor(rep(c("homo_zebu","homo_banteng","mixed_ancestry","allRegions"),each=nrow(myhet)),levels=c("homo_zebu","homo_banteng","mixed_ancestry","allRegions")),het=myhet)

mydat=data.frame(anc_reg=factor(rep(c("zebu","banteng","mixed_ancestry","all"),each=nrow(df)),levels=c("zebu","banteng","mixed_ancestry","all")),het=myhet)



means<-aggregate(het~ anc_reg,mydat,median)

print(paste("heterozygosity_byAnc_PSMC_correctApproach",group,".png",sep=""))
png(paste("heterozygosity_byAnc_PSMC_correctApproach",group,".png",sep=""),width=700,height=450)
boxplot(het ~ anc_reg,data=mydat,col="darkblue",medcol="grey", xlab="Ancestry segments", ylab="Heterozygosity")
points(1:4,means$het,col="white")
#text(1:4,means$het+0.00025,labels=round(means$het,digit=5), col="white")
text((1:4)+0.53,means$het+0.00001,labels=paste("",round(means$het,digit=4),sep=""), col="darkblue",cex=0.7)
dev.off()

