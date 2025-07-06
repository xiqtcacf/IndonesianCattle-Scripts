args <- commandArgs(trailing=T)

inf <- args[1]
outf <- args[2]

df = read.table(inf, h=T)
abc = df[order(df$king, decreasing=T),]
abc1=subset(abc, king>0.05)
maxobs<-max(.5,max(abc1$king))
bitmap(outf, res=300,height=.15*nrow(abc1), width=6)
par(mar=c(1,12,1,1))
barplot(abc1$king, names.arg=abc1$id, las=1,horiz=T, cex.names=0.75, xlim=c(0,maxobs), xaxt='n',border=NA)
ticks=seq(0,maxobs,.1)
axis(1,at=ticks, labels=ticks,pos=0)
dev.off()
