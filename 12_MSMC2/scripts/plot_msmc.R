require(ggplot2)

a <- 1:10
b <- a+log10(5)

args <- commandArgs(trailing=T)

infile <- args[1]
outfile <- args[2]

df <- read.table(infile)

df$V2[which(df$V2 == 0)] <- df$V2[which(df$V2==0)+1]/2
p <- ggplot(df, aes(V2, V3, col=V1, group=V1, linetype=V1)) + geom_step() + scale_x_log10(breaks = 10^a, minor_breaks=10^b) + scale_y_log10(breaks=10^a, minor_breaks=10^b) + scale_linetype_manual(values=c(rep(1:2, times=50))) + theme_bw() + theme(legend.position=c(0.05, 0.95), legend.justification=c("left", "top"), legend.title=element_blank()) + labs(x="Years ago", y="Effective population size") ; ggsave(outfile, plot = p, dpi=500, width=7, height=5)
