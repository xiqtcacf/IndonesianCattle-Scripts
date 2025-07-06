## nk = read.table("Nubian.Kordofan.coalest.txt")
## ns = read.table("Nubian.SouthAfrica.coalest.txt")
## bitmap("msmc2_coalescence.png", res=400)
## plot(nk$V3, nk$V4, type='s', col='blue', xlim=c(0, 500000), ylim=c(0,1), xlab="Years ago", ylab="Relative Coalescence")
## lines(ns$V3, ns$V4, col='red', type='s', xlim=c(0, 500000), ylim=c(0,1))
## legend("bottomright", legend=c("Nubian-Kordofan", "Nubian-SouthAfrica"), fill=c("blue","red"))
## dev.off()

args <- commandArgs(trailing=T)
##
infiles <- args[1]
outpng <- args[2]

a = scan(infiles, what='character')
b = lapply(a, read.table)
n = gsub(".coalest.txt", "", basename(a))

ltys = rep(1:3, 100)[1:length(n)]
cols = rep(1:5, 100)[1:length(n)]
pchs = rep(1:4, 100)[1:length(n)]

bitmap(outpng, res=400)
with(b[[1]],
     plot(V3,V4, type='s',
          col=cols[1], pch=pchs[1],
          xlim=c(0, 500000), ylim=c(0,1),
          xlab="Years ago", ylab="Relative Coalescence"))
with(b[[1]],
     points(V3,V4,
          col=cols[1], pch=pchs[1],
          xlim=c(0, 500000), ylim=c(0,1),
          xlab="Years ago", ylab="Relative Coalescence"))


for (i in 2:length(n)){
    with(b[[i]],
         lines(V3,V4, type='s',
               col=cols[i], lty=ltys[i],
               xlim=c(0, 500000), ylim=c(0,1),
               xlab="Years ago", ylab="Relative Coalescence"))
    with(b[[i]],
         points(V3,V4,
               col=cols[i], pch=pchs[i],
               xlim=c(0, 500000), ylim=c(0,1),
               xlab="Years ago", ylab="Relative Coalescence"))
}
legend("bottomright", legend=n, col=cols, lty=ltys, pch=pchs, cex=0.7)
dev.off()
