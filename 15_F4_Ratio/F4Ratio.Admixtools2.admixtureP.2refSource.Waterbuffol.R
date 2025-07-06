library(admixtools)
library(tidyverse)
args<-commandArgs(TRUE)

left=c('ZebuRef','BantengRef')
right=c('EuropeanTaurine','waterbuffalo')
target=c(args[1])
out=qpadm("/home/wlk579/2.0Indonesia_Bos_project/Indonesia_Bos_project/5.Analyses_King015/qpadmixture/using_Taurine/merge.261inds.waterbuffaloRefPlusBos",left, right, target, blgsize = 5e6, auto_only = FALSE)
write.table(out$weights, "/home/wlk579/2.0Indonesia_Bos_project/Indonesia_Bos_project/5.Analyses_King015/qpadmixture/using_Taurine/xaa/EUTaurine.weights.tsv", sep="\t", quote=FALSE, col.names = T, row.names = F)
