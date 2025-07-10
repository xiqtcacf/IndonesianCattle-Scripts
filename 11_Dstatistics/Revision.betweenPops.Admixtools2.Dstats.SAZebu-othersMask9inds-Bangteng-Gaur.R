library(admixtools)
#library(tidyverse)

pop1='SouthAsianzebu'
pop2=c('Aceh','EastAsianzebu','Jabres','Madura','Pasundan','Pesisir','SumbaOngole','Bangladesh1','Bangladesh2','Bangladesh3','Dehong1','Dehong2','Jiangcheng1','Jiangcheng2','Kailali','Kapilvastu','Longlin1','Longlin2','Myanmar1','Myanmar2','Myanmar3','Saptari','Shigatse1','Shigatse2','Sri_Lanka1','Sri_Lanka2')
pop3='Captivebanteng'
pop4='Gaur'

out = qpdstat("/home/wlk579/2.0Indonesia_Bos_project/Revision_NC/EmuZebu_SouthernIndian/Dstatistics_LoterMask/141inds.BosTau9Ref.AsianBos.SouthIndia.Mask9inds", pop1=pop1,pop2=pop2,pop3=pop3,pop4=pop4,unique_only = TRUE, auto_only = FALSE, f4mode = FALSE, sure = TRUE,blgsize = 5e6)
#write_tsv(out, "/isdata/hellergrp/wlk579/Bos_project/5.Analyses/Admixtools2/betweenPops.SouthAsia-others-Banteng-waterbuffol.dstats.tsv")
write.table(out, "/home/wlk579/2.0Indonesia_Bos_project/Revision_NC/EmuZebu_SouthernIndian/Dstatistics_LoterMask/141inds.BosTau9Ref.AsianBos.SouthIndia.Mask9inds.dstats.tsv", sep="\t", quote=FALSE, col.names = T, row.names = F)
