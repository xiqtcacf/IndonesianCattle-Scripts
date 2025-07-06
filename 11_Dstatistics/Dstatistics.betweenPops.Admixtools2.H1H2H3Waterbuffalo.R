library(admixtools)
#library(tidyverse)

pop1=c('Aceh','Africantaurine','Africanzebu','Asianadmixed','Australia','Bali','Captivebanteng','EastAsiantaurine','EastAsianzebu','Europeantaurine','Gaur','Indonesia','Jabres','Kupang','Madura','Pasundan','Pesisir','SouthAsianzebu','SumbaOngole','waterbuffaloRef')
pop2=c('Aceh','Africantaurine','Africanzebu','Asianadmixed','Australia','Bali','Captivebanteng','EastAsiantaurine','EastAsianzebu','Europeantaurine','Gaur','Indonesia','Jabres','Kupang','Madura','Pasundan','Pesisir','SouthAsianzebu','SumbaOngole','waterbuffaloRef')
pop3=c('Aceh','Africantaurine','Africanzebu','Asianadmixed','Australia','Bali','Captivebanteng','EastAsiantaurine','EastAsianzebu','Europeantaurine','Gaur','Indonesia','Jabres','Kupang','Madura','Pasundan','Pesisir','SouthAsianzebu','SumbaOngole','waterbuffaloRef')
pop4='waterbuffaloRef'

out = qpdstat("/isdata/hellergrp/wlk579/Indonesia_Bos_project/5.Analyses_King015/Admixtools2/between_pops/pure.merge.222inds.waterbuffaloRefPlusBos", pop1=pop1,pop2=pop2,pop3=pop3,pop4=pop4,unique_only = TRUE, auto_only = FALSE,f4mode = FALSE, sure = TRUE, blgsize = 5e6)
#write_tsv(out, "/isdata/hellergrp/wlk579/Bos_project/5.Analyses/Admixtools2/betweenPops.SouthAsia-others-Banteng-waterbuffol.dstats.tsv")
write.table(out, "/isdata/hellergrp/wlk579/Indonesia_Bos_project/5.Analyses_King015/Admixtools2/between_pops/allpopsH1H2H3.waterbuffolH4.dstats.tsv", sep="\t", quote=FALSE, col.names = T, row.names = F)
