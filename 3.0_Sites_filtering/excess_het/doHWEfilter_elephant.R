library(RcppCNPy)

prefix <- "/artemis-internal/casia16/asianbos/sites_QC/results/internal_noDP_noFD_e2"
allbed <- read.table("/artemis-internal/casia16/asianbos/sites_QC/note/bostau_len_chr.txt",stringsAsFactors=F)

estF <- npyLoad(paste0(prefix,".inbreed.sites.npy"))
site <- scan(paste0(prefix,".sites"),what="df")
lrt <- npyLoad(paste0(prefix,".lrt.sites.npy"))
dfsite<-lapply(site,function(x) unlist(strsplit(x,"[_]")))
chr<-sapply(dfsite,function(x) ifelse(length(x)==3,paste(x[1],x[2],sep="_"),x[1]))
pos<-sapply(dfsite,function(x) ifelse(length(x)==3,as.numeric(x[3]),as.numeric(x[2]) ))
#chr <- sub("_[0123456789]+$","",x=site)
#pos <- as.integer(sub(".*_","",x=site))
#pos<-sapply(dfsite,function(x) as.numeric(x[3]))

#cut -f1 /isdata/hellergrp/nzg134/radseq/ref/elephant/fzx.fa.fai > /isdata/hellergrp/nzg134/radseq/ref/elephant/fzxRawScaff.txt
chrs <- scan("/artemis-internal/casia16/asianbos/sites_QC/note/bostau_chr.txt", what="ds")
#k <- chr %in% chr1mb

#estF <- estF[k]
#lrt <- lrt[k]
#chr <- chr[k]
#pos <- pos[k]
maxPos <- tapply(pos,chr,max)

meanF <- tapply(estF,chr,mean)
lenF <- tapply(estF,chr,length)


badRegions <- function(ch, minF=-0.95, reg=50000){
    k <- chr == ch
    # get bad sites (meanF below minF and significantly deviating hwe)
    badpos <- pos[k][lrt[k] > 24 & estF[k] < minF]
    
    if(length(badpos)==0) return(data.frame(chr=ch, start=1, end=1))
    # get regions reg (default 50000) bp in both directions of bad sites                             
    badreg <- matrix(c(badpos-reg, badpos+reg), ncol=2)
    badreg[badreg<0] <- 1

    if(nrow(badreg)==1) return(data.frame(chr=ch,start= badreg[1,1], end= badreg[1,2]))
    # collapse overlapping regions
    badreg2 <- c()
    start <- badreg[1,1]
    end <- badreg[1,2]
    for(i in 2:nrow(badreg)){

        if(badreg[i,1]<end){
            end <- badreg[i,2]
        }else{
            badreg2 <- c(badreg2, start, end)
            start <- badreg[i,1]
            end <- badreg[i,2]
        }
    }

    badreg2 <- c(badreg2,start,end)
    
    badreg <- t(matrix(badreg2,nrow=2))
    #out <- cbind(rep(ch, nrow(badreg)), badreg)
    out <- data.frame(chr=rep(ch,nrow(badreg)), start = badreg[,1], end = badreg[,2])
    out$start <- out$start - 1
    return(out) 
}

badbedl <- lapply(chrs, badRegions,minF=-0.95, reg=25000)

badbed <- do.call('rbind', badbedl)

## cat /isdata/hellergrp/nzg134/radseq/ref/elephant/fzx.fa.fai|awk 'BEGIN {OFS="\t"}{print $1,0,$2}' >  /isdata/hellergrp/nzg134/radseq/ref/elephant/allElpAutoXScaff.bed
#allbed <- read.table("/isdata/hellergrp/nzg134/radseq/ref/elephant/allElpAutoXScaff.bed")

names(allbed) <- names(badbed)
lost <- sum(as.numeric(badbed$end-badbed$start))/sum(as.numeric(allbed$end))

#names(allbed) <- names(badbed)

lenghtbad <- tapply(badbed$end-badbed$start, badbed$chr, sum)

Nbadreg <- tapply(badbed$chr, badbed$chr, length)* (lenghtbad>0)

meanFdf<-data.frame(chr=names(meanF),meanF=meanF)
lengthbed_df<-data.frame(chr=names(lenghtbad),lenghtbad=lenghtbad)
Nbadreg_df<-data.frame(chr=names(Nbadreg),Nbadreg=Nbadreg)
summ<-merge(allbed,meanFdf,all.y=T,by="chr")
summ<-merge(summ,lengthbed_df,all.x=T,by="chr")
summ<-merge(summ,Nbadreg_df,all.x=T,by="chr")
summ$keep=TRUE
summ$proportionBad=summ$lenghtbad/summ$end
allbed=merge(allbed,meanFdf,by="chr",all.x=T,sort=F)

summarydf <- data.frame(chr=allbed$chr, length=allbed$end,meanF=allbed$meanF,
                        proportionBad=lenghtbad/allbed$end,
                        Nbadregions=Nbadreg,
                        keep = TRUE)
#                        keep= !(lenghtbad/allbed$end > 0.2 | meanF < -0.2))



write.table(summarydf,paste0(prefix,"InbreedSummary.tsv"),
            col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")



splitChrs <- function(c, s=5e6, bed=allbed){

    cat("doing chr",c,"\n")
    i <- which(bed[,1]==c)
    idx <- cut(x=pos[chr==c], breaks=seq(allbed[i,2], allbed[i,3], length.out=ceiling(allbed[i,3]/5e6)),
               labels=FALSE)

    

    return(paste(c, idx,sep="_"))

}
# split chromosomes in 5 mb fragments for better visualization
splitt <- unlist(sapply(chrs, splitChrs))

meanF2 <- tapply(estF, splitt, mean)

w <- names(sort(meanF2))
pdf(paste0(prefix,"byChrSplitPlotLen50k.pdf"))
for(i in 1:length(w)){
    reg <- w[i]

    p<-pos[splitt==reg]
    F<-estF[splitt==reg]
    l <- lrt[splitt==reg]
    
    plot(p,F,pch=4,col=ifelse(F< -0.95 & l>24, "goldenrod","grey"),lwd=2,main=reg,ylab="F",xlab="position")
    
    win <- 200
    fun<-function(x,w) ( cumsum(as.numeric(x))[-c(1:w)]-rev(rev(cumsum(as.numeric(x)))[-c(1:w)])) / w
    lines(fun(p,win),fun(F,win))

    
    c <- sub("_[0123456789]+$","",x=reg)

    start <- badbed[badbed$chr==c,]$start
    end <- badbed[badbed$chr==c,]$end
     abline(v=start, col="darkred",lwd=1)
     abline(v=end, col="darkred",lwd=1)

    segments(y0=-1, y1=1, x0=start, x1=end, col="darkred")
    segments(y0=1, y1=-1, x0=start, x1=end, col="darkred")
    
}
dev.off()


#w <- 1:length(meanF)
#pdf("/home/genis/impala/localSitesQC/inbreedSite/impalaV1AllRegionsGoatMappedInbreedByChrPlot.pdf")
#for(i in 1:length(w)){
#    png(sprintf("/home/genis/impala/localSitesQC/inbreedSite/plots/%sInbreedByChrPlot%s.png",basename(prefix),  names(meanF)[w[i]]), width=2000, height=500)
#     plot(pos[chr==names(meanF)[w[i]]],estF[chr==names(meanF)[w[i]]],pch=4,col="goldenrod",lwd=2,main=names(meanF)[w[i]],ylab="F",xlab="position")#
#	 p<-pos[chr==names(meanF)[w[i]]]
#	 F<-estF[chr==names(meanF)[w[i]]]
#	 win <- 200
#	 fun<-function(x,w)
#	    ( cumsum(as.numeric(x))[-c(1:w)]-rev(rev(cumsum(as.numeric(x)))[-c(1:w)])) / w
 #   lines(fun(p,win),fun(F,win))
  #   abline(v=badbed[badbed$chr==chrs[w[i]],]$start, col="darkred",lwd=1,lty=2)
   #  abline(v=badbed[badbed$chr==chrs[w[i]],]$end, col="darkred",lwd=1,lty=2)
    # dev.off()
#}
#dev.off()


badbed <- badbed[badbed$end > 1,]

finalbadbed <- rbind(badbed, data.frame(chr=summarydf$chr[!summarydf$keep],
                                start=rep(0, sum(!summarydf$keep)),
                                end=summarydf$length[!summarydf$keep]))

write.table(finalbadbed, paste0(dirname(prefix), "/excludelistInbreedSites_remapped_unmerged_unsorted.BED"), col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")

# this create blacklist, which then is merged and used to create whitelist with bedtools
