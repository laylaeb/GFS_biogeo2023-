library(vegan)
library(stringr)
setwd(choose.dir())


ASVs<-read.csv("20240221_NOMIS_rarefied_deblur_table.csv", row.names=1)
ASVs.rel<-decostand(ASVs, "total",MARGIN=2)

# filtered dataset
par(mfrow=c(12,13), mar=c(0,0.5,1,0),oma=c(0.1,0.1,0.1,0.1), cex.main=0.75)

ASVs.rel<-ASVs.rel[,str_sort(colnames(ASVs.rel), numeric = TRUE)]
for(i in 1:ncol(ASVs.rel)){
  ASVs.i<-ASVs.rel[,i]
  rad.ln.i<-rad.lognormal(ASVs.i, family=gaussian)
  rad.null.i<-rad.null(ASVs.i)
  plot(rad.ln.i, col="grey", main=colnames(ASVs.rel)[i], xaxt="n", yaxt="n") 
  lines(rad.ln.i, col="red", lwd=2)
  lines(rad.null.i, col="cyan", lwd=2)
}



#unfiltered dataset
ASVs.unfilt<-read.csv("081623_NOMIS_16S_merged_deblur_table.csv", head=T, row.names = 1)
samp<-read.csv("081623_NOMIS_metadata_GFS.csv", head=T)

ASVs.unfilt<-ASVs.unfilt[rowSums(ASVs.unfilt)>0,colnames(ASVs.unfilt) %in% samp$id]
ASVs.unfilt<-ASVs.unfilt[,str_sort(colnames(ASVs.unfilt), numeric = TRUE)]
colnames(ASVs.unfilt)<-sub("\\_.*", "", colnames(ASVs.unfilt))

nms<-unique(colnames(ASVs.unfilt))
par(mfrow=c(12,13), mar=c(0,0.5,1,0),oma=c(0.1,0.1,0.1,0.1), cex.main=0.75)
for(i in 1:151){
  nms.i<-nms[i]
  ASVs.i<-data.frame(ASVs.unfilt[,colnames(ASVs.unfilt) %in% nms.i])
  ASVs.i<-rowMeans(ASVs.i)
  rad.ln<-rad.lognormal(ASVs.i, family=gaussian)
  plot(rad.ln, col="grey", main=nms.i, xaxt="n", yaxt="n") 
  rad.null<-rad.null(ASVs.i)
  lines(rad.ln, col="red", lwd=2)
  lines(rad.null, col="cyan", lwd=2)
}

