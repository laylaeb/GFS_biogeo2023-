### Indicator ASVs ###

library(foreach)
library(doParallel)
library(adiv)
library(vegan)

## Load data and remove GFS from Uganda
ASVs<- read.csv("20240221_NOMIS_rarefied_deblur_table.csv", head=T, row.names = 1)
samp<- read.table("20240221_NOMIS_metadata_GFS.tsv", sep="\t", head=T)

samp2<- samp[samp$region !="Uganda",]
asv_table<- ASVs[colnames(ASVs) %in% samp2$ID]
asv_table<- asv_table[rowSums(asv_table)>0,]
asv_table<- asv_table[rowSums(asv_table>0)>1,]# remove unique ASVs

## Hellinger-transformation
asv_table<- decostand(asv_table, method="hellinger", MARGIN=2)

n.cores<- detectCores()
registerDoParallel(cl <- makeCluster(n.cores-1))

res.f.all <-foreach(f = levels(factor(samp2$region)), .packages = c("adiv")) %dopar% {
  samp2.f <-samp2
  samp2.f$region[!grepl(f, samp2$region)]<-"other"
  Q.f<- dbMANOVAspecies(t(asv_table),samp2.f$region, nrep=999, global=FALSE,species=TRUE, padj="BH")
  Q.f.pairwise_adj <- dbMANOVAspecies_pairwise(Q.f)
  c(summary(Q.f.pairwise_adj))
}
save.image("db.MANOVA.res.RData")

