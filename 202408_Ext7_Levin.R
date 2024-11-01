library(vegan)
library(ggplot2)
library(ape)
library(picante)
library(EcolUtils)

MG <-read.csv("NOMIS_merged_KEGG_counts_20240104.csv", row.names = 1)
prok.KO <-read.table("prokaryote_KOs.txt", head=T)
MG <-MG[rownames(MG) %in% prok.KO$KO,]
MG.md <-read.csv("MG_metadata_clean.csv")
MG.sed <-MG.md[MG.md$substrate=="sed",]
ASVs <-read.csv("20240221_NOMIS_rarefied_deblur_table.csv", row.names = 1)

MG <-MG[,colnames(MG) %in% MG.sed$sample.id]
MG <-MG[rowSums(MG)>0,]

## averaging UP-DN MGs
MG.tab <-data.frame(table(MG.sed$GFS))
MG.comb.tab <-MG.tab[MG.tab$Freq==1,]
MG.tab1 <-MG.sed[MG.sed$GFS %in% MG.comb.tab$Var1,]
MG.comb <-MG[,colnames(MG) %in% MG.tab1$sample.id]
for(i in MG.tab$Var1[MG.tab$Freq==2]){
  MG.sed.i<-MG.sed[MG.sed$GFS==i,]
  MG.i<-MG[,colnames(MG) %in% MG.sed.i$sample.id]
  MG.comb<-cbind(MG.comb, rowMeans(MG.i))
}

MG.comb <-data.frame(MG.comb)
colnames(MG.comb) <-c(MG.tab1$GFS,as.character(MG.tab$Var1[MG.tab$Freq==2]))
MG.comb2 <-decostand(MG.comb, MARGIN=2, method="total")

## match MG and ASV tables
ASV.sub<-ASVs[,colnames(ASVs) %in% colnames(MG.comb2)]
ASV.sub<-ASV.sub[rowSums(ASV.sub)>0,]
ASV.sub<-ASV.sub[,order(colnames(ASV.sub))]
MG.sub<-MG.comb2
MG.sub<-MG.sub[,order(colnames(MG.sub))]
MG.sub<-MG.sub[,!colnames(MG.sub)=="GL61"]
summary(colnames(ASV.sub)==colnames(MG.sub))
MG.comb2<-decostand(MG.comb, MARGIN=2, method="total")
colSums(MG.comb2)
MG.comb2<-round(MG.comb2/(min(MG.comb2[MG.comb2>0])),0)

MG.comb3<-data.frame(t(MG.comb2[rowSums(MG.comb2>0)<84,]))  #remove KOs present in all samples (otherwise algorithm crashes)

spc_gen<-spec.gen(MG.comb3, niche.width.method = "levins", perm.method = "quasiswap", n = 999, probs = c(0.025, 0.975))
#write.csv(spc_gen, "spec_gen_KOs.csv")  #add back KOs present in all samples (KO generalists)
spc_gen<-read.csv("spec_gen_KOs.csv", row.names=1)  #read result
MG.comb2<-MG.comb2[order(rownames(MG.comb2)),]
summary(rownames(MG.comb2)==rownames(spc_gen))
MG.comb3<-decostand(MG.comb2, method="total", MARGIN=2)

plot(rowSums(MG.comb3>0),rowSums(MG.comb3),log="y", col=as.factor(spc_gen$sign),pch=19, xlab="prevalence [nr GFS]", ylab="relative abundance")
legend("bottomright",levels(as.factor(spc_gen$sign)),col=1:3,pch=19,inset=0.01,cex=0.7, bty="n")


