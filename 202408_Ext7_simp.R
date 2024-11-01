
library(vegan)
library(ggplot2)
library(ape)
library(picante)

ASVs <-read.csv("20240221_NOMIS_rarefied_deblur_table.csv", row.names = 1)
MG <-read.csv("NOMIS_merged_KEGG_counts_20240104.csv", row.names = 1)
prok.KO <-read.table("prokaryote_KOs.txt", head=T)
MG <-MG[rownames(MG) %in% prok.KO$KO,]
MG.md <-read.csv("MG_metadata_clean.csv")
MG.sed <-MG.md[MG.md$substrate=="sed",]
MG <-MG[,colnames(MG) %in% MG.sed$sample.id]
MG <-MG[rowSums(MG)>0,]

#averaging UP-DN MGs
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

## Match MG and ASV tables
colnames(ASVs)
ASV.sub <-ASVs[,colnames(ASVs) %in% colnames(MG.comb2)]
ASV.sub <-ASV.sub[rowSums(ASV.sub)>0,]
ASV.sub <-ASV.sub[,order(colnames(ASV.sub))]
MG.sub <-MG.comb2
MG.sub <-MG.sub[,order(colnames(MG.sub))]
MG.sub <-MG.sub[,!colnames(MG.sub)=="GL61"]
summary(colnames(ASV.sub)==colnames(MG.sub))

## read Levin's niche breadth index (specialist_generalists), see Ext Fig 8d
spc_gen <-read.csv("spec_gen_KOs.csv", row.names = 1)

## SIMPER
KO.sim <-simper(t(MG.sub))
KO.sim.res <-data.frame(summary(KO.sim)$total)
KO.sim.res <-KO.sim.res[order(rownames(KO.sim.res)),]
summary(rownames(KO.sim.res)==rownames(spc_gen))

KO.sim.res$spc.gen <-spc_gen$sign
KO.sim.res2 <-KO.sim.res[KO.sim.res$average>0,]

boxplot(KO.sim.res2$average~KO.sim.res2$spc.gen, xlab="", log="y", ylab="contribution to Bray-Curtis dissimilarity", col="white")




