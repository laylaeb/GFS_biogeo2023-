library(vegan)
library(ggplot2)
library(ape)
library(picante)


MG <-read.csv("NOMIS_merged_KEGG_counts_20240104.csv", row.names = 1)
prok.KO <-read.table("prokaryote_KOs.txt", head=T)
MG <-MG[rownames(MG) %in% prok.KO$KO,]
MG.md <-read.csv("MG_metadata_clean.csv")
MG.sed <-MG.md[MG.md$substrate=="sed",]
ASVs <-read.csv("20240221_NOMIS_rarefied_deblur_table.csv", row.names = 1)
tree <-read.tree("20240221_NOMIS_rarefied_deblur.tree")


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

#match MG and ASV tables
colnames(ASVs)
ASV.sub <-ASVs[,colnames(ASVs) %in% colnames(MG.comb2)]
ASV.sub <-ASV.sub[rowSums(ASV.sub)>0,]
ASV.sub <-ASV.sub[,order(colnames(ASV.sub))]
MG.sub <-MG.comb2
MG.sub <-MG.sub[,order(colnames(MG.sub))]
MG.sub <-MG.sub[,!colnames(MG.sub)=="GL61"]
summary(colnames(ASV.sub)==colnames(MG.sub))

#richness vs functional diversity
ASV.rich <-colSums(ASV.sub>0)
MG.rich <-colSums(MG.sub>0)
ASV.pd <-pd(t(ASV.sub), tree, include.root = FALSE)

MG.md2 <-MG.md[MG.md$GFS %in% rownames(ASV.pd),]
MG.md2 <- MG.md2[!duplicated(MG.md2$GFS),]

ASV.KO <-data.frame(cbind(ASV.pd, MG.rich, MG.md2$mean.chla))
colnames(ASV.KO) <-c("ASV_PD","ASV_richness","KO_richness","mean_chla")
ggplot(ASV.KO, aes(x = ASV_richness, y = KO_richness, color = log(mean_chla))) +
  geom_point(size=4)+theme_bw()

ggplot(ASV.KO, aes(x = ASV_PD, y =KO_richness, color = log(mean_chla))) +
  geom_point(size=4)+theme_bw()

