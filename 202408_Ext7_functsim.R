library(vegan)
library(ggplot2)
library(ape)
library(picante)

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

ASVs <-read.csv("20240221_NOMIS_rarefied_deblur_table.csv", row.names = 1)

## match MG and ASV tables
colnames(ASVs)
ASV.sub <-ASVs[,colnames(ASVs) %in% colnames(MG.comb2)]
ASV.sub <-ASV.sub[rowSums(ASV.sub)>0,]
ASV.sub <-ASV.sub[,order(colnames(ASV.sub))]
MG.sub <-MG.comb2
MG.sub <-MG.sub[,order(colnames(MG.sub))]

MG.sub <-MG.sub[,!colnames(MG.sub)=="GL61"]
summary(colnames(ASV.sub)==colnames(MG.sub))

## compositional vs functional dissimilarity 
BC.ASV <-vegdist(t(ASV.sub), dist="Bray")
BC.MG <-vegdist(t(MG.sub), dist="Bray")
BC.ASV.long <-as.data.frame(as.table(t(as.matrix(BC.ASV))))
BC.MG.long <-as.data.frame(as.table(t(as.matrix(BC.MG))))
summary(BC.ASV.long$Var1==BC.MG.long$Var1)
summary(BC.ASV.long$Var2==BC.MG.long$Var2)
BC.comb <-BC.ASV.long
BC.comb$BC.MG <-BC.MG.long$Freq
colnames(BC.comb) <-c("S1","S2","BC.ASV","BC.MG")
BC.comb <-BC.comb[BC.comb$BC.ASV>0,]

ggplot(BC.comb, aes(x=BC.ASV, y=BC.MG) ) +
  geom_hex(bins = 45) +
  scale_fill_continuous(type = "viridis") +
  xlab("taxonomic dissimilarity")+
  ylab("functional dissimilarity")+
  theme_bw()


