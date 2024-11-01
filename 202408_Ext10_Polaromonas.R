### POLAROMONAS ###
library(vegan)
library(picante)
library(ggtree)
library(ggtreeExtra)

ASVs <-read.csv("20240221_NOMIS_rarefied_deblur_table.csv", row.names = 1)
tax <-read.csv("20240221_NOMIS_rarefied_deblur_taxonomy.csv", head=T, row.names = 1)
tree <-read.tree('20240221_NOMIS_rarefied_deblur.tree')
samp <-read.table("20240221_NOMIS_metadata_GFS.tsv", sep="\t", head=T)
ASVs <-ASVs[,colnames(ASVs) %in% samp$sample]

Polaromonas <-tax[tax$Genus==" g__Polaromonas",]

nom.rel <-decostand(ASVs,"total",2)
Pol.ASVs <-ASVs[rownames(Polaromonas),]

Pol.regs<-c()
for(j in levels(factor(samp$site_c))){
  reg.j<-samp[samp$site_c==j,]
  n.j<-rowMeans(Pol.ASVs[reg.j$sample])
  Pol.regs<-cbind(Pol.regs,n.j)
}
Pol.regs<-data.frame(Pol.regs)
colnames(Pol.regs)<-levels(factor(samp$site_c))

Pol.tree <-prune.sample(t(Pol.ASVs),tree)
Pol.tree$edge.length[Pol.tree$edge.length==0] <-quantile(Pol.tree$edge.length,0.1)*0.1


Pol.max.reg <-data.frame(colnames(Pol.regs)[apply(Pol.regs,1,which.max)])
Pol.max.reg$abund <-apply(Pol.regs,1,max)
colnames(Pol.max.reg) <-c("region","abundance")
Pol.max.reg$ASV <-rownames(Pol.regs)
Pol.max.reg$nr.regions <-rowSums(Pol.regs>0)
Pol.max.reg[Pol.max.reg$nr.regions>3,1] <-"widespread"
Pol.max.reg <- Pol.max.reg[order(match(Pol.max.reg$ASV, Pol.tree$tip.label)),]
tiplabs_df <- data.frame(ASV = Pol.tree$tip.label, Colour = Pol.max.reg$region)


Pol.regs2 <-Pol.regs
Pol.regs2[Pol.regs2==0] <-NA

p1 <-ggtree(Pol.tree, layout="circular")
p1 %<+% tiplabs_df + geom_tippoint(aes(color=Colour), size=2)
p2 <-gheatmap(p1, log1p(Pol.regs2), offset=0, width=0.5, 
             colnames=TRUE, legend_title="regional\nabundance") +
  scale_x_ggtree() +
  scale_y_continuous(expand=c(0, 4))
p2 %<+% tiplabs_df + geom_tippoint(aes(color=Colour), size=2)


