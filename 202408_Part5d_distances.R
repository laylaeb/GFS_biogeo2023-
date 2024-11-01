library(picante)
library(progress)
library(ggplot2)
library(ggridges)
library(forcats)


ASVs <-read.csv("20240221_NOMIS_rarefied_deblur_table.csv", row.names=1)
tree <-read.tree("20240221_NOMIS_rarefied_deblur.tree")
tax <-read.csv("20240221_NOMIS_rarefied_deblur_taxonomy.csv", row.names = 1)

levels(factor(tax$Genus))

goi <-data.frame(table(tax$Genus))
goi <-goi[order(goi$Freq, decreasing=TRUE), ]
goi2 <-goi[3:32,]
goi2 <-goi2[!goi2$Var1 %in% c(" g__OM190"," g__vadinHA49"," g__A21b"," g__0319-6G20", " g__TRA3-20"," g__AKYH767", " g__Anaeromyxobacter"," g__Blfdi19"," g__Candidatus_Udaeobacter"),]
goi3 <-goi[3:102,]

pb <- progress_bar$new(total = nrow(goi3))
res1 <-c()
for(g in levels(factor(goi3$Var1))){
  tax.g<-tax[tax$Genus==g,]
  com.g<-ASVs[rownames(ASVs) %in% rownames(tax.g),]
  com.g<-com.g[,colSums(com.g>0)>2]
  tree.g<-prune.sample(t(com.g),tree)
  coph.g<-mean(cophenetic(tree.g))
  res1<-rbind(res1,c(coph.g,g))
  pb$tick()
  }

#write.csv(res1,"overall_phy_dists_genera.csv")

res <-c()
for(g in levels(factor(goi2$Var1))){
  tax.g<-tax[tax$Genus==g,]
  com.g<-ASVs[rownames(ASVs) %in% rownames(tax.g),]
  com.g<-com.g[,colSums(com.g>0)>2]
  tree.g<-prune.sample(t(com.g),tree)
  for(k in 1:ncol(com.g)){
    com.k<-com.g[,k, drop=FALSE]
    com.k<-com.k[com.k>0,,drop=FALSE]
    tree.k<-prune.sample(t(com.k),tree.g)
    coph.k<-mean(cophenetic(tree.k))
    res<-rbind(res, c(g, coph.k, colnames(com.k)))}
   }

res <-data.frame(res)
colnames(res) <-c("genus","phy.dist","sample")

write.csv(res, "phylogenetic_distances_genera1.csv")  # add "order" of genera to be displayed
res2 <-read.csv("phylogenetic_distances_genera.csv")

ggplot(res2,aes(x=phy.dist, 
                y=fct_reorder(genus,order), 
                fill=genus))+
  geom_density_ridges()+
  theme_ridges() + 
  theme(legend.position = "none")




