library(picante)
library(adespatial)
library(speedyseq)
library(phyloseq)
library(ape)
library(progress)

ASVs <- read.csv("20240221_NOMIS_rarefied_deblur_table.csv", head=T, row.names = 1)
samp <- read.table("20240221_NOMIS_metadata_GFS.tsv", sep="\t", head=T)
tree <- read.tree("20240221_NOMIS_rarefied_deblur.tree")

## phylogenetic structure of beta-diversity 
# mean ASV abundance per region 
regs <- c()
for(j in levels(factor(samp$site_c))){
  reg.j<-samp[samp$site_c==j,]
  n.j<-rowMeans(ASVs[reg.j$sample])
  regs<-cbind(regs,n.j)
}
regs <- data.frame(regs)
colnames(regs) <- levels(factor(samp$site_c))

regs <- regs[,1:10] #drop Uganda

MNTDs <- c()
for(k in levels(factor(samp$site_c))){
  reg.k<-samp[samp$site_c==k,]
  n.k<-ASVs[reg.k$sample]
  n.k<-n.k[rowSums(n.k)>0,]
  tree.k<-prune.sample(t(n.k), tree)
  mntd.k<-mntd(t(n.k), cophenetic(tree.k))
  mntd.k<-data.frame(mntd.k)
  mntd.k$reg<-k
  MNTDs<-rbind(MNTDs, mntd.k)
}
boxplot(MNTDs$mntd.k~MNTDs$reg, las=2, xlab="", ylab="MNTD")
abline(h=mean(na.omit(as.numeric(MNTDs$mntd.k))), col="red", lwd=2)


# beta diversity profiles across regions 
across.part <- c()
b.part.across.all <- beta.div.comp(t(regs), coef="S", quant=T)
across.part <- rbind(across.part, c(b.part.across.all$part, nrow(regs),0))

pb <- progress_bar$new(total = length(seq(0.005,0.2,0.005)))
regs.phy <- phyloseq(otu_table(regs, taxa_are_rows = TRUE), phy_tree(tree))
for(h in seq(0.005,0.2,0.005)){
  glom.h<-speedyseq::tree_glom(regs.phy, h)
  OTU.h<-as.data.frame(otu_table(glom.h), "matrix")
  b.part.across<-beta.div.comp(t(OTU.h), coef="S", quant=T)
  across.part<-rbind(across.part, c(b.part.across$part, nrow(OTU.h), h))
  pb$tick()
}
across.part <- data.frame(across.part)
colnames(across.part)[6:7] <- c("tips","phy.agglom")

plot(across.part$phy.agglom, across.part$BDtotal,type="b",lwd=1.75, cex=1, xlim=c(0,0.2), ylim=c(0,0.5), xlab="phylogenetic agglomeration depth", ylab="beta diversity")
abline(v=mean(na.omit(MNTDs$mntd.k)), lty=1, lwd=2)

## Beta diversity profiles within regions 
samp2 <- samp[samp$site_c!="Uganda",]
regs.parts <- c()
for(a in levels(factor(samp2$site_c))){
  samp.a<-samp2[samp2$site_c==a,]
  reg.a<-ASVs[samp.a$sample]
  reg.a<-reg.a[rowSums(reg.a)>0,]
  b.part.a<-beta.div.comp(t(reg.a), coef="S", quant=T)
  regs.parts<-rbind(regs.parts, c(b.part.a$part, nrow(reg.a), nrow(reg.a)/nrow(reg.a), h=0, a))}
for(a in levels(factor(samp2$site_c))){
  samp.a<-samp2[samp2$site_c==a,]
  reg.a<-ASVs[samp.a$sample]
  reg.a<-reg.a[rowSums(reg.a)>0,]
  reg.a.tree<-prune.sample(t(reg.a),tree)
  reg.a.phy<-phyloseq(otu_table(reg.a, taxa_are_rows = TRUE), phy_tree(reg.a.tree))
  for(h in seq(0.005,0.2,0.005)){
    glom.h<-speedyseq::tree_glom(reg.a.phy, h)
    OTU.h<-as.data.frame(otu_table(glom.h), "matrix")
    b.part<-beta.div.comp(t(OTU.h), coef="S", quant=T)
    regs.parts<-rbind(regs.parts, c(b.part$part, nrow(OTU.h), nrow(OTU.h)/nrow(reg.a), h, a))
    print(c(a,h))}
}

regs.parts <- data.frame(regs.parts)
colnames(regs.parts)[6:9] <- c("n_tips","perc_tips","phy_glom","region")

#write.csv(regs.parts,"within.regions.betapart.csv")

plot(regs.parts$phy_glom,regs.parts$BDtotal, type="n", ylim=c(0,6), xlab="phylogenetic agglomeration depth", ylab="beta diversity change (%)")
for(i in levels(factor(regs.parts$region))){
  reg.i<-regs.parts[regs.parts$region==i,]
  points(reg.i$phy_glom[1:40], abs(diff(as.numeric(reg.i$BDtotal)))/as.numeric(reg.i$BDtotal)[1]*100, pch=16, cex=1, col="grey")
}
points(across.part$phy.agglom[1:40], abs(diff(across.part$BDtotal)/across.part$BDtotal[1])*100, pch=16, col="red", type="b")
abline(v=0.08498672, lty=5, col="darkgrey")


x <-as.numeric(across.part$phy.agglom)[1:40]
y <-abs(diff(across.part$BDtotal)/across.part$BDtotal[1])*100

m2 <-nls(y~a*exp(b*x),start=list(a=8,b=-50))
lines(x,predict(m2),col="red",lwd=3)


## Phylogenetic structure of uniqueness 

#across regions

across.unique <-c(sum(rowSums(regs>0)==1), nrow(regs),1,0)
regs.phy<-phyloseq(otu_table(regs, taxa_are_rows = TRUE), phy_tree(tree))
for(h in seq(0.005,0.2,0.005)){
  glom.h<-speedyseq::tree_glom(regs.phy, h)
  OTU.h<-as.data.frame(otu_table(glom.h), "matrix")
  across.unique<-rbind(across.unique, c(sum(rowSums(OTU.h>0)==1), nrow(OTU.h),nrow(OTU.h)/nrow(regs),h))
  print(c(h))
}

across.unique <-data.frame(across.unique)
colnames(across.unique) <-c("uniques","n_tips","n_total","h")

#within regions

regs.uniques <-c()
for(a in levels(factor(samp2$region))){
  samp.a<-samp2[samp2$region==a,]
  reg.a<-ASVs[samp.a$ID]
  reg.a<-reg.a[rowSums(reg.a)>0,]
  reg.a.uni<-sum(rowSums(reg.a>0)==1)
  regs.uniques<-rbind(regs.uniques, c(reg.a.uni, nrow(reg.a), 1, 0, a))}


for(a in levels(factor(samp2$region))){
  samp.a<-samp2[samp2$region==a,]
  reg.a<-ASVs[samp.a$ID]
  reg.a<-reg.a[rowSums(reg.a)>0,]
  reg.a.tree<-prune.sample(t(reg.a),tree)
  reg.a.phy<-phyloseq(otu_table(reg.a, taxa_are_rows = TRUE), phy_tree(reg.a.tree))
  for(h in seq(0.005,0.2,0.005)){
    glom.h<-speedyseq::tree_glom(reg.a.phy, h)
    OTU.h<-as.data.frame(otu_table(glom.h), "matrix")
    reg.uni<-sum(rowSums(OTU.h>0)==1)
    regs.uniques<-rbind(regs.uniques, c(reg.uni, nrow(OTU.h), nrow(OTU.h)/nrow(reg.a), h, a))
    print(c(a,h))
    }
}

regs.uniques <- data.frame(regs.uniques)

colnames(regs.uniques) <-c("uniq_tips","all_tips", "perc_tips","h","region")
plot(regs.uniques$h, regs.uniques$perc_tips, type="n",  xlab="agglomeration level",ylab="unique taxa (%)", ylim=c(0,100))
for(z in levels(factor(regs.uniques$region))){
  regs.uniques.z<-regs.uniques[regs.uniques$region==z,]
  points(regs.uniques.z$h, as.numeric(regs.uniques.z$perc_tips)*100,pch=16, col="grey", cex=1)
}
points(across.unique$h, across.unique$uniques/max(across.unique$uniques)*100, type="b", col="red", pch=16)


