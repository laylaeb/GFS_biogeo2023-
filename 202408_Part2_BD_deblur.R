### Beta diversity across mountain ranges ###
library(phyloseq)
library(phyloseqCompanion)
library(vegan)
library(ggplot2)
library(betadisper)
library(adonis)
library(ecodist)

asv<-read.csv(file="20240221_NOMIS_rarefied_deblur_table.csv.gz",sep=",",header=TRUE,row.names=1)
tax<-read.csv(file="20240221_NOMIS_rarefied_deblur_taxonomy.csv.gz",sep=",",header=TRUE,row.names=1)

## Load metadata file 
metadata_NOMIS="202402_NOMIS_metadata_GFS.tsv"
metadata_NOMIS<-import_qiime_sample_data(metadata_NOMIS)

## Create ASV table 
OTU_NOMIS <- otu_table(asv, taxa_are_rows=TRUE)
tax_NOMIS <- tax_table(as.matrix(tax))

## Create phyloseq object
merged_NOMIS_DEBLUR <- merge_phyloseq(OTU_NOMIS, tax_NOMIS, metadata_NOMIS)

## Be sure to compute this with the correct object - either prune_Uganda for statistical purposes, or merged_NOMIS_DEBLUR to show the NMDS including all samples.
prune_Uganda <- subset_samples(merged_NOMIS_DEBLUR , !site_c %in% "Uganda")
prune_Uganda_df <- as.matrix(otu_table(prune_Uganda, taxa_are_rows=T))

metadata_nmds <- sample.data.frame(prune_Uganda)
vegan_otu <- function(physeq){
  OTU <- otu_table(physeq)
  if(taxa_are_rows(OTU)){
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}

metadata_nmds$lat_attribute <- "A"
metadata_nmds$lat_attribute <- ifelse((metadata_nmds$lat_sp) < 1, metadata_nmds$lat_attribute == "B", metadata_nmds$lat_attribute == "A")
metadata_nmds$lat_attribute <- as.factor(metadata_nmds$lat_attribute)

## NMDS plot
asv_table_nmds <- as.matrix((otu_table(prune_Uganda, taxa_are_rows=T)))
asv_table_nmds_f <- asv_table_nmds[rowSums(asv_table_nmds[])>0,]
nmds_bc_nomis <- metaMDS(t(log1p(asv_table_nmds_f)), distance = "bray", k = 2, trymax=999)

stressplot(nmds_bc_nomis)
nmds_bc_nomis$stress
# nmds_bc_nomis$stress
data.scores = as.data.frame(scores(nmds_bc_nomis)$sites)
#add columns to data frame 
data.scores$Sample = metadata_nmds$Sample
data.scores$Site = metadata_nmds$site_c
head(data.scores)

nmds_bc_GFS_plot = ggplot(data.scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 6, aes(colour = Site))+ 
 #stat_ellipse(aes(x=NMDS1, y=NMDS2,color=Site),type = "norm")+
  theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 12), 
        legend.text = element_text(size = 12, face ="bold", colour ="black"), 
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 14), 
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
        legend.title = element_text(size = 14, colour = "black", face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank()) + 
  labs(x = "NMDS1", colour = "Region", y = "NMDS2", shape = "Type") 

nmds_bc_GFS_plot

## Include Smooth line to show the latitude
colors<-c("#2E2A2BFF", "#CF4E9CFF","#8C57A2FF",
          "#3EBCB6","#82581FFF","#2F509EFF",
          "#E5614CFF","#97A1A7FF","#bee183","#DC9445FF","#EDD03E")

plot(x=data.scores$NMDS1, y=data.scores$NMDS2, type="n", xlim = c(-2.5, 2.9), ylim = c(-1.6, 1.8))
points(nmds_bc_nomis, display = "sites", cex = 2.3, pch=19, col=alpha(colors[factor(data.scores$Site)], 0.8))
ordisurf(nmds_bc_nomis, metadata_nmds$lat_sp, add = TRUE, col="blue", labcex=1)

merge_data <- merge_data[merge_data$sample != "GL140", ]

## Statistical analyses - how mountain ranges shape community composition
vegan_matrix_GFS<- vegan_otu(prune_Uganda)
GFS_bray <- vegdist(log1p(vegan_matrix_GFS), method="bray")

## PERMDISP region -- if significant, compute manyglm
bdisp_nomis<- betadisper(GFS_bray, metadata_nmds$site_c, type=c("centroid"))
bdisp_nomis
aov_bdisp <-anova(bdisp_nomis)
permutest(bdisp_nomis, pairwise=T)

## PERMDISP latitude --- if significant, compute manyglm
bdisp_nomis_lat<- betadisper(GFS_bray, metadata_nmds$lat, type=c("centroid"))
bdisp_nomis_lat
aov_bdisp_lat <-anova(bdisp_nomis_lat)
permutest(bdisp_nomis_lat, pairwise=T)

# if not significant go with adonis and pairwise adonis and report these different values in the manuscript
GFS_ado <- adonis2(GFS_bray ~ site_c, permutations = 999, method = "bray", data=metadata_nmds)
GFS_ado_latitude <- adonis2(GFS_bray ~ lat_attribute, permutations = 999, method = "bray", data=metadata_nmds)

library(pairwiseAdonis)
pairwise_GFS <- pairwise.adonis(GFS_bray, metadata_nmds$site_c, p.adjust.m="holm")

## Test the effects of latitude on community composition
pairwise_GFS_lat <- pairwise.adonis(GFS_bray, metadata_nmds$lat_attribute, p.adjust.m="holm")

# asv_pu <- t(otu_table(prune_Uganda, taxa_are_rows=T))
# 
# ab <- mvabund(asv_pu)
# asv_nb <- manyglm(ab ~ lat_attribute,
#                   data = metadata_nmds, family = 'negative binomial')
# 
# nomis_avo <- anova(asv_nb, p.uni= "adjusted", nBoot = 99, pairwise.comp=metadata_nmds$lat_attribute, show.time=T)
# 
# nomis_manyglm_res <- nomis_avo_site$uni.p
#write.csv(nomis_manyglm_res, "results_manyglm_nomis.csv")

