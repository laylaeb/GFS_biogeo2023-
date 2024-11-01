library(ggtree)
library(ggtreeExtra)
library(ape)
library(phytools)
library(ggplot2)
library(ggnewscale)
library(tidyverse)
library(castor)
library(viridis)
library(phangorn)

ASVs <- read.csv("20240221_NOMIS_rarefied_deblur_table.csv", row.names = 1)
tree <- read.tree('20240221_NOMIS_rarefied_deblur.tree')
tree <- midpoint(tree)

## Mean relative abundance (outer ring)
mean.RA <- data.frame(rowMeans(ASVs))
mean.RA$ASVs <- rownames(ASVs)
mean.RA$log.RA <- log1p(mean.RA$rowMeans.ASVs.)
colnames(mean.RA) <- c("mean.RA","ASVs", "log.RA")

## Heatmap showing HoS, DL and DR for each region
all.ASVs <- read.csv("iCAMP_all_combined.csv")
all.ASVs <- all.ASVs[all.ASVs$process =="DL" | all.ASVs$process == "HoS" | all.ASVs$process == "DR",]

## Taxonomy
tax <- read.csv("20240221_NOMIS_rarefied_deblur_taxonomy.csv", head=T, row.names = 1)

## Selected phyla (HoS and/or DL important across many regions)
p_acid_asvs <- tax[tax$Phylum==" p__Acidobacteriota",]
p_actino_asvs <- tax[tax$Phylum==" p__Actinobacteriota",]
p_bact_asvs <- tax[tax$Phylum==" p__Bacteroidota",]
p_bdello_asvs <- tax[tax$Phylum==" p__Bdellovibrionota",]
p_gemma_asvs <- tax[tax$Phylum==" p__Gemmatimonadota",]
p_latesc_asvs <- tax[tax$Phylum==" p__Latescibacterota",]
p_myxo_asvs <- tax[tax$Phylum==" p__Myxococcota",]
p_plan_asvs <- tax[tax$Phylum==" p__Planctomycetota",]
p_prot_asvs <- tax[tax$Phylum==" p__Proteobacteria",]
p_verr_asvs <- tax[tax$Phylum==" p__Verrucomicrobiota",]

taxa_df <- data.frame(ASV=tree$tip.label, Family=rep('Others', length(tree$tip.label)), Phylum=rep('Others', length(tree$tip.label)))
taxa_df$Phylum[taxa_df$ASV %in% rownames(p_acid_asvs)] <- 'Acidobacteriota'
taxa_df$Phylum[taxa_df$ASV %in% rownames(p_actino_asvs)] <- 'Actinobacteriota'
taxa_df$Phylum[taxa_df$ASV %in% rownames(p_bact_asvs)] <- 'Bacteroidota'
taxa_df$Phylum[taxa_df$ASV %in% rownames(p_bdello_asvs)] <- 'Bdellovibrionota'
taxa_df$Phylum[taxa_df$ASV %in% rownames(p_gemma_asvs)] <- 'Gemmatimonadota'
taxa_df$Phylum[taxa_df$ASV %in% rownames(p_latesc_asvs)] <- 'Latescibacterota'
taxa_df$Phylum[taxa_df$ASV %in% rownames(p_myxo_asvs)] <- 'Myxococcota'
taxa_df$Phylum[taxa_df$ASV %in% rownames(p_plan_asvs)] <- 'Planctomycetota'
taxa_df$Phylum[taxa_df$ASV %in% rownames(p_prot_asvs)] <- 'Proteobacteria'
taxa_df$Phylum[taxa_df$ASV %in% rownames(p_verr_asvs)] <- 'Verrucomicrobiota'

## Plot tree with taxonomy, heatmap and mean relative abundance 
p = ggtree(tree, layout="fan", size=0.25, open.angle=10)
p2<-p + 
  geom_fruit(data=taxa_df, geom=geom_tile,
             mapping=aes(y=ASV, fill=Phylum), width = 0.5,
             offset = 0.02) + scale_fill_manual(values = c('#FA907B','#7bbbfa', '#C42D50','#46CF6B','#529C7E','#84CCFA', '#3234D9', 'white','#fafa7b','#bbfa7b','#cf46aa')) + new_scale_fill() +
  geom_fruit(
    data=all.ASVs,
    geom=geom_tile,
    mapping=aes(y=ID, x=region, fill=process),
    offset=0.08,   # The distance between external layers, default is 0.03 times of x range of tree.
    pwidth=0.25 # width of the external layer, default is 0.2 times of x range of tree.
  ) + 
  scale_fill_manual(
    values=c("#3303df", "#dfac03","#df3e03"),
    guide=guide_legend(keywidth=0.5, keyheight=0.5, order=3)
  ) 


p3<-p2 + 
  new_scale_fill()+
  geom_fruit(data=mean.RA, 
             geom=geom_bar, 
             mapping=aes(y=ASVs, x=log.RA),
             offset = 0.05, orientation='y',
             stat="identity") 


ggsave(p3, filename = 'Figure_phy_tree.pdf', width = 15, height = 15)



