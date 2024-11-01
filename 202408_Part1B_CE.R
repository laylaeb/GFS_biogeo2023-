### Comparisons with other samples from Cryospheric Ecosystems - samples retrieved from Bourquin et al.2022 ###

library(phyloseq)
library(phyloseqCompanion)
library(data.table)
library(tidyverse)
library(Rmisc)
library(vegan)
library(wesanderson)
library(paletteer)
library(FSA)
library(speedyseq)
library(RColorBrewer)
library(R.filesets)

## Read cryobiome table 
asv_table_NOMISCRYO = (fread("cryobiome_table.txt", header=T))

## Load metadata CRYO
metadata_cryo="PP1_metadata_date.tsv"
metadata_cryo<-import_qiime_sample_data(metadata_cryo)

asv_unrarefied<-read.csv(file="20240221_NOMIS_filtered_deblur_table.csv.gz",sep=",",header=TRUE,row.names=1)
tax_unrarefied<-read.csv(file="20240221_NOMIS_filtered_deblur_taxonomy.csv.gz",sep=",",header=TRUE,row.names=1)

## Create ASV table 
OTU_NOMIS <- otu_table(asv_unrarefied, taxa_are_rows=TRUE)
tax_NOMIS <- tax_table(as.matrix(tax_unrarefied))

## Load metadata file
metadata_NOMIS="202402_NOMIS_metadata_GFS.tsv"
metadata_NOMIS<-import_qiime_sample_data(metadata_NOMIS)

## Create phyloseq object
merged_NOMIS_DEBLUR_unrarefied <- merge_phyloseq(OTU_NOMIS, tax_NOMIS, metadata_NOMIS)
asv_unrarefied_nomis <- otu_table(merged_NOMIS_DEBLUR_unrarefied, taxa_are_rows=T)

merge_nomis_cryo <- merge(asv_unrarefied_nomis, asv_table_NOMISCRYO, by.x="row.names", by.y="Feature_ID", all=T)

## Only keep samples from cryospheric ecosystems
## Load metadata CRYO + NOMIS
metadata_merge_cryonom="202403_metadata_merge_nomis_cryo_deblur.txt"
metadata_cryo_nomis <-import_qiime_sample_data(metadata_merge_cryonom)

merge_nomis_cryo[is.na(merge_nomis_cryo)] <- 0

#saveRDS(merge_nomis_cryo,"20240301_merge_nomis_cryo.RDS")
#merge_nomis_cryo <- readRDS("merge_nomis_cryo.RDS")

asv_table_NOMISCRYO = as.data.frame(merge_nomis_cryo)

unAsIs <- function(X) {
   if("AsIs" %in% class(X)) {
     class(X) <- class(X)[-match("AsIs", class(X))]
   }
  X
}

asv_table_NOMISCRYO$Row.names<-unAsIs(asv_table_NOMISCRYO$Row.names)
rownames(asv_table_NOMISCRYO) <- asv_table_NOMISCRYO$Row.names
asv_table_NOMISCRYO$Row.names <- NULL
asv_table_NOMISCRYO_m <- as.matrix(asv_table_NOMISCRYO)
asv_tableNOMISCRYO_final <- otu_table(asv_table_NOMISCRYO_m, taxa_are_rows=T)

merge_asvnomiscryo_metadata <- merge_phyloseq(asv_tableNOMISCRYO_final, metadata_cryo_nomis)

eco=c("Yes")
Cryosphere_sample <- subset_samples(merge_asvnomiscryo_metadata, Cryosphere %in% eco)
#2024
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 79358 taxa and 431 samples ]
# sample_data() Sample Data:       [ 431 samples by 8 sample variables ]
#saveRDS(Cryosphere_sample,"20240301_Cryosphere_sample.RDS")

Cryosphere_sample <- readRDS("20240301_Cryosphere_sample.RDS")
prune_nomiscryo <- prune_samples(sample_sums(Cryosphere_sample) >20000,Cryosphere_sample)

NOMIS_CRYOS<- rarefy_even_depth(prune_nomiscryo, sample.size=min(sample_sums(prune_nomiscryo)), rngseed=678, replace=F, trimOTUs=TRUE, verbose=TRUE)

#saveRDS(NOMIS_CRYOS,file="20240301_NOMISCRYOS.RData")
NOMIS_CRYOS<-readRDS("20240301_NOMISCRYOS.RData")
meta <- sample.data.frame(NOMIS_CRYOS)

NOMIS_CRYOS_sub <- subset_samples(NOMIS_CRYOS, Habitat != "Rock")
NOMIS_CRYOS_sub <- subset_samples(NOMIS_CRYOS_sub, Habitat != "Ice/Snow")
NOMIS_CRYOS_sub <- subset_samples(NOMIS_CRYOS_sub, Habitat !="Water_lake")
NOMIS_CRYOS_sub <- subset_samples(NOMIS_CRYOS_sub, Habitat !="Water")
NOMIS_CRYOS_sub <- subset_samples(NOMIS_CRYOS_sub, Habitat !="Sediment")

asv_nomis_cryo <- as.data.frame(otu_table(NOMIS_CRYOS_sub,taxa_are_rows=T))
asv_nomis_cryo <- asv_nomis_cryo[rowSums(asv_nomis_cryo)>0,]

#write.csv(asv_nomis_cryo, "202408_asv_nomis_cryo.csv")

saveRDS(NOMIS_CRYOS_sub, file="20240301_NOMISCRYOS_sub.RData")

## Calculate ASV richness
diversity_cryosphere <- plot_richness(NOMIS_CRYOS_sub, measures=c("Observed","Shannon"), color="Habitat", shape="Ecosystem")
alphadt_cryosphere <-data.table(diversity_cryosphere$data)

OR_cryo <- subset(alphadt_cryosphere, variable == "Observed")
ASVrichness_cryo <- OR_cryo %>% 
  group_by(Habitat) %>% 
  summarise(average=mean(value), std=sd(value))

median_OR_cryosphere <- OR_cryo %>%
  group_by(Habitat) %>% 
  summarize(median_x=median(value)) 

## Filter out the GFS habitat and its corresponding median_x value
GFS_median <- median_OR_cryosphere %>%
  filter(Habitat == "GFS") %>%
  pull(median_x)

## Calculate the ratios of median_x of GFS to median_x of other habitats
ratios <- median_OR_cryosphere %>%
  filter(Habitat != "GFS") %>%
  mutate(ratio_GFS_to_other = GFS_median / median_x)

Shannon_cryo <- subset(alphadt_cryosphere, variable == "Shannon")
Shannondiv_cryo <- Shannon_cryo %>% 
  group_by(Habitat) %>% 
  summarize(median=median(value))

## Filter out the GFS habitat and its corresponding median Shannon value
GFS_median_shannon <- Shannondiv_cryo %>%
  filter(Habitat == "GFS") %>%
  pull(median)

## Calculate the ratios of median Shannon values of GFS to median Shannon values of other habitats
ratios_shannon <- Shannondiv_cryo %>%
  filter(Habitat != "GFS") %>%
  mutate(ratio_GFS_to_other_shannon = GFS_median_shannon / median)

## Hill numbers 
hill_q1_cryo<-as.data.frame(hill_taxa(t(asv_nomis_cryo), q = 1))
hill_q1_indices_cryo <- phyloseq::sample_data(hill_q1_cryo)

nomis_cryo_hill <- phyloseq::merge_phyloseq(NOMIS_CRYOS_sub,hill_q1_indices_cryo)
meta_diversity_cryo_hill <- sample.data.frame(nomis_cryo_hill)

median_shannon_hill_cryo <- meta_diversity_cryo_hill[,-10]

median_shannon_hill_cryo_med<-median_shannon_hill_cryo %>% 
  summarize(median=median(hill_taxa.t.asv_nomis_cryo...q...1.), x = quantile(hill_taxa.t.asv_nomis_cryo...q...1., c(0.25, 0.5, 0.75)))

# Filter out the GFS habitat and its corresponding median Shannon H value
GFS_median_shannon_entropy <- median_shannon_hill_cryo %>%
  filter(Habitat == "GFS") %>%
  summarize(med_gfs=median(hill_taxa.t.asv_nomis_cryo...q...1.))

# Calculate the ratios of median Shannon H values of GFS to median Shannon H values of other habitats
ratios_shannon <- meta_diversity_cryo_hill  %>%
  filter(Habitat != "GFS") %>%
  group_by(Habitat)%>%
  summarize(med_cryo=median(hill_taxa.t.asv_nomis_cryo...q...1.)) %>%
  mutate(ratio_GFS_to_other_shannon = GFS_median_shannon_entropy$med_gfs / med_cryo)

## NMDS cryospheric ecosystems
metadata_nmds_cryo <- sample.data.frame(sample_data(NOMIS_CRYOS_sub))
#write.csv(metadata_nmds_cryo, "202408_metadata_nmds_cryo.csv")
ord <- ordinate(NOMIS_CRYOS_sub, method = "NMDS", distance = "bray", trymax = 999, k=2)
stressplot(ord)
ord$stress

data.scores = as.data.frame(scores(ord)$sites)

#add columns to data frame 
data.scores$Sample = metadata_nmds_cryo$Sample
data.scores$Site = metadata_nmds_cryo$Habitat
head(data.scores)
#write.csv(data.scores,"2024_datascores_nmds_cryo.csv")

##plot NMDS
##k=2
nmds_bc_plot = ggplot(data.scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 6, aes(color = Site), alpha = 0.6) + 
  theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 12), 
        legend.text = element_text(size = 12, face = "bold", colour = "black"), 
        legend.position = "right", 
        axis.title.y = element_text(face = "bold", size = 14), 
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
        legend.title = element_text(size = 14, colour = "black", face = "bold"), 
        panel.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.2), # Use linewidth instead of size
        legend.key = element_blank()) + 
  labs(x = "NMDS1", y = "NMDS2", colour = "Biomes", shape = "Type") 

## Testing the effects of the habitat on community structure
vegan_otu <- function(physeq){
  OTU <- otu_table(physeq)
  if(taxa_are_rows(OTU)){
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}

vegan_matrix_cryosphere <- vegan_otu(NOMIS_CRYOS_sub)
NOMISCRYO_bray <- vegdist(vegan_matrix_cryosphere, method="bray")
betadisp_cryo <- betadisper(NOMISCRYO_bray, metadata_nmds_cryo$Habitat, type="centroid")
betadisp_cryo_perm <- permutest(betadisp_cryo, permutations=999)

ado_habitat <- adonis2(NOMISCRYO_bray ~ Habitat, data=metadata_nmds_cryo, method="bray", permutations=999)

## Pairwise adonis
library(pairwiseAdonis)
pairwise_habitat <- pairwise.adonis(NOMISCRYO_bray, metadata_nmds_cryo$Habitat, p.adjust.m="holm")

## Taxonomy
asv_table_NOMISCRYO = as.data.frame(fread("cryobiome_table.txt", header=T))
rownames(asv_table_NOMISCRYO)<-asv_table_NOMISCRYO$Feature_ID

asv_table_NOMISCRYO$Feature_ID <- NULL
asv_table_NOMISCRYO <- as.matrix(asv_table_NOMISCRYO)

## Load metadata CRYO
metadata_cryo=(("PP1_metadata_date.tsv"))
metadata_cryo<-import_qiime_sample_data(metadata_cryo)
metadata_cryo$Habitat[metadata_cryo$Habitat %in% "Mat"] <- "Microbial mat"

tax_cryo <- as.data.frame(read_tsv("cryobiome_taxonomy.tsv"),sep="\t")
rownames(tax_cryo)<-tax_cryo$Feature_ID
tax_cryo$Feature_ID <- NULL
tax_cryo_m <- tax_table(as.matrix(tax_cryo))

asv_cryo <- otu_table(asv_table_NOMISCRYO, taxa_are_rows=T)
merge_cryo_taxo <- merge_phyloseq(asv_cryo, tax_cryo_m, metadata_cryo)

merge_cryo_taxo_sub <- subset_samples(merge_cryo_taxo, Habitat != "Rock")
merge_cryo_taxo_sub <- subset_samples(merge_cryo_taxo_sub, Habitat !="Ice/Snow")
merge_cryo_taxo_sub <- subset_samples(merge_cryo_taxo_sub, Habitat !="Water_lake")
merge_cryo_taxo_sub <- subset_samples(merge_cryo_taxo_sub, Habitat !="Water")
merge_cryo_taxo_sub <- subset_samples(merge_cryo_taxo_sub, Habitat !="Sediment")

cryo_taxglom <- tax_glom(merge_cryo_taxo_sub, taxrank=rank_names(merge_cryo_taxo_sub)[6], NArm=F)
transf_cryo <- transform_sample_counts(cryo_taxglom, function(x) x / sum(x))
sample_merge_habitat <- merge_samples(transf_cryo, "Habitat")
region_cryo <- transform_sample_counts(sample_merge_habitat, function(x) x / sum(x))

TopASV_f <- names(sort(taxa_sums(region_cryo), TRUE)[1:25])
top15_NOMIS_f <- prune_species(TopASV_f, region_cryo)
top15_NOMIS_f <- prune_taxa(taxa_sums(top15_NOMIS_f)>0, top15_NOMIS_f)
top_genus<-as.data.frame(tax_table(top15_NOMIS_f))

## Turn all ASVs into Genus counts
cryo_df <- psmelt(region_cryo) # create dataframe from phyloseq object
cryo_df$Genus<- as.character(cryo_df$Genus) # convert to character

## Only Genera that are the most abundant
cryo_df$Genus[!(cryo_df$Genus %in% top_genus$Genus)] <- "Other"
cryo_df$Genus[(cryo_df$Genus == "")] <- "Other"
cryo_df$Genus <- ifelse(is.na(cryo_df$Genus), "Other", cryo_df$Genus)
cryo_df$Genus[(cryo_df$Genus == "g__uncultured")] <- "Other"

## Generating n colors via rbrewer
n <- 20
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

barplot_biogeo <- ggplot(data=cryo_df, aes(x=Sample, y=Abundance, fill=Genus))
barplot_biogeo + geom_bar(aes(), stat="identity", position="stack") +
  scale_fill_manual(values = c("#7FC97F", "#BEAE44", "#F08783", "#FFFF99", "#386CB0", "#F0397F", "#BF5B17", "#1B9E77", "#D95F02",  "#283983","#66A61E",
                               "#E6AB02","#A6761D", "#666666", "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#7579E8", "#FDBF6F", "#FF7F00", "#CAB2D6"
                              )) +
  theme(legend.position="bottom") + guides(fill=guide_legend(ncol=3)) + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                                                                         panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
#write.csv(cryo_df, "202408_cryo_df_taxo.csv")


