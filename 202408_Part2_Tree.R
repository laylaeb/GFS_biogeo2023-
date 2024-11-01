### Taxonomic Tree to represent GFS microbiome ###

library(phyloseq)
library(phyloseqCompanion)
library(vegan)
library(metacoder)
library(ggplot2)
library(dplyr)
library(readr)
library(stringr)
library(agricolae)
library(ape)


## Load ASV and taxonomy tables
asv<-read.csv(file="20240221_NOMIS_rarefied_deblur_table.csv.gz",sep=",",header=TRUE,row.names=1)
tax<-read.csv(file="20240221_NOMIS_rarefied_deblur_taxonomy.csv.gz",sep=",",header=TRUE,row.names=1)

## Load metadata file
metadata_NOMIS="202402_NOMIS_metadata_GFS.tsv"
metadata_NOMIS<-import_qiime_sample_data(metadata_NOMIS)

## Create ASV table
OTU_NOMIS <- otu_table(asv, taxa_are_rows=TRUE)
tax_NOMIS <- tax_table(as.matrix(tax))

## Create phyloseq object
merged_NOMIS_DEBLUR_taxo_tree <- merge_phyloseq(OTU_NOMIS, tax_NOMIS, metadata_NOMIS)

asv_nomis <- otu_table(merged_NOMIS_DEBLUR_taxo_tree, taxa_are_rows=T)
tax_nomis <- tax_table(merged_NOMIS_DEBLUR_taxo_tree)

#write.csv(asv_nomis, "asv_nomis_metacoder.csv")
#write.csv(tax_nomis, "tax_nomis_metacoder.csv")

## Build the object that could then be used by metacoder including sample from Uganda!
asv_taxo_metacoder <- read.csv("20240303_metacoder.csv",header=T)
asv_metacoder_tibble <- as.tibble(asv_taxo_metacoder)

## Load metadata
meta_glaciers_tibble <- as.tibble(sample.data.frame(merged_NOMIS_DEBLUR_taxo_tree))

obj_tree_nomis <- parse_tax_data(asv_metacoder_tibble,
                      class_cols = "lineage",
                      class_sep = "; ", ##careful here because there is a space also!
                      class_regex = "^(.+)__(.+)$", 
                      class_key = c("tax_rank" = "taxon_rank", "name" = "taxon_name"))

no_reads <- rowSums(obj_tree_nomis$data$tax_data[, meta_glaciers_tibble$sample]) == 0
sum(no_reads)

## filter some specific taxa
obj_tree_nomis <- filter_taxa(obj_tree_nomis, taxon_names != "")
obj_tree_nomis <- filter_taxa(obj_tree_nomis, taxon_names != "uncultured")
obj_tree_nomis <- filter_taxa(obj_tree_nomis, taxon_names != "uncultured_bacterium")
obj_tree_nomis <- filter_taxa(obj_tree_nomis, taxon_names != "Bacteria;;;;;")


obj_tree_noms <- obj_tree_nomis %>% 
filter_taxa(taxon_ranks == "g", supertaxa=T)# subset to the order rank
  
## Getting per-taxon information
## "e can sum the abundance per-taxon and add the results to the tax-map
obj_tree_nomis$data$tax_abund <- calc_taxon_abund(obj_tree_nomis, "tax_data",
                                       cols = meta_glaciers_tibble$sample)


## Calculate the nb of samples that have reads for each taxon
obj_tree_nomis$data$tax_occ <- calc_n_samples(obj_tree_nomis, "tax_abund", cols = meta_glaciers_tibble$sample)


## Plot taxonomic tree
set.seed(3) 
heat_tree_nomis <- heat_tree(obj_tree_nomis, 
                            node_color=n_obs,
                            node_size=n_obs,
                            node_label = taxon_names,
                            edge_color_range = c("#CC8394FF", "#AC563BFF", "#CDA97CFF","#7C8EC5FF","#2B3C51FF"),
                            node_color_range = c("#CC8394FF", "#AC563BFF", "#CDA97CFF","#7C8EC5FF","#2B3C51FF"),
                            edge_color=n_samples,
                            initial_layout = "re", layout = "da")
heat_tree_nomis
