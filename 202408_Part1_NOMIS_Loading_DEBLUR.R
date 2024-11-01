## Loading data & Processing  - Filtering, averaging and rarefaction

library(phyloseq)
library(vegan)
library(data.table)
library(speedyseq)
library(phyloseqCompanion)
library(dplyr)
library(tidyverse)

asv<-fread(file="081623_NOMIS_16S_merged_deblur_table.csv",header=TRUE, sep=",")
newnames <- colnames(asv)
asv_df<-as.data.frame(asv)
colnames(asv_df) <- newnames
rownames(asv_df) <- asv_df$Feature_ID
asv_df$Feature_ID <- NULL

tax<-read.csv(file="081623_NOMIS_16S_merged_deblur_taxonomy.csv",sep=",",header=TRUE,row.names=1)

## Import QIIME pipeline map-file
metadata_NOMIS="20230508_NOMIS_metadata_deblur.txt"
metadata_NOMIS<-import_qiime_sample_data(metadata_NOMIS)

OTU_NOMIS <- otu_table((asv_df), taxa_are_rows=TRUE)
tax_NOMIS <- tax_table(as.matrix(tax))

## Create phyloseq object
merged_NOMIS_DEBLUR <- merge_phyloseq(OTU_NOMIS, tax_NOMIS, metadata_NOMIS)

nomis_meta <- sample.data.frame(merged_NOMIS_DEBLUR)

## Subsetting taxa
merged_NOMIS_sub <- subset_taxa(merged_NOMIS_DEBLUR, (Kingdom!="d__Eukaryota") | is.na(Kingdom)) 
merged_NOMIS_sub <- subset_taxa(merged_NOMIS_sub, (Kingdom!="d__Archaea") | is.na(Kingdom))
merged_NOMIS_sub <- subset_taxa(merged_NOMIS_sub, (Order!=" o__Chloroplast") )
merged_NOMIS_sub <- subset_taxa(merged_NOMIS_sub, (Family!=" f__Mitochondria"))

## Subset blank 
merged_NOMIS_blank <- subset_samples(merged_NOMIS_sub, sample_type %in% "Blank")
asv_blank <- otu_table(merged_NOMIS_blank)
asv_blank <- as.data.frame(asv_blank)

asv_blank <- asv_blank[rowSums(asv_blank[])>0,]

pop_taxa = function(physeq, badTaxa){
  allTaxa = taxa_names(physeq)
  allTaxa <- allTaxa[!(allTaxa %in% badTaxa)]
  return(prune_taxa(allTaxa, physeq))
}

## Dataset with no blank reads
merged_NOMIS_noblank <- pop_taxa(merged_NOMIS_sub, rownames(asv_blank))

## Select Substrate 
site = c("sed")
sediment <- subset_samples(merged_NOMIS_noblank, sample_type %in% site)

## ASV table
asv_sedup_unfiltered <- otu_table(sediment)
asv_sedup_unfiltered <- as.data.frame(asv_sedup_unfiltered)
asv_sedup_unfiltered <- asv_sedup_unfiltered[rowSums(asv_sedup_unfiltered[])>0,]

## Filtering step - removes ASVs that were not found in at least 2 of biological replicates
set.seed(132)
filtering_br <- function(phyloseq_obj){
  meta_sediment <- as.data.frame(sample_data(phyloseq_obj))
  asv_sediment <- as.data.frame(otu_table(phyloseq_obj, taxa_are_rows = T))
  
  new_asv <- c()
  for (Sample in meta_sediment$gl_code) {
    sub_asvsed <- asv_sediment[,startsWith(colnames(asv_sediment), paste0(Sample,'_',collapse = ''))]
    if (is.null(sub_asvsed)){}
    else{
      if (!(is.vector(sub_asvsed))){
        sub_asv_preval <- rowSums(sub_asvsed > 0)
        new_asv = c(new_asv, rownames(sub_asvsed)[sub_asv_preval >= 2])} 
      else {
        new_asv = c(new_asv, names(sub_asvsed)[sub_asvsed > 0])}
    }}
  new_asv = unique(new_asv)
  return(new_asv)
}

## Apply filtering on sediment
filtered_asvs = filtering_br(sediment)

## Subset ASVs that were filtered from ASV table.
filtered_table = asv_sedup_unfiltered[rownames(asv_sedup_unfiltered) %in% filtered_asvs, ]
filtered_asv_table<- otu_table(filtered_table, taxa_are_rows=T)

########################################### How many sample(s) per patch?
## Get column names
col_names <- colnames(filtered_asv_table)
count_prefixes <- function(prefix, col_names) {
  count <- sum(grepl(paste0("^", prefix, "_"), col_names))
  return(count)
}

## List of unique prefixes
unique_prefixes <- unique(sub("_.*", "", col_names))

## Create data frame to store results
result_df <- data.frame(Prefix = unique_prefixes, Count = numeric(length(unique_prefixes)))

## Iterate through prefixes and count
for (i in 1:length(unique_prefixes)) {
  prefix <- unique_prefixes[i]
  count <- count_prefixes(prefix, col_names)
  result_df[i, "Count"] <- count
}

## Remove samples with less than 3 ecological replicates
result_df_patch <- result_df %>%
  filter(Count < 3)

######################################
## melt ASV_table
melt_asv <- reshape2::melt(filtered_asv_table)
code_glaciers <- sample.data.frame(nomis_meta[,c("sample","gl_code")])
rownames(code_glaciers)<-NULL

## merge ASV_table and gl_code
melt_merge_as <- melt_asv %>%
  left_join(code_glaciers, join_by(Var2 == sample))%>%
  filter(!(gl_code %in% result_df_patch$Prefix))

## take average
melt_merge_asv <- melt_merge_as %>%
  group_by(gl_code, Var1) %>% 
  summarise(average=mean(value))

## Re-build ASV table
dcast_asv <- reshape2::dcast(melt_merge_asv, formula= Var1 ~ gl_code)
rownames(dcast_asv) <- dcast_asv$Var1
dcast_asv$Var1 <- NULL

min_non_zero = min(dcast_asv[dcast_asv > 0])
dcast_asv = as.data.frame(dcast_asv)*1/min_non_zero
dcast_asv_n = as.data.frame(sapply(dcast_asv, round))
dcast_asv_n = as.data.frame(sapply(dcast_asv_n, as.integer))
rownames(dcast_asv_n) <- rownames(dcast_asv)
asv_table_dcast <-(otu_table(dcast_asv_n, taxa_are_rows=T))

## Re-build phyloseq object --- metadata with glcode
metadata_glaciers="202402_NOMIS_metadata_GFS.tsv"
metadata_glaciers<-import_qiime_sample_data(metadata_glaciers)
tax_sed <- as.matrix(tax_table(sediment))
merge_GFS_full <- merge_phyloseq(asv_table_dcast, tax_sed, metadata_glaciers)

prune_GFS<- prune_samples(sample_sums(merge_GFS_full)>120000,merge_GFS_full)

## save filtered tables
asv_filtered <- otu_table(prune_GFS, taxa_are_rows=T)

## rarefaction
NOMIS_rarefied <- rarefy_even_depth(prune_GFS, sample.size=min(sample_sums(prune_GFS)), rngseed=678, replace=F, trimOTUs=TRUE, verbose=TRUE)

## save rarefied tables
asv_rarefied <- otu_table(NOMIS_rarefied, taxa_are_rows=T)
tax_table_rarefied <- tax_table(NOMIS_rarefied)
#write.csv(asv_rarefied, "20240221_NOMIS_rarefied_deblur_table.csv")
#write.csv(tax_table_rarefied, "20240221_NOMIS_rarefied_deblur_taxonomy.csv")

saveRDS(NOMIS_rarefied,"20240221_NOMIS_rarefied.RData")
