library(tidyverse)
library(purrr)
library(ape)
library(vegan)
library(ggpubr)
library(foreach)
library(doMC)
library(reshape2)

# Functions

filtering_br <- function(asv_tab, meta_tab){
  new_asv <- c()
  for (sample in meta_tab$sample) {
    sub_asvsed <- asv_tab[,colnames(asv_tab) == sample]
    if (!(is.vector(sub_asvsed))){
        sub_asv_preval <- rowMeans(sub_asvsed > 0)
        new_asv = c(new_asv, rownames(sub_asvsed)[sub_asv_preval >= 0.5])}
    else {
        new_asv = c(new_asv, names(sub_asvsed)[sub_asvsed > 0])}}
  new_asv = unique(new_asv)
  return(new_asv)}

rarefy_tab <- function(x, n){
  # Wrapper around the vegan rrarefy function (removes ASVs with no observations after rarefaction)
  rarefied_tab <- as.data.frame(as.data.frame(t(rrarefy(t(x), n))))
  rarefied_tab = rarefied_tab[rowSums(rarefied_tab) > 0,]
  return(rarefied_tab)}

select_one_replicate <- function(x){
  # Function that creates an ASV table with one randomly drawn replicate for each sample
  x = x[,sample(ncol(x))] # randomly shuffle to randomise the replicates kept
  colnames(x) = map_chr(colnames(x), function(x) paste0(strsplit(x, '\\.')[[1]][1])) # the last line adds suffixes to the colnames, this one removes them
  new_x = data.frame(row.names = rownames(x))
  # main for loop to create new data
  for (sample in unique(colnames(x))){
    sample_data = as.data.frame(x[,colnames(x) == sample])
    if (ncol(sample_data) > 1){new_x = cbind(new_x, sample_data[,sample(ncol(sample_data), 1)])}
    else {new_x = cbind(new_x, sample_data)}}
  colnames(new_x) = unique(colnames(x))
  # remove ASVs with no observations after rarefaction
  new_x = new_x[rowSums(new_x) > 0,]
  return(new_x)}

take_sum <- function(x){
  # Function that creates an ASV table with summed replicates for each sample
  new_x = data.frame(row.names = rownames(x))
  # main for loop to create new data
  for (sample in unique(colnames(x))){
    sample_data = as.data.frame(x[,colnames(x) == sample])
    if (ncol(sample_data) > 1){new_x = cbind(new_x, rowSums(sample_data))}
    else {new_x = cbind(new_x, sample_data)}}
  colnames(new_x) = unique(colnames(x))
  # remove ASVs with no observations after rarefaction
  new_x = new_x[rowSums(new_x) > 0,]
  return(new_x)}

take_avg <- function(x){
  # Function that creates an ASV table with summed replicates for each sample
  new_x = data.frame(row.names = rownames(x))
  # main for loop to create new data
  for (sample in unique(colnames(x))){
    sample_data = as.data.frame(x[,colnames(x) == sample])
    if (ncol(sample_data) > 1){new_x = cbind(new_x, rowMeans(sample_data))}
    else {new_x = cbind(new_x, sample_data)}}
  colnames(new_x) = unique(colnames(x))
  # remove ASVs with no observations after rarefaction
  new_x = new_x[rowSums(new_x) > 0,]
  # Set lowest count to 1
  #min_non_zero = min(new_x[new_x > 0])
  #new_x = min_non_zero * as.data.frame(new_x)
  #new_x = as.data.frame(sapply(new_x, as.integer))
  #melt_x = melt(as.matrix(x))
  #avg_x <- melt_x %>%
  #  group_by(Var1, Var2) %>% 
  #  summarise(average=mean(value)) %>% filter(average > 0)
  #dcast_x <- as.data.frame(dcast(avg_x, formula= Var1 ~ Var2, value.var = 'average', fill = 0))
  #
  #rownames(dcast_x) = dcast_x$Var1
  #dcast_x$Var1 = NULL
  
  min_non_zero = min(new_x[new_x > 0])
  new_x = as.data.frame(new_x) / min_non_zero
  new_x = as.data.frame(sapply(new_x, round))
  new_x = as.data.frame(sapply(new_x, as.integer))
  return(new_x)}

filter_taxa_samples <- function(tab, tax){
  # filter taxa
  if ('Taxon' %in% colnames(tax)){
       ASVs_to_keep = tax %>% filter(!grepl("d__Eukaryota", Taxon) , !grepl("d__;", Taxon), !grepl("d__Archaea", Taxon),
                                     !grepl(" o__Chloroplast", Taxon) , !grepl(" f__Mitochondria", Taxon)) %>% select(Feature_ID)}
  else{ASVs_to_keep = tax %>% filter(Kingdom != "d__Eukaryota", Kingdom != '', Kingdom != "d__Archaea", 
                                     Order != " o__Chloroplast", Family != " f__Mitochondria") %>% select(Feature_ID)}
  ASVs_to_remove = rownames(tab)[rowSums(tab[,colnames(tab) %in% c("NC_16S_KAUST2_rock","BLEX1","BLEX2","oldBL1","oldBL2")] > 0)]
  new_tab = tab[!(rownames(tab) %in% ASVs_to_remove),]
  new_tab = tab[rownames(tab) %in% ASVs_to_keep$Feature_ID,]
  
  # filter samples
  samples_to_keep = colnames(new_tab)[(grepl('_UP_', colnames(new_tab)) & grepl('16S_sed', colnames(new_tab)))]
  new_tab = new_tab[,colnames(new_tab) %in% samples_to_keep]
  
  # remove empty ASVs
  new_tab = new_tab[rowSums(new_tab) > 0,]
  return(new_tab)}

load_asv_data <- function(min_reads_cutoff){
  setwd('~/Desktop/leila_data')
  asv_tab_dada2 = read.csv('data/20230508_NOMIS_16S_merged_dada2_table.csv', sep=',')
  asv_tax_dada2 = read.csv('data/20230508_NOMIS_16S_merged_dada2_taxonomy.csv', sep=',')
  rownames(asv_tab_dada2) = asv_tab_dada2$Feature_ID
  asv_tab_dada2 = asv_tab_dada2 %>% select(-GL140_UP_1_16S_sed, -GL140_UP_2_16S_sed, -GL140_UP_3_16S_sed, -Feature_ID)
  asv_tab_dada2 = filter_taxa_samples(asv_tab_dada2, asv_tax_dada2)
  
  asv_tab_unoise = read.csv('data/unoise3_zotu_table_raw.txt', sep='\t')
  asv_tax_unoise = read.csv('data/unoise3_taxonomy_raw.tsv', sep='\t')
  rownames(asv_tab_unoise) = asv_tab_unoise$X.OTU.ID
  asv_tab_unoise = asv_tab_unoise %>% select(-GL140_UP_1_16S_sed, -GL140_UP_2_16S_sed, -GL140_UP_3_16S_sed, -X.OTU.ID)
  asv_tab_unoise = filter_taxa_samples(asv_tab_unoise, asv_tax_unoise)
  
  asv_tab_deblur = read.csv('data/081623_NOMIS_16S_merged_deblur_table.csv', sep=',')
  asv_tax_deblur = read.csv('data/081623_NOMIS_16S_merged_deblur_taxonomy.csv', sep=',')
  rownames(asv_tab_deblur) = asv_tab_deblur$Feature_ID
  asv_tab_deblur = asv_tab_deblur %>% select(-GL140_UP_1_16S_sed, -GL140_UP_2_16S_sed, -GL140_UP_3_16S_sed, -Feature_ID)
  asv_tab_deblur = filter_taxa_samples(asv_tab_deblur, asv_tax_deblur)
  
  asv_tab_dada2 = asv_tab_dada2[names(asv_tab_unoise)]
  asv_tab_deblur = asv_tab_deblur[names(asv_tab_unoise)]
  samples_to_keep = colnames(asv_tab_unoise)[((colSums(asv_tab_unoise) > min_reads_cutoff) & (colSums(asv_tab_dada2) > min_reads_cutoff) & (colSums(asv_tab_deblur) > min_reads_cutoff))]
  asv_tab_dada2 = asv_tab_dada2 %>% select(any_of(samples_to_keep))
  asv_tab_unoise = asv_tab_unoise %>% select(any_of(samples_to_keep))
  asv_tab_deblur = asv_tab_deblur %>% select(any_of(samples_to_keep))
  
  colnames(asv_tab_dada2) = map_chr(colnames(asv_tab_dada2), function(x) paste0(strsplit(x, '_')[[1]][1:2], collapse = '_'))
  colnames(asv_tab_unoise) = map_chr(colnames(asv_tab_unoise), function(x) paste0(strsplit(x, '_')[[1]][1:2], collapse = '_'))
  colnames(asv_tab_deblur) = map_chr(colnames(asv_tab_deblur), function(x) paste0(strsplit(x, '_')[[1]][1:2], collapse = '_'))
  return(list(dada2=asv_tab_dada2, unoise=asv_tab_unoise, deblur=asv_tab_deblur))}

load_meta_data <- function(asv_tab){
  meta_tab = read.csv('data/nomis-20230320-0855-db.csv')
  meta_tab$sample = map_chr(meta_tab$patch, function(x) paste0(strsplit(x, '_')[[1]][1:2], collapse = '_'))
  meta_tab$rep_n = map_int(meta_tab$sample, function(x) sum(colnames(asv_tab) == x))
  meta_tab = meta_tab %>% select(sample, mountain_range, lat_sn..DD., lon_sp..DD., rep_n) %>% 
    filter(mountain_range != 'Rwenzori Mountains', grepl('_UP', sample)) %>% distinct()
  
  meta_tab = meta_tab[meta_tab$sample %in% colnames(asv_tab),]
  return(meta_tab)}

endemic_analysis <- function(x, meta_tab, n_sample, rarefaction_depth){
  colnames(x) = map_chr(colnames(x), function(y) unique(meta_tab$mountain_range[meta_tab$sample == y]))
  
  # Group by mountain range
  new_x = data.frame(row.names = rownames(x))
  for (mountain_range in unique(colnames(x))){
    samples = sample(meta_tab$sample[meta_tab$mountain_range == mountain_range], n_sample)
    mr_data = x[,colnames(x) == mountain_range]
    new_x = cbind(new_x, mr_data[,sample(ncol(mr_data), n_sample)])}
  colnames(new_x) = map_chr(colnames(new_x), function(x) paste0(strsplit(x, '\\.')[[1]][1]))
  new_x = as.data.frame(do.call(cbind, by(t(new_x),INDICES=names(new_x),FUN=colSums)))
  new_x = new_x[rowSums(new_x) > 0,]
  
  res = data.frame()
  for (mountain_range in colnames(new_x)){
    endemic_mask = (new_x[,mountain_range] == rowSums(new_x))
    tot_asv_n = sum(new_x[,mountain_range] > 0)
    end_asv_n = sum((new_x[,mountain_range] > 0) & endemic_mask)
    end_asv_ab = sum(new_x[endemic_mask, mountain_range]) / sum(new_x[,mountain_range])
    res = rbind(res, data.frame(mountain_range=mountain_range, rarefaction_depth=rarefaction_depth,
                                tot_asv_n=tot_asv_n, end_asv_n=end_asv_n, end_asv_ab=end_asv_ab, sample_n=n_sample))}
  return(res)}

analyse_endemics <- function(asv_tab_dada2, asv_tab_unoise, asv_tab_deblur, asv_filtered_tab_dada2, asv_filtered_tab_unoise, asv_filtered_tab_deblur, meta_data, rd){
  all_res_endemics = data.frame()
  
  for (i in 1:5){
    set.seed(i)
    
    # Create summed table, and random replicate table
    tab_sum_dada2 = take_sum(asv_tab_dada2)
    tab_sample_dada2 = select_one_replicate(asv_tab_dada2)
    tab_filtering_dada2 = take_avg(asv_filtered_tab_dada2)
    
    tab_sum_unoise = take_sum(asv_tab_unoise)
    tab_sample_unoise = select_one_replicate(asv_tab_unoise)
    tab_filtering_unoise = take_avg(asv_filtered_tab_unoise)
    
    tab_sum_deblur = take_sum(asv_tab_deblur )
    tab_sample_deblur = select_one_replicate(asv_tab_deblur)
    tab_filtering_deblur = take_avg(asv_filtered_tab_deblur)
    
    # Rarefy the tables
    rarefied_tab_sum_dada2 = rarefy_tab(tab_sum_dada2, rd)
    rarefied_tab_sample_dada2 = rarefy_tab(tab_sample_dada2, rd)
    rarefied_tab_filtered_dada2 = rarefy_tab(tab_filtering_dada2, rd)
    
    rarefied_tab_sum_unoise = rarefy_tab(tab_sum_unoise, rd)
    rarefied_tab_sample_unoise = rarefy_tab(tab_sample_unoise, rd)
    rarefied_tab_filtered_unoise = rarefy_tab(tab_filtering_unoise, rd)
    
    rarefied_tab_sum_deblur = rarefy_tab(tab_sum_deblur, rd)
    rarefied_tab_sample_deblur = rarefy_tab(tab_sample_deblur, rd)
    rarefied_tab_filtered_deblur = rarefy_tab(tab_filtering_deblur, rd)
    
    # Endemic analysis
    res_sum_dada2 = endemic_analysis(rarefied_tab_sum_dada2, meta_data, 6, rd)
    res_sample_dada2 = endemic_analysis(rarefied_tab_sample_dada2, meta_data, 6, rd)
    res_filtered_dada2 = endemic_analysis(rarefied_tab_filtered_dada2, meta_data, 6, rd)
    res_sum_dada2$Method = 'dada2'
    res_sample_dada2$Method = 'dada2'
    res_filtered_dada2$Method = 'dada2'
    res_sum_dada2$Filtering = 'sum'
    res_sample_dada2$Filtering = 'sample'
    res_filtered_dada2$Filtering = 'filtering'
    
    res_sum_unoise = endemic_analysis(rarefied_tab_sum_unoise, meta_data, 6, rd)
    res_sample_unoise = endemic_analysis(rarefied_tab_sample_unoise, meta_data, 6, rd)
    res_filtered_unoise = endemic_analysis(rarefied_tab_filtered_unoise, meta_data, 6, rd)
    res_sum_unoise$Method = 'unoise'
    res_sample_unoise$Method = 'unoise'
    res_filtered_unoise$Method = 'unoise'
    res_sum_unoise$Filtering = 'sum'
    res_sample_unoise$Filtering = 'sample'
    res_filtered_unoise$Filtering = 'filtering'
    
    res_sum_deblur = endemic_analysis(rarefied_tab_sum_deblur, meta_data, 6, rd)
    res_sample_deblur = endemic_analysis(rarefied_tab_sample_deblur, meta_data, 6, rd)
    res_filtered_deblur = endemic_analysis(rarefied_tab_filtered_deblur, meta_data, 6, rd)
    res_sum_deblur$Method = 'deblur'
    res_sample_deblur$Method = 'deblur'
    res_filtered_deblur$Method = 'deblur'
    res_sum_deblur$Filtering = 'sum'
    res_sample_deblur$Filtering = 'sample'
    res_filtered_deblur$Filtering = 'filtering'
    
    all_res_endemics = rbind(all_res_endemics,
                             res_sum_dada2, res_sample_dada2, res_filtered_dada2, 
                             res_sum_unoise, res_sample_unoise, res_filtered_unoise,
                             res_sum_deblur, res_sample_deblur, res_filtered_deblur)}
  return(all_res_endemics)}

analyse_diversity <- function(asv_tab_dada2, asv_tab_unoise, asv_tab_deblur, asv_filtered_tab_dada2, asv_filtered_tab_unoise, asv_filtered_tab_deblur, meta_data, rd){
  all_res_diversity = data.frame()
  
  for (i in 1:5){
    set.seed(i)
    
    # Create summed table, and random replicate table
    tab_sum_dada2 = take_sum(asv_tab_dada2)
    tab_sample_dada2 = select_one_replicate(asv_tab_dada2)
    tab_filtering_dada2 = take_avg(asv_filtered_tab_dada2)
    
    tab_sum_unoise = take_sum(asv_tab_unoise)
    tab_sample_unoise = select_one_replicate(asv_tab_unoise)
    tab_filtering_unoise = take_avg(asv_filtered_tab_unoise)
    
    tab_sum_deblur = take_sum(asv_tab_deblur)
    tab_sample_deblur = select_one_replicate(asv_tab_deblur)
    tab_filtering_deblur = take_avg(asv_filtered_tab_deblur)
    
    # Rarefy the tables
    rarefied_tab_sum_dada2 = rarefy_tab(tab_sum_dada2, rd)
    rarefied_tab_sample_dada2 = rarefy_tab(tab_sample_dada2, rd)
    rarefied_tab_filtered_dada2 = rarefy_tab(tab_filtering_dada2, rd)
    
    rarefied_tab_sum_unoise = rarefy_tab(tab_sum_unoise, rd)
    rarefied_tab_sample_unoise = rarefy_tab(tab_sample_unoise, rd)
    rarefied_tab_filtered_unoise = rarefy_tab(tab_filtering_unoise, rd)
    
    rarefied_tab_sum_deblur = rarefy_tab(tab_sum_deblur, rd)
    rarefied_tab_sample_deblur = rarefy_tab(tab_sample_deblur, rd)
    rarefied_tab_filtered_deblur = rarefy_tab(tab_filtering_deblur, rd)
    
    # Diversity analysis
    sh_div_sum_dada2 = diversity(t(rarefied_tab_sum_dada2))
    sh_div_sample_dada2 = diversity(t(rarefied_tab_sample_dada2))
    sh_div_filtering_dada2 = diversity(t(rarefied_tab_filtered_dada2))
    
    sh_div_sum_unoise = diversity(t(rarefied_tab_sum_unoise))
    sh_div_sample_unoise = diversity(t(rarefied_tab_sample_unoise))
    sh_div_filtering_unoise = diversity(t(rarefied_tab_filtered_unoise))
    
    sh_div_sum_deblur = diversity(t(rarefied_tab_sum_deblur))
    sh_div_sample_deblur = diversity(t(rarefied_tab_sample_deblur))
    sh_div_filtering_deblur = diversity(t(rarefied_tab_filtered_deblur))
    
    ob_asv_sum_dada2 = colSums(rarefied_tab_sum_dada2 > 0)
    ob_asv_sample_dada2 = colSums(rarefied_tab_sample_dada2 > 0)
    ob_asv_filtering_dada2 = colSums(rarefied_tab_filtered_dada2 > 0)
    
    ob_asv_sum_unoise = colSums(rarefied_tab_sum_unoise > 0)
    ob_asv_sample_unoise = colSums(rarefied_tab_sample_unoise > 0)
    ob_asv_filtering_unoise = colSums(rarefied_tab_filtered_unoise > 0)
    
    ob_asv_sum_deblur = colSums(rarefied_tab_sum_deblur > 0)
    ob_asv_sample_deblur = colSums(rarefied_tab_sample_deblur > 0)
    ob_asv_filtering_deblur = colSums(rarefied_tab_filtered_deblur > 0)
    
    p_uniq_sum_dada2 = colSums(rarefied_tab_sum_dada2[(rowSums(rarefied_tab_sum_dada2 > 0) == 1), ] > 0) / colSums(rarefied_tab_sum_dada2 > 0)
    p_uniq_sample_dada2 = colSums(rarefied_tab_sample_dada2[(rowSums(rarefied_tab_sample_dada2 > 0) == 1), ] > 0) / colSums(rarefied_tab_sample_dada2 > 0)
    p_uniq_filtering_dada2 = colSums(rarefied_tab_filtered_dada2[(rowSums(rarefied_tab_filtered_dada2 > 0) == 1), ] > 0) / colSums(rarefied_tab_filtered_dada2 > 0)
    
    p_uniq_sum_unoise = colSums(rarefied_tab_sum_unoise[(rowSums(rarefied_tab_sum_unoise > 0) == 1), ] > 0) / colSums(rarefied_tab_sum_unoise > 0)
    p_uniq_sample_unoise = colSums(rarefied_tab_sample_unoise[(rowSums(rarefied_tab_sample_unoise > 0) == 1), ] > 0) / colSums(rarefied_tab_sample_unoise > 0)
    p_uniq_filtering_unoise = colSums(rarefied_tab_filtered_unoise[(rowSums(rarefied_tab_filtered_unoise > 0) == 1), ] > 0) / colSums(rarefied_tab_filtered_unoise > 0)
    
    p_uniq_sum_deblur = colSums(rarefied_tab_sum_deblur[(rowSums(rarefied_tab_sum_deblur > 0) == 1), ] > 0) / colSums(rarefied_tab_sum_deblur > 0)
    p_uniq_sample_deblur = colSums(rarefied_tab_sample_deblur[(rowSums(rarefied_tab_sample_deblur > 0) == 1), ] > 0) / colSums(rarefied_tab_sample_deblur > 0)
    p_uniq_filtering_deblur = colSums(rarefied_tab_filtered_deblur[(rowSums(rarefied_tab_filtered_deblur > 0) == 1), ] > 0) / colSums(rarefied_tab_filtered_deblur > 0)
 
    ab_uniq_sum_dada2 = colSums(rarefied_tab_sum_dada2[(rowSums(rarefied_tab_sum_dada2 > 0) == 1), ]) / colSums(rarefied_tab_sum_dada2)
    ab_uniq_sample_dada2 = colSums(rarefied_tab_sample_dada2[(rowSums(rarefied_tab_sample_dada2 > 0) == 1), ]) / colSums(rarefied_tab_sample_dada2)
    ab_uniq_filtering_dada2 = colSums(rarefied_tab_filtered_dada2[(rowSums(rarefied_tab_filtered_dada2 > 0) == 1), ]) / colSums(rarefied_tab_filtered_dada2)
    
    ab_uniq_sum_unoise = colSums(rarefied_tab_sum_unoise[(rowSums(rarefied_tab_sum_unoise > 0) == 1), ]) / colSums(rarefied_tab_sum_unoise)
    ab_uniq_sample_unoise = colSums(rarefied_tab_sample_unoise[(rowSums(rarefied_tab_sample_unoise > 0) == 1), ]) / colSums(rarefied_tab_sample_unoise)
    ab_uniq_filtering_unoise = colSums(rarefied_tab_filtered_unoise[(rowSums(rarefied_tab_filtered_unoise > 0) == 1), ]) / colSums(rarefied_tab_filtered_unoise)
    
    ab_uniq_sum_deblur = colSums(rarefied_tab_sum_deblur[(rowSums(rarefied_tab_sum_deblur > 0) == 1), ]) / colSums(rarefied_tab_sum_deblur)
    ab_uniq_sample_deblur = colSums(rarefied_tab_sample_deblur[(rowSums(rarefied_tab_sample_deblur > 0) == 1), ]) / colSums(rarefied_tab_sample_deblur)
    ab_uniq_filtering_deblur = colSums(rarefied_tab_filtered_deblur[(rowSums(rarefied_tab_filtered_deblur > 0) == 1), ]) / colSums(rarefied_tab_filtered_deblur)
    
    all_res_diversity = rbind(all_res_diversity,
                              # Shannon
                              data.frame(Method = rep('dada2', length(sh_div_sum_dada2)),
                                         Filtering = rep('sum', length(sh_div_sum_dada2)),
                                         Metric = rep('Shannon', length(sh_div_sum_dada2)),
                                         Value = sh_div_sum_dada2),
                              data.frame(Method = rep('dada2', length(sh_div_sample_dada2)),
                                         Filtering = rep('sample', length(sh_div_sample_dada2)),
                                         Metric = rep('Shannon', length(sh_div_sample_dada2)),
                                         Value = sh_div_sample_dada2),
                              data.frame(Method = rep('dada2', length(sh_div_filtering_dada2)),
                                         Filtering = rep('filtering', length(sh_div_filtering_dada2)),
                                         Metric = rep('Shannon', length(sh_div_filtering_dada2)),
                                         Value = sh_div_filtering_dada2),
                              
                              data.frame(Method = rep('unoise', length(sh_div_sum_unoise)),
                                         Filtering = rep('sum', length(sh_div_sum_unoise)),
                                         Metric = rep('Shannon', length(sh_div_sum_unoise)),
                                         Value = sh_div_sum_unoise),
                              data.frame(Method = rep('unoise', length(sh_div_sample_unoise)),
                                         Filtering = rep('sample', length(sh_div_sample_unoise)),
                                         Metric = rep('Shannon', length(sh_div_sample_unoise)),
                                         Value = sh_div_sample_unoise),
                              data.frame(Method = rep('unoise', length(sh_div_filtering_unoise)),
                                         Filtering = rep('filtering', length(sh_div_filtering_unoise)),
                                         Metric = rep('Shannon', length(sh_div_filtering_unoise)),
                                         Value = sh_div_filtering_unoise),
                              
                              data.frame(Method = rep('deblur', length(sh_div_sum_deblur)),
                                         Filtering = rep('sum', length(sh_div_sum_deblur)),
                                         Metric = rep('Shannon', length(sh_div_sum_deblur)),
                                         Value = sh_div_sum_deblur),
                              data.frame(Method = rep('deblur', length(sh_div_sample_deblur)),
                                         Filtering = rep('sample', length(sh_div_sample_deblur)),
                                         Metric = rep('Shannon', length(sh_div_sample_deblur)),
                                         Value = sh_div_sample_deblur),
                              data.frame(Method = rep('deblur', length(sh_div_filtering_deblur)),
                                         Filtering = rep('filtering', length(sh_div_filtering_deblur)),
                                         Metric = rep('Shannon', length(sh_div_filtering_deblur)),
                                         Value = sh_div_filtering_deblur),
                              
                              # Obs ASVs
                              data.frame(Method = rep('dada2', length(ob_asv_sum_dada2)),
                                         Filtering = rep('sum', length(ob_asv_sum_dada2)),
                                         Metric = rep('Observed', length(ob_asv_sum_dada2)),
                                         Value = ob_asv_sum_dada2),
                              data.frame(Method = rep('dada2', length(ob_asv_sample_dada2)),
                                         Filtering = rep('sample', length(ob_asv_sample_dada2)),
                                         Metric = rep('Observed', length(ob_asv_sample_dada2)),
                                         Value = ob_asv_sample_dada2),
                              data.frame(Method = rep('dada2', length(ob_asv_filtering_dada2)),
                                         Filtering = rep('filtering', length(ob_asv_filtering_dada2)),
                                         Metric = rep('Observed', length(ob_asv_filtering_dada2)),
                                         Value = ob_asv_filtering_dada2),
                              
                              data.frame(Method = rep('unoise', length(ob_asv_sum_unoise)),
                                         Filtering = rep('sum', length(ob_asv_sum_unoise)),
                                         Metric = rep('Observed', length(ob_asv_sum_unoise)),
                                         Value = ob_asv_sum_unoise),
                              data.frame(Method = rep('unoise', length(ob_asv_sample_unoise)),
                                         Filtering = rep('sample', length(ob_asv_sample_unoise)),
                                         Metric = rep('Observed', length(ob_asv_sample_unoise)),
                                         Value = ob_asv_sample_unoise),
                              data.frame(Method = rep('unoise', length(ob_asv_filtering_unoise)),
                                         Filtering = rep('filtering', length(ob_asv_filtering_unoise)),
                                         Metric = rep('Observed', length(ob_asv_filtering_unoise)),
                                         Value = ob_asv_filtering_unoise),
                              
                              data.frame(Method = rep('deblur', length(ob_asv_sum_deblur)),
                                         Filtering = rep('sum', length(ob_asv_sum_deblur)),
                                         Metric = rep('Observed', length(ob_asv_sum_deblur)),
                                         Value = ob_asv_sum_deblur),
                              data.frame(Method = rep('deblur', length(ob_asv_sample_deblur)),
                                         Filtering = rep('sample', length(ob_asv_sample_deblur)),
                                         Metric = rep('Observed', length(ob_asv_sample_deblur)),
                                         Value = ob_asv_sample_deblur),
                              data.frame(Method = rep('deblur', length(ob_asv_filtering_deblur)),
                                         Filtering = rep('filtering', length(ob_asv_filtering_deblur)),
                                         Metric = rep('Observed', length(ob_asv_filtering_deblur)),
                                         Value = ob_asv_filtering_deblur),

                              # P uniques
                              data.frame(Method = rep('dada2', length(p_uniq_sum_dada2)),
                                         Filtering = rep('sum', length(p_uniq_sum_dada2)),
                                         Metric = rep('Unique #', length(p_uniq_sum_dada2)),
                                         Value = p_uniq_sum_dada2),
                              data.frame(Method = rep('dada2', length(p_uniq_sample_dada2)),
                                         Filtering = rep('sample', length(p_uniq_sample_dada2)),
                                         Metric = rep('Unique #', length(p_uniq_sample_dada2)),
                                         Value = p_uniq_sample_dada2),
                              data.frame(Method = rep('dada2', length(p_uniq_filtering_dada2)),
                                         Filtering = rep('filtering', length(p_uniq_filtering_dada2)),
                                         Metric = rep('Unique #', length(p_uniq_filtering_dada2)),
                                         Value = p_uniq_filtering_dada2),
                              
                              data.frame(Method = rep('unoise', length(p_uniq_sum_unoise)),
                                         Filtering = rep('sum', length(p_uniq_sum_unoise)),
                                         Metric = rep('Unique #', length(p_uniq_sum_unoise)),
                                         Value = p_uniq_sum_unoise),
                              data.frame(Method = rep('unoise', length(p_uniq_sample_unoise)),
                                         Filtering = rep('sample', length(p_uniq_sample_unoise)),
                                         Metric = rep('Unique #', length(p_uniq_sample_unoise)),
                                         Value = p_uniq_sample_unoise),
                              data.frame(Method = rep('unoise', length(p_uniq_filtering_unoise)),
                                         Filtering = rep('filtering', length(p_uniq_filtering_unoise)),
                                         Metric = rep('Unique #', length(p_uniq_filtering_unoise)),
                                         Value = p_uniq_filtering_unoise),
                              
                              data.frame(Method = rep('deblur', length(p_uniq_sum_deblur)),
                                         Filtering = rep('sum', length(p_uniq_sum_deblur)),
                                         Metric = rep('Unique #', length(p_uniq_sum_deblur)),
                                         Value = p_uniq_sum_deblur),
                              data.frame(Method = rep('deblur', length(p_uniq_sample_deblur)),
                                         Filtering = rep('sample', length(p_uniq_sample_deblur)),
                                         Metric = rep('Unique #', length(p_uniq_sample_deblur)),
                                         Value = p_uniq_sample_deblur),
                              data.frame(Method = rep('deblur', length(p_uniq_filtering_deblur)),
                                         Filtering = rep('filtering', length(p_uniq_filtering_deblur)),
                                         Metric = rep('Unique #', length(p_uniq_filtering_deblur)),
                                         Value = p_uniq_filtering_deblur),
                              
                              # Ab uniques
                              data.frame(Method = rep('dada2', length(ab_uniq_sum_dada2)),
                                         Filtering = rep('sum', length(ab_uniq_sum_dada2)),
                                         Metric = rep('Unique rel. ab', length(ab_uniq_sum_dada2)),
                                         Value = ab_uniq_sum_dada2),
                              data.frame(Method = rep('dada2', length(ab_uniq_sample_dada2)),
                                         Filtering = rep('sample', length(ab_uniq_sample_dada2)),
                                         Metric = rep('Unique rel. ab', length(ab_uniq_sample_dada2)),
                                         Value = ab_uniq_sample_dada2),
                              data.frame(Method = rep('dada2', length(ab_uniq_filtering_dada2)),
                                         Filtering = rep('filtering', length(ab_uniq_filtering_dada2)),
                                         Metric = rep('Unique rel. ab', length(ab_uniq_filtering_dada2)),
                                         Value = ab_uniq_filtering_dada2),
                              
                              data.frame(Method = rep('unoise', length(ab_uniq_sum_unoise)),
                                         Filtering = rep('sum', length(ab_uniq_sum_unoise)),
                                         Metric = rep('Unique rel. ab', length(ab_uniq_sum_unoise)),
                                         Value = ab_uniq_sum_unoise),
                              data.frame(Method = rep('unoise', length(ab_uniq_sample_unoise)),
                                         Filtering = rep('sample', length(ab_uniq_sample_unoise)),
                                         Metric = rep('Unique rel. ab', length(ab_uniq_sample_unoise)),
                                         Value = ab_uniq_sample_unoise),
                              data.frame(Method = rep('unoise', length(ab_uniq_filtering_unoise)),
                                         Filtering = rep('filtering', length(ab_uniq_filtering_unoise)),
                                         Metric = rep('Unique rel. ab', length(ab_uniq_filtering_unoise)),
                                         Value = ab_uniq_filtering_unoise),
                              
                              data.frame(Method = rep('deblur', length(ab_uniq_sum_deblur)),
                                         Filtering = rep('sum', length(ab_uniq_sum_deblur)),
                                         Metric = rep('Unique rel. ab', length(ab_uniq_sum_deblur)),
                                         Value = ab_uniq_sum_deblur),
                              data.frame(Method = rep('deblur', length(ab_uniq_sample_deblur)),
                                         Filtering = rep('sample', length(ab_uniq_sample_deblur)),
                                         Metric = rep('Unique rel. ab', length(ab_uniq_sample_deblur)),
                                         Value = ab_uniq_sample_deblur),
                              data.frame(Method = rep('deblur', length(ab_uniq_filtering_deblur)),
                                         Filtering = rep('filtering', length(ab_uniq_filtering_deblur)),
                                         Metric = rep('Unique rel. ab', length(ab_uniq_filtering_deblur)),
                                         Value = ab_uniq_filtering_deblur))}
  return(all_res_diversity)}

main <- function(){
  asv_data = load_asv_data(25000)
  asv_tab_dada2 = asv_data$dada2
  asv_tab_unoise = asv_data$unoise
  asv_tab_deblur = asv_data$deblur
  
  meta_data = load_meta_data(asv_tab_unoise)
  
  # Custom filtering
  filtered_asvs_dada2 = filtering_br(asv_tab_dada2, meta_data)
  asv_filtered_tab_dada2 = asv_tab_dada2[rownames(asv_tab_dada2) %in% filtered_asvs_dada2, ]

  filtered_asvs_unoise = filtering_br(asv_tab_unoise, meta_data)
  asv_filtered_tab_unoise = asv_tab_unoise[rownames(asv_tab_unoise) %in% filtered_asvs_unoise, ]

  filtered_asvs_deblur = filtering_br(asv_tab_deblur, meta_data)
  asv_filtered_tab_deblur = asv_tab_deblur[rownames(asv_tab_deblur) %in% filtered_asvs_deblur, ]
  
  all_res_endemics = analyse_endemics(asv_tab_dada2, asv_tab_unoise, asv_tab_deblur,
                                      asv_filtered_tab_dada2, asv_filtered_tab_unoise, asv_filtered_tab_deblur, meta_data, 76580)

  all_res_diversity = analyse_diversity(asv_tab_dada2, asv_tab_unoise, asv_tab_deblur,
                                        asv_filtered_tab_dada2, asv_filtered_tab_unoise, asv_filtered_tab_deblur, meta_data, 76580)
  
  write.csv(all_res_endemics, 'method_comp_endemics.csv', quote = F, row.names = F)
  write.csv(all_res_diversity, 'method_comp_diversity.csv', quote = F, row.names = F)
  
  # plot 1 : endemics  
  p1 = all_res_endemics %>% group_by(Filtering, Method, mountain_range) %>% summarise(end_prop = mean(end_asv_n / tot_asv_n),
                                                                                 end_prop_sd = sd(end_asv_n / tot_asv_n)) %>%
      ggplot(aes(y=Filtering, x=end_prop, xmin=end_prop-end_prop_sd, xmax=end_prop+end_prop_sd, fill=Method)) + geom_errorbarh(position='dodge') + 
      geom_bar(stat='identity', position = "dodge") + facet_grid(.~mountain_range) + xlab('Proportion of endemic ASVs') + ylab('') + theme_minimal()
    
  p2 = all_res_endemics %>% group_by(Filtering, Method, mountain_range) %>% summarise(end_ab = mean(end_asv_ab),
                                                                                        end_ab_sd = sd(end_asv_ab)) %>%
      ggplot(aes(y=Filtering, x=end_ab, xmin=end_ab-end_ab_sd, xmax=end_ab+end_ab_sd, fill=Method)) + geom_errorbarh(position='dodge') + 
      geom_bar(stat='identity', position = "dodge") + facet_grid(.~mountain_range) + xlab('Relative abundance of endemic ASVs') + ylab('') + theme_minimal()

  p = ggarrange(p1, p2, ncol = 1, nrow = 2)
  ggsave(p, filename = 'comparison_methods.pdf', width = 18, height = 10)
  
  # plot 2,3 : diversity  
  p = all_res_diversity %>% group_by(Filtering, Method, Metric) %>% summarise(diversity = mean(Value), sd = sd(Value)) %>%
      ggplot(aes(y=Filtering, x=diversity, xmin=diversity-sd, xmax=diversity+sd, fill=Method)) + geom_errorbarh(position='dodge') + 
      geom_bar(stat='identity', position = "dodge") + ylab('') + theme_minimal() + facet_grid(.~Metric, scales = 'free')
    
  ggsave(p, filename = 'diversity_methods.pdf', width = 8, height = 5)
    
  all_res_diversity$sample = map_chr(rownames(all_res_diversity), function(x) paste0(c(strsplit(x, 'UP')[[1]][1], 'UP'), collapse = ''))
  all_res_diversity$Metric[all_res_diversity$Metric == 'Unique #'] = 'Unique prop.'
  all_res_diversity$mountain_range = map_chr(all_res_diversity$sample, function(x) ifelse(x %in% meta_data$sample, meta_data$mountain_range[meta_data$sample == x], 'NA'))
  p = all_res_diversity %>% filter(Filtering=='filtering', mountain_range != 'NA') %>% group_by(Filtering, Method, Metric, mountain_range) %>% summarise(diversity = mean(Value), sd = sd(Value)) %>%
      ggplot(aes(x=Filtering, y=diversity, ymin=diversity-sd, ymax=diversity+sd, fill=Method)) + geom_errorbar(position='dodge') + 
      geom_bar(stat='identity', position = "dodge") + ylab('Value') + xlab('') + theme_minimal() + facet_grid(Metric~mountain_range, scales = 'free')
    
  ggsave(p, filename = 'diversity_methods_ranges.pdf', width = 14, height = 6)}

main()


