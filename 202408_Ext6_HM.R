### Heatmap ###

library(circlize)
library(ComplexHeatmap)
library(reshape2)
library(phyloseq)
library(phyloseqCompanion)

## Load indicators
indicator_taxa <- read.csv("2024_0307_indicator_taxa.csv",header=T)

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
merged_NOMIS_DEBLUR <- merge_phyloseq(OTU_NOMIS, tax_NOMIS, metadata_NOMIS)

## Prune Uganda samples
prune_Uganda <- subset_samples(merged_NOMIS_DEBLUR, !site_c %in% "Uganda")
metadata_nomis_inext <- sample.data.frame(prune_Uganda)

## take the relative abundance of prune_uganda
ra_prune_uganda<- transform_sample_counts(prune_Uganda, function(x) x / sum(x))
ra_asv_prune_uganda <- otu_table(ra_prune_uganda, taxa_are_rows=T)
ra_asv_prune_uganda_filtered <- ra_asv_prune_uganda[rowSums(ra_asv_prune_uganda[])>0,]

## Relative abundance of the indicator taxa
## asv_indicator -> table d'abondance relative 
ra_indicator <- ra_asv_prune_uganda_filtered[row.names(ra_asv_prune_uganda_filtered) %in% (indicator_taxa$X)]
asv_ra_indicator <- otu_table(ra_indicator, taxa_are_rows=T)

## Merge everything 
merged_NOMIS_indicator_ab<- merge_phyloseq(asv_ra_indicator, tax_NOMIS, metadata_NOMIS)
metadata_indicator <- sample.data.frame(merged_NOMIS_indicator_ab)
#taxglom_biogeo_barplot <- tax_glom(physeq=NOMIS_FR, taxrank=rank_names(NOMIS_FR)[5], NArm=F)

merged_NOMIS_indicator_ra_asv <- (as.matrix(otu_table(merged_NOMIS_indicator_ab, taxa_are_rows=T)))

melt_asv_indicator <- melt(merged_NOMIS_indicator_ra_asv)
merge_asv_indicator <- merge(as.data.frame(melt_asv_indicator),as.matrix(metadata_NOMIS), by.x="X2",by.y="sample")

## Here we would need to divide by the number of glaciers per mountain ranges. 
## this is what is done: 1) group by site, 2) calculate the sum of the RA within each site and rename it to summmrr
## 3) calculate the number of distinct values of GL code within each region and then rename it to n
## 4) calculate the average summrr divided by the nb of distinct values of GL code and keep site-c
## this gives us the average relative abundance of indicator taxa within each region! greatest in NZ and lowest in the Alps!
## the nb of GL differ between sites, that's why we have to divide by the correct nb of sites to compare 

sum_mr_indicator <- merge_asv_indicator %>% group_by(site_c)%>% summarize(summrr=sum(value), n=n_distinct(X2))%>% summarize(ar_mr=summrr/n, site_c)

median_indicator_ab <- sum_mr_indicator %>% 
  summarise(median=median(ar_mr), x = quantile(ar_mr, c(0.25, 0.5, 0.75)))

## Lets get the most abundant ASVs - les premiers 500
TopASV_f <- names(sort(taxa_sums(merged_NOMIS_indicator_ab), TRUE)[1:500])
top500_NOMIS_f <- prune_species(TopASV_f, merged_NOMIS_indicator_ab)
top500_NOMIS_f <- prune_taxa(taxa_sums(top500_NOMIS_f)>0, top500_NOMIS_f)
top_ASV<-as.data.frame(tax_table(top500_NOMIS_f))

merge_indicator_ab_sub <- merged_NOMIS_indicator_ra_asv [rownames(merged_NOMIS_indicator_ra_asv ) %in% rownames(top_ASV),]
hm_mat = (scale(log1p(as.matrix(merge_indicator_ab_sub))))
col_fun = colorRamp2(c(0, 1, 2), c("steelblue", "white", "darkred"))

#write.csv(hm_mat, "202408_hm_mat.csv")

ht <- ComplexHeatmap::Heatmap((hm_mat), use_raster = T, col = col_fun, column_split = metadata_indicator$site_c,
                              show_column_names = F, show_row_names = F, name = 'Normalised abundance',
                              row_dend_reorder = TRUE, show_row_dend = F, raster_quality = 5)

