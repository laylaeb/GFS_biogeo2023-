### Computing within and across mountain range similarity ###
library(vegan)
library(adespatial)
library(fishualize)
library(phyloseq)
library(RColorBrewer)
library(pheatmap)
library(tidyverse)

## Read in filtered, merged and rarefied deblur ASV table
merged_NOMIS_DEBLUR <- readRDS("20240221_NOMIS_rarefied.RData")

## remove Uganda because only 1 GFS
uganda=c("Uganda")
prune_Uganda <- subset_samples(merged_NOMIS_DEBLUR , !site_c %in% uganda)
prune_Uganda_df <- as.matrix(otu_table(prune_Uganda, taxa_are_rows=T))
merge_glaciers_data <- as.data.frame(sample_data(prune_Uganda))

vegan_otu <- function(physeq){
  OTU <- otu_table(physeq)
  if(taxa_are_rows(OTU)){
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}

## Create Bray-Curtis and Sorensen dissimilarity matrices
vegan_matrix_nomis<-vegan_otu(prune_Uganda)
vegan_matrix_nomis<-(vegan_matrix_nomis)
asv_region_bray<-vegdist(log1p(vegan_matrix_nomis), method="bray") 

asv_region_sorensen<-vegdist(vegan_matrix_nomis, binary=T) ## Sorensen
asv_region_jaccard<-vegdist((vegan_matrix_nomis), method="jaccard") ## Jaccard

## Mean Bray-Curtis Dissimilarity values 
nomis_meandist_bray <- vegan::meandist(asv_region_bray, merge_glaciers_data$site_c)
hist(nomis_meandist_bray)
summary(nomis_meandist_bray)

df_bray <- as.data.frame(nomis_meandist_bray)%>%
  rownames_to_column("key1")%>%
  gather(key2, value, 2:ncol(.))

# within
df_bray%>%
  filter(key1 == key2)%>%
  summarise(mean = mean(value),
            media = median(value),
            IQR1 = quantile(value, c(0.25)),
            IQR2 = quantile(value, c(0.75))
              )

between
df_bray%>%
  filter(key1 != key2)%>%
  distinct(value)%>%
  summarise(mean = mean(value),
            media = median(value),
            IQR1 = quantile(value, c(0.25)),
            IQR2 = quantile(value, c(0.75))
              )

## Mean Sorensen Dissimilarity values 
nomis_meandist_sorensen <- vegan::meandist(asv_region_sorensen, merge_glaciers_data$site_c)
hist(nomis_meandist_sorensen)
summary(nomis_meandist_sorensen)

df_sorensen <- as.data.frame(nomis_meandist_sorensen)%>%
  rownames_to_column("key1")%>%
  gather(key2, value, 2:ncol(.))

# within
df_sorensen%>%
  filter(key1 == key2)%>%
  summarise(mean = mean(value),
            media = median(value),
            IQR1 = quantile(value, c(0.25)),
            IQR2 = quantile(value, c(0.75))
  )

# between
df_sorensen%>%
  filter(key1 != key2)%>%
  distinct(value)%>%
  summarise(mean = mean(value),
            media = median(value),
            IQR1 = quantile(value, c(0.25)),
            IQR2 = quantile(value, c(0.75))
  )


## Plotting Dissimilarity
nomis_braycurtis <- as.matrix(nomis_meandist_bray)
nomis_sorensen <- as.matrix(nomis_meandist_sorensen)

## Heatmap for the paper
fish_color<- fish(n=5, option = "Zebrasoma_velifer", end=1,  begin=0.3)
fish_color_bray<- fish(n=5, option = "Coris_gaimard",direction=1,begin=0.2)

habitat_labeller <- c("Alaska" = "Alaska\nRange", Alps = "European\nAlps", Caucasus = "Caucasus\nMountains",
                      Chile = "Chilean\nAndes", Ecuador = "Ecuadoran\nAndes", Greenland = "Southwest\nGreenland",
                      Kirghizistan = "Pamir &\nTien Shan", Nepal = "Himalayas", New_Zealand = "Southern\nAlps", 
                      Norway = "Scandinavian\nMountains")

## Bray-Curtis 
pdf(file = "images/20230305_Bray_Curtis_heatmap.pdf", width = 6.5, height = 6)
heatmap_bray_curtis <- pheatmap(nomis_braycurtis,
         display_numbers = T,
         color =fish_color_bray, 
         fontsize_number = 8,
         treeheight_row=0,
         treeheight_col=0,cluster_rows=F, cluster_cols=F, na_col="white",
         labels_row = habitat_labeller, labels_col = habitat_labeller  )
heatmap_bray_curtis
dev.off()
saveRDS(heatmap_bray_curtis, "images/20230305_Bray_Curtis_heatmap.rds")

## Sorensen
pdf(file = "images/20230305_Sorensen_heatmap.pdf", width = 6.5, height = 6)
heatmap_sorensen <- pheatmap(nomis_sorensen,
                             display_numbers = T,
                             color =fish_color, 
                             fontsize_number = 8,
                             treeheight_row=0,
                             treeheight_col=0, cluster_rows=F, cluster_cols=F, na_col="white",
                             labels_row = habitat_labeller, labels_col = habitat_labeller  )
heatmap_sorensen
dev.off()
saveRDS(heatmap_sorensen, "images/20230305_Sorensen_heatmap.rds")

## Relationships between Bray and Sorensen
asv_region_bray_mat <- as.matrix(asv_region_bray)
asv_region_sor_mat <- as.matrix(asv_region_sorensen)
dissimilarity <- data.frame(b=as.vector(asv_region_bray_mat[upper.tri(asv_region_bray_mat)]),
                  s=as.vector(asv_region_sor_mat[upper.tri(asv_region_sor_mat)]))

ggplot(dissimilarity, aes(x=s, y=b))+
  geom_point()+
  geom_smooth(method="lm")+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

fit_dissimilarity=lm(b~s,data=dissimilarity)
summary(fit_dissimilarity)





