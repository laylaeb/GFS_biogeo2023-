### Distance Decay ###

library(speedyseq)
library(phyloseq)
library(phyloseqCompanion)
library(geosphere)
library(vegan)
library(rbiom)
library(scales)
library(ggplot2)
library(dplyr)
library(fishualize)
library(ggpubr)
library(reshape2)
library(performance)
library(fANCOVA)

## Load ASV and taxonomy tables
asv<-read.csv(file="20240221_NOMIS_rarefied_deblur_table.csv.gz",sep=",",header=TRUE,row.names=1)
tax<-read.csv(file="20240221_NOMIS_rarefied_deblur_taxonomy.csv.gz",sep=",",header=TRUE,row.names=1)

## Load metadata file ## kept lat and lon_sp UP values
metadata_NOMIS="202402_NOMIS_metadata_GFS.tsv"
metadata_NOMIS<-import_qiime_sample_data(metadata_NOMIS)

## Create ASV table 
OTU_NOMIS <- otu_table(asv, taxa_are_rows=TRUE)
tax_NOMIS <- tax_table(as.matrix(tax))
tree_NOMIS <- read_tree("20240221_NOMIS_rarefied_deblur.tree")

merged_NOMIS_DEBLUR_DD <- merge_phyloseq(OTU_NOMIS, tax_NOMIS, metadata_NOMIS, tree_NOMIS)

## Prune Uganda
uganda=c("Uganda")
prune_Uganda <- subset_samples(merged_NOMIS_DEBLUR_DD,!site_c %in% uganda)
nomis_metadata <- sample.data.frame(prune_Uganda)

vegan_otu <- function(physeq){
  OTU <- otu_table(physeq)
  if(taxa_are_rows(OTU)){
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}

## Distance-Decay patterns on Bray-Curtis
## Full distance decay with all samples
## building dataframe encompassing latitude and longitude from all glaciers
NOMIS_df <- data.frame(nomis_metadata$lon_sp, nomis_metadata$lat_sp)
NOMIS_df$nomis_metadata.lon_sp <- as.numeric(NOMIS_df$nomis_metadata.lon_sp)
NOMIS_df$nomis_metadata.lat_sp <- as.numeric(NOMIS_df$nomis_metadata.lat_sp)

dist_geo_all <- distm(NOMIS_df, NOMIS_df, fun=distGeo)
dist_geo_all <- as.matrix(dist_geo_all)
diag(dist_geo_all)=NA
dist_geo_all_diag <- t(matrix(t(dist_geo_all)[which(!is.na(dist_geo_all))],nrow=150,ncol=151))

min(dist_geo_all_diag)
#write.csv(dist_geo_all_diag, "dist_geo_nomis.csv")

vegan_matrix_all <-vegan_otu(prune_Uganda)
allregion_bray <-vegdist(log1p(vegan_matrix_all), method="bray")
allregion_m <- as.matrix(allregion_bray)
diag(allregion_m)=NA
allregion_diag <-t(matrix(t(allregion_m)[which(!is.na(allregion_m))],nrow=150,ncol=151))
min(allregion_diag)

mantel(allregion_diag, dist_geo_all_diag, method = "pearson", permutations = 999, na.rm = TRUE)
##2024
# Mantel statistic based on Pearson's product-moment correlation 
# 
# Call:
# mantel(xdis = allregion_diag, ydis = dist_geo_all_diag, method = "pearson",      permutations = 999, na.rm = TRUE) 
# 
# Mantel statistic r: 0.5809 
#       Significance: 0.001 
# 
# Upper quantiles of permutations (null model):
#    90%    95%  97.5%    99% 
# 0.0382 0.0501 0.0612 0.0706 
# Permutation: free
# Number of permutations: 999

dist_all_bray <- data.frame(BC_dist_bc=as.vector(allregion_diag), BC_sim_bc=as.vector(1-allregion_diag), geo_dist=as.vector(dist_geo_all_diag), Method="All glaciers")

ggplot(dist_all_bray, aes(x=geo_dist/1000, y=BC_sim_bc)) +
  geom_point(size=1.2) + ylim(0,1) + xlim (0,20000)+
 # facet_wrap(~Method) +
  stat_smooth(method="lm", formula=y ~ (x), size=1.2, se=FALSE, linetype="solid") +
  ylab("Community Similarity - Bray-Curtis") +
  xlab("Geographic distance (km)") + ggtitle("") + theme_bw() +
  theme(strip.text.x=element_text(size=10, color="black", face="bold.italic")) +
  theme(strip.background=element_rect(colour="black", fill="white")) +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=14),
        legend.text=element_text(size=14)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_x_log10()

## Getting the slope coefficients
  dist_all_bray %>% 
  do({
    mod = lm(BC_sim_bc ~ geo_dist, data = dist_all_bray)
    data.frame(Intercept = coef(mod)[1],
               Slope = coef(mod)[2])
  })
  
## 2024
  # Intercept         Slope
  # (Intercept) 0.2655723 -1.173702e-08

## Distance Decay Patterns on Weighted Unifrac
## Full distance decay with all samples
## building dataframe encompassing latitude and longitude from all glaciers
asv_table_unif<-otu_table(prune_Uganda, taxa_are_rows=T)
unifrac_essai<-phyloseq::UniFrac(prune_Uganda, weighted=T, normalized=TRUE, parallel=FALSE, fast=TRUE)
allregion_unif<-as.matrix(unifrac_essai)

diag(allregion_unif)=NA
allregion_diag_unif<-t(matrix(t(allregion_unif)[which(!is.na(allregion_unif))],nrow=150,ncol=151))

mantel(allregion_diag_unif, dist_geo_all_diag, method = "pearson", permutations = 999, na.rm = TRUE)

## 2024

# Mantel statistic based on Pearson's product-moment correlation 
# 
# Call:
# mantel(xdis = allregion_diag_unif, ydis = dist_geo_all_diag,      method = "pearson", permutations = 9999, na.rm = TRUE) 
# 
# Mantel statistic r: 0.1103 
#       Significance: 8e-04 
# 
# Upper quantiles of permutations (null model):
#    90%    95%  97.5%    99% 
# 0.0402 0.0517 0.0627 0.0753 
# Permutation: free
# Number of permutations: 9999

dist_all_wunif <- data.frame(dis_wunif=as.vector(allregion_diag_unif), sim_wunif=as.vector(1-allregion_diag_unif), geo_dist=as.vector(dist_geo_all_diag), Method="All glaciers")

ggplot(dist_all_wunif, aes(x=geo_dist/1000, y=sim_wunif)) +
  geom_point(size=1.2) + ylim(0,1) + xlim (0,20000)+
 # facet_wrap(~Method) +
  stat_smooth(method="lm", formula= y ~ (x), size=1.2, se=FALSE, linetype="solid") +
  ylab("Community Similarity - Weighted Unifrac") +
  xlab("Geographic distance (km)") + ggtitle("") + theme_bw() +
  theme(strip.text.x=element_text(size=10, color="black", face="bold.italic")) +
  theme(strip.background=element_rect(colour="black", fill="white")) +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=14),
        legend.text=element_text(size=14)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_x_log10()

dist_all_wunif %>% 
  do({
    mod = lm(sim_wunif ~ geo_dist, data = .)
    data.frame(Intercept = coef(mod)[1],
               Slope = coef(mod)[2])
  })
# Intercept         Slope
# (Intercept) 0.8525405 -1.394385e-09

## Full distance decay using Unweighted unifrac distance
asv_table_unif <-otu_table(prune_Uganda, taxa_are_rows=T)
uwunifrac_essai<-phyloseq::UniFrac(prune_Uganda, weighted=F, normalized=TRUE, parallel=FALSE, fast=TRUE)

allregion_uwnif<-as.matrix(uwunifrac_essai)

diag(allregion_uwnif)=NA
allregion_diag_uwnif<-t(matrix(t(allregion_uwnif)[which(!is.na(allregion_uwnif))],nrow=150,ncol=151))

mantel(allregion_diag_uwnif, dist_geo_all_diag, method = "pearson", permutations = 999, na.rm = TRUE)

# Mantel statistic based on Pearson's product-moment correlation 
# 
# Call:
# mantel(xdis = allregion_diag_uwnif, ydis = dist_geo_all_diag,      method = "pearson", permutations = 999, na.rm = TRUE) 
# 
# Mantel statistic r: 0.3601 
#       Significance: 1e-04 
# 
# Upper quantiles of permutations (null model):
#    90%    95%  97.5%    99% 
# 0.0439 0.0584 0.0700 0.0845 
# Permutation: free
# Number of permutations: 9999

dist_all_uw <-data.frame(dis_uw=as.vector(allregion_diag_uwnif), sim_uw=as.vector(1-allregion_diag_uwnif), geo_dist=as.vector(dist_geo_all_diag), Method="All glaciers")

ggplot(dist_all_uw, aes(x=geo_dist/1000, y=sim_uw)) +
  geom_point(size=1.2) + ylim(0,1) + xlim (0,20000)+
  # facet_wrap(~Method) +
  stat_smooth(method="lm", formula= y ~ (x), size=1.2, se=FALSE, linetype="solid") +
  ylab("Community Similarity - UnWeighted Unifrac") +
  xlab("Geographic distance (km)") + ggtitle("") + theme_bw() +
  theme(strip.text.x=element_text(size=10, color="black", face="bold.italic")) +
  theme(strip.background=element_rect(colour="black", fill="white")) +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=14),
        legend.text=element_text(size=14)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_x_log10()


dist_all_uw %>% 
  do({
    mod = lm(sim_uw ~ geo_dist, data = .)
    data.frame(Intercept = coef(mod)[1],
               Slope = coef(mod)[2])
  })

# Intercept         Slope
# (Intercept) 0.5028242 -5.418735e-09

## Full distance decay using SORENSEN
vegan_matrix_all <-vegan_otu(prune_Uganda)
allregion_sor <-vegdist(vegan_matrix_all, binary=T)
allregion_sor_m <-as.matrix(allregion_sor)
diag(allregion_sor_m)=NA
allregion_sor_diag <-t(matrix(t(allregion_sor_m)[which(!is.na(allregion_sor_m))],nrow=150,ncol=151))
min(allregion_sor_diag)

mantel(allregion_sor_diag, dist_geo_all_diag, method = "pearson", permutations = 999, na.rm = TRUE)

# Mantel statistic based on Pearson's product-moment correlation 
# 
# Call:
# mantel(xdis = allregion_sor_diag, ydis = dist_geo_all_diag, method = "pearson",      permutations = 999, na.rm = TRUE) 
# 
# Mantel statistic r: 0.6008 
#       Significance: 0.001 
# 
# Upper quantiles of permutations (null model):
#    90%    95%  97.5%    99% 
# 0.0376 0.0517 0.0595 0.0727 
# Permutation: free
# Number of permutations: 999


dist_all_sor <-data.frame(dis_sor=as.vector(allregion_sor_diag), sim_sor=as.vector(1-allregion_sor_diag), geo_dist=as.vector(dist_geo_all_diag), Method="All glaciers")
ggplot(dist_all_sor, aes(x=geo_dist/1000, y=sim_sor)) +
  geom_point(size=1.2) + ylim(0,1) + xlim (0,20000)+
  stat_smooth(method="lm", formula=y ~ (x), size=1.2, se=FALSE, linetype="solid") +
  ylab("Community Similarity - Sorensen") +
  xlab("Geographic distance (km)") + ggtitle("") + theme_bw() +
  theme(strip.text.x=element_text(size=10, color="black", face="bold.italic")) +
  theme(strip.background=element_rect(colour="black", fill="white")) +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=14),
        legend.text=element_text(size=14)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_x_log10()

dist_all_sor %>% 
  do({
    mod = lm(sim_sor ~ geo_dist, data = .)
    data.frame(Intercept = coef(mod)[1],
               Slope = coef(mod)[2])
  })

# Intercept         Slope
# (Intercept) 0.2995777 -1.308125e-08

## Combine different decay patterns together
melt_bc <- melt(allregion_diag)
melt_bc$dis <- "BC"
melt_sor <- melt(allregion_sor_diag)
melt_sor$dis <- "SOR"

melt_geo <- melt(dist_geo_all_diag)
colnames(melt_geo) <- c("Var1","Var2","dist_geo")
bind_bcsor <- rbind(melt_bc, melt_sor)
merge_dis_geo <- merge(bind_bcsor, melt_geo)

## Plot everything!
bcsor_plot<- ggplot(merge_dis_geo) + 
  geom_point(aes(y=1-value, x = dist_geo, color= dis)) + 
  stat_smooth(aes(y=1-value, x = dist_geo, color=dis), method="lm", formula=y ~ (x), size=1.2, se=FALSE, linetype="solid") +
  scale_color_fish(option="Scarus_quoyi",discrete=T)+
  ylab("Community Similarity") +
  xlab("Geographic distance (km)") + ggtitle("") + theme_bw() +
  scale_x_log10() +theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
bcsor_plot

## Lets combine unweighted and weighted unifrac
melt_unifrac <- melt(allregion_diag_unif)
melt_unifrac$dis <- "wunifrac"
melt_uwunifrac <- melt(allregion_diag_uwnif)
melt_uwunifrac$dis <- "unwunifrac"
bind_unifrac<- rbind(melt_unifrac, melt_uwunifrac)
merge_dis_unif <- merge(bind_unifrac, melt_geo)

unifrac_plot<- ggplot(merge_dis_unif) + 
  geom_point(aes(y=1-value, x = dist_geo, color= dis)) + 
  stat_smooth(aes(y=1-value, x = dist_geo, color=dis), method="lm", formula=y ~ (x), size=1.2, se=FALSE, linetype="solid") +
  scale_color_fish(option="Trimma_lantana",discrete=T)+
  ylab("Community Similarity") +
  xlab("Geographic distance (km)") + ggtitle("") + theme_bw() +
  scale_x_log10() +theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
unifrac_plot

## To include in the MS, we will have to subset the matrices because the figure is too heavy
## Create a plot with all the distance matrices
bind_all_dist <- rbind(melt_bc, melt_sor, melt_uwunifrac, melt_unifrac)
merge_alldist_geo <- merge(bind_all_dist, melt_geo)

#write.csv(merge_alldist_geo,"202408_merge_alldist_geo.csv")

# For the figure included within the manuscript, we include a subset
get_fishcol <- fish(4, option="Trimma_lantana")

alldist_plot<- ggplot(merge_alldist_geo ) + 
  geom_point(aes(y=1-value, x = dist_geo/1000, color= dis), alpha=0.1) + 
  stat_smooth(aes(y=1-value, x = dist_geo/1000, color=dis), method="lm", formula=y ~ (x), size=1.2, se=FALSE, linetype="solid") +
  scale_color_fish(option="Scarus_quoyi",discrete=T,direction=-1)+
  ylab("Community Similarity") +
  xlab("Geographic distance (km)") + ggtitle("") + theme_bw() +
  scale_x_log10() +theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

## Ensure GeographicDistance is numeric
merge_alldist_geo$dist_geo <- as.numeric(merge_alldist_geo$dist_geo)

## Convert GeographicDistance to kilometers
merge_alldist_geo$dist_geo_km <- merge_alldist_geo$dist_geo / 1000

## Create a function to calculate R-squared
calculate_r_squared <- function(model) {
  return(format(summary(model)$r.squared, digits = 3))
}

## Fit linear regression models using stat_smooth
lm_models_tax <- by(merge_alldist_geo, merge_alldist_geo$dis, function(sub_df_tax) {
  lm(1 - sub_df_tax$value ~ sub_df_tax$dist_geo_km, data = sub_df_tax)
})

## Extract coefficients from the linear regression models and calculate R-squared
intercepts <- sapply(lm_models_tax, function(model) coef(model)[1])
slopes <- sapply(lm_models_tax, function(model) coef(model)[2])
r_squared <- sapply(lm_models_tax, calculate_r_squared)

subset_size <- 35000
set.seed(42)  # You can change the seed for randomness
subset_merge_alldist <- merge_alldist_geo[sample(nrow(merge_alldist_geo), subset_size), , drop = FALSE]

# ggsave(alldist_plot, path='dissimilarity.pdf', device='pdf', dpi=100)

### DDPs comparisons
## Bray-Curtis vs Sorensen
## Lets try another way
library(tidyverse)
library(ggpubr)
library(rstatix)
library(broom)
  
## Homogeneity of regression slopes  
  
merge_bcsor <- merge(dist_all_bray, dist_all_sor, by="row.names")
anco_bcsor <- merge_bcsor[c("BC_sim_bc","sim_sor","geo_dist_x")]
manco <- melt(anco_bcsor,  id=c("geo_dist_x")) 

ggscatter(
  manco, x = "geo_dist_x", y = "value",
  color = "variable", add = "reg.line"
)+
  stat_regline_equation(
    aes(label =  paste(..eq.label.., ..rr.label.., sep = "~~~~"), color = variable)
  )

manco$geo_dist_log <- log((manco$geo_dist_x)/1000)

## Weighted vs Unweighted unifrac
merge_unifrac <- merge(dist_all_wunif, dist_all_uw, by="row.names")
anco_unifrac <- merge_unifrac[c("sim_wunif","sim_uw","geo_dist_x")]
manco_unifrac <- melt(anco_unifrac, id=c("geo_dist_x")) 
manco_unifrac$geo_dist_log <- log((manco_unifrac$geo_dist_x)/1000)

ggscatter(
  manco_unifrac, x = "geo_dist_x", y = "value",
  color = "variable", add = "reg.line"
)+
  stat_regline_equation(
    aes(label =  paste(..eq.label.., ..rr.label.., sep = "~~~~"), color = variable)
  )
manco_unifrac$geo_dist_log <- log((manco_unifrac$geo_dist_x)/1000)
manco %>% anova_test(sqrt(value) ~ variable*geo_dist_log)  

# Normality of residuals
# Fit the model, the covariate goes first
model_bcsor <- lm((value) ~ geo_dist_log + variable, data = manco)

# Inspect the model diagnostic metrics
model_metrics_bcsor <- augment(model_bcsor) %>%
  select(-.hat, -.sigma, -.fitted) # Remove details
head(model_metrics_bcsor, 3)

## Anderson-Darling normality test
library(nortest)
ad.test(model_bcsor$residuals)

# Anderson-Darling normality test
# 
# data:  model$residuals
# A = 555.58, p-value < 2.2e-16

## Homogeneity of regression slopes
manco_unifrac %>% anova_test(log10(value) ~variable*geo_dist_x)

# Violation of the assumption
# Center the geo_dist_log variable
manco$geo_dist_log_centered <- scale(manco$geo_dist_log, center = TRUE, scale = FALSE)

# Fit the model with the centered variable
glm_dd_bcsor <- lm(sqrt(value) ~ geo_dist_log_centered * variable, data = manco)
check_model(glm_dd_bcsor)
summary(glm_dd_bcsor)

manco_unifrac$geo_dist_log <- log((manco_unifrac$geo_dist_x)/1000)
manco_unifrac$geo_dist_log_centered <- scale(manco_unifrac$geo_dist_log, center = TRUE, scale = FALSE)

glm_dd_wufrac <- lm(manco_unifrac, formula = value ~ geo_dist_log_centered*variable)
check_model(glm_dd_wufrac)
summary(glm_dd_wufrac)
