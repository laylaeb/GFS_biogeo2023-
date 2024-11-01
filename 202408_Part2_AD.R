### Alpha diversity ###
library(dplyr)
library(tidyverse)
library(breakaway)
library(data.table)
library(phyloseq)
library(microViz)
library(viridis)
library(hrbrthemes)
library(paletteer)
library(microbiome)
library(mgcv)

## Remove Uganda because only 1 GFS
uganda=c("Uganda")
prune_Uganda <- subset_samples(merged_NOMIS_DEBLUR,!site_c %in% uganda)
prune_Uganda_df <- as.matrix(otu_table(prune_Uganda, taxa_are_rows=T))
prune_Uganda_df <- prune_Uganda_df[rowSums(prune_Uganda_df[])>0,]

## Calculate diversity metrics
diversity_nomis <- phyloseq::sample_data(estimate_richness(prune_Uganda,measures=c("Observed","Shannon")))
alphadiv_nomis <- merge_phyloseq(prune_Uganda, diversity_nomis)
alphadt_nomis_df<- sample_data(alphadiv_nomis)

## Compute the average observed richness per mountain range
OR <- alphadt_nomis_df[, c("sample", "gl_name", "site_c", "Observed")]
#hist(OR$Observed) ##normally distributed
ASVrichness_GFS <- OR %>% 
  group_by(site_c) %>% 
  summarise(average=mean(Observed), std=sd(Observed))
ASVrichness_GFS

average_richness <- OR %>% 
  summarise(average=mean(Observed), std=sd(Observed))
average_richness

median_richness <- OR %>% 
  summarise(median=median(Observed), x = quantile(Observed, c(0.25, 0.5, 0.75)))
median_richness

## Plot ASV richness per region using violin plot and include jitter  -->for figure 1
habitat_labeller <- c("Alaska" = "Alaska\nRange", Alps = "European\nAlps", Caucasus = "Caucasus\nMountains",
                      Chile = "Chilean\nAndes", Ecuador = "Ecuadoran\nAndes", Greenland = "Southwest\nGreenland",
                      Kirghizistan = "Pamir &\nTien Shan", Nepal = "Himalayas", New_Zealand = "Southern\nAlps", 
                      Norway = "Scandinavian\nMountains")

plot_OR <- ggplot(OR,aes(x=site_c,y=Observed, color=site_c)) + 
  geom_violin(width=1.4,alpha=0.5) +
  geom_boxplot(width=0.1, color="black", alpha=1, outlier.shape=NA) +
  geom_jitter(position=position_jitter(0.2), alpha=0.9) +  # Set alpha value for transparency
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.position = "none") +
  labs(y = "ASV Richness", x = "") +
  scale_x_discrete(labels = habitat_labeller)+
  scale_colour_manual(values=c("#2E2A2BFF","#CF4E9CFF","#8C57A2FF",
                               "#3EBCB6","#82581FFF","#2F509EFF",
                               "#E5614CFF","#97A1A7FF","#bee183","#DC9445FF","#EDD03E"))
plot_OR
# ggsave("images/20240226_ASV_richness_violins.pdf", plot_OR, height = 5, width = 8)
# saveRDS(plot_OR, file = "images/20240226_ASV_richness_violins.rds")

## Shannon diversity
Shannon <- alphadt_nomis_df[, c("sample", "gl_name", "site_c", "Shannon")]
Shannon_GFS <- Shannon %>% 
  group_by(site_c) %>% 
  summarise(average=mean(Shannon), std=sd(Shannon))
Shannon_GFS

median_shannon <- Shannon %>% 
  reframe(median=median(Shannon), x = quantile(Shannon, c(0.25, 0.5, 0.75)))
median_shannon

## Plot Shannon diversity per region using violin plot and include jitter 
plot_Shannon <- ggplot(Shannon,aes(x=site_c,y=Shannon, color=site_c)) + 
  geom_violin(width=1.4,alpha=0.5) +
  geom_boxplot(width=0.1, color="black", alpha=1, outlier.shape=NA) +
  geom_jitter(position=position_jitter(0.2), alpha=0.9) +  # Set alpha value fShannon transparency
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90,vjust = 0.5, hjust = 1),
        legend.position = "none") +
  scale_x_discrete(labels = habitat_labeller)+
  labs(x = "")+
  scale_colour_manual(values=c("#2E2A2BFF","#CF4E9CFF","#8C57A2FF",
                               "#3EBCB6","#82581FFF","#2F509EFF",
                               "#E5614CFF","#97A1A7FF","#bee183","#DC9445FF","#EDD03E"))
plot_Shannon
# ggsave("images/20240226_Shannon.pdf", plot_Shannon, height = 5, width = 8)
# saveRDS(plot_Shannon, file = "images/20240226_Shannon_violins.rds")

## Hill number - q=1
hill_q1_nomis<-as.data.frame(hill_taxa(t(prune_Uganda_df), q = 1))
hill_q1_indices_nomis <- phyloseq::sample_data(hill_q1_nomis)

prune_Uganda_hill_nomis <- phyloseq::merge_phyloseq(prune_Uganda,hill_q1_indices_nomis)
meta_diversity_nomis_hill <- sample.data.frame(prune_Uganda_hill_nomis)

median_shannon_hill<- meta_diversity_nomis_hill %>% 
  summarise(median=median(hill_taxa.t.prune_Uganda_df...q...1. ), x = quantile(hill_taxa.t.prune_Uganda_df...q...1. , c(0.25, 0.5, 0.75)))
# median        x
# 1 327.0075 200.8401
# 2 327.0075 327.0075
# 3 327.0075 464.0519

#Evenness with Bulla index
bulla_estimate <- phyloseq::sample_data(microbiome::evenness(prune_Uganda, index="all"))
alphadiv_nomis <- merge_phyloseq(alphadiv_nomis, bulla_estimate)
alphadt_nomis_df<- sample_data(alphadiv_nomis)
#hist(alphadt_nomis_df$bulla)

evenness_GFS <- alphadt_nomis_df %>% 
  group_by(site_c) %>% 
  summarise(average=mean(bulla), std=sd(bulla))
evenness_GFS

evenness_GFS_median <- alphadt_nomis_df %>% 
  reframe(median=median(bulla), x = quantile(bulla, c(0.25, 0.5, 0.75)))
evenness_GFS_median

evenness_GFS_total <- alphadt_nomis_df %>% 
  summarise(average=mean(bulla), std=sd(bulla))

## Effect of altitude on species richness ##

# First we check how latitude and altitude are linked

## read in glaciological metadata 

metadata <- fread("20240226_nomis_glaciological_metadata.csv")%>%
  separate(location, into = c("sample","site"), sep = "_")%>%
  summarise(.by = sample,
            lat = `lat_sp [DD]`[site == "UP"],
            lon = `lon_sp [DD]`[site == "UP"],
            gl_cov = mean(`gl_cov [%]`),
            ele_sp = `ele_sp [m]`[site == "UP"])%>%
  drop_na()

alphadt_nomis_df<- sample_data(alphadiv_nomis)%>%
  data.frame()%>%
  left_join(metadata, by = "sample")

## Effect of latitude on altitude
lat_ele <- ggplot(data=alphadt_nomis_df,aes(x=lat,y=ele_sp))+
  geom_point()+
  geom_smooth(method="gam", formula = y ~ s(x, bs = 'tp')) +
  theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
# lat_ele

mod_lat <- gam(data=alphadt_nomis_df, formula = ele_sp ~ s(lat, bs='tp'))
summary(mod_lat)

# mod_lat$residuals

## Get the residuals
alphadt_nomis_df$ele_resids = mod_lat$residuals

## We retrieve the residuals of the models and investigate the correlations between altitude or latitude and species richness

pred <- predict(lm(data=alphadt_nomis_df, formula = Observed ~ ele_resids),
                se.fit = TRUE, interval = "confidence")
limits <- as.data.frame(pred$fit)
spe_ele <- ggplot(alphadt_nomis_df, aes(x=ele_resids, y=Observed)) +
  geom_point(size=3, alpha=0.4,color="#3F459BFF") +
  geom_smooth(method='lm', se=T, formula = y ~ x, fill="blue", color="black", alpha=0.2, span=0.3)  +
  geom_line(aes(x = ele_resids, y = limits$lwr),linetype = 2) +
  geom_line(aes(x = ele_resids, y = limits$upr),linetype = 2) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  labs(x = "Elevation model residuals", y = "ASV richness")
spe_ele

#ggsave("images/20240226_richness_elevation_model_residuals.pdf", spe_ele, height = 5, width = 7)
#saveRDS(spe_ele, file = "images/20240226_richness_elevation_model_residuals.rds")

model_elevation<-lm(data=alphadt_nomis_df, formula = Observed ~ ele_resids)
summary(model_elevation)

#write.csv(alphadt_nomis_df, "202408_alphadt_nomis.csv")

# There is no relationship between the %of glacier coverage and the altitude
cov_ele <- ggplot(alphadt_nomis_df,aes(x=(gl_cov),y=ele_sp))+
  geom_point()+
  geom_smooth(method="gam", formula = y ~ s(x, bs = 'tp'))+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
cov_ele

## We test the effect % of glacier coverage on species diversity
pred_cov <- predict(lm(alphadt_nomis_df, formula = Observed ~ gl_cov),
                    se.fit = TRUE, interval = "confidence")
limits_cov<- as.data.frame(pred_cov$fit)

spe_cov <- ggplot(alphadt_nomis_df, aes(x=gl_cov, y=Observed)) +
  geom_point(size=3, alpha=0.4,color="#3F459BFF") +
  geom_smooth(method='lm', se=T, formula = y ~ x, fill="blue", color="black", alpha=0.2, span=0.3)  +
  geom_line(aes(x = gl_cov, y = limits_cov$lwr),linetype = 2) +
  geom_line(aes(x = gl_cov, y = limits_cov$upr),linetype = 2) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  labs(x = "Percent glacier coverage (%)", y = "ASV richness")
spe_cov
#ggsave("images/20240226_richness_coverage.pdf", spe_cov, height = 5, width = 7)
#saveRDS(spe_cov, file = "images/20240226_richness_coverage.rds")

model_cov<-lm(data=alphadt_nomis_df, formula = Observed ~ gl_cov)
summary(model_cov)

