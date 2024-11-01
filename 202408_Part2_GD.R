### Gamma diversity ###

library(phyloseq)
library(phyloseqCompanion)
library(biomformat)
library(ggplot2)
library(vegan)
library(knitr)
library(dplyr)
library(data.table)
library(DescTools)
library(devtools)
library(iNEXT)

## Load ASV and taxonomy tables
asv<-read.csv(file="20240221_NOMIS_rarefied_deblur_table.csv.gz",sep=",",header=TRUE,row.names=1)
tax<-read.csv(file="20240221_NOMIS_rarefied_deblur_taxonomy.csv.gz",sep=",",header=TRUE,row.names=1)

## Load metadata file
metadata_NOMIS="202402_NOMIS_metadata_GFS.tsv"
metadata_NOMIS<-import_qiime_sample_data(metadata_NOMIS)

## Create OTU table 
OTU_NOMIS <- otu_table(asv, taxa_are_rows=TRUE)
tax_NOMIS <- tax_table(as.matrix(tax))

## Create phyloseq object
merged_NOMIS_DEBLUR <- merge_phyloseq(OTU_NOMIS, tax_NOMIS, metadata_NOMIS)

## Prune Uganda samples
prune_Uganda <- subset_samples(merged_NOMIS_DEBLUR, !site_c %in% "Uganda")
metadata_nomis_inext <- sample.data.frame(prune_Uganda)

sites <- c("Greenland", "Norway", "Alps", "New_Zealand", "Ecuador", "Caucasus", "Kirghizistan", "Chile", "Alaska", "Nepal")

sample_lists <- list()
asv_lists <- list()
vec_lists <- list()

# Loop through each site and perform the subset and merge operations
for (site in sites) {
   subset_result <- subset_samples(prune_Uganda, site_c == site)
   merged_result <- merge_samples(subset_result,"sample")

   # ## create asv table
   asv_table_site <- otu_table(merged_result, taxa_are_rows=T)
   asv_table_site_t <- t(asv_table_site)
   asv_lists[[site]] <- asv_table_site_t
   #
   # ## calculate rowSums and
   asv_table_df <- asv_table_site_t
   sumrow_site <- unname(rowSums(asv_table_df>0))
   sort_site<- sort(sumrow_site, decreasing=T)
   vec_site <- sort_site[sort_site >0]
   vec_lists[[site]] <- vec_site
}

list_exped_all <- list(alps=c(ncol(asv_lists$Alps),vec_lists$Alps),nz=c(ncol(asv_lists$New_Zealand),vec_lists$New_Zealand),
                       caucasus=c(ncol(asv_lists$Caucasus),vec_lists$Caucasus),
                       kh=c(ncol(asv_lists$Kirghizistan),vec_lists$Kirghizistan),
                       greenland=c(ncol(asv_lists$Greenland),vec_lists$Greenland),
                       norway=c(ncol(asv_lists$Norway),vec_lists$Norway), 
                       ecuador=c(ncol(asv_lists$Ecuador),vec_lists$Ecuador),
                       nepal=c(ncol(asv_lists$Nepal),vec_lists$Nepal),
                       alaska=c(ncol(asv_lists$Alaska),vec_lists$Alaska),
                       chile=c(ncol(asv_lists$Chile),vec_lists$Chile))
                       

out_all_exped <- iNEXT(list_exped_all, q=0, datatype="incidence_freq", se=T, conf=0.95, nboot=99)

df <- fortify(out_all_exped, type =1)

df.point <- df[which(df$Method=="Observed"),]
df.line <- df[which(df$Method!="Observed"),]
df.line$Method <- factor(df.line$Method, 
                         c("Rarefaction", "Extrapolation"),
                       )

df.asympote <- data.frame(y = c(24,10),
                          Asymptote = c("alps","nz","caucasus","kh","greenland","norway","ecuador","nepal","alaska","chile"))


ggplot(df, aes(x=x, y=y, colour=Assemblage)) + 
  #geom_point(aes(shape=Assemblage), size=5, data=df.point) +
  geom_line(aes(linetype= Method), lwd=1.5, data=df.line) +
  geom_ribbon(aes(ymin=y.lwr, ymax=y.upr,
                  fill=Assemblage, colour=NULL), alpha=0.2) +
  labs(x="Number of GFS", y="Species diversity") +
scale_fill_manual(values=c("#2E2A2BFF", "#CF4E9CFF","#8C57A2FF",
                              "#3EBCB6","#82581FFF","#2F509EFF",
                              "#E5614CFF","#97A1A7FF","#DC9445FF","#bee183")
)+
  scale_color_manual(values=c("#2E2A2BFF", "#CF4E9CFF","#8C57A2FF",
                             "#3EBCB6","#82581FFF","#2F509EFF",
                             "#E5614CFF","#97A1A7FF","#DC9445FF","#bee183")
  )+
  scale_linetype_discrete(name ="Method")+
theme_bw() + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                 panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

inext_freq_results <- out_all_exped$AsyEst  
inext_freq_results$prop <- inext_freq_results$Observed/inext_freq_results$Estimator
inext_freq_results<-inext_freq_results[inext_freq_results$Diversity == 'Species richness',]
median_GD_freq <- inext_freq_results %>% 
  summarise(med = median(prop), 
            lower_quartile = quantile(prop, 0.25),
            median = quantile(prop, 0.5),
            upper_quartile = quantile(prop, 0.75))

