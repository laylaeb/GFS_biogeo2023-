### Rank Abundance curves across mountain ranges ###
library(phyloseq)
library(phyloseqCompanion)
library(tidyverse)
library(ggplot2)

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
prune_Uganda <- subset_samples(merged_NOMIS_DEBLUR,!site_c %in% "Uganda")
prune_Uganda_df <- as.matrix(otu_table(prune_Uganda,taxa_are_rows=T))

## Rank abundance Curves, all Glaciers! ## Check 151 GFSs (-1 Uganda)
comm<-as(otu_table(prune_Uganda), "matrix")
env<-as(sample.data.frame(prune_Uganda), "matrix")

comm<-as.data.frame(comm)
env<-as.data.frame(env)

## Applying this to the full dataset
comm_melt <- reshape2::melt(as.matrix(comm), na.rm=TRUE)
colnames(comm_melt)<-c("Var1","sample","Abundance")
head(comm_melt)
comm_melt<-as.data.frame(comm_melt)
comm_melt_scaling<-comm_melt[comm_melt$Abundance!=0,] ## This removes taxa with 0 ab.
head(comm_melt_scaling)

comm_melt_RA<-comm_melt_scaling %>%
  group_by(sample) %>%
  mutate(Scale_RA=Abundance/sum(Abundance)) ## relative abundance
head(comm_melt_RA)

## Sanity check
comm_melt_RA %>%
  group_by(sample)%>% summarize(somme=sum(Scale_RA))

comm_melt_RANK<-comm_melt_RA %>%
  group_by(sample) %>%
  mutate(my_ranks = order(order(Scale_RA,decreasing=T))) #Ranks

## Get everything in order
Rank_merge<-merge(comm_melt_RANK, env, by="sample")

Rank_merge_mod <- expand.grid(site_c = unique(Rank_merge$site_c),
                              Rank = 1:max(Rank_merge$my_ranks))

## Adding the mean value for each mountain range
Rank_merge_mod$Scale_RA = vapply(1:nrow(Rank_merge_mod), 
                                 function(i) mean(Rank_merge$Scale_RA[(Rank_merge$my_ranks == Rank_merge_mod$Rank[i]) & 
                                                                        (Rank_merge$site_c == Rank_merge_mod$site_c[i])]),FUN.VALUE = numeric(1))
colors<-c("#2E2A2BFF", "#CF4E9CFF", "#8C57A2FF",
          "#3EBCB6","#82581FFF","#2F509EFF",
          "#E5614CFF","#97A1A7FF","#bee183","#DC9445FF","#EDD03E")
rare_curve <- ggplot() +
  geom_line(data=Rank_merge, aes(x=my_ranks, y=log(Scale_RA), group=sample),size=0.2,colour="black",alpha=0.2)+
  geom_line(data=Rank_merge_mod, aes(x=Rank, y=log(Scale_RA), colour=site_c))+
  scale_color_manual(values= c("Alaska"="#2E2A2BFF","Alps"="#CF4E9CFF","Caucasus"="#8C57A2FF",
                               "Chile"="#3EBCB6","Ecuador"="#82581FFF","Greenland"="#2F509EFF",
                               "Kirghizistan"="#E5614CFF","Nepal"="#97A1A7FF","New_Zealand"="#bee183","Norway"="#DC9445FF"))+
  xlab("Rank")+
  facet_wrap(~site_c, nrow=3)+
  theme(panel.background = element_rect(fill = 'white', color="grey60"))
