### Specific ASVs ###

library(speedyseq)
library(phyloseq)
library(tidyverse)
library(phyloseqCompanion)

## Load ASV and taxonomy tables
asv<-read.csv(file="20240221_NOMIS_rarefied_deblur_table.csv.gz",sep=",",header=TRUE,row.names=1)
tax<-read.csv(file="20240221_NOMIS_rarefied_deblur_taxonomy.csv.gz",sep=",",header=TRUE,row.names=1)

## Load metadata file
metadata_NOMIS="202402_NOMIS_metadata_GFS.tsv"
metadata_NOMIS<-import_qiime_sample_data(metadata_NOMIS)

## Create ASV table 
OTU_NOMIS <- otu_table(asv, taxa_are_rows=TRUE)
tax_NOMIS <- tax_table(as.matrix(tax))

merged_NOMIS_DEBLUR <- merge_phyloseq(OTU_NOMIS, tax_NOMIS, metadata_NOMIS)

uganda=c("Uganda")
prune_Uganda <- subset_samples(merged_NOMIS_DEBLUR , !site_c %in% uganda)
prune_Uganda_df <- as.matrix(otu_table(prune_Uganda, taxa_are_rows=T))
row_sums <- rowSums(prune_Uganda_df)

## Subset the dataframe to remove rows where the sum is zero
prune_Uganda_df <- prune_Uganda_df[row_sums > 0, ] ### 54019 ASVs and 151 GFSs

metadata_nomis <- sample.data.frame(prune_Uganda)
asv_df <- as.data.frame(t(otu_table(prune_Uganda_df, taxa_are_rows=T)))

## Here we are investigating the prevalence of ASV across the dataset
howmanyasv<- as.data.frame(colSums(asv_df != 0))
colnames(howmanyasv) <- c("Count_nb")
howmanyasv$ASV <- rownames(howmanyasv)
rownames(howmanyasv) <- NULL

## Melt asv table
asvdfmelt <- melt(as.matrix(asv_df))
## Keep only values that are > 0
asvdfmelt <- asvdfmelt[asvdfmelt$value >0,]

##### Specific ASVs
#up + down - Liste des ASVs qui ont une prévalence de 1: 13594 soit 13594/54019=25.2%
endemic_oneGFS <- howmanyasv %>% 
  filter(Count_nb == 1) 
#write.csv(endemic_oneGFS, "2024_unique_list.csv")

endemic_oneGFS_n <- howmanyasv %>% 
  filter(Count_nb == 1) %>%
  summarise(n = n())

## total number of specific ASVs found in one single GFS
endemic_oneGFS_n/length(unique(row.names(prune_Uganda_df)))

# value is the nb of count and count_nb nb of samples 
merge_equalone <- merge(endemic_oneGFS, asvdfmelt, by.x="ASV",by.y="Var2")
## assign the name of the mountain range
merge_equalone$MR <- vapply((merge_equalone$Var1), function(x) metadata_nomis$site_c[metadata_nomis$sample == x], FUN.VALUE = character(1))
## rename column
colnames(merge_equalone) <- c("ASV","prev","glname","nb_count","MR")

## Nb of ASVs that are range-specific 
endemic_MR_prev <- merge_equalone %>% group_by(MR) %>% summarize(prev=sum(prev))
endemic_oneGFS_plot <- merge_equalone %>% group_by(MR) 
endemic_oneGFS_plot<- endemic_oneGFS_plot[c("MR","ASV","prev")]
endemic_oneGFS_plot$Color <- "C_uniquetoone"

## Calculate unique ASVs per mountain range. This means that these ASVs are found in only one GFS within their mountain range
prop_unique_MR<- endemic_MR_prev%>%
  group_by(MR)%>% 
  summarize(prop=prev/endemic_oneGFS_n)
colnames(prop_unique_MR) <- c("mountain_range","prop_unique")

## Filter the table of prevalence >1 & <10 GFS
twoandnine<- howmanyasv %>% filter(Count_nb > 1 & Count_nb < 10) ## 20468 ASVs 52.5% des ASVs //2024 -> 29445/54019= 54.5% of ASVs found in 1 to 9 GFSs
twoandnine_n <- twoandnine %>% summarize(n())
twoandnine_n/length(unique(row.names(prune_Uganda_df)))

# merging twoandnine with ASV table to get sample id and for each sample we attribute a MR
merge_twoandnine <- merge(twoandnine, asvdfmelt, by.x="ASV",by.y="Var2")
merge_twoandnine$MR <- vapply((merge_twoandnine$Var1), function(x) metadata_nomis$site_c[metadata_nomis$sample == x], FUN.VALUE = character(1))

colnames(merge_twoandnine) <- c("ASV","prev_GFS","glname","nb_count","MR")

# How many times an ASV is identified in a mountain range, sum by gl_code
# we count the nb of lines by ASV and gl_code. 
twoandnine_end<- merge_twoandnine %>% group_by(MR,ASV) %>% summarize(prev=n()) 
twoandnine_end$Color <- "B_twoandnine"

# now we would like it to be unique in this MR! we mutate to count the nb of times this ASVs is observed (in how many MR),then we filter ASVs found in one MR
# then we group by MR and we sum by N to know how many ASVs per MR
endemism_twoandnine <- twoandnine_end %>% group_by(ASV) %>% mutate(n=n()) %>% filter(n==1)%>%
  ungroup()%>% group_by(MR) %>%summarize(number=n()) 

endemism_twoandnine_plot <- twoandnine_end %>% group_by(ASV) %>% mutate(n=n()) %>% filter(n==1)

# Number of ASVs that are specific to a mountain range (present in 2-9 GFSs)
sum_endemic_twoandnine <- endemism_twoandnine %>% summarize(sum=sum(number))
sum_endemic_twoandnine/length(unique(row.names(prune_Uganda_df)))

## Proportion by MR
prop_unique_twonine<- endemism_twoandnine%>%
  group_by(MR)%>% 
  summarize(prop=number/twoandnine_n)

colnames(prop_unique_twonine)<-c("mountain_range", "prop_unique")

tenandmore<- howmanyasv %>% filter(Count_nb >= 10) ## 7284 -> 18.7% of all ASVs //// 2024 10980 ASVs soit 20.3% of ASVs that are found in more than 10 GFSs
tenandmore_n <- tenandmore %>% summarize(n())
tenandmore_n/length(unique(row.names(prune_Uganda_df)))

# merge tenandmore with asv table to get the sample id and then for each sample we attribute a MR 
merge_tenandmore <- merge(tenandmore, asvdfmelt, by.x="ASV",by.y="Var2")
merge_tenandmore$MR <- vapply((merge_tenandmore$Var1), function(x) metadata_nomis$site_c[metadata_nomis$sample == x], FUN.VALUE = character(1))

colnames(merge_tenandmore) <- c("ASV","prev_GFS","glname","nb_count","MR")

# sum ASvs
tenandmore_end<- merge_tenandmore %>% group_by(MR,ASV) %>% summarize(prev=sum(nb_count>0))

tenandmore_end$Color <-"A_tenandemore"

# We would like it to be unique from this MR! mutate to count the nb of times when ASV has been identified, then filter only ASVs that are observed 1X in the MR
# group by MR and sum by N to know how many ASVs per MR

endemism_tenandmore <- tenandmore_end %>% group_by(ASV) %>% mutate(n=n()) %>% filter(n==1)%>%
  ungroup()%>% group_by(MR) %>%summarize(number=n()) ## 329/38983=0.8% (>=10 GFSs).
###2024 -- 648/54019=1.2% of total ASVs are specific to this section (>=10GFs)

endemism_tenandmore_plot <- tenandmore_end %>% group_by(ASV) %>% mutate(n=n()) %>% filter(n==1)

sum(endemism_tenandmore$number)

## Control Sanity check! ##
control_asv <- asvdfmelt
control_asv$MR <- vapply((control_asv$Var1), function(x) metadata_nomis$site_c[metadata_nomis$sample == x], FUN.VALUE = character(1))
## We want to know how many ASVs are specific, and found in only one mountain range! We start from the original ASV table
## First we group by MR and ASV, and we count the nb of times the ASV appears for a given mountain range (in how many GFSs are they present?)
## Then we ungroup and we want to know how many ASVs are specific so we just group by ASV and count the nb of times an ASV is observed within the dataset 
## meaning in how many MR does it appear? Since we would like it to be only in 1 MR, we apply the filter(n==1)..
## and then we count the nb of lines -- meaning the nb of ASVs that are specific (unique to one MR)

controlasv_end<- control_asv %>% group_by(MR,Var2) %>% summarize(prev=n()) %>%
  ungroup()%>% group_by(Var2) %>% mutate(n=n())%>% filter(n==1)%>%
  ungroup()%>% group_by(MR) %>%summarize(number=n())

#2024 dataset
# # A tibble: 10 × 2
# MR           number
# <chr>         <int>
# 1 Alaska         2703
# 2 Alps           3180
# 3 Caucasus       2273
# 4 Chile          4071
# 5 Ecuador        3874
# 6 Greenland       734
# 7 Kirghizistan   3788
# 8 Nepal          4392
# 9 New_Zealand    7085
# 10 Norway         1517

#33617 that's the sum. 33617/54019=62.2%

#write.csv(controlasv_end,"202403_endemic_list.csv")

## Number of Observed ASVs per mountain range
controlasv_total<- control_asv %>% group_by(MR,Var2) %>% summarize(sumi=sum(n()))%>%
  ungroup()%>% group_by(MR)%>%summarize(sumii=sum(n()))

##Data 2024
# MR number sumii
# 1        Alaska   2703 10110
# 2          Alps   3180 14345
# 3      Caucasus   2273 13679
# 4         Chile   4071  9865
# 5       Ecuador   3874  7847
# 6     Greenland    734  8594
# 7  Kirghizistan   3788 11138
# 8         Nepal   4392 13174
# 9   New_Zealand   7085 12839
# 10       Norway   1517  7146


## Proportion of specific ASVs per mountain range
prop_endemic_MR <- merge(controlasv_end, controlasv_total, by="MR")
prop_endemic_MR$prop_endemic <- prop_endemic_MR$number/prop_endemic_MR$sumii

## Data 2024
# MR number sumii prop_endemic
# 1        Alaska   2703 10110   0.26735905
# 2          Alps   3180 14345   0.22168003
# 3      Caucasus   2273 13679   0.16616712
# 4         Chile   4071  9865   0.41267106
# 5       Ecuador   3874  7847   0.49369186
# 6     Greenland    734  8594   0.08540842
# 7  Kirghizistan   3788 11138   0.34009697
# 8         Nepal   4392 13174   0.33338394
# 9   New_Zealand   7085 12839   0.55183426
# 10       Norway   1517  7146   0.21228659

## Proportion of unique ASVs -- 49.3% of specific ASVs are unique to one GFS!
endemic_oneGFS_n/sum(controlasv_end$number)
## Data 2024
# > endemic_oneGFS_n/sum(controlasv_end$number)
# n
# 1 0.4043787

## Graphs
# We rbind the 3 different dataframes to plot the ASVs that are specific to the different mountain ranges
endemic_oneGFS_plot$n <- 1
df_full <- rbind(endemic_oneGFS_plot, endemism_tenandmore_plot, endemism_twoandnine_plot)
niveaux <- c("New_Zealand","Nepal","Alps","Ecuador","Chile","Kirghizistan","Caucasus","Norway","Alaska","Greenland")
df_full$MR<- factor(df_full$MR, levels = niveaux)
df_full$MR <- fct_rev(df_full$MR)

ggplot(df_full, aes(fill=Color, y=MR)) + 
  geom_bar(position="stack", stat="count")+theme_minimal()

## We would like to create a stacked plot with the taxonomy of the three different fractions (unique, 2-9 and >10)
## For all the taxa but also for the endemic only! =)
merge_asv_endemic<- merge(t(asv_df), df_full, by.x="row.names",by.y="ASV")

## remove the columns that we dont need anymore
merge_asv_endemic<- as.data.frame(merge_asv_endemic[,-c(149:152)])
## create phyloseq object
row.names(merge_asv_endemic)<-merge_asv_endemic$Row.names
merge_asv_endemic$Row.names <- NULL

endemic_table <- as.matrix(otu_table(merge_asv_endemic, taxa_are_rows=T))
endemic_table <- endemic_table[rowSums(endemic_table[])>0,] ## remove rows containing 0 values 

merge_endemic_phylo <- merge_phyloseq(endemic_table, metadata_nomis, tax_NOMIS)
#saveRDS(merge_endemic_phylo, "merge_endemic_phylo2023.RDS")
             
endemic_taxglom <- tax_glom(merge_endemic_phylo, taxrank=rank_names(merge_endemic_phylo)[5], NArm=F)
transf_endemic <- transform_sample_counts(endemic_taxglom, function(x) x / sum(x))
sample_merge_region <- merge_samples(transf_endemic, "site_c")
region_endemic <- transform_sample_counts(sample_merge_region, function(x) x / sum(x))

TopASV_f <- names(sort(taxa_sums(region_endemic), TRUE)[1:19])
top15_NOMIS_f <- prune_species(TopASV_f, region_endemic)
top15_NOMIS_f <- prune_taxa(taxa_sums(top15_NOMIS_f)>0, top15_NOMIS_f)
top_family<-as.data.frame(tax_table(top15_NOMIS_f))

## Turn all ASVs into Family counts
endemic_df <- psmelt(region_endemic) # create dataframe from phyloseq object
endemic_df$Family<- as.character(endemic_df$Family) #convert to character

## We put it in "other" only ASVs that are not part of the most abundant families
endemic_df$Family[!(endemic_df$Family %in% top_family$Family)] <- "Other"
endemic_df$Family[(endemic_df$Family == "")] <- "Other"
endemic_df$Family[(endemic_df$Family == "g__uncultured")] <- "Other"

n <- 20
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

barplot_biogeo <- ggplot(data=endemic_df, aes(x=Sample, y=Abundance, fill=Family))
barplot_biogeo + geom_bar(aes(), stat="identity", position="stack") +
  scale_fill_manual(values = c("#7FC97F", "#BEAED4", "#FDC086", "#FFFF99", "#386CB0", "#F0027F", "#BF5B17", "#666666", "#1B9E77", "#D95F02", "#7570B3", "#E7298A","#66A61E",
                               "#E6AB02","#A6761D", "#666666", "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A",
                               "#FFFF99")) +
  theme(legend.position="bottom") + guides(fill=guide_legend(nrow=5))

## % relative abundance of ASVs
prune_uganda_abondance= transform_sample_counts(prune_Uganda, function(x) x / sum(x))
prune_uganda_asv <- otu_table(prune_uganda_abondance, taxa_are_rows=T)
prune_uganda_asv <- prune_uganda_asv[rowSums(prune_uganda_asv[])>0,]
merge_endemic_abondance <- merge(df_full, prune_uganda_asv, by.x="ASV", by.y="row.names")

## remove the columns that we dont need anymore
merge_endemic_abondance<- as.data.frame(merge_endemic_abondance[,-c(2:5)])

## create phyloseq object
row.names(merge_endemic_abondance)<-merge_endemic_abondance$ASV
merge_endemic_abondance$ASV <- NULL
endemic_table_abondance <- otu_table(merge_endemic_abondance, taxa_are_rows=T)
merge_endemic_abondance_phylo <- merge_phyloseq(endemic_table_abondance, tax_NOMIS,metadata_nomis)

end_ab_phylo_table <- (as.matrix(otu_table(merge_endemic_abondance_phylo, taxa_are_rows=T)))

melt_asv <- melt(end_ab_phylo_table)
merge_asv_data <- merge(as.data.frame(melt_asv), as.matrix(metadata_nomis), by.x="Var2",by.y="sample")

## Here we would need to divise by the number of glaciers per mountain ranges 
sum_mr <- merge_asv_data %>% group_by(site_c)%>% summarize(summrr=sum(value), n=n_distinct(Var2))%>% summarize(ar_mr=summrr/n, site_c)

###Data 2024
# # A tibble: 10 × 2
# ar_mr site_c      
# <dbl> <chr>       
# 1 0.0745 Alaska      
# 2 0.0370 Alps        
# 3 0.0287 Caucasus    
# 4 0.255  Chile       
# 5 0.300  Ecuador     
# 6 0.0303 Greenland   
# 7 0.120  Kirghizistan
# 8 0.104  Nepal       
# 9 0.291  New_Zealand 
# 10 0.0370 Norway    

mean(sum_mr$ar_mr)
# [1] 0.1347194
##2024
# > mean(sum_mr$ar_mr)
# [1] 0.1278457

sd(sum_mr$ar_mr)
# [1] 0.1157538
##2024
# > sd(sum_mr$ar_mr)
# [1] 0.111434

median_endemism_ab <- sum_mr %>% 
  summarise(median=median(ar_mr), x = quantile(ar_mr, c(0.25, 0.5, 0.75)))
median_endemism_ab

## This represents the median and other IQR of the mean of the contribution of ASV in term of relative abundance per region
## Data 2024
# # A tibble: 3 × 2
# median      x
# <dbl>  <dbl>
# 1 0.0895 0.0370
# 2 0.0895 0.0895
# 3 0.0895 0.222 

## filter only ASVs > 0 and we are interested in the median of relative abundance 
filtered_merge_asv = merge_asv_data[merge_asv_data$value > 0,]
#write.csv(filtered_merge_asv,"2024_ra_specific.csv")
abasv_spe <- filtered_merge_asv %>% group_by(site_c)%>%summarize(med_asv=median(value), n=n_distinct(Var2))

## Median of relative abundance of ASVs for a given mountain range
# # A tibble: 10 × 3
# site_c         med_asv     n
# <chr>            <dbl> <int>
# 1 Alaska       0.0000663    15
# 2 Alps         0.0000369    26
# 3 Caucasus     0.0000221    19
# 4 Chile        0.0000663    10
# 5 Ecuador      0.0000590    10
# 6 Greenland    0.0000590     7
# 7 Kirghizistan 0.0000369    18
# 8 Nepal        0.0000442    17
# 9 New_Zealand  0.0000295    20
# 10 Norway     0.0000295     9

## The most abundant phyla for the specific ASVs
endemic_taxglom_phyla <- tax_glom(merge_endemic_phylo, taxrank=rank_names(merge_endemic_phylo)[2], NArm=F)
tax_table_end_phyla <- tax_table(endemic_taxglom_phyla)
transf_endemic_phyla = transform_sample_counts(endemic_taxglom_phyla, function(x) x / sum(x))
endemic_region = merge_samples(transf_endemic_phyla, "Site_c")
trans_ra_end_phy = transform_sample_counts(endemic_region, function(x) x / sum(x))

asv_endemic_phy <- otu_table(trans_ra_end_phy, taxa_are_rows=T)
melt_endemic_phy <- psmelt(asv_endemic_phy)
merge_taxo_endemic <- merge(melt_endemic_phy, tax_table_end_phyla, by.x="OTU", by.y="row.names")

## Now we would like to know what are the most abundant phyla across the dataset
sumtot_endo_phylum <- merge_taxo_endemic %>% group_by(Phylum) %>% summarize(sum=sum(Abundance/10))%>%
  filter(!(Phylum %in% c(""," g__uncultured"))) %>%
  slice_max(n=15, order_by=sum)

## Most abundant genera
endemic_taxglom_genera <- tax_glom(merge_endemic_phylo, taxrank=rank_names(merge_endemic_phylo)[6], NArm=F)
tax_table_end_genera <- tax_table(endemic_taxglom_genera)
transf_endemic_genera = transform_sample_counts(endemic_taxglom_genera, function(x) x / sum(x))
endemic_region_genera = merge_samples(transf_endemic_genera, "Site_c")
trans_ra_end_gen = transform_sample_counts(endemic_region_genera, function(x) x / sum(x))

asv_endemic_gen <- otu_table(trans_ra_end_gen, taxa_are_rows=T)
melt_endemic_gen <- psmelt(asv_endemic_gen)
merge_taxo_endemic_gen <- merge(melt_endemic_gen, tax_table_end_genera, by.x="OTU", by.y="row.names")

sumtot_endo_genera <- merge_taxo_endemic_gen %>% group_by(Genus) %>% summarize(sum=sum(Abundance/10))%>%
  filter(!(Genus %in% c(""," g__uncultured"))) %>%
  slice_max(n=15, order_by=sum)

phylumfac = factor(tax_table(prune_Uganda)[, "Phylum"])
classfac = factor(tax_table(prune_Uganda)[, "Class"])
orderfac = factor(tax_table(prune_Uganda)[, "Order"])
familyfac = factor(tax_table(prune_Uganda)[, "Family"])
genusfac = factor(tax_table(prune_Uganda)[, "Genus"])


