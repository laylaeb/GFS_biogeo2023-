### Core microbiome across the different regions ###

library("reshape2")
library("speedyseq")
library("phyloseq")
library("phyloseqCompanion")

core_data <- as.data.frame(sample_data(prune_Uganda))
asv_count_level <- otu_table(prune_Uganda, taxa_are_rows=T)
asv_count_levelf <- asv_count_level[rowSums(asv_count_level[])>0,]

asv_core_ab  = transform_sample_counts(prune_Uganda, function(x) x / sum(x))
asv_core <- otu_table(asv_core_ab, taxa_are_rows=T)
asv_coref <- asv_core[rowSums(asv_core[])>0,]

tax_core <-as.matrix(tax_table(prune_Uganda))
#write.csv(asv_core, "asv_core.csv")

metadata_nomis <- sample.data.frame(prune_Uganda)

## Calculating the binomial for the proba of occurrence
prob_occ <- function(asv_table, metadata, asv, ab_thr){
  subset = as.data.frame(t(as.matrix(asv_table[rownames(asv_table) == asv,])))
  
  subset[subset < ab_thr] = 0
  subset[subset > 0] = 1
  subset$site_c = vapply(1:nrow(subset), function(x){as.character(metadata$site_c[metadata$sample == rownames(subset)[x]])}, FUN.VALUE = character(1))
  colnames(subset)[colnames(subset) == asv] = 'asv'
  fit = glm(asv ~ site_c, data=subset, family = binomial())
  # fit = gam(asv ~ s(latitude, longitude, type='sos'), data=subset, family=binomial())
  # pred_df = metadata[,c('latitude','longitude')]
  pred_df = expand.grid(site_c=unique(metadata$site_c))
  pred_df$pred = predict(fit, newdata=pred_df, type='response')
  
  pred_out = data.frame(ASV = asv, Ab_thr = ab_thr)
  pred_out$average_pred = mean(pred_df$pred)
  for (site_c in unique(metadata$site_c)){
    pred_out[site_c] = pred_df$pred[pred_df$site_c == site_c]}
  return(pred_out)}

seq_log <- function(v1, v2, n) {exp(seq(from = log(v1), to = log(v2), length.out = n))}

library(foreach)
library(doMC)
registerDoMC(4)

proba_df = foreach(abundance=seq_log(1/5000,1,30), .combine='rbind') %:% foreach(asv=rownames(asv_coref), .combine='rbind') %dopar% {prob_occ(asv_coref, core_data, asv, abundance)}
write.csv(proba_df, "20240301_proba_df.csv")

ggplot(data=proba_df, aes(x=(Ab_thr), y=(average_pred), group=ASV))+
  geom_line() + scale_x_log10() 

core_size <- expand.grid(Abundance = seq_log(1/5000, 1, 30),
                         Prevalence = seq_log(1/5000, 1, 30))

core_size$N = vapply(1:nrow(core_size), function(x){length(unique(proba_df$ASV[(proba_df$Ab_thr > core_size$Abundance[x]) & (proba_df$average_pred > core_size$Prevalence[x])]))}, FUN.VALUE = numeric(1))

ggplot(core_size, aes(x=Abundance,y=Prevalence)) + geom_tile(aes(fill=N)) +
  geom_line(data=data.frame(x=c(0,1),y=c(0,1)), aes(x=x,y=y), color='darkgrey', linetype = "dashed") + geom_point(aes(x=0.0008685272, y=0.2), color='red',size=5) +
  scale_x_log10() + scale_y_log10() + scale_fill_gradient(name = 'Core size', trans = "log", breaks = c(1,10,100,1000),na.value = 'transparent') + theme_linedraw()

#ggsave('core_microbiome_biogeo.pdf', width=6, height=5)

## We can choose a threshold -- we would want the ASVs to be present in at least 6 out of 10 regions
test_fwrite_uncomp <- function(large_table) {
  data.table::fwrite(large_table, "test_dt.csv")
  return(TRUE)
}
test_fwrite_uncomp(proba_df)

# proba_df_core <- fread(file="proba_df_202308.csv",header=TRUE, sep=",")
# format(proba_df_core$Ab_thr, scientific=F)

## Look at the abundance threshold 
proba_df_thre <- 
  filter(proba_df, Ab_thr=="0.00116502202129508")

## Then we want to keep taxa that are present at least in 6 regions out of 10! 
## In the Alps, there are 26 glaciers, so we mainly divided 1/26 to get the threshold 
binary_transform <- sapply(proba_df_thre[,4:13], function(x) ifelse(x > 0.03846154, TRUE, FALSE),
                           USE.NAMES = F)
binary_transform_merge <- merge(x=binary_transform, y=proba_df_thre, by="row.names")

binary_sum <- binary_transform_merge %>% 
  rowwise()%>% mutate(sum=sum(c_across(New_Zealand.x:Norway.x)))%>%filter(sum>5)

binary_tax_merge <- merge(x=binary_sum, y=tax_core, by.x="ASV", by.y=0)

#write.csv(binary_tax_merge,"20240301_binary_tax_merge.csv")
#save(binary_tax_merge, file="binary_tax_merge202308.RData")
#load("binary_tax_merge.RData")

## Relative abundance of the core microbiome
merge_core_abundance <- merge(asv_coref, binary_tax_merge, by.x="row.names", by.y="ASV" )
merge_core_abundance <- merge_core_abundance %>% select(-c(Row.names.y, New_Zealand.x, Alps.x, Nepal.x,
                                                           Kirghizistan.x, Chile.x, Alaska.x, Greenland.x, Caucasus.x, Ecuador.x, Norway.x, Ab_thr, average_pred, New_Zealand.y, Alps.y, Nepal.y, Kirghizistan.y, Chile.y, Alaska.y, Greenland.y, Caucasus.y, Ecuador.y, Norway.y, sum, Kingdom, Phylum, Class, Order, Family, Genus, Species))

unAsIs <- function(X) {
  if("AsIs" %in% class(X)) {
    class(X) <- class(X)[-match("AsIs", class(X))]
  }
  X
}

merge_core_abundance$Row.names<-unAsIs(merge_core_abundance$Row.names)
rownames(merge_core_abundance) <- merge_core_abundance$Row.names
merge_core_abundance$Row.names <- NULL
merge_core_abundance_m <- as.matrix(merge_core_abundance)
merge_core_abundance_final <- otu_table(merge_core_abundance_m, taxa_are_rows=T)

tax_NOMIS <- tax_table(tax_core)

## Merge everything (mapping file, tree and OTUs) 
merged_NOMIS_core_ab<- merge_phyloseq(merge_core_abundance_final, tax_NOMIS, metadata_nomis)
merged_NOMIS_core_ab_otu <- (as.matrix(otu_table(merged_NOMIS_core_ab, taxa_are_rows=T)))

melt_asv <- melt(merged_NOMIS_core_ab_otu)
merge_asv_data <- merge(as.data.frame(melt_asv),as.matrix(metadata_nomis), by.x="Var2",by.y="sample")

## Here we would need to divide by the number of glaciers per mountain ranges. 
sum_mr <- merge_asv_data %>% group_by(site_c)%>% summarize(summrr=sum(value), n=n_distinct(Var2))%>% summarize(ar_mr=summrr/n, site_c)
median(sum_mr$ar_mr)
quantile(sum_mr$ar_mr, prob=c(.25,.5,.75), type=2)

