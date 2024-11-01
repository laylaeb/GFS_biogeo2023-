## Multivariate analysis and Variation Partitioning ###

require(ggplot2)
require(reshape)
require(scales)
library(scales)
library(ggthemes)
library(ggpubr)
library(phyloseq)
library(phyloseqCompanion)
library(fmsb)
library(ggvanced)
library(fishualize)

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

uganda=c("Uganda")
prune_Uganda <- subset_samples(merged_NOMIS_DEBLUR , !site_c %in% uganda)
prune_Uganda_df <- as.matrix(otu_table(prune_Uganda, taxa_are_rows=T))

prune_Uganda_df <- (otu_table(prune_Uganda, taxa_are_rows=T))
metadata_multi <- sample.data.frame(prune_Uganda)

## Some adjustments to our comm data matrix for later
comm_table <- as.data.frame(otu_table(prune_Uganda, taxa_are_rows=T))
comm_table_t <- as.data.frame(t(comm_table))
comm_table_t$sample<- rownames(comm_table_t)

## Merging the different environmental datasets 
metadata_nomis_env <- read.csv("2024_nomis_environmental_data_v2.csv",header=T)
metadata_nomis_env$gl_code <- sub("_.*", "", metadata_nomis_env$sample_ID)

## Load the minerals
minerals_nomis <- read.csv("minerals_nomis_2023.csv")
minerals_nomis <- minerals_nomis[,c(1:5)]

## Load the climatic variable as well
climatic_nomis <- as.data.frame(read.csv("climate_nomis_2023.csv"))

metadata_multi_db <- merge(metadata_nomis_env, minerals_nomis, by.x="gl_code", by.y="sample", no.dups = T)
metadata_multi_db <- merge(metadata_multi_db, metadata_multi, by.x="gl_code", by.y="sample", no.dups = T)
metadata_multi_db <- merge(metadata_multi_db, climatic_nomis, by.x="gl_code", by.y="Sample", no.dups = T)

## function to add small constant so that we do not encounter any troubles when logtransforming later
add_const <- function(x) {
  min_nonzero <- min(x[which(x > 0)]) 
  return((x + (min_nonzero/2)))
}

## Specify the column names you want to modify
columns_to_modify <- c("chla","scd","water_temp","turb")

## Apply add_const only to the specified columns
metadata_multi_db <- metadata_multi_db %>%
  mutate_at(vars(all_of(columns_to_modify)), add_const)

## Correct values for up and down sites for nutrients
metadata_multi_corrected_updown <- metadata_multi_db %>%
  group_by(gl_code) %>%
  mutate(across(c(srp, NH4, NO2, NO3, doc), ~ ifelse(is.na(.), na.omit(.), .)))

## Include position
metadata_multi_corrected_updown$position <- ifelse(grepl("UP", metadata_multi_corrected_updown$sample_ID), "UP", "DOWN")

## Distance to glacier snout for up and down samples
snspdist_up <- metadata_multi_corrected_updown %>%
  filter(position == "DOWN") %>%
  group_by(position) %>%
  summarize(
    median_value = median(sn_sp_dist),
    q1 = quantile(sn_sp_dist, 0.25),
    q3 = quantile(sn_sp_dist, 0.75),
    iqr = q3 - q1
  )

## Keep the data you want 
multi_data_subset <- metadata_multi_corrected_updown[c("gl_code","site_c","water_temp", "pH", "cond", "turb","doc","srp", "NH4","NO3", "chla","gl_sa","gl_cov","lat_sp.x","lon_sp.x","ele_sp","pr","scd","Clays","Quartz","Calcite","Feldspar","position")]

## Remove lines with NAs
multi_data_subset_na <- na.omit(multi_data_subset)

## Sum inorganic nitrogen into DIN
rowsum_nut <- rowSums(multi_data_subset_na[, c("NO3", "NH4")])
multi_data_subset_na$DIN <- rowsum_nut 

  group_by(gl_code) %>%
  mutate(lat_sp.x = ifelse(all(position == "DOWN"), lat_sp.x, 
                           ifelse(position == "DOWN", first(lat_sp.x[position == "UP"]), lat_sp.x)),
         lon_sp.x = ifelse(all(position == "DOWN"), lon_sp.x, 
                           ifelse(position == "DOWN", first(lon_sp.x[position == "UP"]), lon_sp.x)),
         ele_sp = ifelse(all(position == "DOWN"), ele_sp, 
                         ifelse(position == "DOWN", first(ele_sp[position == "UP"]), ele_sp))) %>%
  ungroup()

## Average the values for each patch
## Select only the numeric columns (excluding lat_sp and lon_sp, gl_code, site_c and position)
numeric_columns <- setdiff(names(multi_data_subset_na_coordinates), c("gl_code","site_c","position","lat_sp.x","lon_sp.x","ele_sp"))

## Group by gl_code
grouped_data <- multi_data_subset_na_coordinates %>%
  group_by(gl_code)

## Extract lon,lat, ele_sp and site_c 
extracted_data <- grouped_data %>%
  group_by(gl_code) %>%
  summarize(lat_sp.x = first(lat_sp.x),
            lon_sp.x = first(lon_sp.x),
            ele_sp = first(ele_sp),
            site_c = first(site_c))  

## Calculate the average for each numeric column
metadata_multi_average <- grouped_data %>%
  group_by(gl_code) %>%
  summarize(across(all_of(numeric_columns), mean, na.rm = F))

## Merge with the original data based on gl_code
metadata_multi_sub_u_with_site <- merge(metadata_multi_average , extracted_data, by = "gl_code", all.x = TRUE)
numeric_cols <- sapply(metadata_multi_sub_u_with_site, is.numeric)

## Plot histograms for each numeric variable
par(mfrow = c(6, 4))  # Set up a 6x4 grid for the plots
for (i in 1:length(numeric_cols)) {
  if (numeric_cols[i]) {  # Check if the column is numeric
    hist(metadata_multi_sub_u_with_site[[i]], 
         main = colnames(metadata_multi_sub_u_with_site)[i],
         xlab = "",
         col = "skyblue", 
         border = "white",
         xlim = range(metadata_multi_sub_u_with_site[[i]])) 
  }
}
## Log-transform data prior to computing dbRDA

## Select columns to log transform
cols_to_log_transform <- c("water_temp", "pH", "cond", "turb", "doc", "srp", "NH4", "NO3", "chla", "pr", "scd", "DIN","Quartz","Calcite","Clays","Feldspar")

metadata_multi_sub_u_with_site_logged <- metadata_multi_sub_u_with_site

## Apply log transformation
metadata_multi_sub_u_with_site_logged[cols_to_log_transform] <- lapply(metadata_multi_sub_u_with_site_logged[cols_to_log_transform], function(x) log(x))
metadata_multi_logged_sub <- metadata_multi_sub_u_with_site_logged[c("gl_code","site_c","water_temp", "pH", "cond", "turb","doc","srp","DIN", "chla","gl_sa","gl_cov","lat_sp.x","lon_sp.x","ele_sp","pr","scd","Feldspar","Calcite","Quartz","Clays")]

## Geographic Distance Matrix
metadata_distance <- as.data.frame(metadata_multi_logged_sub  %>% select(lon_sp.x, lat_sp.x))
metadata_dist_df <- as.data.frame(metadata_multi_logged_sub  %>% select(lon_sp.x, lat_sp.x, ele_sp))
metadata_dist_y <- metadata_dist_df
metadata_dist_df <- as.data.frame(sapply(metadata_dist_df, as.numeric))
library(geosphere)
dist_geo<-distm(metadata_distance, fun=distGeo)

## Environmental data without geographic coordinates and altitude
multi_data_trans_ming <-subset(metadata_multi_logged_sub, select=-c(lat_sp.x, lon_sp.x, ele_sp))

## Community matrix with filtered samples 
asv_community_trim  <- comm_table_t %>% 
  filter(sample %in% multi_data_trans_ming$gl_code)

## We have to remove the last column, which contains a "sample"
row.names(asv_community_trim) == multi_data_trans_ming$sample
colnames(asv_community_trim)[54838] ## sample
Y.com <- asv_community_trim[,1:(dim(asv_community_trim)[2]-1)]

## Convert into numeric values
Y.com <- sapply(Y.com, as.numeric)

## Saving Files 
save(Y.com,file="20240408_commY.Rdata")
save(multi_data_trans_ming, file="20240408_Envdata.Rdata")
save(dist_geo,file= "20240408_geo.Rdata")
save(metadata_dist_df,file="20240408_geocoordinates.Rdata")

