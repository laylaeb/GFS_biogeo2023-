### KEGG Alpha diversity ###

library(tidyverse)
library(vegan)
library(data.table)
library(microeco)
library(mgcv)

clean_sample_name <- function(raw_name){
  raw_name = gsub('DownB', 'DN', raw_name)
  raw_name = gsub('UpB', 'UP', raw_name)
  raw_name = gsub('Down', 'DN', raw_name)
  clean_name = gsub('Up', 'UP', raw_name)
  return(clean_name)}

kegg <- read.csv(file = "NOMIS_merged_KEGG_counts_20240104.csv", check.names = FALSE)

#Identify unique prefixes
prefixes <- unique(sub("_.*", "", names(kegg)))

#Calculate averages for each prefix
averages_list <- map(prefixes, ~{
  prefix <- .x
  
  # Select columns with "_UP" or "_DN" for the current prefix
  prefix_columns <- select(kegg, starts_with(paste0(prefix, "_UP")), starts_with(paste0(prefix, "_DN")))
  
  # Calculate row-wise averages for the selected columns
  prefix_average <- rowMeans(prefix_columns, na.rm = TRUE)
  
  # Round the averages to integers
  prefix_average <- round(prefix_average)
  
  # Create a dataframe with the prefix and its corresponding average
  tibble(!!paste0(prefix, "_Average") := prefix_average)
})

## Combine the average dataframes into a single dataframe
averages_kegg <- bind_cols(averages_list)

## Remove columns where all values are NaN
averages_kegg <- averages_kegg %>%
  select(which(!apply(is.na(averages_kegg), 2, function(x) all(x))))

## Remove the string "_Average" from the column names
names(averages_kegg) <- gsub("_Average", "", names(averages_kegg))

## Combine the averages with non-UP and DN columns
result <- cbind(select(kegg, -ends_with("_UP"), -ends_with("_DN")), averages_kegg)

#write.csv(result, "averaged_NOMIS_KEGG_counts.csv", row.names = FALSE)

## Read in metadata
metadata <- fread("20240226_nomis_glaciological_metadata.csv")%>%
  separate(location, into = c("sample","site"), sep = "_")%>%
  summarise(.by = sample,
            lat = `lat_sp [DD]`[site == "UP"],
            lon = `lon_sp [DD]`[site == "UP"],
            gl_cov = mean(`gl_cov [%]`),
            ele_sp = `ele_sp [m]`[site == "UP"])%>%
  drop_na()


## read in averaged KEGG data 
kegg_average <- fread("averaged_NOMIS_KEGG_counts.csv")

prokaryote_kegg <- fread("prokaryoteKOs.txt", header = F, col.names = c("Geneid"))

kegg_average <- kegg_average%>%
  semi_join(prokaryote_kegg)

## Remove Uganda because only 1 GFS
prune_Uganda_df <- kegg_average%>%
  select(-starts_with("3"), -starts_with("GL140"))%>% # this removes Uganda samples
  select(-starts_with("GLR"))%>% # this removes the remaining rock samples
  column_to_rownames("Geneid")

prune_Uganda_df <- prune_Uganda_df[rowSums(prune_Uganda_df[])>0,]

hist(colSums(prune_Uganda_df))

## make log transformation

log_df <- log1p(t(prune_Uganda_df))

## Calculate diversity metrics
microtable_object <- microtable$new(as.data.frame(t(log_df)))

trans_alpha$new(microtable_object)

alphadt_nomis_df<- microtable_object$alpha_diversity%>%
  rownames_to_column("sample")%>%
  left_join(metadata, by = "sample")

ggplot(data=alphadt_nomis_df,aes(x=lat,y=ele_sp))+
  geom_point()+
  geom_smooth(method="gam", formula = y ~ s(x, bs = 'tp')) +
  theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

mod_lat <- gam(data=alphadt_nomis_df, formula = ele_sp ~ s(lat, bs='tp'))
summary(mod_lat)

alphadt_nomis_df$ele_resids = mod_lat$residuals

pred <- predict(lm(data=alphadt_nomis_df, formula = Shannon ~ ele_resids),
                se.fit = TRUE, interval = "confidence")
limits <- as.data.frame(pred$fit)
spe_ele <- ggplot(alphadt_nomis_df, aes(x=ele_resids, y=Shannon)) +
  geom_point(size=3, alpha=0.4,color="#3F459BFF") +
  geom_smooth(method='lm', se=T, formula = y ~ x, fill="blue", color="black", alpha=0.2, span=0.3)  +
  geom_line(aes(x = ele_resids, y = limits$lwr),linetype = 2) +
  geom_line(aes(x = ele_resids, y = limits$upr),linetype = 2) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  labs(x = "Elevation model residuals", y = "KEGG gene Shannon diversity")
spe_ele

#ggsave("images/20240415_prokaryotic_kegg_shannon_elevation_model_residuals.pdf", spe_ele, height = 5, width = 7)
#ggsave("images/20240415_prokaryotic_kegg_shannon_elevation_model_residuals.png", spe_ele, height = 5, width = 7)

model_elevation<-lm(data=alphadt_nomis_df, formula = Shannon ~ ele_resids)
summary(model_elevation)
