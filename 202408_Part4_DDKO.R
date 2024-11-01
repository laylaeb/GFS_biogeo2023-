### Distance Decay for KOs ###

library(tidyverse)
library(mgcv)
library(data.table)
library(gt)
library(vegan)
library(randomForest)
library(smacof)
library(ampvis2)
library(Maaslin2)
library(readxl)
library(ggsci)
library(ggalluvial)
library(pheatmap)
library(ggpubr)
library(here)
library(permuco)
library(effectsize)
library(permutes)
library(buildmer)
library(geosphere)

# functions
clean_sample_name <- function(raw_name){
  raw_name = gsub('DownB', 'DN', raw_name)
  raw_name = gsub('UpB', 'UP', raw_name)
  raw_name = gsub('Down', 'DN', raw_name)
  clean_name = gsub('Up', 'UP', raw_name)
  return(clean_name)}

# load data
kegg <- read.csv(file = "NOMIS_merged_kegg_counts_20240104.csv", check.names = FALSE)
merged_metadata <- read.csv(file="merged_metadata_Susheel.csv")
proka_filter <- read.csv("2024_KO_Proka.csv")

kegg_filtered <- kegg[kegg$Geneid %in% proka_filter$KO_proka,]

# Identify unique prefixes
prefixes <- unique(sub("_.*", "", names(kegg_filtered)))

# Calculate averages for each prefix
averages_list <- map(prefixes, ~{
  prefix <- .x
  
  # Select columns with "_UP" or "_DN" for the current prefix
  prefix_columns <- select(kegg_filtered, starts_with(paste0(prefix, "_UP")), starts_with(paste0(prefix, "_DN")))
  
  # Calculate row-wise averages for the selected columns
  prefix_average <- rowMeans(prefix_columns, na.rm = TRUE)
  
  # Round the averages to integers
  prefix_average <- round(prefix_average)
  
  # Create a dataframe with the prefix and its corresponding average
  tibble(!!paste0(prefix, "_Average") := prefix_average)
})

# Combine the average dataframes into a single dataframe
averages_kegg_filtered <- bind_cols(averages_list)

# Remove columns where all values are NaN
averages_kegg_filtered <- averages_kegg_filtered %>%
  select(which(!apply(is.na(averages_kegg_filtered), 2, function(x) all(x))))

# Remove the string "_Average" from the column names
names(averages_kegg_filtered) <- gsub("_Average", "", names(averages_kegg_filtered))

# Combine the averages with non-UP and DN columns
result <- as.data.frame(cbind(select(kegg_filtered, -ends_with("_UP"), -ends_with("_DN")), averages_kegg_filtered))

rownames(result) <- result$Geneid
result$Geneid <- NULL

# Remove all rock data
otu_mat <- result %>% 
  select(-matches("R"))

# Preparing data
# creating df with colnames and keeping relevant metadata only
otu_mat_df <- data.frame(Sample = colnames(otu_mat))
otu_mat_df <- otu_mat_df %>% mutate(edited_Sample = str_replace(Sample, "_.*$", ""))

le_meta <- read_tsv("metadata_NOMIS_Leila.txt") 

otu_mat_df <- otu_mat_df %>% left_join(le_meta, by = c("edited_Sample"="sample")) %>% 
  select(-edited_Sample) 

# There are two GFS samples that were not in the taxonomy but we can include them so that we increase the nb of samples
otu_mat_df <- otu_mat_df %>%
  mutate(site_c = ifelse(Sample == 'GL61', 'Ecuador', site_c),
         site_c = ifelse(Sample == 'GL126', 'Kirghizistan', site_c))

otu_mat_df[otu_mat_df$Sample == "GL61", c("lat", "lon")] <- c(0.0402, -77.9969)
otu_mat_df[otu_mat_df$Sample == "GL126", c("lat", "lon")] <- c(39.4713, 72.8967)

# Remove Uganda
otu_mat_df <- otu_mat_df%>%
  filter(site_c!="Uganda")

# Adjust otu_mat with the metadata from otu_df
otu_mat <- otu_mat[colnames(otu_mat) %in% otu_mat_df$Sample]

# Calculate Bray-Curtis dissimilarity matrix
bray_curtis_matrix <- vegdist(t(otu_mat), method = "bray") 
bray_curtis_matrix_mat <- as.matrix(bray_curtis_matrix)

# Calculate Sørensen dissimilarity matrix
sorensen_matrix <- vegdist(t(otu_mat), method = "bray", binary=T) 
sorensen_matrix_mat <- as.matrix(sorensen_matrix)

# Calculate geographic distances between points using geosphere package
# Assuming you have latitude and longitude information in 'matrix2'
otu_mat_dist <- otu_mat_df %>% select(Sample, lon, lat) %>% column_to_rownames(., var = "Sample") %>% as.matrix()
geographic_distances <- distm(otu_mat_dist, fun=distGeo)
rownames(geographic_distances) <- rownames(otu_mat_dist)
colnames(geographic_distances) <- rownames(otu_mat_dist)

# Create a dataframe for Bray-Curtis dissimilarity vs. Geographic distances
bray_curtis_df <- data.frame(
  Dissimilarity  = as.vector(bray_curtis_matrix_mat),
  GeographicDistance = as.vector(geographic_distances),
  Type="Bray-Curtis"
)

braycurtis_df_filtered <- bray_curtis_df %>%
  filter(Dissimilarity != 0)

allbray_m<- as.matrix(bray_curtis_matrix)
diag(allbray_m)=NA
allbray_diag<-t(matrix(t(allbray_m)[which(!is.na(allbray_m))],nrow=83,ncol=84))
min(allbray_diag)

geographic_distances <- distm(otu_mat_dist, fun=distGeo)
dist_geo_all <- as.matrix(geographic_distances)
diag(dist_geo_all)= NA
dist_geo_all_diag<-t(matrix(t(dist_geo_all)[which(!is.na(dist_geo_all))],nrow=83,ncol=84))

mantel(allbray_diag, dist_geo_all_diag, method = "pearson", permutations = 999, na.rm = TRUE)

# # Create the scatter plot for Bray-Curtis dissimilarity vs. Geographic distance
# bray_curtis_plot <- ggplot(bray_curtis_df, aes(x = log10(GeographicDistance), y = 1-Dissimilarity )) +
#   geom_point() +
#   geom_smooth(method = "lm", se = FALSE, color = "blue") +  # Add linear trendline
#   stat_cor(method = "pearson", label.x = 0.8, label.y = 0.6, label.sep = 0.1, size = 4) +  # Add correlation coefficient
#   labs(
#     x = "Geographic Distance",
#     y = "Bray-Curtis Dissimilarity",
#     title = "Bray-Curtis Dissimilarity vs. Geographic Distance",
#   )
# # Print the plot
# bray_curtis_plot

# Create a dataframe for Sørensen dissimilarity vs. Geographic distances
# Check for complete cases (no NA values)
sorensen_df <- data.frame(
  Dissimilarity  = as.vector(sorensen_matrix_mat),
  GeographicDistance = as.vector(geographic_distances),
  Type = "Sorensen"
)

sorensen_df_filtered <- sorensen_df %>%
  filter(Dissimilarity != 0)

# mantel test 
diag(sorensen_matrix_mat)=NA
sorensen_diag<-t(matrix(t(sorensen_matrix_mat)[which(!is.na(sorensen_matrix_mat))],nrow=83,ncol=84))
mantel(sorensen_diag, dist_geo_all_diag, method = "pearson", permutations = 9999, na.rm = TRUE)

# Combine data frames
combined_df <- rbind(braycurtis_df_filtered, sorensen_df_filtered)
#write.csv(combined_df, "202408_combined_dd_function.csv")

get_fishcol <- fish(4, option="Scarus_quoyi")
custom_colors <- c("Bray-Curtis" = "#3F459BFF", "Sorenson" = "#009E9EFF")
ggplot(combined_df) + 
  geom_point(aes(y=1-Dissimilarity, x = GeographicDistance/1000, color= Type), alpha=0.1) + 

  scale_color_manual(values = custom_colors)+
  #scale_color_fish(option="Scarus_quoyi",discrete=T,direction=1)+
  ylab("Community Similarity") +
  xlab("Geographic distance (km)") + ggtitle("") + theme_bw() +
  scale_x_log10() +theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))


# Ensure GeographicDistance is numeric
combined_df$GeographicDistance <- as.numeric(combined_df$GeographicDistance)

# Convert GeographicDistance to kilometers
combined_df$GeographicDistance_km <- combined_df$GeographicDistance / 1000

# Create a function to calculate R-squared
calculate_r_squared <- function(model) {
  return(format(summary(model)$r.squared, digits = 3))
}

# Fit linear regression models using stat_smooth
lm_models <- by(combined_df, combined_df$Type, function(sub_df) {
  lm(1 - sub_df$Dissimilarity ~ sub_df$GeographicDistance_km, data = sub_df)
})

# Extract coefficients from the linear regression models and calculate R-squared
intercepts <- sapply(lm_models, function(model) coef(model)[1])
slopes <- sapply(lm_models, function(model) coef(model)[2])
r_squared <- sapply(lm_models, calculate_r_squared)

get_fishcol <- fish(4, option="Scarus_quoyi")
custom_colors <- c("Bray-Curtis" = "#3F459BFF", "Sorensen" = "#009E9EFF")

ggplot(combined_df, aes(x = GeographicDistance/1000, y = 1 - Dissimilarity, color = Type)) +
  geom_point(alpha = 0.1) +
  scale_color_manual(values = custom_colors) +
  stat_smooth(method = "lm", formula = y ~ x, size = 1.2, se = FALSE, linetype = "solid") +
  ylab("Functional similarity") +
  xlab("Geographic Distance (km)") +
  ggtitle("Scatter Plot of Similarity vs. Geographic Distance in km") +
  theme_bw() +
  scale_x_log10() +
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black")
  )

# Test difference in slopes
dist_bray <- data.frame(sim_bc=as.vector(1-allbray_diag), geo.dist=as.vector(dist_geo_all_diag), Method="All glaciers")
dist_sor <- data.frame(sim_sor=as.vector(1-sorensen_diag), geo.dist=as.vector(dist_geo_all_diag), Method="All glaciers")

merge_bcsor <- merge(dist_bray, dist_sor, by="row.names")
anco_bcsor <- merge_bcsor[c("sim_bc","sim_sor","geo.dist.x")]
manco <- melt(anco_bcsor, id=c("geo.dist.x")) 

ggscatter(
  manco, x = "geo.dist.x", y = "value",
  color = "variable", add = "reg.line"
)+
  stat_regline_equation(
    aes(label =  paste(..eq.label.., ..rr.label.., sep = "~~~~"), color = variable)
  )

manco$geo.dist.log <- log((manco$geo.dist.x)/1000)

glm_dd_bcsor <- (lm(manco, formula = value ~ geo.dist.log + geo.dist.log:variable + variable))# violation of assumptions!

## Non-parametric ANCOVA
formula <- value ~ geo.dist.log + geo.dist.log:variable + variable
model <- lm(formula, data = manco)
## Perform permutation test for ANCOVA
perm_test <- aovperm(model, nperm = 999, method="freedman_lane")
## Extract the resampled p-value for the interaction term
interaction_p_value <- perm_test$results$`geo.dist.log:variable`$p.value


