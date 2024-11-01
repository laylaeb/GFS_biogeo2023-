### Radar mountain range ###
library(phyloseq)
library(phyloseqCompanion)
require(ggplot2)
require(reshape)
require(scales)
library(ggthemes)
library(ggpubr)
library(fmsb)
library(tidyverse)
library(ggvanced)
library(fishualize)
library(gt)

## Load ASV and taxonomy tables
asv <-read.csv(file="20240221_NOMIS_rarefied_deblur_table.csv.gz",sep=",",header=TRUE,row.names=1)
tax <-read.csv(file="20240221_NOMIS_rarefied_deblur_taxonomy.csv.gz",sep=",",header=TRUE,row.names=1)

## Load metadata file
metadata_NOMIS="202402_NOMIS_metadata_GFS.tsv"
metadata_NOMIS <-import_qiime_sample_data(metadata_NOMIS)

## Create ASV table 
OTU_NOMIS <- otu_table(asv, taxa_are_rows=TRUE)
tax_NOMIS <- tax_table(as.matrix(tax))

## Create phyloseq object
merged_NOMIS_DEBLUR <- merge_phyloseq(OTU_NOMIS, tax_NOMIS, metadata_NOMIS)

## Remove Uganda because only 1 GFS
uganda=c("Uganda")
prune_Uganda <- subset_samples(merged_NOMIS_DEBLUR , !site_c %in% uganda)
prune_Uganda_df <- as.matrix(otu_table(prune_Uganda, taxa_are_rows=T))
prune_Uganda_df <- (otu_table(prune_Uganda, taxa_are_rows=T))
metadata_radar_site <- sample.data.frame(prune_Uganda)

## Merging the different environmental datasets - nomis db, minerals, chl and climatic variables
## Load nomis_metadata from nomis_portal
metadata_nomis_env <- read.csv("2024_nomis_environmental_data_v2.csv",header=T)
metadata_nomis_env$gl_code <- sub("_.*", "", metadata_nomis_env$sample_ID)

## Load the minerals
minerals_nomis <- read.csv("minerals_nomis_2023.csv")
minerals_nomis <- minerals_nomis[,c(1:5)]

## Load the climatic variables
climatic_nomis <- as.data.frame(read.csv("climate_nomis_2023.csv"))
metadata_radar_db <- merge(metadata_nomis_env, minerals_nomis, by.x="gl_code", by.y="sample", no.dups = T)
metadata_radar_db <- merge(metadata_radar_db, metadata_radar_site, by.x="gl_code", by.y="sample", no.dups = T)
metadata_radar_db <- merge(metadata_radar_db, climatic_nomis, by.x="gl_code", by.y="Sample", no.dups = T)

## Add small constant to chla values
add_const <- function(x) {
  min_nonzero <- min(x[which(x > 0)])  
  return((x + (min_nonzero/2)))
}

## Specify the column names you want to modify
columns_to_modify <- c("chla")

## Apply add_const only to the specified columns
metadata_radar_db <- metadata_radar_db %>%
  mutate_at(vars(all_of(columns_to_modify)), add_const)

## Keep Nutrient values for similar patch
metadata_radar_corrected_updown <- metadata_radar_db %>%
  group_by(gl_code) %>%
  mutate(across(c(srp, NH4, NO2, NO3, doc), ~ ifelse(is.na(.), na.omit(.), .)))

## Keep the data you want 
radar_data_subset <- metadata_radar_corrected_updown[c("gl_code","site_c","water_temp","pH", "cond", "turb","doc","srp", "NH4","NO3", "sba", "chla","gl_sa","gl_cov","lat_sp.x","lon_sp.x","ele_sp","Feldspar","Clays","Calcite","Quartz","pr","sn_sp_dist","scd")]

## Remove NA values
radar_data_subset_na <- na.omit(radar_data_subset)

## Sum inorganic nitrogen into DIN
rowsum_nut <- rowSums(radar_data_subset_na[, c("NO3", "NH4")])
radar_data_subset_na$DIN <- rowsum_nut 

## Select only the numeric columns 
numeric_columns <- setdiff(names(radar_data_subset_na), c("site_c","gl_code","lat_sp.x","lon_sp.x")) # we have to average values for each patch

## Group by gl_code
grouped_data <- radar_data_subset_na %>%
  group_by(gl_code)

## Extract site_c
extracted_data <- grouped_data %>%
  group_by(gl_code) %>%
  summarize(
            site_c = first(site_c))  # Assuming site_c is the same for each gl_code

## Calculate the average for each numeric column except lat_sp and lon_sp and site_c
metadata_radar_average <- grouped_data %>%
  group_by(gl_code) %>%
  summarize(across(all_of(numeric_columns), mean, na.rm = F))

## Merge with the original data based on gl_code
metadata_radar_sub_u_with_site <- merge(metadata_radar_average, extracted_data, by = "gl_code", all.x = TRUE)

sba_quantiles <- metadata_radar_sub_u_with_site %>%
  summarize(
    q25 = quantile(sba, 0.25, type = 2),
    q50 = quantile(sba, 0.50, type = 2),
    q75 = quantile(sba, 0.75, type = 2)
  )

#write.csv(metadata_radar_sub_u_with_site, "metadata_radar_Hannes.csv")

## Subset for radar plots
metadata_radar_subset_table <- metadata_radar_sub_u_with_site[c('site_c','water_temp', 'pH', 'cond', 'turb', 'gl_sa', 'gl_cov', 'srp', 'DIN','doc','Feldspar','Calcite','Quartz','Clays','chla','sn_sp_dist','sba', 'pr','scd')]
metadata_radar_subset_plot <- metadata_radar_sub_u_with_site[c('site_c','water_temp', 'pH', 'cond', 'turb', 'gl_sa', 'gl_cov', 'srp', 'DIN','doc')]

## Compute quantiles
metadata_median <- as.data.frame(metadata_radar_subset_plot %>% group_by(site_c) %>% summarize_if(is.numeric, list("q25" = ~quantile(., 0.25,type=2),"q75" = ~quantile(., 0.75,type=2),"q50" = ~quantile(., 0.50,type=2))))
median_only <- as.data.frame(metadata_radar_subset_plot %>% group_by(site_c) %>% summarize_if(is.numeric, list("q50" = ~quantile(., 0.50,type=2))))

metadata_median_table<- as.data.frame(metadata_radar_subset_table %>% group_by(site_c) %>% summarize_if(is.numeric, list("q50" = ~quantile(., 0.50,type=2))))

process_site_data <- function(data, site_name) {
  # Choose a region
  metadata_site <- subset(data, site_c == site_name)
  melt_site <- melt(metadata_site)
  
  # Filter data for q25, q50, and q75
  q25_data <- subset(melt_site, grepl("_q25$", variable))
  q50_data <- subset(melt_site, grepl("_q50$", variable))
  q75_data <- subset(melt_site, grepl("_q75$", variable))
  
  # Pivot the data
  q25_pivot <- pivot_wider(q25_data, names_from = variable, values_from = value)
  q50_pivot <- pivot_wider(q50_data, names_from = variable, values_from = value)
  q75_pivot <- pivot_wider(q75_data, names_from = variable, values_from = value)
  
  # rename column
  colnames(q25_pivot) <- c("group","water_temp","pH","cond","turb","gl_sa","gl_cov","SRP","DIN","doc")
  colnames(q50_pivot) <- c("group","water_temp","pH","cond","turb","gl_sa","gl_cov","SRP","DIN","doc")
  colnames(q75_pivot) <- c("group","water_temp","pH","cond","turb","gl_sa","gl_cov","SRP","DIN","doc")
  
  # Combine q25, q50, and q75 data
  min <- data.frame(matrix(rep(c(0), 10), nrow = 1))
  colnames(min)<-c("group","water_temp","pH","cond","turb","gl_sa","gl_cov","SRP","DIN","doc")

  combined_data <- rbind(min, q75_pivot, q50_pivot)
  
  # Add a new column "group1"
  combined_data$group1 <- c("min", "max", "Q50")
  combined_data$group <- NULL
  
  combined_data <- combined_data %>%
    select(group1, everything())
  
  return(combined_data)
}

generate_and_save_spider_plot <- function(data, site_name) {
  gg <- ggspider(data, axis_name_offset = 0.1, background_color = "white", fill_opacity = 0.4) +
    scale_fill_manual(values = c("Q25" = "black", "Q50" = "#3F459BFF", "Q75" = "white")) +
    scale_color_manual(values = c("Q25" = "black", "Q50" = "#3F459BFF", "Q75" = "#3F459BFF"))
  
  # Start a PDF device
  pdf(paste0(site_name, "_spider_plot.pdf"), width = 8, height = 8)
  
  # Print the plot
  print(gg)
  
  # Close the PDF device
  dev.off()
}

## Compute dataframe
alaska_data <- process_site_data(metadata_median, "Alaska")
caucasus_data <- process_site_data(metadata_median, "Caucasus")
norway_data <- process_site_data(metadata_median, "Norway")
nepal_data <- process_site_data(metadata_median, "Nepal")
alps_data <- process_site_data(metadata_median, "Alps")
greenland_data <- process_site_data(metadata_median, "Greenland")
ecuador_data <- process_site_data(metadata_median, "Ecuador")
chile_data <- process_site_data(metadata_median, "Chile")
kh_data <- process_site_data(metadata_median, "Kirghizistan")
nz_data <- process_site_data(metadata_median, "New_Zealand")


## Generate and save the plots
generate_and_save_spider_plot(alaska_data, "Alaska")
generate_and_save_spider_plot(caucasus_data, "Caucasus")
generate_and_save_spider_plot(norway_data, "Norway")
generate_and_save_spider_plot(nepal_data, "Nepal")
generate_and_save_spider_plot(alps_data, "Alps")
generate_and_save_spider_plot(greenland_data, "Greenland")
generate_and_save_spider_plot(ecuador_data, "Ecuador")
generate_and_save_spider_plot(chile_data, "Chile")
generate_and_save_spider_plot(nz_data, "New_Zealand")
generate_and_save_spider_plot(kh_data, "Kirghizistan")

scale_column_to_0_100 <- function(x) { # Function to scale values between 0 and 100 for each column 
  rescaled <- scale(x, center = FALSE, scale = max(x))
  scaled <- rescaled * 100
  return(scaled)
}

## Apply scaling to each numeric column
for (col_name in colnames(median_only)[-1]) {
  median_only[, col_name] <- scale_column_to_0_100(median_only[, col_name])
}

median_only_rescaled <-median_only

## Now we need to include two additional rows max and min to each of the variable
min <- data.frame( 
  site_c = "Min",
  water_temp_q50 = 0,
  pH_q50 = 0,
  cond_q50 = 0,
  turb_q50 = 0,
  gl_sa_q50 = 0,
  gl_cov_q50 =0,
  srp_q50 = 0,
  DIN_q50 = 0,
  doc_q50=0
)

max <- data.frame(
  site_c = "max",
  water_temp_q50 = 100,
  pH_q50 = 100,
  cond_q50 = 100,
  turb_q50 = 100,
  gl_sa_q50 = 100,
  gl_cov_q50 =100,
  srp_q50 = 100,
  DIN_q50 = 100,
  doc_q50=100
)

## Replace spr into SRP
colnames(median_only_rescaled)[colnames(median_only_rescaled) == "SRP_q50"] <- "srp_q50"

## Now bind the different columns together
bind_columns <- rbind(max,min,median_only_rescaled)
rownames(bind_columns) <- bind_columns$site_c
bind_columns$site_c <- NULL

colors <- c("#2E2A2BFF", "#CF4E9CFF","#8C57A2FF",
            "#3EBCB6","#82581FFF","#2F509EFF",
            "#E5614CFF","#97A1A7FF","#DC9445FF","#bee183")
titles <- c("Alaska", "Alps", "Caucasus","Chile","Ecuador","Greenland","Kirghizistan","Nepal","New_Zealand","Norway")

op <- par(mar = c(1, 1, 1, 1))# Adjust the margins as needed
par(mfrow = c(1, 2))# Create a single row with 2 columns

## Create the radar chart 2 by 2 
create_beautiful_radarchart <- function(data, color = "#00AFBB", 
                                        vlabels = colnames(data), vlcex = 0.7,
                                        caxislabels = NULL, title = NULL, ...){
  radarchart(
    data, axistype = 1,
    # Customize the polygon
    pcol = color, pfcol = scales::alpha(color, 0.5), plwd = 2, plty = 1,
    # Customize the grid
    cglcol = "grey", cglty = 1, cglwd = 0.8,
    # Customize the axis
    axislabcol = "grey", 
    # Variable labels
    vlcex = vlcex, vlabels = vlabels,
    caxislabels = caxislabels, title = title, ...
  )
}

op <- par(mar = c(1, 1, 1, 1))
par(mfrow = c(1, 2))

for(i in 9:10){
  create_beautiful_radarchart(
    data = bind_columns[c(1, 2, i+2), ], caxislabels =  c(0, 25, 50, 75, 100),
    color = colors[i], title = titles[i]
  )
}
par(op)

## Export table 

tbl_q50 <- gt(metadata_median_table) %>%
  tab_spanner(
    label = "Median",
    columns = pH_q50:sn_sp_dist_q50
  ) %>%
  cols_label(
    site_c = "Mountain range",
    water_temp_q50 = "Water Temp",
    pH_q50 = "pH",
    cond_q50 = "Cond",
    turb_q50 = "Turb",
    gl_sa_q50 = "Gl_Sa",
    gl_cov_q50 = "Gl_Cov",
    srp_q50 = "SRP",
    DIN_q50 = "DIN",
    doc_q50 = "DOC",
    chla_q50 = "Chla",
    Clays_q50 = "Clays",
    Feldspar_q50 = "Feldspar",
    Calcite_q50 = "Calcite",
    Quartz_q50 = "Quartz",
    pr_q50 = "Pr",
    scd_q50 = "Scd",
    sba_q50 = "Sba",
    sn_sp_dist_q50 = "Sn_Sp_Dist"
  ) %>%
  tab_header(
    title = "Summary Statistics",
    subtitle = "Water Quality Parameters"
  )
tbl_q50 |> gtsave("tab_4.docx")


tbl_q25 <- gt(metadata_median_table) %>%
  tab_spanner(
    label = "q25",
    columns = pH_q25:sn_sp_dist_q25
  ) %>%
  cols_label(
    site_c = "Mountain range",
    water_temp_q25 = "Water Temp",
    pH_q25 = "pH",
    cond_q25 = "Cond",
    turb_q25 = "Turb",
    gl_sa_q25 = "Gl_Sa",
    gl_cov_q25 = "Gl_Cov",
    srp_q25 = "SRP",
    DIN_q25 = "DIN",
    doc_q25 = "DOC",
    chla_q25 = "Chla",
    Clays_q25 = "Clays",
    Feldspar_q25 = "Feldspar",
    Calcite_q25 = "Calcite",
    Quartz_q25 = "Quartz",
    pr_q25 = "Pr",
    scd_q25 = "Scd",
    sba_q25 = "Sba",
    sn_sp_dist_q25 = "Sn_Sp_Dist"
  ) %>%
  tab_header(
    title = "Summary Statistics",
    subtitle = "Water Quality Parameters"
  )

tbl_q25 |> gtsave("tab_1.docx")

tbl_q75 <- gt(metadata_median_table) %>%
  tab_spanner(
    label = "q75",
    columns = pH_q75:sn_sp_dist_q75
  ) %>%
  cols_label(
    site_c = "Mountain range",
    water_temp_q75 = "Water Temp",
    pH_q75 = "pH",
    cond_q75 = "Cond",
    turb_q75 = "Turb",
    gl_sa_q75 = "Gl_Sa",
    gl_cov_q75 = "Gl_Cov",
    srp_q75 = "SRP",
    DIN_q75 = "DIN",
    doc_q75 = "DOC",
    chla_q75 = "Chla",
    Clays_q75 = "Clays",
    Feldspar_q75 = "Feldspar",
    Calcite_q75 = "Calcite",
    Quartz_q75 = "Quartz",
    pr_q75 = "Pr",
    scd_q75 = "Scd",
    sba_q75 = "Sba",
    sn_sp_dist_q75 = "Sn_Sp_Dist"
  ) %>%
  tab_header(
    title = "Summary Statistics",
    subtitle = "Water Quality Parameters"
  )

tbl_q75 |> gtsave("tab_3.docx")
#write.csv(metadata_radar_sub_u_with_site,"202408_metadata_radar_graphs.csv")

######## Graphs of environmental variables ########
### Water temp
## Calculate IQR
iqr_value <- quantile(metadata_radar_sub_u_with_site$water_temp, probs = c(0.25, 0.75))
iqr <- iqr_value[2] - iqr_value[1]
median_value <- median(metadata_radar_sub_u_with_site$water_temp)

gg_water <- ggplot(metadata_radar_sub_u_with_site, aes(x = site_c, y = water_temp, color = site_c)) +
  geom_violin(width = 0.9, fill = "grey", alpha = 0.5) +
  geom_boxplot(width = 0.2, color = "black", alpha = 0.9, outlier.shape = NA) +
  geom_jitter(position = position_jitter(0.2), alpha = 0.7) +
  theme_bw() +
  scale_colour_manual(values = c("#2E2A2BFF", "#CF4E9CFF", "#8C57A2FF",
                                 "#3EBCB6", "#82581FFF", "#2F509EFF",
                                 "#E5614CFF", "#97A1A7FF", "#bee183", "#DC9445FF")) +
  
  
  geom_hline(yintercept = median_value, color = "black", size = 2) +
  geom_hline(yintercept = median_value - iqr / 2, linetype = "dashed", color = "black", size = 1) +
  geom_hline(yintercept = median_value + iqr / 2, linetype = "dashed", color = "black", size = 1) +
  
 
  ylim(c(min(metadata_radar_sub_u_with_site$water_temp) - 1, max(metadata_radar_sub_u_with_site$water_temp) + 1)) +
  
  
  theme(legend.position = "none") +
  
  
  ggtitle("Water Temperature Distribution") +
  xlab("Site") +
  ylab("Water Temperature")

mean_value_wt <- mean(metadata_radar_sub_u_with_site$water_temp)
sd_value_wt <- sd(metadata_radar_sub_u_with_site$water_temp)

## pH
iqr_value <- quantile(metadata_radar_sub_u_with_site$pH, probs = c(0.25, 0.75))
iqr <- iqr_value[2] - iqr_value[1]
median_value <- median(metadata_radar_sub_u_with_site$pH)

gg_pH <- ggplot(metadata_radar_sub_u_with_site, aes(x = site_c, y = pH, color = site_c)) +
  geom_violin(width = 0.9, fill = "grey", alpha = 0.5) +
  geom_boxplot(width = 0.2, color = "black", alpha = 0.9, outlier.shape = NA) +
  geom_jitter(position = position_jitter(0.2), alpha = 0.7) +
  theme_bw() +
  scale_colour_manual(values = c("#2E2A2BFF", "#CF4E9CFF", "#8C57A2FF",
                                 "#3EBCB6", "#82581FFF", "#2F509EFF",
                                 "#E5614CFF", "#97A1A7FF", "#bee183", "#DC9445FF")) +
  
  
  geom_hline(yintercept = median_value, color = "black", size = 2) +
  geom_hline(yintercept = median_value - iqr / 2, linetype = "dashed", color = "black", size = 1) +
  geom_hline(yintercept = median_value + iqr / 2, linetype = "dashed", color = "black", size = 1) +
  
  
  ylim(c(min(metadata_radar_sub_u_with_site$pH) - 1, max(metadata_radar_sub_u_with_site$pH) + 1)) +
  
  
  theme(legend.position = "none") +
  
 
  ggtitle("pH") +
  xlab("Site") +
  ylab("pH")

mean_value_ph <- mean(metadata_radar_sub_u_with_site$pH)
sd_value_ph <- sd(metadata_radar_sub_u_with_site$pH)

## Conductivity
Q <- quantile(metadata_radar_sub_u_with_site$cond, probs=c(.25, .75), na.rm = FALSE) 

iqr <- IQR(metadata_radar_sub_u_with_site$cond)
up <-  Q[2]+1.5*iqr # Upper Range  
low<- Q[1]-1.5*iqr # Lower Range

outliers_cond<- subset(metadata_radar_sub_u_with_site,metadata_radar_sub_u_with_site$cond > (Q[1] - 1.5*iqr) & metadata_radar_sub_u_with_site$cond < (Q[2]+1.5*iqr))

iqr_value <- quantile(outliers_cond$cond, probs = c(0.25, 0.75))
iqr <- iqr_value[2] - iqr_value[1]
median_value <- median(outliers_cond$cond)

gg_cond <- ggplot(outliers_cond, aes(x = site_c, y = cond, color = site_c)) +
  geom_violin(width = 1.5, fill = "grey", alpha = 0.5) +
  geom_boxplot(width = 0.2, color = "black", alpha = 0.9, outlier.shape = NA) +
  geom_jitter(position = position_jitter(0.2), alpha = 0.7) +
  theme_bw() +
  scale_colour_manual(values = c("#2E2A2BFF", "#CF4E9CFF", "#8C57A2FF",
                                 "#3EBCB6", "#82581FFF", "#2F509EFF",
                                 "#E5614CFF", "#97A1A7FF", "#bee183", "#DC9445FF")) +
  

  geom_hline(yintercept = median_value, color = "black", size = 2) +
  geom_hline(yintercept = median_value - iqr / 2, linetype = "dashed", color = "black", size = 1) +
  geom_hline(yintercept = median_value + iqr / 2, linetype = "dashed", color = "black", size = 1) +
  
  
  ylim(c(min(outliers_cond$cond) - 1, max(outliers_cond$cond) + 1)) +
  
 
  theme(legend.position = "none") +
  
  # Labels and title
  ggtitle("Conductivity") +
  xlab("Site") +
  ylab("Conductivity")

mean_value_cond <- mean(outliers_cond$cond)
sd_value_cond <- sd(outliers_cond$cond)

## Turbidity
df <- subset(metadata_radar_sub_u_with_site, !((metadata_radar_sub_u_with_site$gl_code) %in% "GL27"))
Q <- quantile(df$turb, probs=c(.25, .50, .75), na.rm = FALSE)

# Calculate IQR and median
iqr_value <- quantile(df$turb, probs = c(0.25, 0.75))
iqr <- iqr_value[2] - iqr_value[1]
median_value <- median(df$turb)

# Create the plot
gg_turbidity <- ggplot(df, aes(x = site_c, y = turb, color = site_c)) +
  geom_violin(width = 2, fill = "grey", alpha = 0.5) + # Adjusted width
  geom_boxplot(width = 0.2, color = "black", alpha = 0.9, outlier.shape = NA) +
  geom_jitter(position = position_jitter(0.2), alpha = 0.7) +
  theme_bw() +
  scale_colour_manual(values = c("#2E2A2BFF", "#CF4E9CFF", "#8C57A2FF",
                                 "#3EBCB6", "#82581FFF", "#2F509EFF",
                                 "#E5614CFF", "#97A1A7FF", "#bee183", "#DC9445FF")) +
  
  geom_hline(yintercept = median_value, color = "black", size = 2) +
  geom_hline(yintercept = median_value - iqr / 2, linetype = "dashed", color = "black", size = 1) +
  geom_hline(yintercept = median_value + iqr / 2, linetype = "dashed", color = "black", size = 1) +
  
  ylim(c(min(df$turb) - 1.5 * iqr, max(df$turb) + 1.5 * iqr)) +
 
  theme(legend.position = "none") +
  
  ggtitle("Turbidity") +
  xlab("Site") +
  ylab("Turbidity")

print(gg_turbidity)

## SRP

df_srp <- subset(metadata_radar_sub_u_with_site, !((metadata_radar_sub_u_with_site$gl_code) %in% c("GL142","GL160")))
Q <- quantile(df_srp$Quartz, probs=c(.25, .50, .75), na.rm = FALSE)

iqr_value <- quantile(df_srp$srp, probs = c(0.25, 0.75))
iqr <- iqr_value[2] - iqr_value[1]
median_value <- median(df_srp$srp)

# Create the plot
gg_srp <- ggplot(df_srp, aes(x = site_c, y = srp, color = site_c)) +
  geom_violin(width = 2, fill = "grey", alpha = 0.5) + # Adjusted width
  geom_boxplot(width = 0.2, color = "black", alpha = 0.9, outlier.shape = NA) +
  geom_jitter(position = position_jitter(0.2), alpha = 0.7) +
  theme_bw() +
  scale_colour_manual(values = c("#2E2A2BFF", "#CF4E9CFF", "#8C57A2FF",
                                 "#3EBCB6", "#82581FFF", "#2F509EFF",
                                 "#E5614CFF", "#97A1A7FF", "#bee183", "#DC9445FF")) +
  
  geom_hline(yintercept = median_value, color = "black", size = 2) +
  geom_hline(yintercept = median_value - iqr / 2, linetype = "dashed", color = "black", size = 1) +
  geom_hline(yintercept = median_value + iqr / 2, linetype = "dashed", color = "black", size = 1) +

  ylim(c(min(df_srp$srp) - 1.5 * iqr, max(df_srp$srp) + 1.5 * iqr)) +
  
  theme(legend.position = "none") +
  
  ggtitle("SRP") +
  xlab("Site") +
  ylab("SRP")

print(gg_srp)

mean_value_srp <- mean(metadata_radar_sub_u_with_site$srp)
sd_value_srp <- sd(metadata_radar_sub_u_with_site$srp)

## DIN
Q <- quantile(metadata_radar_sub_u_with_site$DIN, probs=c(.25, .50,.75), na.rm = FALSE)

iqr_value <- quantile(metadata_radar_sub_u_with_site$DIN, probs = c(0.25, 0.75))
iqr <- iqr_value[2] - iqr_value[1]
median_value <- median(metadata_radar_sub_u_with_site$DIN)

# Create the plot
gg_DIN <- ggplot(metadata_radar_sub_u_with_site, aes(x = site_c, y = DIN, color = site_c)) +
  geom_violin(width = 2, fill = "grey", alpha = 0.5) + # Adjusted width
  geom_boxplot(width = 0.2, color = "black", alpha = 0.9, outlier.shape = NA) +
  geom_jitter(position = position_jitter(0.2), alpha = 0.7) +
  theme_bw() +
  scale_colour_manual(values = c("#2E2A2BFF", "#CF4E9CFF", "#8C57A2FF",
                                 "#3EBCB6", "#82581FFF", "#2F509EFF",
                                 "#E5614CFF", "#97A1A7FF", "#bee183", "#DC9445FF")) +
  
  geom_hline(yintercept = median_value, color = "black", size = 2) +
  geom_hline(yintercept = median_value - iqr / 2, linetype = "dashed", color = "black", size = 1) +
  geom_hline(yintercept = median_value + iqr / 2, linetype = "dashed", color = "black", size = 1) +
  
  ylim(c(min(metadata_radar_sub_u_with_site$DIN) - 1.5 * iqr, max(metadata_radar_sub_u_with_site$DIN) + 1.5 * iqr)) +
  
  theme(legend.position = "none") +
  
  ggtitle("DIN") +
  xlab("Site") +
  ylab("DIN")

## DOC
df_doc <- subset(metadata_radar_sub_u_with_site, !((metadata_radar_sub_u_with_site$gl_code) %in% c("GL142","GL160")))

Q <- quantile(metadata_radar_sub_u_with_site$doc, probs=c(.25, .75), na.rm = FALSE)
# how to find outliers  - calculate Interquartile Range
iqr <- IQR(metadata_radar_sub_u_with_site$doc)
up <-  Q[2]+1.5*iqr # Upper Range  
low<- Q[1]-1.5*iqr # Lower Range

#outliers_doc<- subset(metadata_radar_sub_u_with_site,metadata_radar_sub_u_with_site$doc > (Q[1] - 1.5*iqr) & metadata_radar_sub_u_with_site$doc < (Q[2]+1.5*iqr))

iqr_value <- quantile(metadata_radar_sub_u_with_site$doc, probs = c(0.25, 0.75))
iqr <- iqr_value[2] - iqr_value[1]
median_value <- median(metadata_radar_sub_u_with_site$doc)

# Create the plot
gg_DOC <- ggplot(metadata_radar_sub_u_with_site, aes(x = site_c, y = doc, color = site_c)) +
  geom_violin(width = 2, fill = "grey", alpha = 0.5) + # Adjusted width
  geom_boxplot(width = 0.2, color = "black", alpha = 0.9, outlier.shape = NA) +
  geom_jitter(position = position_jitter(0.2), alpha = 0.7) +
  theme_bw() +
  scale_colour_manual(values = c("#2E2A2BFF", "#CF4E9CFF", "#8C57A2FF",
                                 "#3EBCB6", "#82581FFF", "#2F509EFF",
                                 "#E5614CFF", "#97A1A7FF", "#bee183", "#DC9445FF")) +
  geom_hline(yintercept = median_value, color = "black", size = 2) +
  geom_hline(yintercept = median_value - iqr / 2, linetype = "dashed", color = "black", size = 1) +
  geom_hline(yintercept = median_value + iqr / 2, linetype = "dashed", color = "black", size = 1) +
  ylim(c(min(metadata_radar_sub_u_with_site$doc) - 1.5 * iqr, max(metadata_radar_sub_u_with_site$doc) + 1.5 * iqr)) +
  theme(legend.position = "none") +
  ggtitle("DOC") +
  xlab("Site") +
  ylab("DOC")

## Chlorophyll a
Q <- quantile(metadata_radar_sub_u_with_site$chla, probs=c(.25, .75), na.rm = FALSE)
# how to find outliers  - calculate Interquartile Range
iqr <- IQR(metadata_radar_sub_u_with_site$chla)
up <-  Q[2]+1.5*iqr # Upper Range  
low<- Q[1]-1.5*iqr # Lower Range

outliers_chla<- subset(metadata_radar_sub_u_with_site, metadata_radar_sub_u_with_site$chla > (Q[1] - 1.5*iqr) & metadata_radar_sub_u_with_site$chla < (Q[2]+1.5*iqr))
iqr_value <- quantile(metadata_radar_sub_u_with_site$chla, probs = c(0.25, 0.75))
iqr <- iqr_value[2] - iqr_value[1]
median_value <- median(metadata_radar_sub_u_with_site$chla)

# Create the plot
gg_chla <- ggplot(metadata_radar_sub_u_with_site, aes(x = site_c, y = chla, color = site_c)) +
  geom_violin(width = 2, fill = "grey", alpha = 0.5) + # Adjusted width
  geom_boxplot(width = 0.2, color = "black", alpha = 0.9, outlier.shape = NA) +
  geom_jitter(position = position_jitter(0.2), alpha = 0.7) +
  theme_bw() +
  scale_colour_manual(values = c("#2E2A2BFF", "#CF4E9CFF", "#8C57A2FF",
                                 "#3EBCB6", "#82581FFF", "#2F509EFF",
                                 "#E5614CFF", "#97A1A7FF", "#bee183", "#DC9445FF")) +
  geom_hline(yintercept = median_value, color = "black", size = 2) +
  geom_hline(yintercept = median_value - iqr / 2, linetype = "dashed", color = "black", size = 1) +
  geom_hline(yintercept = median_value + iqr / 2, linetype = "dashed", color = "black", size = 1) +
  ylim(c(min(metadata_radar_sub_u_with_site$chla) - 1.5 * iqr, max(metadata_radar_sub_u_with_site$chla) + 1.5 * iqr)) +
  theme(legend.position = "none") +
  ggtitle("chla") +
  xlab("Site") +
  ylab("chla")

## Clays
Q <- quantile(metadata_radar_sub_u_with_site$Clays, probs=c(.25, .50,.75), na.rm = FALSE)

iqr_value <- quantile(metadata_radar_sub_u_with_site$Clays, probs = c(0.25, 0.75))
iqr <- iqr_value[2] - iqr_value[1]
median_value <- median(metadata_radar_sub_u_with_site$Clays)

# Create the plot
gg_Clays <- ggplot(metadata_radar_sub_u_with_site, aes(x = site_c, y = Clays, color = site_c)) +
  geom_violin(width = 2, fill = "grey", alpha = 0.5) + # Adjusted width
  geom_boxplot(width = 0.2, color = "black", alpha = 0.9, outlier.shape = NA) +
  geom_jitter(position = position_jitter(0.2), alpha = 0.7) +
  theme_bw() +
  scale_colour_manual(values = c("#2E2A2BFF", "#CF4E9CFF", "#8C57A2FF",
                                 "#3EBCB6", "#82581FFF", "#2F509EFF",
                                 "#E5614CFF", "#97A1A7FF", "#bee183", "#DC9445FF")) +
  geom_hline(yintercept = median_value, color = "black", size = 2) +
  geom_hline(yintercept = median_value - iqr / 2, linetype = "dashed", color = "black", size = 1) +
  geom_hline(yintercept = median_value + iqr / 2, linetype = "dashed", color = "black", size = 1) +
  ylim(c(min(metadata_radar_sub_u_with_site$Clays) - 1.5 * iqr, max(metadata_radar_sub_u_with_site$Clays) + 1.5 * iqr)) +
  theme(legend.position = "none") +
  ggtitle("Clays") +
  xlab("Site") +
  ylab("Clays")

## Calcite
Q_Calcite <- quantile(metadata_radar_sub_u_with_site$Calcite, probs=c(.25, .50, .75), na.rm = FALSE)

iqr_value <- quantile(metadata_radar_sub_u_with_site$Calcite, probs = c(0.25, 0.75))
iqr <- iqr_value[2] - iqr_value[1]
median_value <- median(metadata_radar_sub_u_with_site$Calcite)

# Create the plot
gg_Calcite <- ggplot(metadata_radar_sub_u_with_site, aes(x = site_c, y = Calcite, color = site_c)) +
  geom_violin(width = 2, fill = "grey", alpha = 0.5) + # Adjusted width
  geom_boxplot(width = 0.2, color = "black", alpha = 0.9, outlier.shape = NA) +
  geom_jitter(position = position_jitter(0.2), alpha = 0.7) +
  theme_bw() +
  scale_colour_manual(values = c("#2E2A2BFF", "#CF4E9CFF", "#8C57A2FF",
                                 "#3EBCB6", "#82581FFF", "#2F509EFF",
                                 "#E5614CFF", "#97A1A7FF", "#bee183", "#DC9445FF")) +
  geom_hline(yintercept = median_value, color = "black", size = 2) +
  geom_hline(yintercept = median_value - iqr / 2, linetype = "dashed", color = "black", size = 1) +
  geom_hline(yintercept = median_value + iqr / 2, linetype = "dashed", color = "black", size = 1) +
  ylim(c(min(metadata_radar_sub_u_with_site$Calcite) - 1.5 * iqr, max(metadata_radar_sub_u_with_site$Calcite) + 1.5 * iqr)) +
  theme(legend.position = "none") +
  ggtitle("Calcite") +
  xlab("Site") +
  ylab("Calcite")


## Quartz
Q <- quantile(metadata_radar_sub_u_with_site$Quartz, probs=c(.25, .50, .75), na.rm = FALSE)

iqr_value <- quantile(metadata_radar_sub_u_with_site$Quartz, probs = c(0.25, 0.75))
iqr <- iqr_value[2] - iqr_value[1]
median_value <- median(metadata_radar_sub_u_with_site$Quartz)

# Create the plot
gg_Quartz <- ggplot(metadata_radar_sub_u_with_site, aes(x = site_c, y = Quartz, color = site_c)) +
  geom_violin(width = 2, fill = "grey", alpha = 0.5) + # Adjusted width
  geom_boxplot(width = 0.2, color = "black", alpha = 0.9, outlier.shape = NA) +
  geom_jitter(position = position_jitter(0.2), alpha = 0.7) +
  theme_bw() +
  scale_colour_manual(values = c("#2E2A2BFF", "#CF4E9CFF", "#8C57A2FF",
                                 "#3EBCB6", "#82581FFF", "#2F509EFF",
                                 "#E5614CFF", "#97A1A7FF", "#bee183", "#DC9445FF")) +
  geom_hline(yintercept = median_value, color = "black", size = 2) +
  geom_hline(yintercept = median_value - iqr / 2, linetype = "dashed", color = "black", size = 1) +
  geom_hline(yintercept = median_value + iqr / 2, linetype = "dashed", color = "black", size = 1) +
  ylim(c(min(metadata_radar_sub_u_with_site$Quartz) - 1.5 * iqr, max(metadata_radar_sub_u_with_site$Quartz) + 1.5 * iqr)) +
  theme(legend.position = "none") +
  ggtitle("Quartz") +
  xlab("Site") +
  ylab("Quartz")

## Feldspar
Q <- quantile(metadata_radar_sub_u_with_site$Feldspar, probs=c(.25, .50,.75), na.rm = FALSE)
iqr_value <- quantile(metadata_radar_sub_u_with_site$Feldspar, probs = c(0.25, 0.75))
iqr <- iqr_value[2] - iqr_value[1]
median_value <- median(metadata_radar_sub_u_with_site$Feldspar)

# Create the plot
gg_Feldspar <- ggplot(metadata_radar_sub_u_with_site, aes(x = site_c, y = Feldspar, color = site_c)) +
  geom_violin(width = 2, fill = "grey", alpha = 0.5) + # Adjusted width
  geom_boxplot(width = 0.2, color = "black", alpha = 0.9, outlier.shape = NA) +
  geom_jitter(position = position_jitter(0.2), alpha = 0.7) +
  theme_bw() +
  scale_colour_manual(values = c("#2E2A2BFF", "#CF4E9CFF", "#8C57A2FF",
                                 "#3EBCB6", "#82581FFF", "#2F509EFF",
                                 "#E5614CFF", "#97A1A7FF", "#bee183", "#DC9445FF")) +
  geom_hline(yintercept = median_value, color = "black", size = 2) +
  geom_hline(yintercept = median_value - iqr / 2, linetype = "dashed", color = "black", size = 1) +
  geom_hline(yintercept = median_value + iqr / 2, linetype = "dashed", color = "black", size = 1) +
  ylim(c(min(metadata_radar_sub_u_with_site$Feldspar) - 1.5 * iqr, max(metadata_radar_sub_u_with_site$Feldspar) + 1.5 * iqr)) +
  theme(legend.position = "none") +
  ggtitle("Feldspar") +
  xlab("Site") +
  ylab("Feldspar")

## Glacier surface area
Q <- quantile(metadata_radar_sub_u_with_site$gl_sa, probs=c(.25, .50,.75), na.rm = FALSE)

iqr_value <- quantile(metadata_radar_sub_u_with_site$gl_sa, probs = c(0.25, 0.75))
iqr <- iqr_value[2] - iqr_value[1]
median_value <- median(metadata_radar_sub_u_with_site$gl_sa)

# Create the plot
gg_gl_sa <- ggplot(metadata_radar_sub_u_with_site, aes(x = site_c, y = gl_sa, color = site_c)) +
  geom_violin(width = 2, fill = "grey", alpha = 0.5) + # Adjusted width
  geom_boxplot(width = 0.2, color = "black", alpha = 0.9, outlier.shape = NA) +
  geom_jitter(position = position_jitter(0.2), alpha = 0.7) +
  theme_bw() +
  scale_colour_manual(values = c("#2E2A2BFF", "#CF4E9CFF", "#8C57A2FF",
                                 "#3EBCB6", "#82581FFF", "#2F509EFF",
                                 "#E5614CFF", "#97A1A7FF", "#bee183", "#DC9445FF")) +
  geom_hline(yintercept = median_value, color = "black", size = 2) +
  geom_hline(yintercept = median_value - iqr / 2, linetype = "dashed", color = "black", size = 1) +
  geom_hline(yintercept = median_value + iqr / 2, linetype = "dashed", color = "black", size = 1) +
  ylim(c(min(metadata_radar_sub_u_with_site$gl_sa) - 1.5 * iqr, max(metadata_radar_sub_u_with_site$gl_sa) + 1.5 * iqr)) +
  theme(legend.position = "none") +
  ggtitle("gl_sa") +
  xlab("Site") +
  ylab("gl_sa")

## % Glacier coverage
Q <- quantile(metadata_radar_sub_u_with_site$gl_cov, probs=c(.25, .50,.75), na.rm = FALSE)

iqr_value <- quantile(metadata_radar_sub_u_with_site$gl_cov, probs = c(0.25, 0.75))
iqr <- iqr_value[2] - iqr_value[1]
median_value <- median(metadata_radar_sub_u_with_site$gl_cov)

# Create the plot
gg_gl_cov <- ggplot(metadata_radar_sub_u_with_site, aes(x = site_c, y = gl_cov, color = site_c)) +
  geom_violin(width = 2, fill = "grey", alpha = 0.5) + # Adjusted width
  geom_boxplot(width = 0.2, color = "black", alpha = 0.9, outlier.shape = NA) +
  geom_jitter(position = position_jitter(0.2), alpha = 0.7) +
  theme_bw() +
  scale_colour_manual(values = c("#2E2A2BFF", "#CF4E9CFF", "#8C57A2FF",
                                 "#3EBCB6", "#82581FFF", "#2F509EFF",
                                 "#E5614CFF", "#97A1A7FF", "#bee183", "#DC9445FF")) +
  geom_hline(yintercept = median_value, color = "black", size = 2) +
  geom_hline(yintercept = median_value - iqr / 2, linetype = "dashed", color = "black", size = 1) +
  geom_hline(yintercept = median_value + iqr / 2, linetype = "dashed", color = "black", size = 1) +
  ylim(c(min(metadata_radar_sub_u_with_site$gl_cov) - 1.5 * iqr, max(metadata_radar_sub_u_with_site$gl_cov) + 1.5 * iqr)) +
  theme(legend.position = "none") +
  ggtitle("gl_cov") +
  xlab("Site") +
  ylab("gl_cov")


ggarrange(gg_water,gg_pH,gg_cond,gg_turbidity,gg_DIN,gg_srp,gg_DOC,gg_chla,gg_Clays,gg_Calcite,gg_Quartz,gg_Feldspar,nrow=3, ncol=4,align = "hv")
ggarrange(gg_DIN,gg_srp,gg_DOC,gg_chla)
ggarrange(gg_Clays,gg_Calcite,gg_Quartz,gg_Feldspar)
ggarrange(gg_gl_sa,gg_gl_cov)

