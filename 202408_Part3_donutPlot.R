### Donut plots to highlight the different fractions of microbiome - specific, indicator and core taxa ###

library("RColorBrewer")
library("ggplot2")
library(tidyverse)
library(reshape2)
library(ggpubr)

## any reference to the term "endemic" is attributed to "specific"

## Load data tables
dat	<- as.data.frame(otu_table(prune_Uganda_df, taxa_are_rows=T))
cores <- read_csv("20240301_core.csv")
endemic <- read_csv("202403_endemic_list.csv")
indicator <- read_csv("2024_0307_indicator_taxa.csv")
uniq <- read_csv("2024_unique_list.csv")

ASV <- rownames(dat)
ASV <- as.data.frame(ASV)

cores <- cores %>%
  select(ASV)%>%
  mutate(type = "core")

endemic <- endemic %>%
  select(ASV)%>%
  mutate(type = "endemic")

colnames(endemic) <- c("ASV", "type")

indicator$type <- "indicator"
indicator <- indicator[,c(1,12)]
colnames(indicator) <- c("ASV","type")

uniq <- uniq %>%
  select(ASV)%>%
  mutate(type = "uniq")

dat_numb <- rbind(cores, endemic, uniq, indicator)

ASV_other <- ASV %>%
  filter(!(ASV %in% dat_numb$ASV))%>%
  mutate(type = "Other")

datToPlot <- dat_numb %>%
  arrange(ASV, type) %>%
  group_by(ASV)%>%
  summarise(type = str_c(type, collapse="_"))%>%
  ungroup()%>%
  bind_rows(ASV_other)%>%
  group_by(type)%>%
  summarise(sum = n())

datToPlot$type <- factor(datToPlot$type, levels = c("core","core_indicator", "indicator",
                                                    "endemic_indicator", "endemic_uniq", "endemic", "Other"))

colors <- c("#034e7b", "#a6bddb", "#addd8e",
            "#feb24c", "#fc4e2a", "#b10026", "#737373")

names(colors) <- c("core","core_indicator", "indicator",
                   "endemic_indicator", "endemic_uniq", "endemic", "Other")

datToPlot <- datToPlot[order(datToPlot$type),]
# Compute percentages
datToPlot$fraction = datToPlot$sum / sum(datToPlot$sum)

# Compute the cumulative percentages (top of each rectangle)
datToPlot$ymax = cumsum(datToPlot$fraction)

# Compute the bottom of each rectangle
datToPlot$ymin = c(0, head(datToPlot$ymax, n=-1))

datToPlot$labelPosition <- (datToPlot$ymax + datToPlot$ymin) / 2

p1 <- ggplot(datToPlot, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=type)) +
  geom_rect() +
  geom_label( x=4, aes(y=labelPosition, label=type), size=6) +
  scale_fill_manual(values = colors) +
  coord_polar(theta="y") +
  #xlim(c(2, 4)) +
  xlim(c(1, 4)) +
  theme_void() +
  theme(legend.position = "none")
p1

dat_sum <- rowSums(dat[1:151])

dat_sum <- as.data.frame(dat_sum)
dat_sum$ASV <- rownames(dat_sum)

datToPlotCov <- dat_numb %>%
  arrange(ASV, type) %>%
  group_by(ASV)%>%
  summarise(type = str_c(type, collapse="_"))%>%
  ungroup()%>%
  bind_rows(ASV_other)%>%
  left_join(dat_sum)%>%
  filter(!(dat_sum == 0))%>%
  mutate(perc = dat_sum / sum(dat_sum))%>%
  group_by(type)%>%
  summarise(fraction = sum(perc))

sum(datToPlotCov$fraction
    )

datToPlotCov$type <- factor(datToPlotCov$type, levels = c("core","core_indicator", "indicator",
                                                    "endemic_indicator", "endemic_uniq", "endemic", "Other"))
datToPlotCov <- datToPlotCov[order(datToPlotCov$type),]

## Compute percentages
# datToPlotCov$fraction = datToPlotCov$sum / sum(datToPlotCov$sum)

## Compute the cumulative percentages (top of each rectangle)
datToPlotCov$ymax = cumsum(datToPlotCov$fraction)

## Compute the bottom of each rectangle
datToPlotCov$ymin = c(0, head(datToPlotCov$ymax, n=-1))

datToPlotCov$labelPosition <- (datToPlotCov$ymax + datToPlotCov$ymin) / 2

p2 <- ggplot(datToPlotCov, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=type)) +
  geom_rect() +
  geom_label(x=4, aes(y=labelPosition, label=type), size=6) +
  scale_fill_manual(values = colors) +
  coord_polar(theta="y") +
  #xlim(c(2, 4)) +
  xlim(c(1, 4)) +
  theme_void() +
  theme(legend.position = "none")
p2

p3 <- ggarrange(p1, p2)
