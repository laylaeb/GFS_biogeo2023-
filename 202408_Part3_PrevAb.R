### Prevalence, Relative Abundance ###

library(tidyverse)
library(ggrepel)
library(RColorBrewer)

## Load dataset
dat <- read_csv("20240221_NOMIS_rarefied_deblur_table.csv.gz")
tax <- read_csv("20240221_NOMIS_rarefied_deblur_taxonomy.csv.gz")
meta <- read_tsv("202402_NOMIS_metadata_GFS.tsv")

dat_m <- dat %>%
  left_join(tax)%>%
  group_by(Genus)%>%
  filter(!(Genus == "g__uncultured"))%>%
  summarise_if(is.numeric, sum)%>%
  na.omit()%>%
  pivot_longer(cols = !Genus)%>%
  left_join(meta, by = c("name" = "sample"))%>%
  select(-c(gl_name,bassin))

sum_all <- sum(dat_m$value)
num_all <- length(unique(dat_m$name))

dat_abun <- dat_m %>%
  group_by(Genus)%>%
  summarise(Abundance = sum(value)/sum_all)

colnames(dat_abun) <- c("Genus", "Abundance")

dat_prev <- dat_m %>%
  group_by(Genus, name)%>%
  filter(value > 0)%>%
  ungroup()%>%
  group_by(Genus)%>%
  summarise(Prevalence = n()/num_all)

tax_sel <- tax%>%
  select(-c(...1,Species))%>%
  distinct()%>%
  filter(Genus %in% dat_abun$Genus)


dat_final <- dat_abun%>%
  left_join(dat_prev)%>%
  left_join(tax_sel)

classToPlot <- dat_final %>%
  group_by(Class)%>%
  summarise(sum= sum(Abundance))%>%
  arrange(desc(sum))%>%
  top_n(9)

dat_final %>%
  dplyr::arrange(desc(Abundance)) %>%
  .[1:8, ] -> text_size

dat_final <- dat_final%>%
  mutate(Class = if_else(Class %in% classToPlot$Class, Class, "Other"))

colors <- c(brewer.pal(9, "Set1"), "black")
names(colors) <- c(classToPlot$Class, "Other")

p1 <- ggplot(dat_final, aes(x = Prevalence, y = log10(Abundance), color = Class))+
  geom_point()+
  scale_y_continuous(limits = c(-7,0))+#, breaks = c(1e-5, 1e-4, 1e-3, 1e-2,1e-1, 1))+
  scale_color_manual(values = colors)+
  geom_text_repel(data = text_size, aes(label = Genus))+
  theme_minimal()

p1
ggsave_fitmax("PrevalenceAbundanceNOMIS_Genus.pdf", maxwidth = 10, p1)

dat_sum <- dat_final %>%
  select(Abundance, Prevalence)%>%
  group_by(Prevalence)%>%
  reframe(sum = sum(Abundance))


p2 <- ggplot(dat_sum, aes(x = Prevalence, y = log10(sum)))+
  geom_line()+
  scale_y_continuous(limits = c(-7,0))+
  theme_minimal()+
  geom_smooth(method = "gam")
  
p2

breaks <- (0:10)/10

dat_cut <- dat_final%>%
  mutate( ints = cut(Prevalence ,breaks = 40)) %>% 
  group_by(ints) %>% 
  summarise(sum =  sum(Abundance))

dat_cut$Prevalence <- as.numeric(dat_cut$ints)
  
p3 <- ggplot(dat_cut, aes(x = as.numeric(ints), y = log10(sum)))+
  geom_line()+
  scale_y_continuous(limits = c(-7,0))+
  theme_bw()
p3

# ggsave_fitmax("PrevalenceAbundanceNOMISLine.pdf",p3)

temp <- dat_final %>%
  group_by(Prevalence)%>%
  summarise(sum = sum(Abundance))

temp$cumsum <- cumsum(temp$sum)

p4 <- ggplot(dat_final,aes(x = Prevalence, y = log10(Abundance)))+
  geom_point(aes(color = Class))+
  scale_y_continuous(limits = c(-7,0))+#, breaks = c(1e-5, 1e-4, 1e-3, 1e-2,1e-1, 1))+
  scale_color_manual(values = colors)+
  geom_text_repel(data = text_size, aes(label = Genus))+
  geom_line(data = temp, aes(x = Prevalence, y = log10(cumsum)))+
  # geom_smooth(data = dat_sum, aes(x = Prevalence, y = log10(sum)), method = "gam", se = F)+
  theme_minimal()
  
p4

ggsave_fitmax("PrevalenceAbundanceNOMISLineGenus.pdf",maxwidth = 10, p4)
