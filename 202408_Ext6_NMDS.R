### dbMANOVA NMDS ###

library(vegan)
library(RColorBrewer)
library(reshape2)

ASVs <-read.csv("20240221_NOMIS_rarefied_deblur_table.csv.gz", head=T, row.names = 1)
ASVs_uganda <-ASVs %>% select(-GL140)

## Load metadata file
metadata_NOMIS="/Users/ezzat/Leilas_Drive/Projet SBER/NOMIS/NOMIS_16S/NOMIS_August2023/DEBLUR/revision_2/202402_NOMIS_metadata_GFS.tsv"
metadata_NOMIS <-import_qiime_sample_data(metadata_NOMIS)
metadata_glaciers <-sample.data.frame(metadata_NOMIS)
metadata_glaciers_sub <-metadata_glaciers[metadata_glaciers$sample %in% colnames(ASVs_uganda),]

ind.ASVs <-read.csv("2024_0307_indicator_taxa.csv", head=T)

# indicator NMDS
ind.ASVs.long <-reshape2::melt(ind.ASVs, id.vars="X")
# Replace "0" with NA
ind.ASVs.long[ind.ASVs.long == 0] <- NA
ind.ASVs.long <-na.omit(ind.ASVs.long)
ind.ASVs.long <-ind.ASVs.long[order(ind.ASVs.long$X),]

cols=c("#2E2A2BFF","#CF4E9CFF","#8C57A2FF",
       "#3EBCB6","#82581FFF","#2F509EFF",
       "#E5614CFF","#97A1A7FF","#bee183","#DC9445FF","#EDD03E","#000000")

nMDS<-metaMDS(t(ASVs_uganda))
# [1] 0.1541989
plot(nMDS, type="n")
ordispider(nMDS, groups=metadata_glaciers_sub$site_c, col=cols, lwd=2)
points(nMDS, disp="sites", pch=21, bg=cols[factor(metadata_glaciers_sub$site_c)], cex=1.5)

spec.scores <-scores(nMDS, display="species")
spec.scores.db.MANOVA <-data.frame(spec.scores[rownames(spec.scores) %in% ind.ASVs$X,])
spec.scores.db.MANOVA <-spec.scores.db.MANOVA[order(row.names(spec.scores.db.MANOVA)),]
summary(ind.ASVs.long$X==rownames(spec.scores.db.MANOVA))

spec.scores.db.MANOVA <-data.frame(spec.scores[rownames(spec.scores) %in% X$ASV,])
points(spec.scores.db.MANOVA$NMDS1,spec.scores.db.MANOVA$NMDS2, pch=16, cex=0.2, col=cols[factor(ind.ASVs.long$variable)])
#write.csv(spec.scores, "202408_specscores_dbMAN.csv")
