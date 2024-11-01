library(vegan)

ASVs <-read.csv("20240221_NOMIS_rarefied_deblur_table.csv", row.names = 1)
tax <-read.csv("20240221_NOMIS_rarefied_deblur_taxonomy.csv", head=T, row.names = 1)
samp <-read.table("202402_NOMIS_metadata_GFS.tsv", sep="\t", head=T)

ASVs <-ASVs[,colnames(ASVs) %in% samp$sample]

## NMDS of microdiverse genera
pol <-tax[tax$Genus==" g__Polaromonas",]
pols <-ASVs[rownames(ASVs) %in% rownames(pol),]
pol.MDS <-metaMDS(t(pols))
pol.anosim <-anosim(t(pols),samp$site_c)

rhod <-tax[tax$Genus==" g__Rhodoferax",]
rhods <-ASVs[rownames(ASVs) %in% rownames(rhod),]
rhod.MDS <-metaMDS(t(rhods))
rhod.anosim <-anosim(t(rhods),samp$site_c)

meth <-tax[tax$Genus==" g__Methylotenera",]
meths <-ASVs[rownames(ASVs) %in% rownames(meth),]
meth.MDS <-metaMDS(t(meths))
meth.anosim <-anosim(t(meths),samp$site_c)

rhizo <-tax[tax$Genus==" g__Rhizobacter",]
rhizo <-ASVs[rownames(ASVs) %in% rownames(rhizo),]
rhizo <-rhizo[,!(names(rhizo) %in% "GL64")]
rhizo.samp <-samp[samp$sample!="GL64",]
rhizo.MDS <-metaMDS(t(rhizo))
rhizo.anosim <-anosim(t(rhizo),rhizo.samp$site_c)

par(mfrow=c(2,2), mar=c(3,3,2,1))
plot(pol.MDS, main="Polaromonas")
ordispider(pol.MDS, samp$site_c)
legend("topright",paste0("Anosim R: ",round(pol.anosim$statistic,2),"\n", "p-value: ",round(pol.anosim$signif,3)), bty="n")
plot(rhod.MDS, main="Rhodoferax")
ordispider(rhod.MDS, samp$site_c)
legend("topright",paste0("Anosim R: ",round(rhod.anosim$statistic,2),"\n", "p-value: ",round(rhod.anosim$signif,3)), bty="n")
plot(meth.MDS, main="Methylotenera")
ordispider(meth.MDS, samp$site_c)
legend("topright",paste0("Anosim R: ",round(meth.anosim$statistic,2),"\n", "p-value: ",round(meth.anosim$signif,3)), bty="n")
plot(rhizo.MDS,  main="Rhizobacter")
ordispider(rhizo.MDS, rhizo.samp$site_c)
legend("topright",paste0("Anosim R: ",round(rhizo.anosim$statistic,2),"\n", "p-value: ",round(rhizo.anosim$signif,3)), bty="n")




