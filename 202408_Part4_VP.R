### Variance partitioning ###
library(here)
library(ggplot2)
library(tidyverse)
library(geosphere)   
library(vegan)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(svd)
library(propack.svd)

# load functions (source Legendre 2016) - not currently used (need to install package AEM)
'PCNM' <-
  function(matdist, thresh=NULL, dbMEM=FALSE, moran=NULL, all=FALSE, include.zero=FALSE, silent=FALSE)
    #
    # Compute the PCNM or dbMEM eigenfunctions corresponding to
    # all eigenvalues (+, 0, -).
    #    In PCNM computation, the diagonal of D = 0.
    #    In dbMEM, the diagonal of D = 4*threshh.
    #    Distance-based MEM are described in Dray et al. 2006.
    #    The name was abbreviated to db-MEM by PPN & PL (subm.)
    # Input file: distance matrix produced by the function "dist".
    # Computation of the threshold requires a function of the library "ape".
    #
    # Original PCNM function: Stephane Dray, November 11, 2004
# The present version: Pierre Legendre, August 2007, January and March 2009
  {
    require(vegan)
    epsilon <- sqrt(.Machine$double.eps)
    a <- system.time({
      if(is.null(moran)) {
        if(dbMEM) { moran=FALSE } else { moran=TRUE }
      }
      single <- FALSE
      if(moran) {
        # cat("The site coordinates were computed from 'matdist'.",'\n')
        pcoa.xy <- pcoa.all(matdist)
        
        if(is.na(pcoa.xy$values[2]) | (pcoa.xy$values[2] < epsilon)) {
          if(!silent) cat("The sites form a straight line on the map.",'\n')
          xy <- pcoa.xy$vectors
          single <- TRUE
        } else {
          xy <- pcoa.xy$vectors[,1:2]
        }
      }
      
      matdist <- as.matrix(matdist)
      n <- nrow(matdist)
      
      # Truncation of distance matrix
      if(is.null(thresh)) {
        spanning <- vegan::spantree(as.dist(matdist))
        threshh <- max(spanning$dist)
        if(!silent) cat("Truncation level =",threshh+0.000001,'\n')
      } else {
        threshh = thresh
        if(!silent) cat("User-provided truncation threshold =",thresh,'\n')
      }
      matdist[matdist > threshh] <- 4*threshh
      
      if(dbMEM==FALSE) { diagonal <- 0 } else { diagonal <- 4*threshh }
      
      mypcnm.all <- pcoa.all(matdist, diagonal=diagonal, all=all, include.zero=include.zero, rn=rownames(matdist))
      
      # Compute Moran's I
      if(moran) {
        require(AEM)
        if(single) {
          nb <- dnearneigh(matrix(c(xy,rep(0,n)),n,2), 0, (threshh + epsilon))
        } else {
          nb <- dnearneigh(xy, 0, (threshh + epsilon))
        }
        fr.to.pcnm2 <- as.matrix(listw2sn(nb2listw(nb))[,1:2])
        weight.dist.coord.mat <- as.matrix(1-(as.dist(matdist)/(4*threshh))^2)
        weight <- weight.dist.coord.mat[fr.to.pcnm2]
        res <- moran.I.multi(mypcnm.all$vectors, link=fr.to.pcnm2, weight=weight)
        Moran <- res$res.mat[,1:2]
        positive <- rep(FALSE,length(mypcnm.all$values))
        positive[which(Moran[,1] > res$expected)] <- TRUE
        Moran <- cbind(as.data.frame(Moran), positive)
        colnames(Moran) <- c("Moran","p.value","Positive")
      }
    })
    a[3] <- sprintf("%2f",a[3])
    if(!silent) cat("Time to compute PCNMs =",a[3]," sec",'\n')
    if(is.null(thresh)) {
      if(moran) {
        res <- list(values=mypcnm.all$values, vectors=mypcnm.all$vectors, Moran_I=Moran, expected_Moran=res$expected, spanning=spanning, thresh=threshh+0.000001)
      } else {
        res <- list(values=mypcnm.all$values, vectors=mypcnm.all$vectors, spanning=spanning, thresh=threshh+0.000001)
      }
    } else {
      if(moran) {
        res <- list(values=mypcnm.all$values, vectors=mypcnm.all$vectors, Moran_I=Moran, expected_Moran=res$expected, thresh=thresh)
      } else {
        res <- list(values=mypcnm.all$values, vectors=mypcnm.all$vectors, thresh=threshh+0.000001)
      }
    }
    res
  }
'pcoa.all' <- function(D, diagonal=0, all=FALSE, include.zero=FALSE, rn=NULL)
  # Principal coordinate decomposition of a square distance matrix D
  # Get the eigenvectors corresponding to all eigenvalues, positive and negative
  # Pierre Legendre, 2005, 2007
  #
  # D : A distance matrix of class 'dist' or 'matrix'.
  # all : If TRUE, the eigenvectors corresponding to all eigenvalues, positive and negative, are shown in the output list.
  # include.zero : If FALSE (default value), the zero eigenvalues as well as their eigenvectors are excluded from the output list.
  # rn : An optional vector of row names, of length n, for the objects.
{
  epsilon <- sqrt(.Machine$double.eps)
  # replace by:     epsilon <- .Machine$double.eps * 10^2
  D <- as.matrix(D)
  n <- nrow(D)
  D <- D + diag(rep(diagonal,n))
  
  # Gower centring, matrix formula
  One <- matrix(1,n,n)
  mat <- diag(n) - One/n
  Dpr2 <- -0.5 * mat %*% (D^2) %*% mat
  trace <- sum(diag(Dpr2))
  
  # Eigenvalue decomposition
  D.eig <- eigen(Dpr2, symmetric=TRUE)
  rel.values <- D.eig$values/trace
  rel.cum <- cumsum(rel.values)
  if(length(rn)!=0) {
    rownames(D.eig$vectors) <- rn
  } else {
    rownames(D.eig$vectors) <- rownames(D)
  }
  
  # Output the results: k eigenvalues and eigenvectors
  if(all) {
    select <- 1:n
    if(!include.zero) {
      exclude <- which(abs(D.eig$values) < epsilon)
      select <- select[-exclude]
    }
    k <- length(select)
    res <- list(values=D.eig$values[select], rel.values=rel.values[select], rel.cum.values=rel.cum[select], vectors=D.eig$vectors[,select], trace=trace)
    # cat("k =",k,"Select =",select,'\n')
    
  } else {
    
    k <- length(which(D.eig$values > epsilon))        
    weight <- sqrt(D.eig$values[1:k])
    if(k == 1) {
      vectors <- D.eig$vectors[,1]*sqrt(D.eig$values[1])
    } else {
      vectors <- D.eig$vectors[,1:k]%*%diag(weight)
    }
    res <- list(values=D.eig$values[1:k], rel.values=rel.values[1:k], rel.cum.values=rel.cum[1:k], vectors=vectors, trace=trace)
  }
  res
}

## load data
load("comY_20240408.Rdata")# asv_community_trim
load("Env_20240408.Rdata")# multi_data_na_filter
load("Geo_2024048.Rdata")# dist_geo
load("geocoordinates_20240408.Rdata")#metadata_dist_df

## log1p the community matrix
bio.b <- log1p(Y.com)

## Rename environmental data for ease of use
## We checked several times and we had to remove Clays + Feldspar for multicollinearity
env.data <- multi_data_trans_ming
env.data <- multi_data_trans_ming %>%
  select(-gl_sa.x,-Feldspar) ## We remove them for multicollinearity !

env.names <- c("water_temp", "pH","cond", "turb","doc", "srp", "mean_chla", "gl_cov.x","Calcite","Quartz","pr","scd","DIN")
env.mat <- multi_data_trans_ming[,env.names]
env.mat <- as.data.frame(sapply(env.mat, as.numeric))

## rename geographic coordinates for ease of us
xyz.dat <- metadata_dist_df

## check data structure ## nb of glaciers
with(env.data, table(site_c))

## 1) Linear trend in community structure 
rda.linear <- dbrda(bio.b~., data=xyz.dat, distance="bray")
anova(rda.linear, step=1000)

RsquareAdj(rda.linear)$adj.r.squared 

## Nested spatial model 
## Analyse of microbial spatial variation among and within regions by means of a two-level spatial model
## The among-region component is modeled by a set of dummy variables (N-1 variables with N the number of regions)
## The within-region component is modeled by a set of db-MEM variables for each region
## The db-MEM variables were arranged in blocks corresponding to each region
## within each block, all sites belonging to other regions received the value 0
## See Borcard et al. 2011 (Num ecol with R), Declerck et al. (2011) and function create.MEM.model()

# creating data.frame to store the dbMEMs
var.data.hier <- env.data
var.data.hier <- as.tibble(var.data.hier)
library(geosphere)
##create matrix of geographic coordinates
xy.dat = xyz.dat[,1:2]

## ****************************** ##
## SUBSET = Alps                  ##
## ****************************** ##
i = "Alps"
var.data.temp <- subset(env.data, site_c == i)
xy.dat.temp <- subset(xy.dat, env.data$site_c == i)
xyz.dat.temp <- subset(xyz.dat, env.data$site_c == i)
dist_geo.temp <- distm(xy.dat.temp, fun=distGeo)
bio.b.temp <- subset(bio.b, env.data$site_c == i)
dbMEM.temp <- pcnm(dist_geo.temp)
sel.pos <- which(dbMEM.temp$values > 0)
## Eigenfunctions with positive spatial correlation
dbMEM.region <- as.matrix(dbMEM.temp$vectors[,sel.pos])
names(dbMEM.region)[2:(length(sel.pos)+1)] <- paste(i, names(dbMEM.region)[2:(length(sel.pos)+1)], sep=".")
dbMEM.region <- cbind(data.frame(sample=var.data.temp$sample), dbMEM.region)
var.data.hier <- merge(var.data.hier, dbMEM.region, by="sample", sort=F, all.x=T, nomatch=0)

## ****************************** ##
## SUBSET = Greenland             ##
## ****************************** ##
i = "Greenland"
var.data.temp <- subset(env.data, site_c == i)
xy.dat.temp <- subset(xy.dat, env.data$site_c == i)
xyz.dat.temp <- subset(xyz.dat, env.data$site_c == i)
dist_geo.temp <- distm(xy.dat.temp, fun=distGeo)
bio.b.temp <- subset(bio.b, env.data$site_c == i)
dbMEM.temp <- pcnm(dist_geo.temp)
sel.pos <- which(dbMEM.temp$values > 0)

## Eigenfunctions with positive spatial correlation
dbMEM.region <- as.matrix(dbMEM.temp$vectors[,sel.pos])
dbMEM.region <- cbind(data.frame(sample=var.data.temp$sample), dbMEM.region)
names(dbMEM.region)[2:(length(sel.pos)+1)] <- paste(i, names(dbMEM.region)[2:(length(sel.pos)+1)], sep=".")
var.data.hier <- merge(var.data.hier, dbMEM.region, by="sample", sort=F, all.x=T, nomatch=0)

# dbMEM.region <- as.matrix(dbMEM.temp$vectors[,sel.pos])
# rda.region <- dbrda(bio.b.temp~. + Condition(as.matrix(xy.dat.temp)), data=as.data.frame(dbMEM.region),distance="bray")
# R2a <- RsquareAdj(rda.region)$adj.r.squared  
# mod0 <- rda(bio.b.temp ~ 1 + Condition(as.matrix(xy.dat.temp)) , data=as.data.frame(dbMEM.region))
# rda.region.fwd <- ordistep(mod0, scope=formula(rda.region), direction="forward", perm.max=200)
# rda.region.fwd$call
# 
# R2a <- RsquareAdj(rda.region)$adj.r.squared   # R2a = 0.169
# eigenvals(rda.region)[1:2]/sum(eigenvals(rda.region))*RsquareAdj(rda.region)$adj.r.squared/RsquareAdj(rda.region)$r.squared
# axes.region <- scores(rda.region, choices=c(1:2), display="lc", scaling=1)
# pdf(here("Figures", "Figure_Greeland.pdf"), width=12, height=12)
# PCNM.region(axes.region, 1, "Axe1 (16,9%%)")
# dev.off()

## ****************************** ##
## SUBSET = Caucasus              ##
## ****************************** ##
i = "Caucasus"
var.data.temp <- subset(env.data, site_c == i)
xy.dat.temp <- subset(xy.dat, env.data$site_c == i)
xyz.dat.temp <- subset(xyz.dat, env.data$site_c == i)
dist_geo.temp <- distm(xy.dat.temp, fun=distGeo)
bio.b.temp <- subset(bio.b, env.data$site_c == i)
dbMEM.temp <- pcnm(dist_geo.temp)
sel.pos <- which(dbMEM.temp$values > 0)

## Eigenfunctions with positive spatial correlation
dbMEM.region <- as.matrix(dbMEM.temp$vectors[,sel.pos])
dbMEM.region <- cbind(data.frame(sample=var.data.temp$sample), dbMEM.region)
names(dbMEM.region)[2:(length(sel.pos)+1)] <- paste(i, names(dbMEM.region)[2:(length(sel.pos)+1)], sep=".")
var.data.hier <- merge(var.data.hier, dbMEM.region, by="sample", sort=F, all.x=T, nomatch=0)


# dbMEM.region <- as.matrix(dbMEM.temp$vectors[,sel.pos])
# rda.region <- dbrda(bio.b.temp~. + Condition(as.matrix(xy.dat.temp)), data=as.data.frame(dbMEM.region),distance="bray")
# R2a <- RsquareAdj(rda.region)$adj.r.squared  
# mod0 <- rda(bio.b.temp ~ 1 + Condition(as.matrix(xy.dat.temp)) , data=as.data.frame(dbMEM.region))
# rda.region.fwd <- ordistep(mod0, scope=formula(rda.region), direction="forward", perm.max=200)
# rda.region.fwd$call
# 
# eigenvals(rda.region)[1:2]/sum(eigenvals(rda.region))*RsquareAdj(rda.region)$adj.r.squared/RsquareAdj(rda.region)$r.squared
# axes.region <- scores(rda.region, choices=c(1:2), display="lc", scaling=1)
# pdf(here("Figures", "Figure_Caucasus"), width=12, height=12)
# PCNM.region(axes.region, 1, "Axe1 (5.1%)")
# dev.off()


## ****************************** ##
## SUBSET = Chile                 ##
## ****************************** ##
i = "Chile"
var.data.temp <- subset(env.data, site_c == i)
xy.dat.temp <- subset(xy.dat, env.data$site_c == i)
xyz.dat.temp <- subset(xyz.dat, env.data$site_c == i)
dist_geo.temp <- distm(xy.dat.temp, fun=distGeo)
bio.b.temp <- subset(bio.b, env.data$site_c == i)
dbMEM.temp <- pcnm(dist_geo.temp)
sel.pos <- which(dbMEM.temp$values > 0)

## Eigenfunctions with positive spatial correlation
dbMEM.region <- as.matrix(dbMEM.temp$vectors[,sel.pos])
dbMEM.region <- cbind(data.frame(sample=var.data.temp$sample), dbMEM.region)
names(dbMEM.region)[2:(length(sel.pos)+1)] <- paste(i, names(dbMEM.region)[2:(length(sel.pos)+1)], sep=".")
var.data.hier <- merge(var.data.hier, dbMEM.region, by="sample", sort=F, all.x=T, nomatch=0)


# dbMEM.region <- as.matrix(dbMEM.temp$vectors[,sel.pos])
# rda.region <- dbrda(bio.b.temp~. + Condition(as.matrix(xy.dat.temp)), data=as.data.frame(dbMEM.region), distance="bray")
# R2a <- RsquareAdj(rda.region)$adj.r.squared   # R2a = [1] -0.1062698
# mod0 <- rda(bio.b.temp ~ 1 + Condition(as.matrix(xy.dat.temp)) , data=as.data.frame(dbMEM.region))
# rda.region.fwd <- ordistep(mod0, scope=formula(rda.region), direction="forward", perm.max=200)
# rda.region.fwd$call
# 
# 
# eigenvals(rda.region)[1:2]/sum(eigenvals(rda.region))*RsquareAdj(rda.region)$adj.r.squared/RsquareAdj(rda.region)$r.squared
# axes.region <- scores(rda.region, choices=c(1:2), display="lc", scaling=1)
# pdf(here("Figures", "Figure_Caucasus"), width=12, height=12)
# PCNM.region(axes.region, 1, "Axe1 (5.1%)")
# dev.off()
# 

## ****************************** ##
## SUBSET = Norway                ##
## ****************************** ##
i = "Norway"
var.data.temp <- subset(env.data, site_c == i)
xy.dat.temp <- subset(xy.dat, env.data$site_c == i)
xyz.dat.temp <- subset(xyz.dat, env.data$site_c == i)
dist_geo.temp <- distm(xy.dat.temp, fun=distGeo)
bio.b.temp <- subset(bio.b, env.data$site_c == i)
dbMEM.temp <- pcnm(dist_geo.temp)
(sel.pos <- which(dbMEM.temp$values > 0))

## Eigenfunctions with positive spatial correlation
dbMEM.region <- as.matrix(dbMEM.temp$vectors[,sel.pos])
dbMEM.region <- cbind(data.frame(sample=var.data.temp$sample), dbMEM.region)
names(dbMEM.region)[2:(length(sel.pos)+1)] <- paste(i, names(dbMEM.region)[2:(length(sel.pos)+1)], sep=".")
var.data.hier <- merge(var.data.hier, dbMEM.region, by="sample", sort=F, all.x=T, nomatch=0)

# dbMEM.region <- as.matrix(dbMEM.temp$vectors[,sel.pos])
# rda.region <- dbrda(bio.b.temp~. + Condition(as.matrix(xy.dat.temp)), data=as.data.frame(dbMEM.region), distance="bray")
# R2a <- RsquareAdj(rda.region)$adj.r.squared   # R2a = [1] -0.1062698
# mod0 <- rda(bio.b.temp ~ 1 + Condition(as.matrix(xy.dat.temp)) , data=as.data.frame(dbMEM.region))
# rda.region.fwd <- ordistep(mod0, scope=formula(rda.region), direction="forward", perm.max=200)
# rda.region.fwd$call
# 


## ****************************** ##
## SUBSET = Nepal                 ##
## ****************************** ##
i = "Nepal"
var.data.temp <- subset(env.data, site_c == i)
xy.dat.temp <- subset(xy.dat, env.data$site_c == i)
xyz.dat.temp <- subset(xyz.dat, env.data$site_c == i)
dist_geo.temp <- distm(xy.dat.temp, fun=distGeo)
bio.b.temp <- subset(bio.b, env.data$site_c == i)
dbMEM.temp <- pcnm(dist_geo.temp)
sel.pos <- which(dbMEM.temp$values > 0)

## Eigenfunctions with positive spatial correlation
dbMEM.region <- as.matrix(dbMEM.temp$vectors[,sel.pos])
dbMEM.region <- cbind(data.frame(sample=var.data.temp$sample), dbMEM.region)
names(dbMEM.region)[2:(length(sel.pos)+1)] <- paste(i, names(dbMEM.region)[2:(length(sel.pos)+1)], sep=".")
var.data.hier <- merge(var.data.hier, dbMEM.region, by="sample", sort=F, all.x=T, nomatch=0)

# dbMEM.region <- as.matrix(dbMEM.temp$vectors[,sel.pos])
# rda.region <- dbrda(bio.b.temp~. + Condition(as.matrix(xy.dat.temp)), data=as.data.frame(dbMEM.region), distance="bray")
# R2a <- RsquareAdj(rda.region)$adj.r.squared  
# mod0 <- rda(bio.b.temp ~ 1 + Condition(as.matrix(xy.dat.temp)) , data=as.data.frame(dbMEM.region))
# rda.region.fwd <- ordistep(mod0, scope=formula(rda.region), direction="forward", perm.max=200)
# rda.region.fwd$call


## ****************************** ##
## SUBSET = Pamir                 ##
## ****************************** ##
i = "Kirghizistan"
var.data.temp <- subset(env.data, site_c == i)
xy.dat.temp <- subset(xy.dat, env.data$site_c == i)
xyz.dat.temp <- subset(xyz.dat, env.data$site_c == i)
dist_geo.temp <- distm(xy.dat.temp, fun=distGeo)
bio.b.temp <- subset(bio.b, env.data$site_c == i)
dbMEM.temp <- pcnm(dist_geo.temp)
sel.pos <- which(dbMEM.temp$values > 0)

## Eigenfunctions with positive spatial correlation
dbMEM.region <- as.matrix(dbMEM.temp$vectors[,sel.pos])
dbMEM.region <- cbind(data.frame(sample=var.data.temp$sample), dbMEM.region)
names(dbMEM.region)[2:(length(sel.pos)+1)] <- paste(i, names(dbMEM.region)[2:(length(sel.pos)+1)], sep=".")
var.data.hier <- merge(var.data.hier, dbMEM.region, by="sample", sort=F, all.x=T, nomatch=0)

# dbMEM.region <- as.matrix(dbMEM.temp$vectors[,sel.pos])
# rda.region <- dbrda(bio.b.temp~. + Condition(as.matrix(xy.dat.temp)), data=as.data.frame(dbMEM.region), distance="bray")
# R2a <- RsquareAdj(rda.region)$adj.r.squared 
# mod0 <- rda(bio.b.temp ~ 1 + Condition(as.matrix(xy.dat.temp)) , data=as.data.frame(dbMEM.region))
# rda.region.fwd <- ordistep(mod0, scope=formula(rda.region), direction="forward", perm.max=200)
# rda.region.fwd$call
# 


## ****************************** ##
## SUBSET = Ecuador               ##
## ****************************** ##
i = "Ecuador"
var.data.temp <- subset(env.data, site_c == i)
xy.dat.temp <- subset(xy.dat, env.data$site_c == i)
xyz.dat.temp <- subset(xyz.dat, env.data$site_c == i)
dist_geo.temp <- distm(xy.dat.temp, fun=distGeo)
bio.b.temp <- subset(bio.b, env.data$site_c == i)
dbMEM.temp <- pcnm(dist_geo.temp)
sel.pos <- which(dbMEM.temp$values > 0)

## Eigenfunctions with positive spatial correlation
dbMEM.region <- as.matrix(dbMEM.temp$vectors[,sel.pos])
dbMEM.region <- cbind(data.frame(sample=var.data.temp$sample), dbMEM.region)
names(dbMEM.region)[2:(length(sel.pos)+1)] <- paste(i, names(dbMEM.region)[2:(length(sel.pos)+1)], sep=".")
var.data.hier <- merge(var.data.hier, dbMEM.region, by="sample", sort=F, all.x=T, nomatch=0)


# dbMEM.region <- as.matrix(dbMEM.temp$vectors[,sel.pos])
# rda.region <- dbrda(bio.b.temp~. + Condition(as.matrix(xy.dat.temp)), data=as.data.frame(dbMEM.region), distance="bray")
# R2a <- RsquareAdj(rda.region)$adj.r.squared 
# mod0 <- rda(bio.b.temp ~ 1 + Condition(as.matrix(xy.dat.temp)) , data=as.data.frame(dbMEM.region))
# rda.region.fwd <- ordistep(mod0, scope=formula(rda.region), direction="forward", perm.max=200)
# rda.region.fwd$call


## ****************************** ##
## SUBSET = New_Zealand           ##
## ****************************** ##
i = "New_Zealand"
var.data.temp <- subset(env.data, site_c == i)
xy.dat.temp <- subset(xy.dat, env.data$site_c == i)
xyz.dat.temp <- subset(xyz.dat, env.data$site_c == i)
dist_geo.temp <- distm(xy.dat.temp, fun=distGeo)
bio.b.temp <- subset(bio.b, env.data$site_c == i)
dbMEM.temp <- pcnm(dist_geo.temp)
(sel.pos <- which(dbMEM.temp$values > 0))

## Eigenfunctions with positive spatial correlation
dbMEM.region <- as.matrix(dbMEM.temp$vectors[,sel.pos])
dbMEM.region <- cbind(data.frame(sample=var.data.temp$sample), dbMEM.region)
names(dbMEM.region)[2:(length(sel.pos)+1)] <- paste(i, names(dbMEM.region)[2:(length(sel.pos)+1)], sep=".")
var.data.hier <- merge(var.data.hier, dbMEM.region, by="sample", sort=F, all.x=T, nomatch=0)

# dbMEM.region <- as.matrix(dbMEM.temp$vectors[,sel.pos])
# rda.region <- dbrda(bio.b.temp~. + Condition(as.matrix(xy.dat.temp)), data=as.data.frame(dbMEM.region), distance="bray")
# R2a <- RsquareAdj(rda.region)$adj.r.squared 
# mod0 <- rda(bio.b.temp ~ 1 + Condition(as.matrix(xy.dat.temp)) , data=as.data.frame(dbMEM.region))
# rda.region.fwd <- ordistep(mod0, scope=formula(rda.region), direction="forward", perm.max=200)
# rda.region.fwd$call


## ****************************** ##
## SUBSET = Alaska                ##
## ****************************** ##
i = "Alaska"
var.data.temp <- subset(env.data, site_c == i)
xy.dat.temp <- subset(xy.dat, env.data$site_c == i)
xyz.dat.temp <- subset(xyz.dat, env.data$site_c == i)
dist_geo.temp <- distm(xy.dat.temp, fun=distGeo)
bio.b.temp <- subset(bio.b, env.data$site_c == i)
dbMEM.temp <- pcnm(dist_geo.temp)
(sel.pos <- which(dbMEM.temp$values > 0))

## Eigenfunctions with positive spatial correlation
dbMEM.region <- as.matrix(dbMEM.temp$vectors[,sel.pos])
dbMEM.region <- cbind(data.frame(sample=var.data.temp$sample), dbMEM.region)
names(dbMEM.region)[2:(length(sel.pos)+1)] <- paste(i, names(dbMEM.region)[2:(length(sel.pos)+1)], sep=".")
var.data.hier <- merge(var.data.hier, dbMEM.region, by="sample", sort=F, all.x=T, nomatch=0)

# dbMEM.region <- as.matrix(dbMEM.temp$vectors[,sel.pos])
# rda.region <- dbrda(bio.b.temp~. + Condition(as.matrix(xy.dat.temp)), data=as.data.frame(dbMEM.region), distance="bray")
# R2a <- RsquareAdj(rda.region)$adj.r.squared 
# mod0 <- rda(bio.b.temp ~ 1 + Condition(as.matrix(xy.dat.temp)) , data=as.data.frame(dbMEM.region))
# rda.region.fwd <- ordistep(mod0, scope=formula(rda.region), direction="forward", perm.max=200)
# rda.region.fwd$call


## ****************************** ##
## Editing variables              ##
## ****************************** ##
# Replace all NA in db-MEM by 0
var.data.hier[is.na(var.data.hier)] <- 0
# reorder to match community matrix
var.data.hier$sample == env.data$sample
var.data.hier2 <- merge(env.data, var.data.hier[,!(colnames(var.data.hier)%in%env.names)], by=c("sample","site_c"), sort=F, all.x=T, nomatch=0)
var.data.hier2$sample == env.data$sample
# create a region dummy variable (N-1, with N the number of regions)
region.dum <- as.matrix(model.matrix(~-1+var.data.hier2$site_c))
region.dum <- region.dum[,1:(dim(region.dum)[2]-1)]
# db-MEM only dataset  
region.dbMEM <- var.data.hier2[,16:dim(var.data.hier2)[2]]

## 3) Models ##
# 3a) variation among regions at the global scale
# Spatial model
rda1.S <-dbrda(bio.b ~., data=as.data.frame(region.dum), distance="bray")

R2a <- RsquareAdj(rda1.S)$adj.r.squared #[1] 0.3045347

## Environmental model
rda1.E <- dbrda(bio.b~., data=env.mat, distance="bray")
vif.cca(rda1.E)

eigenvals(rda1.E)[1:2]/sum(eigenvals(rda1.E))*RsquareAdj(rda1.E)$adj.r.squared/RsquareAdj(rda1.E)$r.squared

R2a <- RsquareAdj(rda1.E)$adj.r.squared # [1] 0.2061079
## Compute forward selection
mod0 <- dbrda(bio.b ~ 1, data=env.mat, distance="bray")
rda1.E.fwd <- ordistep(mod0, scope=formula(rda1.E), direction="forward", perm.max=200)
rda1.E.fwd$call

rda2 = dbrda(formula = bio.b ~ pH + pr + DIN + turb + cond + mean_chla + 
               water_temp + scd + Calcite + srp + Quartz, data = env.mat, 
             distance = "bray")

anova(rda2, step=1000)
anova(rda2, by="terms")

anova(rda2, by="axis")

R2a <- RsquareAdj(rda2)$adj.r.squared 
eigenvals(rda2)[1:2]/sum(eigenvals(rda2))*RsquareAdj(rda2)$adj.r.squared/RsquareAdj(rda2)$r.squared

## Subset environmental variables to enter it in the variance partitioning
env.mat.adj <- subset(env.mat, select = -c(doc,gl_cov.x))

## Var Part!
vp1 <- varpart(vegdist(bio.b), xyz.dat, region.dum, region.dbMEM, env.mat.adj)
vp1;plot(vp1)

## Test of individual fractions by mean of RDA
rda1.Sx <- dbrda(bio.b ~. +Condition(as.matrix(xyz.dat)), data=as.data.frame(region.dum), distance="bray")
anova(rda1.Sx, step=1000)
R2a <- RsquareAdj(rda1.Sx)$adj.r.squared 

rda1.Ex <- dbrda(bio.b ~. + Condition(as.matrix(region.dum)) + Condition(as.matrix(region.dbMEM)) +Condition(as.matrix(xyz.dat)), data=as.data.frame(env.mat.adj),distance="bray")
anova(rda1.Ex, step=1000)
(R2a <- RsquareAdj(rda1.Ex)$adj.r.squared) 

R2a <- RsquareAdj(rda2.S)$adj.r.squared 


