rm(list=ls())
###############
## libraries ##
library(terra)
library(spBayes)

##################
## read in data ##
dat <- vect("/projects/my-private-bucket/Data/NFI_data/Sudan/Sudan_NFI.gpkg")
strata <- vect("my-public-bucket/Data/NASA_CMS_2023/SUDAN/GEDI_Strata/Sudan_country.gpkg")

## extract strata at plot locations
dat$Id <- extract(strata, dat)$Id

## reproject data to UTM (can't use lat/long)
dat <- project(dat, "+proj=utm +zone=35 +a=6378249.145 +rf=293.465 +towgs84=-161,-14,205,0,0,0,0 +units=m +no_defs +type=crs")
strata <- project(strata, "+proj=utm +zone=35 +a=6378249.145 +rf=293.465 +towgs84=-161,-14,205,0,0,0,0 +units=m +no_defs +type=crs")

## drop strata that are too sparse
dat <- dat[dat$Id == 3 | dat$Id == 4,]
strata <- strata[strata$Id == 3 | strata$Id == 4,]


## isolate AG_Carbon plot values
y <- dat$`AG_Carbon (Mg/ha)`

y.3 <- y[dat$Id == 3]

y.4 <- y[dat$Id == 4]


################################
## get design-based estimates ##

## y
y.bar.hat <- mean(y)

s2 <- sum((y - y.bar.hat)^2)/length(y)

y.se.hat <- sqrt(s2/length(y))

y.db.ci <- y.bar.hat + c(-1,1)*1.96*y.se.hat

## y.3
y.3.bar.hat <- mean(y.3)

s2 <- sum((y.3 - y.3.bar.hat)^2)/length(y.3)

y.3.se.hat <- sqrt(s2/length(y.3))

y.3.db.ci <- y.3.bar.hat + c(-1,1)*1.96*y.3.se.hat

## y.4
y.4.bar.hat <- mean(y.4)

s2 <- sum((y.4 - y.4.bar.hat)^2)/length(y.4)

y.4.se.hat <- sqrt(s2/length(y.4))

y.4.db.ci <- y.4.bar.hat + c(-1,1)*1.96*y.4.se.hat

#######################
## get gmb estimates ##

## y
load("sp.pred.RData")

y.ppd <- sp.pred$p.y.predictive.samples

y.gmb.ppd <- colMeans(y.ppd)

y.gmb.mn <- mean(y.gmb.ppd)
y.gmb.sd <- sd(y.gmb.ppd)
y.gmb.ci <- quantile(y.gmb.ppd, p = c(.025,.975))

## y.3
load("sp.pred.3.RData")

y.3.ppd <- sp.pred$p.y.predictive.samples

y.3.gmb.ppd <- colMeans(y.3.ppd)

y.3.gmb.mn <- mean(y.3.gmb.ppd)
y.3.gmb.sd <- sd(y.3.gmb.ppd)
y.3.gmb.ci <- quantile(y.3.gmb.ppd, p = c(.025,.975))

## y.4
load("sp.pred.4.RData")

y.4.ppd <- sp.pred$p.y.predictive.samples

y.4.gmb.ppd <- colMeans(y.4.ppd)

y.4.gmb.mn <- mean(y.4.gmb.ppd)
y.4.gmb.sd <- sd(y.4.gmb.ppd)
y.4.gmb.ci <- quantile(y.4.gmb.ppd, p = c(.025,.975))



y.results <- as.data.frame(rbind(c(y.bar.hat,y.se.hat,y.db.ci),c(y.gmb.mn,y.gmb.sd,y.gmb.ci)))
rownames(y.results) <- c("DB", "GMB")
names(y.results) <- c("Estimate", "Standard Error","95% lower bound", "95% upper bound")

y.results

y.3.results <- as.data.frame(rbind(c(y.3.bar.hat,y.3.se.hat,y.3.db.ci),c(y.3.gmb.mn,y.3.gmb.sd,y.3.gmb.ci)))
rownames(y.3.results) <- c("DB", "GMB")
names(y.3.results) <- c("Estimate", "Standard Error","95% lower bound", "95% upper bound")

y.3.results

y.4.results <- as.data.frame(rbind(c(y.4.bar.hat,y.4.se.hat,y.4.db.ci),c(y.4.gmb.mn,y.4.gmb.sd,y.4.gmb.ci)))
rownames(y.4.results) <- c("DB", "GMB")
names(y.4.results) <- c("Estimate", "Standard Error","95% lower bound", "95% upper bound")

y.4.results
