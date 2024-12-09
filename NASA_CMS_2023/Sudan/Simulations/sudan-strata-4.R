rm(list=ls())
###############
## libraries ##
library(terra)   ## to read in spatial data and reproject
library(spBayes) ## running Bayesian hierachical spatial model
library(geoR)    ## for semi-variogram fitting

## read in spatial data
dat <- vect("Sudan_NFI.gpkg")
strata <- vect("Strata/Sudan_Strata_Layer.shp")

## extract strata at plot locations
dat$Id <- extract(strata, dat)$Id

## reproject data to UTM (can't use lat/long)
dat <- project(dat, "+proj=utm +zone=35 +a=6378249.145 +rf=293.465 +towgs84=-161,-14,205,0,0,0,0 +units=m +no_defs +type=crs")
strata <- project(strata, "+proj=utm +zone=35 +a=6378249.145 +rf=293.465 +towgs84=-161,-14,205,0,0,0,0 +units=m +no_defs +type=crs")

## isolate strata 4
dat <- dat[dat$Id == 4,]
strata <- strata[strata$Id == 4,]

## set prediction locations in strata (have to do this since no covariates are involved here)
pred.locs <- spatSample(strata, 2500, method = "regular")
coords.pred <- crds(pred.locs)/1000 ## convert coords to km. It helps for fitting phi.

## set prediction locations in strata for mapping (more intense grid just for mapping)
pred.locs.map <- spatSample(strata, 50000, method = "regular")
coords.pred.map <- crds(pred.locs.map)/1000 ## convert coords to km. It helps for fitting phi.

## set up y and coords
y <- dat$`AG_Carbon (Mg/ha)`
coords <- crds(dat)/1000 ## convert coords to km. It helps for fitting phi.

## fit non-spatial model to check residual spatial autocorrelation
mod <- lm(y ~ 1)

vario <- variog(data = resid(mod), coords = coords, max.dist = 800)

plot(vario)

#######################
## fit spatial model ##
priors <- list("beta.Norm"=list(rep(0,1), diag(1000,1)),
               "phi.Unif"=list(3/1400, 3/20),
               "sigma.sq.IG"=list(2, 13),
               "tau.sq.IG"=c(2, 10))

starting <- list("phi"=3/600, "sigma.sq"=13, "tau.sq"=10)

tuning <- list("phi"=0.02, "sigma.sq"=0.02, "tau.sq"=0.02)

sp.mod <- spSVC(y ~ 1, coords = coords, starting = starting, tuning = tuning, 
             priors = priors, cov.model = "exponential", 
             n.omp.threads = 50, n.samples = 10000)

plot(sp.mod$p.theta.samples)

sp.mod <- spRecover(sp.mod, start = 2501, thin = 5, n.omp.threads = 50)

save(sp.mod, file = "sp.mod.4.RData")

#######################################
## predict for                       ##
## areal estimation (joint sampling) ##
sp.pred <- spPredict(sp.mod, pred.coords = coords.pred,
                     pred.covars = matrix(rep(1, nrow(coords.pred)),ncol = 1), 
                     joint = TRUE, n.omp.threads = 50)

save(sp.pred, file = "sp.pred.4.RData")

###############################################
## predict for mapping (point-wise sampling) ##
sp.pred.map <- spPredict(sp.mod, pred.coords = coords.pred.map,
                         pred.covars = matrix(rep(1, nrow(coords.pred.map)),ncol = 1), 
                         joint = FALSE, n.omp.threads = 50)

save(sp.pred.map, file = "sp.pred.map.4.RData")

###################################################
## generate strata estimates (block predictions) ##
y.ppd <- sp.pred$p.y.predictive.samples

y.gmb.ppd <- colMeans(y.ppd)

y.gmb.mn <- mean(y.gmb.ppd)
y.gmb.sd <- sd(y.gmb.ppd)
y.gmb.ci <- quantile(y.gmb.ppd, p = c(.025,.975))

#################
## db estimate ##
y.bar.hat <- mean(y)

s2 <- sum((y - y.bar.hat)^2)/length(y)

y.se.hat <- sqrt(s2/length(y))

y.db.ci <- y.bar.hat + c(-1,1)*1.96*y.se.hat

##################################
## compare DB and GMB estimates ##
y.gmb.mn
y.bar.hat

y.gmb.sd
y.se.hat

y.gmb.ci
y.db.ci

## make maps of carbon and associated uncertainty ##
y.ppd.map <- sp.pred.map$p.y.predictive.samples

y.pred.map <- rowMeans(y.ppd.map)
y.pred.sd.map <- apply(y.ppd.map,1,sd)

y.pred.rast <- rast(data.frame(cbind(coords.pred.map*1000,y.pred.map)))
y.pred.sd.rast <- rast(data.frame(cbind(coords.pred.map*1000,y.pred.sd.map)))

plot(y.pred.rast)
plot(y.pred.sd.rast)


