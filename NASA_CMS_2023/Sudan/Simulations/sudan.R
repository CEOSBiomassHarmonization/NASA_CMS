rm(list=ls())
###############
## libraries ##
library(terra)   ## to read in spatial data
library(spBayes) ## running Bayesian hierachical spatial model
library(geoR)    ## for semi-variogram fitting

##################
## read in data ##
dat <- vect("/projects/my-public-bucket/Data/NASA_CMS_2023/SUDAN/GMB_Test/Sudan_NFI.gpkg")
strata <- vect("Strata/Sudan_Strata_Layer.shp")

## extract strata at plot locations
dat$Id <- extract(strata, dat)$Id

## reproject data to UTM (can't use lat/long)
dat <- project(dat, "+proj=utm +zone=35 +a=6378249.145 +rf=293.465 +towgs84=-161,-14,205,0,0,0,0 +units=m +no_defs +type=crs")
strata <- project(strata, "+proj=utm +zone=35 +a=6378249.145 +rf=293.465 +towgs84=-161,-14,205,0,0,0,0 +units=m +no_defs +type=crs")

## drop strata that are too sparse
dat <- dat[dat$Id == 3 | dat$Id == 4,]
strata <- strata[strata$Id == 3 | strata$Id == 4,]

## check if reprojection and stratum Id extraction worked out
plot(strata,"Id")
points(dat)

## set prediction locations in strata (have to do this since no covariates are involved here)
pred.locs <- spatSample(strata, 2500, method = "regular")
coords.pred <- crds(pred.locs)/1000 ## convert coords to km. It helps for fitting phi.
x.pred <- pred.locs$Id - 3

## set prediction locations in strata for mapping (more intense grid just for mapping)
pred.locs.map <- spatSample(strata, 50000, method = "regular")
coords.pred.map <- crds(pred.locs.map)/1000 ## convert coords to km. It helps for fitting phi.
x.pred.map <- pred.locs.map$Id - 3


## set up y and coords
y <- dat$`AG_Carbon (Mg/ha)`
x <- dat$Id - 3
coords <- crds(dat)/1000

mod <- lm(y ~ x)

mean(y[x == 0])
coef(mod)[1]

mean(y[x == 1])
coef(mod)[1] + coef(mod)[2]

vario <- variog(data = resid(mod), coords = coords, max.dist = 1000)

plot(vario)

#######################
## fit spatial model ##
priors <- list("beta.Norm"=list(rep(0,2), diag(1000,2)),
               "phi.Unif"=list(3/1300, 3/1),
               "sigma.sq.IG"=list(2, 13),
               "tau.sq.IG"=c(2, 7))

starting <- list("phi"=3/750, "sigma.sq"=13, "tau.sq"=7)

tuning <- list("phi"=0.0175, "sigma.sq"=0.0175, "tau.sq"=0.0175)

sp.mod <- spSVC(y ~ x, coords = coords, starting = starting, tuning = tuning, 
             priors = priors, cov.model = "exponential", 
             n.omp.threads = 50, n.samples = 10000)

plot(sp.mod$p.theta.samples)

sp.mod <- spRecover(sp.mod, start = 2501, thin = 5, n.omp.threads = 50)

save(sp.mod, file = "sp.mod.RData")

#######################################
## predict for                       ##
## areal estimation (joint sampling) ##
sp.pred <- spPredict(sp.mod, pred.coords = coords.pred, pred.covars = cbind(1,x.pred), 
                     joint = TRUE, n.omp.threads = 50)

save(sp.pred, file = "sp.pred.RData")

###############################################
## predict for mapping (point-wise sampling) ##
sp.pred.map <- spPredict(sp.mod, pred.coords = coords.pred.map, pred.covars = cbind(1,x.pred.map), 
                         joint = FALSE, n.omp.threads = 50)

save(sp.pred.map, file = "sp.pred.map.RData")

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
