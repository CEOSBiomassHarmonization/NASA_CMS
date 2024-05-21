## ----label = "settings", include = FALSE, results = 'hide', message = FALSE, warning = FALSE----
# Used to be called {#ch:ngns}
source("R/initial_setup.R") 
opts_chunk$set(
  fig.path = 'figs/rain-'
)
library(lattice) 
library(gridExtra) 

## ----label = "start", results = 'hide'-----------------------------------
data(PRprec) 

## ----label = "headdat"---------------------------------------------------
ii <- c(537, 304, 610, 388)
PRprec[ii, 1:8]

## ----label = "meanjan"---------------------------------------------------
PRprec$precMean <- rowMeans(PRprec[, 3 + 1:31], na.rm = TRUE)
summary(PRprec$precMean)

## ------------------------------------------------------------------------
table(rowSums(is.na(PRprec[, 3 + 1:31])))

## ----label = "prborder"--------------------------------------------------
data(PRborder) 

## ----label = "paranastations", echo = FALSE, fig = TRUE, fig.align = "center", fig.width = 10, fig.height = 6, fig.cap = "Locations of Paraná stations, altitude and average of daily accumulated precipitation (mm) in January 2011. Circles in red denote stations with missing observations and red crossess denote the four stations shown in the example in the main text."----

par(mfrow = c(1, 2), mar = c(0, 0, 2, 0))
 
# Paraná state boundary
plot(PRborder, type = 'l', asp = 1, axes = FALSE, main = 'Altitude') 
# Locations of the stations: circle radius prop. to altitude
points(PRprec[1:2], col = is.na(PRprec$Altitude) + 1, 
  cex = ifelse(is.na(PRprec$Altitude), 1, 0.3 + PRprec$Altitude / 1500)) 
# Atlantic Ocean
lines(PRborder[1034:1078, ], col = 'cyan') 
# Legend
legend('topright', format(0:4 * 350), bty = 'n', pch = 1, 
  pt.cex = 0.3 + 0:4 * 35 / 150) 

# Paraná state boundary
plot(PRborder, type = 'l', asp = 1, axes = FALSE, 
  main = paste('Mean of daily accumulated precipitation (mm)')) 
# Locations of the stations: circle radius prop. to mean precip.
points(PRprec[1:2], cex = 0.3 + PRprec$precMean / 20)  
# Locations of  4 stations in the example
points(PRprec[ii, 1:2], pch = 3, col = 2, cex = 2) 
# Atlantic Ocean
lines(PRborder[1034:1078, ], col = 'cyan') 
# Legend
legend('topright', format(seq(1, 21, 5)), 
  bty = 'n', pch = 1, pt.cex = 0.3 + seq(1, 21, 5) / 20) 

## ----label = "project"---------------------------------------------------
coords <- as.matrix(PRprec[, 1:2]) 
mat.dists <- spDists(coords, PRborder[1034:1078, ],
  longlat = TRUE) 

## ----label = "distoceancalc"---------------------------------------------
PRprec$oceanDist <- apply(mat.dists, 1, min) 

## ----label = "disppred", echo = FALSE, fig = TRUE, fig.align = "center", fig.width = 7.5, fig.height = 5, out.width='90%', fig.cap = "Dispersion plots of average of daily accumulated precipitation by longitude (top left), latitude (top right), altitude (bottom left) and distance to ocean (bottom right)."----

par(mfrow = c(2, 2), mar = c(3, 3, 0.5, 0.5), mgp = c(1.7, 0.7, 0),
  las = 1) 
for (i in c(1:3, ncol(PRprec)))
  plot(PRprec[c(i, ncol(PRprec) - 1)], cex = 0.5) 

## ----label = "pcprec"----------------------------------------------------
pcprec <- list(prior = 'pcprec', param = c(1, 0.01))

## ----label = "prmesh"----------------------------------------------------
pts.bound <- inla.nonconvex.hull(coords, 0.3, 0.3)
mesh <- inla.mesh.2d(coords, boundary = pts.bound, 
  max.edge = c(0.3, 1), offset = c(1e-5, 1.5), cutoff = 0.1)

## ----label = "projA"-----------------------------------------------------
A <- inla.spde.make.A(mesh, loc = coords)

## ----label = "spde"------------------------------------------------------
spde <- inla.spde2.pcmatern(mesh = mesh,
  prior.range = c(0.05, 0.01), # P(practic.range < 0.05) = 0.01
  prior.sigma = c(1, 0.01)) # P(sigma > 1) = 0.01

## ----label = "stackpr"---------------------------------------------------
stk.dat <- inla.stack(
  data = list(y = PRprec$precMean), 
  A = list(A,1),
  effects = list(list(s = 1:spde$n.spde), 
    data.frame(Intercept = 1, 
      gWest = inla.group(coords[, 1]), 
      gOceanDist = inla.group(PRprec$oceanDist), 
      oceanDist = PRprec$oceanDist)),
  tag = 'dat') 

## ----label = "fitwestcoo"------------------------------------------------
f.west <- y ~ 0 + Intercept + # f is short for formula
  f(gWest, model = 'rw1', # first random walk prior 
    scale.model = TRUE, # scaling this random effect
    hyper = list(theta = pcprec)) + # use the PC prior
  f(s, model = spde)

r.west <- inla(f.west, family = 'Gamma', # r is short for result
  control.compute = list(cpo = TRUE),
  data = inla.stack.data(stk.dat), 
  control.predictor = list(A = inla.stack.A(stk.dat), link = 1))

## ------------------------------------------------------------------------
inla.models()$likelihood$gamma$link

## ----label = "fitdistocean"----------------------------------------------
f.oceanD <- y ~ 0 + Intercept + 
  f(gOceanDist, model = 'rw1', scale.model = TRUE, 
    hyper = list(theta = pcprec)) + 
  f(s, model = spde)

r.oceanD <- inla(f.oceanD, family = 'Gamma', 
  control.compute = list(cpo = TRUE),
  data = inla.stack.data(stk.dat), 
  control.predictor = list(A = inla.stack.A(stk.dat), link = 1))

## ----label = "fitdistoceanl"---------------------------------------------
f.oceanD.l <- y ~ 0 + Intercept + oceanDist + 
  f(s, model = spde)
r.oceanD.l <- inla(f.oceanD.l, family = 'Gamma', 
  control.compute = list(cpo = TRUE),
  data = inla.stack.data(stk.dat), 
  control.predictor = list(A = inla.stack.A(stk.dat), link = 1))

## ----label = "slcpo-function"--------------------------------------------
slcpo <- function(m, na.rm = TRUE) {
  - sum(log(m$cpo$cpo), na.rm = na.rm)
}

## ----label = "cpos"------------------------------------------------------
c(long = slcpo(r.west), oceanD = slcpo(r.oceanD),
  oceanD.l = slcpo(r.oceanD.l))

## ----label = "betasumary1"-----------------------------------------------
round(r.oceanD.l$summary.fixed, 2)

## ----label = "disp"------------------------------------------------------
round(r.oceanD.l$summary.hyperpar[1, ], 3)

## ----label = "resfield"--------------------------------------------------
round(r.oceanD.l$summary.hyperpar[-1, ], 3)

## ----label = "doceaneff", echo = FALSE, fig = TRUE, fig.align = "center", out.width = "99%", fig.width = 7.5, fig.height = 4, fig.cap = "Posterior marginal distributions for $\\beta_0$ (top left), the distance to ocean coefficient (top center), the Gamma likelihood precision (bottom left), the practical range (bottom center) and the standard deviation of the spatial field (bottom right). The top-right plot represents the posterior mean (continuous line) and 95\\% credibility interval (dashed lines) for the distance to ocean effect."----

#PMD for $\\beta_0$ (top left), PMD for the distance to ocean coefficient (top mid), posterior mean (continuous line) and 95% credibility interval (dashed lines) for the distance to ocean effect (top right), PMD for the Gamma likelihood precision (bottom left), PMD for the practical range (bottom mid) and PMD for the standard deviation of the spatial field (bottom right)."}

par(mfrow = c(2, 3), mar = c(3, 3.5, 0, 0), mgp = c(1.5, 0.5, 0), las = 0) 

plot(r.oceanD.l$marginals.fix[[1]], type = 'l', 
  xlab = 'Intercept', ylab = 'Density') 
plot(r.oceanD.l$marginals.fix[[2]], type = 'l', 
  xlab = 'Ocean distance coefficient', ylab = 'Density') 
plot(r.oceanD$summary.random[[1]][,1:2], type = 'l', 
  xlab = 'Distance to ocean (Km)', ylab = 'Effect') 
abline(h = 0, lty = 3) 

for (i in c(4, 6)) 
  lines(r.oceanD$summary.random[[1]][, c(1, i)], lty = 2) 
abline(h = 0, lty = 3)

#Fix name for plotting
names(r.oceanD.l$marginals.hyperpar)[1] <- 
  "Precision for the Gamma observations"
for (j in 1:3) 
  plot(r.oceanD.l$marginals.hy[[j]], type = 'l', 
    ylab = 'Density', xlab = names(r.oceanD.l$marginals.hyperpar)[j])

## ----label = "mod0"------------------------------------------------------
r0.oceanD.l <- inla(y ~ 0 + Intercept + oceanDist, 
  family = 'Gamma', control.compute = list(cpo = TRUE),
  data = inla.stack.data(stk.dat), 
  control.predictor = list(A = inla.stack.A(stk.dat), link = 1))

## ----label = "cpo0"------------------------------------------------------
c(oceanD.l = slcpo(r.oceanD.l), oceanD.l0 =  slcpo(r0.oceanD.l))

## ----label = "stepsize"--------------------------------------------------
stepsize <- 4 * 1 / 111

## ----label = "ncoords"---------------------------------------------------
x.range <- diff(range(PRborder[, 1]))
y.range <- diff(range(PRborder[, 2]))
nxy <- round(c(x.range, y.range) / stepsize)

## ------------------------------------------------------------------------
nxy

## ----label = "projgrid"--------------------------------------------------
projgrid <- inla.mesh.projector(mesh, xlim = range(PRborder[, 1]), 
  ylim = range(PRborder[, 2]), dims = nxy)

## ----label = "projpred"--------------------------------------------------
xmean <- inla.mesh.project(projgrid,
  r.oceanD$summary.random$s$mean)
xsd <- inla.mesh.project(projgrid, r.oceanD$summary.random$s$sd)

## ----label = "sp"--------------------------------------------------------
library(splancs)
xy.in <- inout(projgrid$lattice$loc, PRborder)

## ------------------------------------------------------------------------
table(xy.in)

## ------------------------------------------------------------------------
xmean[!xy.in] <- NA
xsd[!xy.in] <- NA

## ------------------------------------------------------------------------
Aprd <- projgrid$proj$A[which(xy.in), ]

## ----label = "prdcoo"----------------------------------------------------
prdcoo <- projgrid$lattice$loc[which(xy.in), ]

## ----label = "oceanDcov0"------------------------------------------------
oceanDist0 <- apply(spDists(PRborder[1034:1078, ], 
  prdcoo, longlat = TRUE), 2, min)

## ----label = "oceanDistk"------------------------------------------------
OcDist.k <- sort(unique(stk.dat$effects$data$gOceanDist))

## ----label = "goceanD"---------------------------------------------------
oceanDist.b <- (OcDist.k[-1] + OcDist.k[length(OcDist.k)]) / 2

## ------------------------------------------------------------------------
i0 <- findInterval(oceanDist0, oceanDist.b) + 1
gOceanDist0 <- OcDist.k[i0]

## ----label = "stkprd"----------------------------------------------------
stk.prd <- inla.stack(
  data = list(y = NA),
  A = list(Aprd, 1), 
  effects = list(s = 1:spde$n.spde, 
    data.frame(Intercept = 1, oceanDist = oceanDist0)),
  tag = 'prd') 
stk.all <- inla.stack(stk.dat, stk.prd)

## ----label = "predfit"---------------------------------------------------
r2.oceanD.l <- inla(f.oceanD.l, family = 'Gamma', 
  data = inla.stack.data(stk.all), 
  control.predictor = list(A = inla.stack.A(stk.all),
    compute = TRUE, link = 1),
  quantiles = NULL, 
  control.inla = list(strategy = 'adaptive'), 
  control.results = list(return.marginals.random = FALSE,
    return.marginals.predictor = FALSE), 
  control.mode = list(theta = r.oceanD.l$mode$theta,
    restart = FALSE))

## ----label = "rprevprep"-------------------------------------------------
id.prd <- inla.stack.index(stk.all, 'prd')$data
sd.prd <- m.prd <- matrix(NA, nxy[1], nxy[2])
m.prd[xy.in] <- r2.oceanD.l$summary.fitted.values$mean[id.prd]
sd.prd[xy.in] <- r2.oceanD.l$summary.fitted.values$sd[id.prd]

## ----label = "xrain1", echo = FALSE, fig.width=12, fig.cap = "Posterior mean and standard deviation of the random field (top left and top right, respectively). Posterior mean and standard deviation for the response (bottom left and bottom right, respectively)."----
par(mfrow = c(2, 2), mar = c(0, 0, 0, 4))
book.plot.field(list(x = projgrid$x, y = projgrid$y, z = xmean))
book.plot.field(list(x = projgrid$x, y = projgrid$y, z = xsd),
  col = book.color.c2())
book.plot.field(list(x = projgrid$x, y = projgrid$y, z = m.prd))
book.plot.field(list(x = projgrid$x, y = projgrid$y, z = sd.prd),
  col = book.color.c2())

## ----label = "dsmesh"----------------------------------------------------
oceanDist.mesh <- apply(
  spDists(PRborder[1034:1078, ], mesh$loc[, 1:2], longlat = TRUE),
  2, min)

## ----label = "stkmesh"---------------------------------------------------
stk.mesh <- inla.stack(
  data = list(y = NA),
  A = list(1, 1),  
  effects = list(s = 1:spde$n.spde, 
    data.frame(Intercept = 1, oceanDist = oceanDist.mesh)),
  tag = 'mesh') 

stk.b <- inla.stack(stk.dat, stk.mesh)

## ----label = "fittmesh"--------------------------------------------------
rm.oceanD.l <- inla(f.oceanD.l, family = 'Gamma', 
  data = inla.stack.data(stk.b), 
  control.predictor = list(A = inla.stack.A(stk.b),
    compute = TRUE, link = 1), 
  quantiles = NULL, 
  control.results = list(return.marginals.random = FALSE,
    return.marginals.predictor = FALSE), 
  control.compute = list(config = TRUE)) # Needed to sample

## ----label = "sampl"-----------------------------------------------------
sampl <- inla.posterior.sample(n = 1000, result = rm.oceanD.l)

## ----label = "idexs"-----------------------------------------------------
id.prd.mesh <- inla.stack.index(stk.b, 'mesh')$data
pred.nodes <- exp(sapply(sampl, function(x) 
  x$latent[id.prd.mesh])) 

## ------------------------------------------------------------------------
dim(pred.nodes)

## ----label = "projsamples"-----------------------------------------------
sd.prd.s <- matrix(NA, nxy[1], nxy[2])
m.prd.s <- matrix(NA, nxy[1], nxy[2])

m.prd.s[xy.in] <- drop(Aprd %*% rowMeans(pred.nodes))
sd.prd.s[xy.in] <- drop(Aprd %*% apply(pred.nodes, 1, sd))

## ----label = "comparepreds"----------------------------------------------
cor(as.vector(m.prd.s), as.vector(m.prd), use = 'p')
cor(log(as.vector(sd.prd.s)), log(as.vector(sd.prd)), use = 'p')

