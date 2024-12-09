rm(list=ls())
###############
## libraries ##
library(MASS)

################################
## generate synthetic dataset ##
coords.all <- as.matrix(expand.grid(seq(0,1,length.out = 100),seq(0,1,length.out = 100)))

## get number of gridcells
N <- nrow(coords.all)

## get intersite distances between gridcells
dmat <- as.matrix(dist(coords.all, diag = TRUE, upper = TRUE))

## simulate an x variable with spatial structure
x.sigma.sq <- 2
x.phi <- 3/.15
x.mu <- rep(5,N)

x.cov <- x.sigma.sq*exp(-x.phi*dmat)

x.all <- mvrnorm(1,x.mu,x.cov)

## simulate spatially structured model error
w.sigma.sq <- 7

w.phi <- 3/.25

w.cov <- w.sigma.sq*exp(-w.phi*dmat)

w <- mvrnorm(1,rep(0,N),w.cov)

## set regression parameters
beta.0 <- 1
beta.1 <- 3

## set non-spatial error variance
tau.sq <- 2

## simulate y
y.mu <- beta.0 + x.all*beta.1 + w
y.all <- rnorm(N,y.mu, sqrt(tau.sq))^2

save(x.all,y.all,coords.all, file = "synthetic-data.RData")

## done simulating dataset ##
#############################

par(mfrow = c(1,2))
plot(coords, pch = 19, cex = .5, col = hcl.colors(50)[cut(x,50)],
     main = "x",xlab = "", ylab = "")
plot(coords, pch = 19, cex = .5, col = hcl.colors(50)[cut(y,50)],
     main = "y",xlab = "", ylab = "")









