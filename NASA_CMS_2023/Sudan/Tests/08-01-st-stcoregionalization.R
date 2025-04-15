## ----sett, echo = FALSE, results = 'hide', warning=FALSE, message=FALSE----
source("R/initial_setup.R")
opts_chunk$set(
  fig.path = 'figs/stcoreg-'
)
set.seed(1)

## ----param---------------------------------------------------------------
alpha <- c(-5, 3, 10) # intercept on reparametrized model
z.sigma = c(0.5, 0.6, 0.7) # random field marginal std
range = c(0.2, 0.3, 0.4) # GRF scales: range parameters
beta <- c(0.7, 0.5, -0.5) # copy par.: reparam. coreg. par.
rho <- c(0.7, 0.8, 0.9) # temporal correlations
n <- 50 # number of spatial locations
k <- 4  # number of time points
e.sigma <- c(0.3, 0.2, 0.15) # The measurement error marginal std

## ----sloc----------------------------------------------------------------
loc <- cbind(runif(n), runif(n)) 

## ----rfs, results = 'hide'-----------------------------------------------
x1 <- book.rMatern(k, loc, range = range[1], sigma = z.sigma[1])
x2 <- book.rMatern(k, loc, range = range[2], sigma = z.sigma[2])
x3 <- book.rMatern(k, loc, range = range[3], sigma = z.sigma[3])

## ----st------------------------------------------------------------------
z1 <- x1
z2 <- x2
z3 <- x3

for (j in 2:k) {
  z1[, j] <- rho[1] * z1[, j - 1] + sqrt(1 - rho[1]^2) * x1[, j]
  z2[, j] <- rho[2] * z2[, j - 1] + sqrt(1 - rho[2]^2) * x2[, j]
  z3[, j] <- rho[3] * z3[, j - 1] + sqrt(1 - rho[3]^2) * x3[, j]
}   

## ----yyy-----------------------------------------------------------------
y1 <- alpha[1] + z1 + rnorm(n, 0, e.sigma[1])
y2 <- alpha[2] + beta[1] * z1 + z2 + rnorm(n, 0, e.sigma[2])
y3 <- alpha[3] + beta[2] * z1 + beta[3] * z2 + z3 +
  rnorm(n, 0, e.sigma[3])

## ----mesh----------------------------------------------------------------
mesh <- inla.mesh.2d(loc, max.edge = 0.2, offset = 0.1,
  cutoff = 0.1)

## ----spde----------------------------------------------------------------
spde <- inla.spde2.pcmatern(mesh = mesh, 
  prior.range = c(0.05, 0.01), # P(range < 0.05) = 0.01
  prior.sigma = c(1, 0.01)) # P(sigma > 1) = 0.01

## ----idx-----------------------------------------------------------------
s1 <- rep(1:spde$n.spde, times = k)
s2 <- s1
s3 <- s1
s12 <- s1
s13 <- s1
s23 <- s1

g1 <- rep(1:k, each = spde$n.spde)
g2 <- g1
g3 <- g1
g12 <- g1
g13 <- g1
g23 <- g1

## ----pbeta---------------------------------------------------------------
rho1p <- list(theta = list(prior = 'pccor1', param = c(0, 0.9))) 
ctr.g <- list(model = 'ar1', hyper = rho1p)

## ----pcopy---------------------------------------------------------------
hc1 <- list(theta = list(prior = 'normal', param = c(0, 10)))

## ----form----------------------------------------------------------------
form <- y ~ 0 + intercept1 + intercept2 + intercept3 + 
  f(s1, model = spde, group = g1, control.group = ctr.g) + 
  f(s2, model = spde, group = g2, control.group = ctr.g) + 
  f(s3, model = spde, group = g3, control.group = ctr.g) + 
  f(s12, copy = "s1", group = g12, fixed = FALSE, hyper = hc1) + 
  f(s13, copy = "s1", group = g13, fixed = FALSE, hyper = hc1) + 
  f(s23, copy = "s2", group = g23, fixed = FALSE, hyper = hc1) 

## ----stlokA--------------------------------------------------------------
stloc <- kronecker(matrix(1, k, 1), loc)# repeat coord. each time
A <- inla.spde.make.A(mesh, stloc, n.group = k,
  group = rep(1:k, each = n))

## ----stack---------------------------------------------------------------
stack1 <- inla.stack(
  data = list(y = cbind(as.vector(y1), NA, NA)),
  A = list(A), 
  effects = list(list(intercept1 = 1, s1 = s1, g1 = g1))) 

stack2 <- inla.stack(
  data = list(y = cbind(NA, as.vector(y2), NA)),
  A = list(A), 
  effects = list(list(intercept2 = 1, s2 = s2, g2 = g2, 
    s12 = s12, g12 = g12))) 

stack3 <- inla.stack(
  data = list(y = cbind(NA, NA, as.vector(y3))),
  A = list(A), 
  effects = list(list(intercept3 = 1, s3 = s3, g3 = g3, 
    s13 = s13, g13 = g13, s23 = s23, g23 = g23))) 

stack <- inla.stack(stack1, stack2, stack3) 

## ----fixnugget-----------------------------------------------------------
eprec <- list(hyper = list(theta = list(prior = 'pc.prec',
  param = c(1, 0.01))))

## ----initheta------------------------------------------------------------
theta.ini <- c(log(1 / e.sigma^2), 
  c(log(range), log(z.sigma), 
  qlogis(rho))[c(1, 4, 7, 2, 5, 8, 3, 6, 9)], beta)

# We jitter the starting values to avoid artificially
# recovering the true values
theta.ini = theta.ini + rnorm(length(theta.ini), 0, 0.1)

## ----result,results='hide'-----------------------------------------------
result <- inla(form, rep('gaussian', 3),
  data = inla.stack.data(stack), 
  control.family = list(eprec, eprec, eprec), 
  control.mode = list(theta = theta.ini, restart = TRUE),
  control.inla = list(int.strategy = 'eb'), 
  control.predictor = list(A = inla.stack.A(stack)))

## ----cpu, echo = FALSE---------------------------------------------------
result$cpu

## ----tabstcoreg0, echo = FALSE-------------------------------------------
## Intercepts
tab.stcoreg1 <- cbind(true = alpha, result$summary.fixed[, 1:6])
## Precision errors
tab.stcoreg2 <- cbind(true = c(e = e.sigma^-2), result$summary.hyperpar[1:3, ])
## Temporal correlations
tab.stcoreg3 <- cbind(true = rho, result$summary.hyperpar[c(6, 9, 12), ])
## Copy parameters
tab.stcoreg4 <- cbind(true = beta, result$summary.hyperpar[13:15, ])
## Range of field
tab.stcoreg5 <- cbind(true = range,
  result$summary.hyperpar[c(4, 7, 10), ])
## St. dev. of field
tab.stcoreg6 <- cbind(true = z.sigma, result$summary.hyperpar[c(5, 8, 11), ])

## ----label = "tabstcoreg", echo = FALSE----------------------------------
tab.stcoreg <- rbind(tab.stcoreg1, tab.stcoreg2, tab.stcoreg3, tab.stcoreg4,
  tab.stcoreg5, tab.stcoreg6)

tab.stcoreg <- cbind(Parameter = rownames(tab.stcoreg), tab.stcoreg)

names(tab.stcoreg) <- c("Parameter", "True", "Mean", "St. Dev.",
  "2.5\\% quant.", "50\\% quant.", "97.5\\% quant.", "Mode")

knitr::kable(tab.stcoreg[, c(1:5, 7)],
  row.names = FALSE,
  caption = "Summary of the posterior distributions of the parameters in the model.",
  format = "pandoc")

## ----stzfit, echo = FALSE, fig.width = 10, fig.height = 3.5, fig.cap = "True and fitted random field values."----

par(mfrow = c(1, 3), mar = c(2, 2, 0.5, 0.5), mgp=c(1.5, 0.5, 0))

plot(drop(A %*% result$summary.ran$s1$mean), as.vector(z1),
  xlab = '', ylab = '', asp = 1)
abline(0:1)

plot(drop(A %*% result$summary.ran$s2$mean), as.vector(z2),
  xlab = '', ylab = '', asp = 1)
abline(0:1)
plot(drop(A %*% result$summary.ran$s3$mean), as.vector(z3),
  xlab = '', ylab = '', asp = 1)
abline(0:1)

