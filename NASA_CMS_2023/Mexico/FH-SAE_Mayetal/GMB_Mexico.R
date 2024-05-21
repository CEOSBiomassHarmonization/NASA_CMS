library(INLA)

library(sf)
library(terra)
library(dplyr)
library(spdep)
library(exactextractr)

library(ggplot2)
library(viridis)

########################################
# Geostatistical Model
########################################

# 'Geostatistical' typically means a spatial model that utilized point-referenced
# data: every data point is associated with a particular location, 
# rather than an area, for instance.

# Our target variable at location s (AGBD, here) is assumed to follow the model

# y(s) = (\alpha + \tilde{\alpha}(s)) + (\beta + \tilde{\beta}(s))*x(s) + \epsilon(s)

# Most of the values should be familiar. \alpha and \beta are constant regression
# coefficients. There spatially-varying versions are now assumed to follow
# point-referenced Gaussian processes, which vary continuously with location s.

# Modeling with Gaussian processes natively is very computationally expensive
# when the number of observations becomes large. We can circumvent this with
# the stochastic partial differential equation (SPDE) approach.
# This will hopefully become more clear as the demonstation goes on, 
# and I have further explanations, but I will give a short overview now.
# We first define an underlying "mesh". This draws a graph of triangles across the
# study area. On the vertices of triangles, we define a neighbor-based spatial model
# (much like the area models before) where two vertices are neighbors if there
# is a triangle edge connecting them. To exted this to continuous locations, the
# value of the Gaussian process at any location s lying within a triangle is
# linearly interpolated from the 3 vertices of the containing triangle!
# The reason why this is more computationally effecient is that the
# process on the vertices is SPARSE: for any given vertex, it only needs
# the information from the connected vertices.

# Say we have a mesh with k vertices. Taking \tilde{\alpha}(s) as an example,
# we give the representation
#     \tilde{\alpha}(s) = \vec{a}(s)^T \vec{w}_\alpha
# where \vec{w}_\alpha is a vector of values at the k vertices. To
# interpolate these into any location s, we multiply them by weights
# \vec{a}(s) and add everything together. Weights \vec{a}(s) has at most
# 3 non-zero entries, corresponding to the three vertices of the triangle
# containing location s.




# Load CCI raster
# I wanted everything to be in meters, so I projected to
# the EASE grid. Takes a while on my machine...
cci.rast = rast("./data/CCI_AOImasked.tif") %>%
  project("epsg:6933")

# Load plot data
plot = st_read("./data/INFYS_Arbolado_2015_2020_AGBD.gpkg") %>%
  st_transform(crs = 6933)

agbd.plot = plot$NFI_AGBD_flt

loc.plot = st_coordinates(plot)

n.plot = nrow(loc.plot)


#################################

# What's this relationship look like at the plot level?

cci.plot = extract(x = cci.rast, y = loc.plot)[,1] 

ggplot() + geom_point(aes(x = cci.plot, y = agbd.plot), alpha = 0.2) 

#pretty grim... I think we need something better.

####################

# I'm loading the areas and merging them into one mexico shapefile
mexico = st_read("./Data/ecort08gw_DESECON1_DISS.gpkg") %>%
  st_union() %>%
  st_transform(crs = 6933)

# Okay, this is the most complicated extra bit of the geostatistical model.
# We are going to model the Gaussian processes with the 
# stochastic partial differential equation (SPDE) approach.
# A Gaussian process is the solution to a stochastic differential equation.
# Instead of using the solution directly, we can use a finite element approximation
# (FEM, very common in engineering) to represent the solution.
# The mesh partitions our study are into triangles. These triangles are the finite elements.
# A neighbor-based spatial model is defined on the triangle vertices.
# Then, for all values within the triangles, the spatial process is 
# linearly interpolated from the vertices.
# The finer the mesh (to an extent), the better the approximation to the full
# Gaussian process. But the finer the mesh, the more computationally expensive
# the model fitting will be. It's important that the max triangle edge length (more later)
# is at most 1/3 of the spatial range that you estimate. Usually this entails a couple
# rounds of model fitting.

# Max edge of 15 kilometers, to start
max.edge = 15*10^3

# We need a buffer on the domain to prevent boundary effects from
# the differential equation boundary conditions...
mexico.buffer = st_buffer(mexico, dist = max.edge*5)

mesh = inla.mesh.2d(boundary = list(as(mexico, "Spatial"), 
                                    as(mexico.buffer, "Spatial")),
                    max.edge = c(max.edge, 3*max.edge),
                    cutoff = 2*max.edge/3,
                    offset = c(max.edge, 5*max.edge))

# This is the resulting number of vertices.
# It's the ultimate determinant of computational expense
k = mesh$n

k

# This is what the mesh looks like
plot(mesh)

# Let's zoom in to a corner, so you can see what the mesh and buffer
# looks like
ggplot() + inlabru::gg(data = mesh) +
  coord_sf(crs = 6933, xlim = c(-9.75*10^6, -9.25*10^6), ylim = c(3*10^6, 3.3*10^6))


# We need to put priors on the Gaussian process sd and range parameters.
# We are using the penalized complexity (PC) priors.
# They operate by putting prior probability bound on the data.

# Priors for the varying intercept
spde.alpha = inla.spde2.pcmatern(mesh, 
                                 prior.range = c(30*10^3, 0.01), # This says the probability that the range is LESS than 30 km is 0.01
                                 prior.sigma = c(60, 0.01)) # This says the probability that the sd is GREATER than 60 Mg/ha is 0.01

# Priors for the varying coefficient
spde.beta = inla.spde2.pcmatern(mesh, 
                                prior.range = c(30*10^3, 0.01), # Same interpretations here
                                prior.sigma = c(1, 0.5))

# This is the mesh projection matrix that interpolates the vertices to the plot locations.
# Each row interpolates each observation. So the ith row is the vector of weights that
# interpolates the vertices to the ith location. Results in an n x k matrix.
A.plot = inla.spde.make.A(mesh = mesh, loc = loc.plot)

# formula
formula = agbd ~
  -1 + intercept + 
  cci +
  f(alpha.spat, model = spde.alpha) +
  f(beta.spat, model = spde.beta)

# stack
stack = inla.stack(data = list(agbd = agbd.plot),
                   A = list(1,
                            1,
                            A.plot,
                            Diagonal(x = cci.plot)%*%A.plot),
                   effects = list(
                     intercept = rep(1, n.plot),
                     cci = cci.plot,
                     alpha.spat = 1:k,
                     beta.spat = 1:k
                   ))

# Fit this thing, its takes a while (an hour on my laptop).
# I saved my model fit so you can just load it and start messing around


# model_fit = inla(formula = formula, # Provide the formula
#                  family = 'gaussian', # We assume our data follows a Gaussian generalized linear model (GLM)
#                  data = inla.stack.data(stack), # These line and the one below are simply how we feed INLA the stack.
#                  control.predictor = list(A = inla.stack.A(stack)),
#                  control.compute = list(config = T, dic = T, waic = T),
#                  control.inla = list(int.strategy = "eb"),
#                  verbose = T)

load("INLA_model_fit.RData")

# returns the model parameters.
summary(model_fit)
# One thing to notice is that our expected ranges are ~20 km.
# Our max triangle edge length was 15 km, so not much smaller....
# We almost certainly want to reiterate with a finer mesh,
# but for now we will carry on.

##################################
# Prediction example
##################################


# Let's get a test region to make a prediction for
MU = st_read("./Data/ecort08gw_DESECON4_DISS.gpkg") %>%
  st_transform(crs = 6933)

currentMU = MU[25,]

currentMU$DESECON4

currentPlots = st_filter(plot, currentMU)

# Pretty far away from a random sample...
ggplot() + geom_sf(data = currentMU) + geom_sf(data = currentPlots) 


# To make predictions for an area, we first make a grid
# of point-level predictions. We can then collapse them
# into a prediction for the whole area.
# In general, it's good for your grid to be dense relative to the
# area you are predicting over, with a minimum resolution of the plot area,
# i.e., you are predicting side by side plotse


# We'll make a 2 km grid
grid = st_make_grid(currentMU, cellsize = c(2000,2000), what = "centers") %>%
  st_as_sf() %>%
  st_filter(currentMU) %>%
  st_coordinates()

plot(grid, pch = '.')

# We'll need our covariate values at the grid locations
cci.pred = extract(cci.rast, grid)[,1]

# The plots over the cci grid values...
ggplot() + 
  geom_raster(aes(x = grid[,1], y = grid[,2], fill = cci.pred)) +
  geom_sf(data = currentMU, fill = NA) + geom_sf(data = currentPlots, aes(col = "Plots")) +
  scale_fill_viridis()
# This makes the plot samples look even worse. It looks like
# the higher bio areas are systematically not sampled...

# Our projection matrix for the grid locations
A.pred = inla.spde.make.A(mesh = mesh, loc = grid)

# Draw posterior samples from the model.
# We can use the 'selection = ' option to
# limit the sampler to the parameters/effects we actually
# need for prediction. You also give the indices for 
# each of the parameters/effects, in case you only need
# some indices of particular parameter/effect
samples = inla.posterior.sample(n = 250, result = model_fit,
                                selection = list(intercept = 1, # Only one intercept term
                                                 cci = 1, # Only one regression coefficient
                                                 alpha.spat = 1:k, # we want to sample at all k mesh nodes
                                                 beta.spat = 1:k)) # same as above



# Our prediction equation. Essentially just the model, but
# we omit the random error, as it will not matter as we aggregate over
# moderate-to-large areas.
pred_fun = function(...){
  drop(intercept + 
    cci.pred*cci + 
    A.pred%*%alpha.spat[1:mesh$n] +
    Diagonal(x = cci.pred)%*%A.pred%*%beta.spat)
}

# Generate prediction samples
pred.samples = inla.posterior.sample.eval(fun = pred_fun,
                                          samples = samples)

# Model expectations and SD's at the grid locations
pred.mu = rowMeans(pred.samples)
pred.sd = apply(pred.samples, 1, sd)


ggplot() + geom_raster(aes(x = grid[,1], y = grid[,2], fill = pred.mu)) +
  scale_fill_viridis() + coord_sf(crs = 6933) +
  labs(fill = "Expected AGBD (Mg/ha)")

ggplot() + geom_raster(aes(x = grid[,1], y = grid[,2], fill = pred.sd)) +
  scale_fill_viridis(option = "turbo") + coord_sf(crs = 6933) +
  labs(fill = "SD (Mg/ha)")

# We can tell the mesh isn't quite fine enough, because the triangle
# vertices are occasionally visible in the map (especially bottom right)
# I'm not actually sure if this will make much of a difference for
# predicting large area averages, but better safe...

# Collapse these samples to get an estimate for the whole MU
MU.posterior = colMeans(pred.samples, na.rm = T)


hist(MU.posterior) 

mean(MU.posterior)
sd(MU.posterior)

###################################
# Model refitting
###################################

# Often it's better to start with a coarser mesh and then refine it as needed.
# When you refit the model with the finer mesh, you don't need to start 
# the model fitting from square one. You can give INLA the final parameter 
# values from the previous fit as a starting value, which should speed things up.

# Below is the same as the above model fitting, but with a finer mesh.

# Max edge of 10 kilometers, this time.
max.edge = 10*10^3


mexico.buffer = st_buffer(mexico, dist = max.edge*5)

mesh = inla.mesh.2d(boundary = list(as(mexico, "Spatial"), 
                                    as(mexico.buffer, "Spatial")),
                    max.edge = c(max.edge, 3*max.edge),
                    cutoff = 2*max.edge/3,
                    offset = c(max.edge, 5*max.edge))


# Many more vertices, this time
k = mesh$n

k

# This is what the mesh looks like
plot(mesh)

# Let's zoom in to a corner, so you can see what the mesh and buffer
# looks like
ggplot() + inlabru::gg(data = mesh) +
  coord_sf(crs = 6933, xlim = c(-9.75*10^6, -9.25*10^6), ylim = c(3*10^6, 3.3*10^6))

# Priors. These are the same
# Priors for the varying intercept
spde.alpha = inla.spde2.pcmatern(mesh, 
                                 prior.range = c(30*10^3, 0.01), # This says the probability that the range is LESS than 30 km is 0.01
                                 prior.sigma = c(60, 0.01)) # This says the probability that the sd is GREATER than 60 Mg/ha is 0.01

# Priors for the varying coefficient
spde.beta = inla.spde2.pcmatern(mesh, 
                                prior.range = c(30*10^3, 0.01), # Same interpretations here
                                prior.sigma = c(1, 0.5))

# mesh projection matrix 
A.plot = inla.spde.make.A(mesh = mesh, loc = loc.plot)

# formula. same as before
formula = agbd ~
  -1 + intercept + 
  cci +
  f(alpha.spat, model = spde.alpha) +
  f(beta.spat, model = spde.beta)

# stack
stack = inla.stack(data = list(agbd = agbd.plot),
                   A = list(1,
                            1,
                            A.plot,
                            Diagonal(x = cci.plot)%*%A.plot),
                   effects = list(
                     intercept = rep(1, n.plot),
                     cci = cci.plot,
                     alpha.spat = 1:k,
                     beta.spat = 1:k
                   ))

# Fit this thing. Notice the arguments in 'control.mode = '
# We tell it the starting parameters in theta = 
# but tell it to restart the optimization from there (restart = T),
# and not just settle at these values
# Still will probably take a while, but less than it would have...

model_fit_v2 = inla(formula = formula, # Provide the formula
                 family = 'gaussian', # We assume our data follows a Gaussian generalized linear model (GLM)
                 data = inla.stack.data(stack), # These line and the one below are simply how we feed INLA the stack.
                 control.predictor = list(A = inla.stack.A(stack)),
                 control.compute = list(config = T, dic = T, waic = T),
                 control.inla = list(int.strategy = "eb"),
                 control.mode = list(theta = model_fit$mode$theta, restart = T),
                 verbose = T)
