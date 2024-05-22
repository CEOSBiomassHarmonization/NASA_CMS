############# OPEN R, THEN RUN THE FOLLOWING COMMAND ###################

# options(repos=c(CRAN="https://cran.r-project.org"))
# install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
# install.packages("fmesher", dependencies = TRUE)
# install.packages("MatrixModels", type = "source")
# install.packages("exactextractr")
# install.packages("sn" ,dependencies = TRUE)
# packages <- c("terra","dplyr","spdep", "exactextractr", "sf","ggplot2","viridis","sn","fmesher","exactextractr","fields","Matrix","inlabru")
# package.check <- lapply(packages, FUN = function(x) {
#     if (!require(x, character.only = TRUE)) {
#         install.packages(x, dependencies = TRUE)
#         library(x, character.only = TRUE, quietly=TRUE)
#     }
# })

############## LOAD PACKAGES ###########

library("fmesher")
library(MatrixModels)
library(Matrix)
library(INLA)
library(inlabru)
library(sf)
library(terra)
library(dplyr)
library(spdep)
library(exactextractr)
library(sn)

########################################
args <- commandArgs(trailingOnly = TRUE)

#### LOAD THE DATA USED TO MAKE THE MODEL ######
mexico = st_read(args[3], quiet = TRUE) %>% st_union() %>% st_transform(crs = 6933) # 
cci.rast = rast(args[4])
hei.rast = rast(args[5])

DATA <- read.csv(args[6]) #Read the saved data
cci.plot <- DATA$cci.plot
hei.plot <- DATA$hei.plot
loc.plot <- data.frame(matrix(ncol = 2, nrow = length(hei.plot)))
loc.plot[,1] <- DATA$X
loc.plot[,2] <- DATA$X.1
agbd.plot <- sqrt(DATA$agbd.plot)
loc.plot <- data.matrix(loc.plot)

###### REMAKE THE MESH IN CASE IT IS NOT STORED IN MEMORY ######
max.edge = 10*10^3
mexico.buffer = st_buffer(mexico, dist = max.edge*5) 
mesh = inla.mesh.2d(boundary = list(as(mexico, "Spatial"), as(mexico.buffer, "Spatial")), max.edge = c(max.edge, 3*max.edge), cutoff = 2*max.edge/3, offset = c(max.edge, 5*max.edge)) 
k = mesh$n 
A.plot = inla.spde.make.A(mesh = mesh, loc = loc.plot)

##### LOAD THE MODEL RESULTS ####
load(args[7])
samples = inla.posterior.sample(n = 250, result = model_fit_v2) 

#### LOAD AN URBAN MASK AND PROBABILITY OF FOREST MASK FOR MEXICO ####
URBAN = st_read(args[8],quiet=TRUE)
FNF = rast(args[9])

##### MAKE PREDICTIONS ###########

PROJECTS = st_read(args[2], quiet = TRUE) %>% st_transform(crs = 6933)
PROJECTS$ID <- seq.int(nrow(PROJECTS))
PROJECTS$GMB_mean_AGBD_dense = NA
PROJECTS$GMB_sd_AGBD_dense = NA

for (f in (1:nrow(PROJECTS)))  { 
    PROJECTS_ID = PROJECTS[f,]$ID
    grid_to_predict = st_make_grid(PROJECTS[f,], cellsize = c(100,100), what = "centers") %>% st_as_sf() %>% st_filter(PROJECTS[f,])
    grid_to_predict = grid_to_predict[URBAN, ,op=st_disjoint] 

    if (nrow(grid_to_predict[URBAN, ,op=st_disjoint])>0){
        grid = grid_to_predict %>% st_coordinates()
        cci.pred = extract(cci.rast, grid)[,1]
        hei.pred = extract(hei.rast, grid)[,1]

        if (length(cci.pred) > 0 && length(hei.pred) > 0){
            A.pred = inla.spde.make.A(mesh = mesh, loc = grid) 
            pred_fun = function(...){
              drop(intercept + 
                cci.pred*cci + 
                hei.pred*hei +
                A.pred%*%alpha.spat[1:k] +   
                Diagonal(x = cci.pred)%*%A.pred%*%beta.spat[1:k] +
                Diagonal(x = hei.pred)%*%A.pred%*%eta.spat[1:k]) +
                rnorm(nrow(A.pred), sd = sqrt(1/theta[1])) 
            }

            # Generate prediction samples
            pred.samples = (inla.posterior.sample.eval(fun = pred_fun, samples = samples))
            pred.samples[pred.samples < 0] = 0

            # Model expectations and SD's at the grid locations
            pred.mu_dense = Matrix::rowMeans(pred.samples^2, na.rm = TRUE) 
            pred.sd_dense = apply(pred.samples^2, 1, sd, na.rm = TRUE) 
            
            # Enter information to CSV file
            PROJECTS$GMB_mean_AGBD_dense[PROJECTS$ID == PROJECTS_ID] <- round(mean(pred.mu_dense,na.rm=TRUE),digits=2)
            PROJECTS$GMB_sd_AGBD_dense[PROJECTS$ID == PROJECTS_ID] <- round(sd(colMeans(pred.samples^2, na.rm = T),na.rm=TRUE),digits=2)
        }
    }
}

#### SAVE OUR RESULTS 
write.csv(st_drop_geometry(PROJECTS),paste0(args[1]))