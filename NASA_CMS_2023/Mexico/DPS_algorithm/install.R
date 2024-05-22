############# Install R packages ###################

cpu <- parallel::detectCores()

options(repos=c(CRAN="https://cloud.r-project.org/"))
install.packages("INLA"
                 ,repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/stable")
                 , dep=FALSE
                 , Ncpus=cpu
                 )
