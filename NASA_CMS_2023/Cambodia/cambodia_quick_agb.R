library(BIOMASS)
library(dplyr)
library(plyr)

treedata <- read.csv("L:/vclgp/minord/Mondolkiri_Srepok_treedata.csv")
allometricDefinitions <- read.csv("L:/vclgp/minord/git/gedicalval/config/allometricDefinitions.csv")
speciesdata <- read.csv("L:/vclgp/minord/Mondolkiri_Srepok_speciesdata.csv")

treedata <- treedata[1:1947,]

treedata <- dplyr::left_join(treedata, speciesdata, by = c("X6..Khmer.Name" = "K_Name"))

#rename DBH and height columns to work with existing columns
colnames(treedata)[11] <- 'd.stem'
colnames(treedata)[17] <- 'h.t'
colnames(treedata)[21] <- 'species'

treedata$d.stem <- treedata$d.stem/100
treedata$allom.key <- 1
treedata$allom.name <- "chave2014a"


fill_wsg_biomass <- function(data) {
 
  if ( any(!is.na(data$species)) ) {
    
   
    species <- unlist( lapply(data$species, function(x) strsplit(x," ")[[1]][2]) )
    genus <- unlist( lapply(data$species, function(x) strsplit(x," ")[[1]][1]) )
    
    invisible(capture.output( dataWD <- BIOMASS::getWoodDensity(genus = genus,
                                      species = species)))
    
    data$wsg <- dataWD$meanWD
    data$wsg.sd <- dataWD$sdWD
    
  }
  
  data
}


fill_height <- function(data, height_model) {
  
  ii <- is.na(data$h.t) & !is.na(data$d.stem) & (data$d.stem > 0)
  if ( any(ii, na.rm=TRUE) ) {
    jj <- !is.na(data$d.stem) & !is.na(data$h.t) & (data$d.stem > 0) & (data$h.t > 0)
    if ( any(jj, na.rm=TRUE) ) {
      if (  !is.na(height_model) ) {
        HDmodel <- BIOMASS::modelHD(D=data$d.stem[jj]*100, H=data$h.t[jj], useWeight=TRUE, 
                                    method=height_model)
        HDest <- BIOMASS::retrieveH(D=data$d.stem*100, model=HDmodel)
        data$h.t.mod <- HDest$H
        errH <- HDest$RSE
      } else {
        errH <- NA
      }
    } else {
      errH <- NA
    }
  } else {
    errH <- NA
  }
  
  list(data=data,HtModelRSE=errH) 
}


predict_biomass <- function(data, allometricDefinitions) {
 
  
  # Only process valid measurements
  a <- join(data[!is.na(data$d.stem),], allometricDefinitions, by="allom.key", type="left")  

  ## Chave et al. (2014). Order of allometric selection.
  ## 1. All height measurements present
  ## 2. Local D-H model heights
  ## 3. Height modelled using E
  ii <- is.na(a$m.agb) & (a$allom.name %in% c("chave2014a","chave2014b")) & ( !is.na(a$h.t) | !is.na(a$h.t.mod) )
  if ( any(ii, na.rm=TRUE) ) {
    h.t.tmp <- a$h.t[ii]
    jj <- is.na(h.t.tmp) | (h.t.tmp <= 0)
    h.t.tmp[jj] <- a$h.t.mod[ii][jj]
    a$m.agb[ii] <- BIOMASS::computeAGB(a$d.stem[ii]*100, a$wsg[ii], H=h.t.tmp) * 1e3
    a$allom.key[ii] <- 2
  }
  ii <- is.na(a$m.agb) & (a$allom.name %in% c("chave2014a","chave2014b")) & is.na(a$h.t) & is.na(a$h.t.mod)
  if ( any(ii, na.rm=TRUE) ) {
    geo.coords <- cbind(a$longitude[ii],a$latitude[ii])
	print(geo.coords)
    a$m.agb[ii] <- BIOMASS::computeAGB(a$d.stem[ii]*100, a$wsg[ii], coord=geo.coords) * 1e3
    a$allom.key[ii] <- 1
  }
  

  
  ## Combine results
  data$m.agb[!is.na(data$d.stem)] <- a$m.agb
  data$allom.key[!is.na(data$d.stem)] <- a$allom.key
  data$d.stem.valid[!is.na(data$d.stem)] <- a$d.stem.valid
  
  
  data
}

## Fill wood specific gravity estimates
data <- fill_wsg_biomass(treedata)

## Fill height using local/regional DBH-Height relationship
height_model <- 'log1'
result <- fill_height(data, height_model)
data <- result$data
data$m.agb <- NA

## Estimate above-ground biomass
data <- predict_biomass(data, allometricDefinitions)

## Some trees were duplicated due to the same Khmer name being matched to multiple scientific names. 
## Since we don't know what's correct, average wsg and m.agb, and remove the duplicates
## To be fixed after we get more info about the dataset

# Columns to use to ID trees
key <- c("X1..Cluster.Number","X2..Plot.Number","X4..Tree.No.","X6..Khmer.Name","X7..Distance","X8..Bearing","d.stem")

# Find duplicates and extract IDs of duplicates
dup <- duplicated(data[,key])
dup_ids <- unique(data[dup,key])

# Format sigfigs of IDs so they will match full dataset
dup_ids$X7..Distance <- as.character(formatC(dup_ids$X7..Distance, digits=2, format="f", flag="#"))
dup_ids$X8..Bearing <- as.character(formatC(dup_ids$X8..Bearing, digits=1, format="f", flag="#"))
dup_ids$d.stem <- as.character(formatC(dup_ids$d.stem, digits=4, format="f", flag="#"))

for(i in 1:nrow(dup_ids)){
	# Concatenate IDs from full dataset and duplicates to find all matches for averaging
	# The number of times a tree was duplicated varies
	data_id_string <- gsub("[[:space:]]", "", apply(data[,key], 1, paste, collapse=','))
	dup_id_string <- gsub("[[:space:]]", "", paste(dup_ids[i,],collapse=','))


	ii <- data_id_string %in% dup_id_string
	print(sum(ii))
	
	# Average wood specific gravity and biomass
	data$wsg[ii] <- mean(data$wsg[ii])
	data$m.agb[ii] <- mean(data$m.agb[ii])
}

# Remove duplicate rows
data <- data[-which(dup),]
		



write.csv(data, "L:/vclgp/minord/Mondolkiri_Srepok_tree_agb_20240612.csv")