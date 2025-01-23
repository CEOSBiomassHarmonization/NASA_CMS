library(terra)

args <- commandArgs(trailingOnly = TRUE)

input_file <- args[1]
split_string <- args[2] #"TMST_"
outfile_string <- args[3] #"_JRC_Transition_Map.tif"

# Generate output file name
outname_extent <- strsplit(basename(input_file), split_string)[[1]][2]
outname_extent <- strsplit(outname_extent, ".tif")[[1]][1]
output_file <- file.path(paste0(outname_extent, outfile_string)) 

tile_string <- outname_extent
parts <- strsplit(tile_string, "_")[[1]]

ymax_str <- substr(parts[1], 1, 2) # Extract the numeric part of ymax
ymax_dir <- substr(parts[1], 3, 3) # Extract the direction (N/S)
xmin_str <- substr(parts[2], 1, 3) # Extract the numeric part of xmin
xmin_dir <- substr(parts[2], 4, 4) # Extract the direction (E/W)

ymax <- as.numeric(ymax_str)
if (ymax_dir == "S") {ymax <- -ymax}

xmin <- as.numeric(xmin_str)
if (xmin_dir == "W") {xmin <- -xmin}

xmax <- xmin + 10
ymin <- ymax - 10
target_extent <- ext(xmin, xmax, ymin, ymax)

# Define resolution and nodata value
resolution <- c(0.00025, 0.00025)
nodata_value <- 0.0 
r <- rast(input_file)

# Create an empty raster with the target extent and resolution
target_raster <- rast(extent = target_extent, resolution = resolution, crs = crs(r))

# Apply warp (resampling and reprojecting)
r_warped <- resample(r, target_raster, method = "near")

# Write the warped raster to the output file
writeRaster(r_warped, output_file, overwrite = TRUE, gdal = c("COMPRESS=DEFLATE", "PREDICTOR=2", "ZLEVEL=9"))