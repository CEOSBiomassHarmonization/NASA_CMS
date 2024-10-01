###################### SUMMARY OF FUNCTIONS USED FOR BIOMASS MAP COMPARISONS ################################################
###################### SUMMARY OF FUNCTIONS USED FOR BIOMASS MAP COMPARISONS ################################################
###################### SUMMARY OF FUNCTIONS USED FOR BIOMASS MAP COMPARISONS ################################################

Apply_E_FMASK <- function(iso3,aoi,maps,maps_AOImasked,Out_folder){
    
    ################################### COUNTRY-SPECIFIC MASKS FIRST NEED TO BE DEFINED AND PREPARED FOR ANALYSIS ################################
    
    if (iso3=='LBR'){
        E_FMASK_type <- "raster"
        EMask_folder <- "/projects/my-public-bucket/Data/CI_Liberia"
        EMask_file <- "Ref_class_9_100m.tif"
        non_forest_value <- 0
        resample_method <- "near"
        RES_name <- file.path(Out_folder,paste0(iso3,'_Product_FMASK'),paste0(iso3,"_External_FMASK.tif"))
        if (!file.exists(RES_name)){
            FMASK_class <- rast(file.path(EMask_folder,EMask_file))
            FMASK_class[!is.finite(FMASK_class)] <- 0
            FMASK_class[FMASK_class <= non_forest_value] <- 0
            if (!identical(crs(FMASK_class), crs("+proj=longlat +datum=WGS84 +no_defs"))){FMASK_class <- projectrast(FMASK_class, crs=crs("+proj=longlat +datum=WGS84 +no_defs"))}
            writerast(FMASK_class, filename=RES_name,datatype="Int16",overwrite=TRUE)
        }
    }
    
    if (iso3=='COL'){
        E_FMASK_type <- "raster"
        EMask_folder <- "/projects/shared-buckets/leitoldv/Country_specific_data/COL/Wgs84_Z18N"
        EMask_file <- "bnb_2020_v8_210702_32618 - copia.img"
        non_forest_value <- 2
        resample_method <- "near"
        RES_name <- file.path(Out_folder,paste0(iso3,'_Product_FMASK'),paste0(iso3,"_External_FMASK.tif"))
        if (!file.exists(RES_name)){
            FMASK_class <- rast(file.path(EMask_folder,EMask_file))
            FMASK_class[FMASK_class >= non_forest_value] <- 0
            if (!identical(crs(FMASK_class), crs("+proj=longlat +datum=WGS84 +no_defs"))){FMASK_class <- projectrast(FMASK_class, crs=crs("+proj=longlat +datum=WGS84 +no_defs"))}
            writerast(FMASK_class, filename=RES_name,datatype="Int16",overwrite=TRUE)
        }
    }
    
    if (iso3=='PRY'){
        E_FMASK_type <- "raster"
        EMask_folder <- "/projects/shared-buckets/minord/data/paraguay/"
        EMask_file <- "1_forest_cog.tif"
        band_name <- 'mode'
        non_forest_value <- 0
        resample_method <- "near"
        RES_name <- file.path(Out_folder,paste0(iso3,'_Product_FMASK'),paste0(iso3,"_External_FMASK.tif"))
        if (!file.exists(RES_name)){gdal_translate(file.path(EMask_folder,EMask_file),RES_name,co=c("COMPRESS=DEFLATE","PREDICTOR=2","ZLEVEL=9"),ot="Int16",a_srs="+proj=longlat +datum=WGS84 +no_defs", of="GTiff",overwrite=TRUE)}
    }
    
    if (iso3=='ECU'){
        E_FMASK_type <- "vector"
        EMask_folder <- "/projects/shared-buckets/leitoldv/Country_specific_data/ECU"
        EMask_file <- "Vegetal_NoVegetal.shp"
        if (!file.exists(file.path(Out_folder,paste0(iso3,'_Product_FMASK'),EMask_file))) {
            Vector <- spTransform(readOGR(file.path(EMask_folder,EMask_file)), CRS('+init=epsg:4326')) 
            writeOGR(obj=Vector, dsn=file.path(Out_folder,paste0(iso3,'_Product_FMASK')), layer=sub('.shp','',EMask_file), driver="ESRI Shapefile")
        }
    }
    
    ################################### MASKING COUNTRY AGB PRODUCTS WITH THE MASKS SPECIFIED ABOVE ################################
    
    maps_FMASKS <- ''
    
    if (E_FMASK_type=='vector'){
        for (i in 1:length(maps_AOImasked)) {
            Map_with_Fexternal_mask_temp <- file.path(Out_folder,paste0(iso3,'_Product_FMASK'),paste0(maps[i],"_EXTERNAL_FMASK_",iso3,"_trans.tif"))
            Map_with_Fexternal_mask <- file.path(Out_folder,paste0(iso3,'_Product_FMASK'),paste0(maps[i],"_EXTERNAL_FMASK_",iso3,".tif"))
            vector_file <- file.path(Out_folder,paste0(iso3,'_Product_FMASK'),EMask_file)
            vector_vals<-readOGR(vector_file)
            if (!file.exists(Map_with_Fexternal_mask)){
                file.copy(from = file.path(Out_folder,paste0(iso3,'_Product_AOI'),paste0(maps[i],"_AOImasked.tif")), to = file.path(Out_folder,paste0(iso3,'_Product_FMASK'),paste0(maps[i],"_AOImasked.tif")))
                file.rename(file.path(Out_folder,paste0(iso3,'_Product_FMASK'),paste0(maps[i],"_AOImasked.tif")),Map_with_Fexternal_mask_temp)
                gdalUtils::gdal_rasterize(vector_file,Map_with_Fexternal_mask_temp,b=1,add=TRUE,at=TRUE,burn=50000,l=sub('.shp','',EMask_file),output_Raster=TRUE)
                input <- rast(Map_with_Fexternal_mask_temp)
                input[input > 49999] <- 0 
                writerast(input, filename=Map_with_Fexternal_mask,datatype="FLT4S", format="GTiff", options="COMPRESS=LZW", overwrite=TRUE)
            }
            maps_FMASKS[[length(maps_FMASKS) + 1]] <- Map_with_Fexternal_mask
        }
    }
    
    if (E_FMASK_type=='raster'){
        for (i in 1:length(maps_AOImasked)) { 
            RESOLUTION_maps_AOImasked <- res(eval(parse(text=maps_AOImasked[i])))
            dims <- dim(eval(parse(text=maps_AOImasked[i])))
            dim1 <- dims[[1]]
            dim2 <- dims[[2]]
            extents <- ext(eval(parse(text=maps_AOImasked[i])))
            ex1 <- extents@xmin
            ex2 <- extents@ymin
            ex3 <- extents@xmax
            ex4 <- extents@ymax
            FMASK_country_file <- RES_name
            FMASK_out_name <- paste0(sub('.tif','_trans_',FMASK_country_file),maps[i],'.tif')   
            AOI_masked_country_file <- rast(file.path(Out_folder,paste0(iso3,'_Product_AOI'),paste0(maps[i],"_AOImasked.tif")))
            Map_with_Fexternal_mask <- file.path(Out_folder,paste0(iso3,'_Product_FMASK'),paste0(maps[i],"_EXTERNAL_FMASK_",iso3,".tif"))
            if (!file.exists(Map_with_Fexternal_mask)){
                gdal_translate(FMASK_country_file,FMASK_out_name,co=c("COMPRESS=DEFLATE","PREDICTOR=2","ZLEVEL=9"),tr=c(as.double(RESOLUTION_maps_AOImasked[[1]]),as.double(RESOLUTION_maps_AOImasked[[2]])),ot="Int16", of="GTiff",overwrite=TRUE)
                gdalwarp(FMASK_out_name,dstfile=sub('_trans','_warps',FMASK_out_name),ts=c(dim2,dim1),te=c(ex1,ex2,ex3,ex4),r=resample_method,t_srs="+proj=longlat +datum=WGS84 +no_defs",ot="Int16", co=c("COMPRESS=DEFLATE","PREDICTOR=2","ZLEVEL=9"), overwrite=TRUE, VERBOSE=FALSE)
                FMASK_ready <- rast(sub('_trans','_warps',FMASK_out_name))
                map_Fmasked <- mask(AOI_masked_country_file, FMASK_ready, filename="", inverse=FALSE, maskvalue=0, updatevalue=NA, updateNA=TRUE)
                writerast(map_Fmasked, filename=Map_with_Fexternal_mask,datatype="FLT4S", format="GTiff", options="COMPRESS=LZW", overwrite=TRUE)
            }
            maps_FMASKS[[length(maps_FMASKS) + 1]] <- Map_with_Fexternal_mask 
        }
    }
    
    maps_FMASKS <- maps_FMASKS[maps_FMASKS != ""]
    mydir <- file.path(Out_folder,paste0(iso3,'_Product_FMASK'))
    delfiles <- dir(path=mydir ,pattern="*.aux.xml")
    file.remove(file.path(mydir, delfiles))
    delfiles <- dir(path=mydir ,pattern="*warps*")
    file.remove(file.path(mydir, delfiles))
    delfiles <- dir(path=mydir ,pattern="*trans*")
    file.remove(file.path(mydir, delfiles))
    return(maps_FMASKS)
}



####################### Calculate AGB totals and summaries  ##################################################################

agbd_table <- function(year=2020, agb_thresholded, agb_masked){
    scale_factor=0.0001
    agb_boundpoly <- as(ext(agb_thresholded), 'SpatExtent')
    agb_boundpoly <- sf::st_as_sf(vect(agb_boundpoly))
    agb_boundpoly_wgs84 <- agb_boundpoly %>% st_transform(crs=4326)
    aoi_wgs84 <- aoi %>% st_transform(crs=4326)
    if (!st_intersects(agb_boundpoly, aoi_wgs84)) {stop('aoi does not intersect supplied Biomass extract')}

    if ((((xmin(agb_thresholded) >=-180) & (xmax(agb_thresholded) <=180)) || ((xmin(agb_thresholded) >=0) & (xmax(agb_thresholded) <=360))) &&
        (ymin(agb_thresholded) >=-90) & (ymax(agb_thresholded) <= 90)) {
      # Use the included calc_pixel_area function to calculate the area of 
      # one cell in each line of the raster, allowing for accurate areal 
      # estimates in square meters even when imagery is in 
      # WGS84.
      message('Data appears to be in latitude/longitude. Calculating cell areas on a sphere.')
      spherical_areas <- TRUE
      # Calculate the area of a single pixel in each line of the image (to 
      # avoid repeating this calculation later on)
      #x <- cci
      x <- agb_thresholded
      xleft <- xmin(x)
      xright <- xmin(x) + xres(x)
      ylower <- seq(from=(ymax(x) - yres(x)), by=(-yres(x)), length.out=nrow(x))
      yupper <- seq(from=ymax(x), by=(-yres(x)), length.out=nrow(x))
      poly_areas <- function(xl, xr, yl, yu) {
      areaPolygon(matrix(c(xl, yl,
                       xr, yl,
                       xr, yu,
                       xl, yu), ncol=2, byrow=TRUE))
                      }
        pixel_areas <- mapply(poly_areas, xleft, xright, ylower, yupper)

    } else {
      spherical_areas <- FALSE
      pixel_areas <- xres(agb_thresholded) * yres(agb_thresholded)
    }
    aoi <- spTransform(aoi, CRS(proj4string(agb_thresholded)))
    if (!('label' %in% names(aoi))) {
      aoi$label <- paste('AOI', seq(1:nrow(aoi@data)))
    }
    uniq_aoi_labels <- unique(aoi$label)
    years <- as.numeric(year)
    agb_table <- data.frame(year=rep(years, length(uniq_aoi_labels)),aoi=rep(uniq_aoi_labels, each=length(years)))
    agb_table$agb <- 0
    agb_table$agbd <- 0
    agb_table$area <- 0
    agb_table$n_pixel <- 0
    n_years <- length(years)
    for (n in 1:nrow(aoi)) {
        agb_table_st_row <- match(aoi[n, ]$label, agb_table$aoi)
        bs <- blockSize(agb_masked)
        for (block_num in 1:bs$n) {
            bl_st_row <- bs$row[block_num]
            bl_nrows <- bs$nrows[block_num]
            agb_bl <- getValuesBlock(agb_masked, bl_st_row, bl_nrows)
            if (spherical_areas) {
              bl_pixel_areas <- rep(pixel_areas[bl_st_row:(bl_st_row + bl_nrows - 1)], each=ncol(agb_thresholded))
            } else {
              bl_pixel_areas <- pixel_areas
            }
            for (i in 0:(n_years-1)) {
              agb_table$agb[agb_table_st_row + i] <- agb_table$agb[agb_table_st_row + i] +
                sum(agb_bl * bl_pixel_areas * scale_factor, na.rm=TRUE) 
              agb_table$area[agb_table_st_row + i] <- agb_table$area[agb_table_st_row + i] +
                sum((!is.na(agb_bl))*bl_pixel_areas * scale_factor, na.rm=TRUE)   
              agb_table$n_pixel[agb_table_st_row + i] <- agb_table$n_pixel[agb_table_st_row + i] +
                sum(!is.na(agb_bl))     
              if(bl_st_row==1){
                  agb_table$agbd[agb_table_st_row + i] <- mean(agb_bl, na.rm=TRUE)
              } else{
                  agb_table$agbd[agb_table_st_row + i] <- weighted.mean(c(agb_table$agbd[agb_table_st_row + i], 
                    mean(agb_bl, na.rm=TRUE)), w=c((bl_st_row-1), bl_nrows), na.rm=TRUE)
              } 
            }
        }
    }
    return(agb_table)
}

####################### Find suitable number of rows and columns for plots ######################
CR <- function(MAPS_READY) {
    if (length(MAPS_READY) < 4) {
            rows <- 1
            cols <- length(MAPS_READY)
    }
    if (length(MAPS_READY) > 3) {
        rows <- 2
        cols <- as.integer(ceiling(length(MAPS_READY)*(length(MAPS_READY) - 1 )/4))
    }
    return(c(rows,cols))
}

####################### Overlay and compare various AGB map products #############################

overlay_products <- function(Out_folder,maps,maps_AOImasked,AOI_biomes,coarsest_res,CI_map_AOI){
    jpeg(file.path(Out_folder,paste0(iso3,'_AGB_maps_Comparison.jpeg')),width=700, height=700*(length(maps)-1)/2)
    par(mfrow=c((length(maps_AOImasked)-1),2), mar=c(3,3,3,8), xpd=TRUE)
    options(repr.plot.width=15, repr.plot.height=15)
    base_vals <- eval(parse(text=paste0(maps_AOImasked[coarsest_res[1]])))
    NAflag(base_vals) <- 0 
    base_map_name <- file.path(Out_folder,paste0(iso3,'_Product_Comparisons'),paste0(maps[coarsest_res[1]],"_AOImasked_warped.tif"))
    writerast(base_vals, filename=base_map_name,datatype="FLT4S",overwrite=TRUE)
    base_map <- rast(base_map_name)
    CI <- CI_map_AOI
    base_map_low <- base_map - CI
    base_map_high <- base_map + CI
    set_res <-  res(base_map)
    set_extent <- ext(base_map)
    set_width_height <- dim(base_map)
    dim1 <- set_width_height[[1]]
    dim2 <- set_width_height[[2]]
    ex1 <- set_extent@xmin
    ex2 <- set_extent@ymin
    ex3 <- set_extent@xmax
    ex4 <- set_extent@ymax
    pal <-colorRampPalette(c("purple","blue","gray","yellow","red"))
    pal2 <-colorRampPalette(c("gray","green","red"))

    for (i in 1:length(maps_AOImasked)){
        if (i != coarsest_res[1]) { 
            input_file <- file.path(Out_folder,paste0(iso3,'_Product_AOI'),paste0(maps[i],"_AOImasked.tif"))
            output_vrt <- file.path(Out_folder,paste0(iso3,'_Product_Comparisons'),paste0(maps[i],"_AOImasked.vrt"))
            output_vrt_ND <- file.path(Out_folder,paste0(iso3,'_Product_Comparisons'),paste0(maps[i],"_AOImasked_ND.tif"))
            output_trans <- file.path(Out_folder,paste0(iso3,'_Product_Comparisons'),paste0(maps[i],"_AOImasked_trans.tif"))
            output_warps <- file.path(Out_folder,paste0(iso3,'_Product_Comparisons'),paste0(maps[i],"_AOImasked_warped.tif"))
            if (!file.exists(output_warps)) {
                gdalbuildvrt(gdalfile=input_file,output.vrt=output_vrt)
                vals <- rast(output_vrt)
                NAflag(vals) <- 0 
                writerast(vals, filename=output_vrt_ND,datatype="FLT4S",overwrite=TRUE)
                gdal_translate(output_vrt_ND,output_trans,r='average',co=c("BIGTIFF=IF_NEEDED","COMPRESSION=LZW","TILED=YES"),tr=set_res,overwrite=TRUE)
                gdalwarp(output_trans,dstfile=output_warps,ts=c(dim2,dim1),te=c(ex1,ex2,ex3,ex4),r='average',t_srs="+proj=longlat +datum=WGS84 +no_defs",overwrite=TRUE,VERBOSE=FALSE)
                if (file.exists(output_warps)){file.remove(output_vrt)}
                if (file.exists(output_warps)){file.remove(output_vrt_ND)}
                if (file.exists(output_warps)){file.remove(output_trans)}
            }
            if (file.exists(paste0(output_vrt,'.aux.xml'))){paste0(output_vrt,'.aux.xml')}
            if (file.exists(paste0(output_vrt_ND,'.aux.xml'))){paste0(output_vrt_ND,'.aux.xml')}
            if (file.exists(paste0(output_trans,'.aux.xml'))){paste0(output_trans,'.aux.xml')}
            difference <- rast(output_warps)-base_map
            plot(difference, main=paste0(maps[i]," AGBD ", eval(parse(text=paste0(maps[i],"_year"))), ' - ', maps[coarsest_res[1]]," AGBD ", eval(parse(text=paste0(maps[coarsest_res[1]],"_year"))), ' (Mg/ha)'), zlim=c((mean(values(difference),na.rm=TRUE)-2*sd(values(difference),na.rm=TRUE)),(mean(values(difference),na.rm=TRUE)+2*sd(values(difference),na.rm=TRUE))),col = pal(50))#;plot(AOI_biomes, add=T)
            difference_filename <- file.path(Out_folder,paste0(iso3,'_Product_Comparisons'),paste0(maps[i],"_AGBD_", eval(parse(text=paste0(maps[i],"_year"))), '-', maps[coarsest_res[1]],"_AGBD_", eval(parse(text=paste0(maps[coarsest_res[1]],"_year"))), '.tif'))
            writerast(difference, filename=difference_filename,datatype="FLT4S",overwrite=TRUE)
            agb_for_CI_comparison <- rast(output_warps)
            agb_for_CI_comparison[(agb_for_CI_comparison >= base_map_low) & (agb_for_CI_comparison <= base_map_high)] <- 50000
            agb_for_CI_comparison[agb_for_CI_comparison < base_map_low ] <- 1
            agb_for_CI_comparison[(agb_for_CI_comparison > base_map_high) &  (agb_for_CI_comparison != 50000)] <- 1
            agb_for_CI_comparison[agb_for_CI_comparison > 49999.999] <- 0
            agb_for_CI_comparison[(agb_for_CI_comparison != 0) &  (agb_for_CI_comparison != 1)] <- -1
            plot(agb_for_CI_comparison,col = pal2(3),main=paste0(maps[i]," AGBD ", eval(parse(text=paste0(maps[i],"_year"))), ' within 95% CI of ', maps[coarsest_res[1]]," AGBD ", eval(parse(text=paste0(maps[coarsest_res[1]],"_year")))),legend=FALSE)
            legend("topright", inset=c(-0.4, 0), xpd=TRUE, legend=c("FALSE", "TRUE", "NoData"), fill=c("red", "green", "gray"),bty="n",border=c("black", "black", "black"), cex=1)
            if (file.exists(output_warps)){file.remove(output_warps)}
        }
    }
    if (file.exists(base_map_name)){file.remove(base_map_name)}
    file.remove(dir(file.path(Out_folder,paste0(iso3,'_Product_Comparisons'))), pattern = "^.xml$", full.names = TRUE)
    mydir <- file.path(Out_folder,paste0(iso3,'_Product_Comparisons'))
    delfiles <- dir(path=mydir ,pattern="*.aux.xml")
    file.remove(file.path(mydir, delfiles))
    dev.off()
}


############## Generate maximum of 99th percentile of all AGB product values ###################################
Percentile_xlim <- function(maps) {
   Save_Perc999 <- 100 # default upper limit of Mg/ha for maps 
    for (i in 1:length(maps)){
        map_values <- values(eval(parse(text=paste0(tolower(maps[i]),"_AOImasked"))))
        map_values <- map_values[map_values > 0]
        Perc999 <- quantile(map_values, probs = c(.999),na.rm = TRUE)
        if (is.na(Perc999) | !exists("Perc999")) {Perc999 <- 5}
        if (Perc999 > Save_Perc999){Save_Perc999 <- Perc999}
    }
    Save_Perc999 <- round_any(Save_Perc999, 50, f = ceiling)  
    return(Save_Perc999)
}

############## Generate a list of resolutions (pixel sizes) of various maps ###################################
Resolution_list <- function(map_list){
    Reses <- ''
    for (i in 1:length(maps)){
        resolution <- res(eval(parse(text=tolower(maps[i]))))
        Reses[[length(Reses) + 1]] <- resolution[[1]]
    }
    Reses <- Reses[Reses != ""]
    return(Reses)
}

############## Find the minimum of resolutions (pixel sizes) of various maps ###################################
F_res <- function(map_list){
    Fine_res <- 0.1
    for (i in 1:length(maps)){
        resolution <- res(eval(parse(text=tolower(maps[i]))))
        if (resolution < Fine_res){Fine_res <- resolution}
    }
    return(Fine_res)
}

############## Find the maximum of resolutions (pixel sizes) of various maps ###################################
coarsest_resolution <- function(maps_AOImasked){
    coarsest_res <- 0.0000000000000001
    for (i in 1:length(maps_AOImasked)){
        resolution <- res(eval(parse(text=maps_AOImasked[i])))
        coarsest_resolution <- max(c(resolution[[1]],resolution[[2]]))
        if (coarsest_resolution > coarsest_res){
            coarsest_res <- coarsest_resolution
            save_i <- i}
    }
    return(c(save_i,coarsest_res))
}

#############################################################################################################################
#### Generate 10x10 degree WGS84 grid that the CCI Biomass and Hansen Forest Cover products are tiled on
#### This needs to be made more generic, so that any grid size can be generated? 
Grid_set <- function(maps,verbose=FALSE){
    if (maps == 'CCI'){
        gt <- GridTopology(c(-180 + 5, -60 + 5), c(10, 10), c(36,14))
        grd <- SpatialGrid(gt, proj4string="+init=epsg:4326")
        spix <- as(grd, "SpatialPixels")
        spol <- as(spix, "SpatialPolygons")
        return(spol)
    }
}

#############################################################################################################################
#############################################################################################################################
#### RESAMPLE and APPLY HANSEN FOREST COVER MAP ##################
Apply_HANSEN_TC2020_MASK <- function(tiles,EMask_folder,maps,maps_AOImasked,HOut_folder,Out_folder,non_forest_value,iso3,aoi){
    maps_AOImasked_Hansen_masked <- ''
    file_prefixs <- ''
    for (n in 1:length(tiles)) {
        tile <- tiles[n]
        min_x <- bbox(tile)[1, 1]
        max_y <- bbox(tile)[2, 2]
        if (min_x < 0) {
          min_x <- paste0(sprintf('%03i', abs(min_x)),'W')
        } else {
          min_x <- paste0(sprintf('%03i', min_x),'E')
        }
        if (max_y < 0) {
          max_y <- paste0(sprintf('%02i', abs(max_y)),'S')
        } else {
          max_y <- paste0(sprintf('%02i', max_y),'N')
        }
        file_prefix <- paste0(EMask_folder,'/','Hansen_GFC-2020-v1.8_treecover2020_', max_y,'_', min_x,'.tif')
        file_prefixs[n] <- file_prefix
    }

    file_prefixs <- file_prefixs[file_prefixs != ""]
    HANSENlist <- file_prefixs

    for (i in 1:length(maps_AOImasked)) { 
        RESOLUTION_maps_AOImasked <- res(eval(parse(text=maps_AOImasked[i])))
        dims <- dim(eval(parse(text=maps_AOImasked[i])))
        dim1 <- dims[[1]]
        dim2 <- dims[[2]]
        extents <- ext(eval(parse(text=maps_AOImasked[i])))
        ex1 <- extents@xmin
        ex2 <- extents@ymin
        ex3 <- extents@xmax
        ex4 <- extents@ymax
        HANSENlist_RES_files <- ''
        for (each_file in HANSENlist) {
                WARP_name <- file.path(HOut_folder,paste0(strsplit(basename(each_file),'.tif')[[1]][1],"_FMASK_",toString(non_forest_value),"_RESAMPLED_",maps[i],".vrt"))
                RES_name <- file.path(HOut_folder,paste0(strsplit(basename(each_file),'.tif')[[1]][1],"_FMASK_",toString(non_forest_value),"_",maps[i],".tif"))
                HANSENlist_RES_files[[length(HANSENlist_RES_files) + 1]] <- RES_name
                if (!file.exists(RES_name)){
                    print("HANSEN layers don't exists...this will take a while...hang in there!")
                    gdalwarp(each_file,dstfile=WARP_name,tr=c(as.double(RESOLUTION_maps_AOImasked[[1]]),as.double(RESOLUTION_maps_AOImasked[[2]])),r='average',t_srs="+proj=longlat +datum=WGS84 +no_defs",overwrite=TRUE,VERBOSE=FALSE)
                    gdal_translate(WARP_name,sub('.vrt','.tif',WARP_name),co=c("BIGTIFF=IF_NEEDED","COMPRESSION=LZW","TILED=YES"),overwrite=TRUE)
                    m <- c(-Inf, non_forest_value, 0,  non_forest_value, Inf, 1)
                    rclmat <- matrix(m, ncol=3, byrow=TRUE)
                    FMASK_class <- rast(sub('.vrt','.tif',WARP_name))
                    FMASK_class <- reclassify(FMASK_class, rclmat, include.lowest=TRUE)
                    RES_file <- writerast(FMASK_class, filename=RES_name,datatype="LOG1S",overwrite=TRUE)
                }
                if (file.exists(RES_name) && file.exists(WARP_name)){file.remove(WARP_name)}
                if (file.exists(RES_name) && file.exists(sub('.vrt','.tif',WARP_name))){file.remove(sub('.vrt','.tif',WARP_name))}
        }
        HANSENlist_RES_files <- HANSENlist_RES_files[HANSENlist_RES_files != ""]

        HANSEN_country_file <- file.path(HOut_folder,paste0("HANSEN_TC2020_",maps[i],"_",iso3,".tif"))
        HANSEN_country_file_AOI <- file.path(HOut_folder,paste0("HANSEN_TC2020_",maps[i],"_AOI_",iso3,".tif"))
        Map_with_Hansen_mask <- file.path(Out_folder,paste0(iso3,'_Product_FMASK'),paste0(maps[i],"_HANSEN_FMASK_",iso3,".tif"))
        AOI_masked_country_file <- rast(file.path(Out_folder,paste0(iso3,'_Product_AOI'),paste0(maps[i],"_AOImasked.tif")))
        if (!file.exists(sub('.tif','_warps.tif',HANSEN_country_file)) || !file.exists(Map_with_Hansen_mask)){
            AOI_masked_country_file <- rast(file.path(Out_folder,paste0(iso3,'_Product_AOI'),paste0(maps[i],"_AOImasked.tif")))
            Hansen_mosaic_CMD <- mosaic_rasters(gdalfile=HANSENlist_RES_files,dst_dataset=HANSEN_country_file,datatype="INT1U", format="GTiff",r='nearest',options=c("BIGTIFF=IF_NEEDED","COMPRESSION=LZW","TILED=YES"), overwrite=TRUE,VERBOSE=FALSE) #,tr=c(as.double(RESOLUTION_maps_AOImasked[[1]]),as.double(RESOLUTION_maps_AOImasked[[2]]))
            mosaic_tif <- rast(HANSEN_country_file)          
                            gdal_translate(HANSEN_country_file,sub('.tif','_trans.tif',HANSEN_country_file),co=c("BIGTIFF=IF_NEEDED","COMPRESSION=LZW","TILED=YES"),tr=c(as.double(RESOLUTION_maps_AOImasked[[1]]),as.double(RESOLUTION_maps_AOImasked[[2]])),overwrite=TRUE)
            mosaic_tif_cropped_to_AGBmap <- crop(rast(sub('.tif','_trans.tif',HANSEN_country_file)),ext(AOI_masked_country_file),snap='near')
            gdalwarp(sub('.tif','_trans.tif',HANSEN_country_file),dstfile=sub('.tif','_warps.tif',HANSEN_country_file),ts=c(dim2,dim1),te=c(ex1,ex2,ex3,ex4),r='near',t_srs="+proj=longlat +datum=WGS84 +no_defs",overwrite=TRUE,VERBOSE=FALSE)
            Hansen_ready <- rast(sub('.tif','_warps.tif',HANSEN_country_file))
            map_Fmasked <- mask(AOI_masked_country_file, Hansen_ready, filename="", inverse=FALSE, maskvalue=0, updatevalue=NA, updateNA=TRUE)
            writerast(map_Fmasked, filename=Map_with_Hansen_mask,datatype="FLT4S", format="GTiff", options="COMPRESS=LZW", overwrite=TRUE)
            if (file.exists(sub('.tif','_warps.tif',HANSEN_country_file))){file.remove(sub('.tif','_trans.tif',HANSEN_country_file))}
            if (file.exists(sub('.tif','_warps.tif',HANSEN_country_file))){file.remove(HANSEN_country_file)}
            if (file.exists(sub('.tif','_warps.tif',HANSEN_country_file))){file.rename(sub('.tif','_warps.tif',HANSEN_country_file),HANSEN_country_file)}
        }
        maps_AOImasked_Hansen_masked[[length(maps_AOImasked_Hansen_masked) + 1]] <- Map_with_Hansen_mask
    }
    maps_AOImasked_Hansen_masked <- maps_AOImasked_Hansen_masked[maps_AOImasked_Hansen_masked != ""]
    return(maps_AOImasked_Hansen_masked)
}


#############################################################################################################################
#############################################################################################################################
#### Read in biomass products, except CCI, which is included here but needs another fuction to actually work (!!)
readAGBmap <- function(AOI_file,aoi,continent,map,In_folder,In_year,cci_version,Out_folder,verbose=FALSE){
    each_map <- map
    Band_name <- c(paste0('agbd_',In_year))
    if (each_map == 'ISB'){
        ISB_S3 <- readOGR(file.path(Out_folder,paste0(iso3,'_icesat_json.shp')))
        ISB_list <- ISB_S3$s3_path
        for (ISB_tile in ISB_list){
            ISB_tile_resampled <- paste0(strsplit(basename(ISB_tile), ".tif")[1],"_resampled.tif")
            if (!file.exists(file.path(Out_folder,paste0(iso3,'_Product_ISB'),ISB_tile_resampled))){
                save_object(ISB_tile,file=file.path(Out_folder,paste0(iso3,'_Product_ISB'),basename(ISB_tile)))
                gdalwarp(file.path(Out_folder,paste0(iso3,'_Product_ISB'),basename(ISB_tile)),file.path(Out_folder,paste0(iso3,'_Product_ISB'),ISB_tile_resampled),t_srs="+proj=longlat +datum=WGS84 +no_defs",ot="Float32", co=c("COMPRESS=DEFLATE","PREDICTOR=2","ZLEVEL=9"), overwrite=TRUE, VERBOSE=FALSE, tr=c(0.00176,0.00176),r='average')
                file.remove(file.path(Out_folder,paste0(iso3,'_Product_ISB'),basename(ISB_tile)))
            }
        }
        tiles <- list.files(path = file.path(Out_folder,paste0(iso3,'_Product_ISB')),pattern="*_resampled.tif",full.names=TRUE)
        mosaic_rasters(gdalfile=tiles,dst_dataset=file.path(Out_folder,paste0(iso3,'_Product_ISB'),"ISB_AOImasked.tif"),of="GTiff", gdalwarp_params = list(r = "average",ot="Float32"), co=c("COMPRESS=DEFLATE","PREDICTOR=2","ZLEVEL=9"), overwrite=TRUE, VERBOSE=FALSE)
        tile_stack <- mask(rast(file.path(Out_folder,paste0(iso3,'_Product_ISB'),"ISB_AOImasked.tif")), aoi)
        return(tile_stack)
    } else {
        if (each_map == 'CCI'){
            cci <- readCCI_v2(data_folder=In_folder, year=In_year, version=cci_version, aoi=aoi)
        } else {
            if (each_map == 'GEDI'){
                file_name <- paste0('GEDI04_B_MW019MW223_02_002_02_R01000M_MU.tif') ###########
                GEDI_NAflag <- -1
                aoi <- aoi %>% st_transform(crs = 6933)
            }
            if (each_map == 'GEDI_SE'){
                file_name <- paste0('GEDI04_B_MW019MW138_02_002_05_R01000M_SE.tif') ###########
                GEDI_SE_NAflag<- -1
                aoi <- aoi %>% st_transform(crs = 6933)
            }
            if (each_map == 'JPL'){
                file_name <- paste0('global_agb_mean_2020_20231024.tiff')
                JPL_NAflag<- -1
            }
            if (each_map == 'NCEO'){
                file_name <- 'AGB_map_2017v0m.tif'
                NCEO_NAflag<- -1
            } 
            if (each_map == 'MEX_AGB'){
                file_name <- 'Biomasa_mosaico_por_formacion_WGS84.tif'
                MEX_AGB_NAflag<- -1
            } 
            skip_to_next <- FALSE
            tryCatch(crop(rast(c(file.path(In_folder,file_name))), aoi, mask=TRUE), error = function(e) { skip_to_next <<- TRUE})
            if(skip_to_next) { 
                tile_stack <- rast(matrix(5:5, 1))
            } else {tile_stack <- crop(rast(c(file.path(In_folder,file_name))), aoi, mask=TRUE)
            if (!(crs(tile_stack)) == "epsg:4326"){tile_stack <- tile_stack %>% project("epsg:4326")}
            names(tile_stack) <- Band_name
            NAflag(tile_stack) <- eval(parse(text=paste0(each_map,'_NAflag')))}
            return(tile_stack)
        }
    }
}

#############################################################################################################################
#############################################################################################################################
#### Read in JPL products which are split over many continents!! 
fixJPLfile <- function(aoi,verbose=FALSE){
    tile_stacks <- c()
    All_links <- c("/projects/my-public-bucket/Data/Biomass_maps/JPL2020/Asia_clip_SE_Oceania.tif","/projects/shared-buckets/alanxuliang/AGB_2020_latest/global_008_06dc_agb_mean_prediction_2020_mosaic_veg_gfccorr_scale1_Africa_cog.tif","/projects/shared-buckets/alanxuliang/AGB_2020_latest/global_008_06dc_agb_mean_prediction_2020_mosaic_veg_gfccorr_scale1_Asia_cog.tif","/projects/shared-buckets/alanxuliang/AGB_2020_latest/global_008_06dc_agb_mean_prediction_2020_mosaic_veg_gfccorr_scale1_Europe_cog.tif","/projects/shared-buckets/alanxuliang/AGB_2020_latest/global_008_06dc_agb_mean_prediction_2020_mosaic_veg_gfccorr_scale1_NAmerica_cog.tif","/projects/shared-buckets/alanxuliang/AGB_2020_latest/global_008_06dc_agb_mean_prediction_2020_mosaic_veg_gfccorr_scale1_Oceania_cog.tif","/projects/shared-buckets/alanxuliang/AGB_2020_latest/global_008_06dc_agb_mean_prediction_2020_mosaic_veg_gfccorr_scale1_SAmerica_cog.tif")
    for (i in 1:length(All_links)){    
        skip_to_next <- FALSE
        tryCatch(crop(stack(All_links[i]), aoi), error = function(e) { skip_to_next <<- TRUE})
        if(skip_to_next) { next }  
        else{
            clip <- crop(stack(All_links[i]), aoi, datatype="FLT4S", format='GTiff', options="COMPRESS=LZW")
            writerast(clip,file.path(Out_folder,paste0(iso3,'_Product_AOI'),paste0("Temp_",i,".tif")) ,datatype="FLT4S", format="GTiff", options=c("BIGTIFF=IF_NEEDED","COMPRESSION=LZW","TILED=YES"), overwrite=TRUE)
        }
    }
    JPL_list <- list.files(path=file.path(Out_folder,paste0(iso3,'_Product_AOI/')), pattern="Temp_",full.names = TRUE)
    mosaic_rasters(gdalfile=JPL_list,dst_dataset=file.path(Out_folder,paste0(iso3,'_Product_AOI'),paste0("JPL_AOImasked.tif")),of="GTiff")
    JPL_file <-rast(file.path(Out_folder,paste0(iso3,'_Product_AOI'),paste0("JPL_AOImasked.tif")))
    JPL_file <-(mask(JPL_file, aoi, updatevalue=NA))
    delfiles <- dir(path=file.path(Out_folder,paste0(iso3,'_Product_AOI/')) ,pattern="Temp*")
    file.remove(file.path(Out_folder,paste0(iso3,'_Product_AOI/'), delfiles))
    return(JPL_file)
}

#############################################################################################################################
#############################################################################################################################
#### What does this function do? ##################################
mosaic_list <- function(x, fun, datatype, format, options, overwrite, 
                        tolerance=0.05, filename="") {
  mosaic_args <- x
  if (!missing(fun)) mosaic_args$fun <- fun
  if (!missing(tolerance)) mosaic_args$tolerance <- tolerance
  if (!missing(datatype)) mosaic_args$datatype <- datatype
  if (!missing(format)) mosaic_args$format <- format
  if (!missing(options)) mosaic_args$options <- options
  if (!missing(overwrite)) mosaic_args$overwrite <- overwrite
  mosaic_args$filename <- filename
  do.call(mosaic, mosaic_args)
}

#############################################################################################################################
#############################################################################################################################
#### Read in CCI biomass within a given area of interest. Multiple years can be provided as a vector. 
#### Combines AGBD and SD into as layers in one raster mosaic
#### modified functions added on by Vero (Oct-27-2021) mainly to correct for CCI values in PER notebook
readCCI_v2 <- function(data_folder, year=CCI_year, version=CCI_version, aoi, verbose=FALSE){
  year <- as.character(year)
  if(version==4){
    #v4.0 file format
    file_root <- paste0('_ESACCI-BIOMASS-L4-AGB-MERGED-100m-', year, '-fv4.0.tif')
    file_root_sd <- paste0('_ESACCI-BIOMASS-L4-AGB_SD-MERGED-100m-', year, '-fv4.0.tif')
  }
  if(version==3){
    #v3.0 file format
    file_root <- paste0('_ESACCI-BIOMASS-L4-AGB-MERGED-100m-', year, '-fv3.0.tif')
    file_root_sd <- paste0('_ESACCI-BIOMASS-L4-AGB_SD-MERGED-100m-', year, '-fv3.0.tif')
  }
  if(version ==2){
    #v2.0 file format
    file_root <- paste0('_ESACCI-BIOMASS-L4-AGB-MERGED-100m-', year, '-fv2.0.tif')
    file_root_sd <- paste0('_ESACCI-BIOMASS-L4-AGB_SD-MERGED-100m-', year, '-fv2.0.tif')
  }  
  filename=""
  tile_stacks <- c()
  for (n in 1:length(tiles)) {
    tile <- tiles[n]
    min_x <- bbox(tile)[1, 1]
    max_y <- bbox(tile)[2, 2]
    if (min_x < 0) {
      min_x <- paste0('W', sprintf('%03i', abs(min_x)))
    } else {
      min_x <- paste0('E', sprintf('%03i', min_x))
    }
    if (max_y < 0) {
      max_y <- paste0('S', sprintf('%02i', abs(max_y)))
    } else {
      max_y <- paste0('N', sprintf('%02i', max_y))
    }
    file_prefix <- paste0(max_y, min_x)
    band_names2 <- c(paste0('agbd_', year))
    filenames2 <- c(file.path(data_folder, paste0(file_prefix, file_root)))
    if(verbose){
        print(filenames2)
    }
    if (!file.exists(filenames2[1])){
      print(paste0(filenames2, " does not exist"))
      next()
    }    
    tile_stack <- terra::crop(rast(filenames2), aoi, datatype="FLT4S")
    names(tile_stack) <- band_names2
    tile_stacks <- c(tile_stacks, list(tile_stack))
  }
  #If there are multiple tiles, create mosaic. Otherwise mosaic only includes 1 tile
  if (length(tile_stacks) > 1) {
    print(class(tile_stacks))
    tile_mosaic <- terra::mosaic(sprc(tile_stacks), fun='mean',filename=filename)
  } else {
    tile_mosaic <- tile_stacks[[1]]
    if (filename != '') {
      tile_mosaic <- writerast(tile_mosaic, filename=filename)
    }
  }
  names(tile_mosaic) <- band_names2
  NAflag(tile_mosaic) <- -1
  cci <- tile_mosaic
  return(cci)
} 
########################################################################################################
########################################################################################################