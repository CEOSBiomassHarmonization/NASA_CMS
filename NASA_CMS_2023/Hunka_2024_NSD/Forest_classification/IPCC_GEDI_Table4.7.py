import os
import numpy as np
import sys
import rasterio
import argparse
from rasterio.io import MemoryFile
from rasterio.rio import options
from rasterio.shutil import copy, delete
from rasterio.vrt import WarpedVRT
from rio_cogeo.cogeo import cog_translate
from rio_cogeo.profiles import cog_profiles
from memory_profiler import profile
import copy

####################### TEST RUN ON TERMINAL ON NASA-ESA MAAP ######################
##### INSTALL PACKAGES #####
# bash /projects/my-public-bucket/_ENV/env/build_above_env.sh

##### DOWNLOAD FILES #####
# wget https://glad.umd.edu/users/Potapov/GLCLUC2020/Forest_height_loss_2000_2020/netloss_10S_060W.tif -O /projects/my-public-bucket/Data/Harris_et_al_PAPER/_TEST/netloss_10S_060W.tif

# wget https://glad.umd.edu/users/Potapov/GLCLUC2020/Forest_height_2020/2020_10S_060W.tif -O /projects/my-public-bucket/Data/Harris_et_al_PAPER/_TEST/2020_10S_060W.tif

# wget https://glad.umd.edu/users/Potapov/GLCLUC2020/Forest_height_2000/2000_10S_060W.tif -O /projects/my-public-bucket/Data/Harris_et_al_PAPER/_TEST/2000_10S_060W.tif

# wget https://glad.umd.edu/Potapov/TCC_2010/treecover2010_10S_060W.tif -O /projects/my-public-bucket/Data/Harris_et_al_PAPER/_TEST/treecover2010_10S_060W.tif

# wget https://storage.googleapis.com/earthenginepartners-hansen/GFC-2022-v1.10/Hansen_GFC-2022-v1.10_lossyear_10S_060W.tif -O /projects/my-public-bucket/Data/Harris_et_al_PAPER/_TEST/Hansen_GFC-2022-v1.10_lossyear_10S_060W.tif

##### RUN SCRIPT FROM COMMAND LINE #####
# python /projects/ADE_biomass_harmonization/biomass_harmonization/country_summaries/IPCC_classes_DPS/IPCC_GEDI_Table4.7.py --FOREST_HEIGHT_2000 /projects/my-public-bucket/Data/Harris_et_al_PAPER/_TEST/2000_10S_060W.tif --FOREST_HEIGHT_2020 /projects/my-public-bucket/Data/Harris_et_al_PAPER/_TEST/2020_10S_060W.tif --FOREST_LOSS /projects/my-public-bucket/Data/Harris_et_al_PAPER/_TEST/netloss_10S_060W.tif --ESA_FOREST_COVER /projects/my-public-bucket/Data/Harris_et_al_PAPER/ESA_WorldCover2021_v2/WC_10S_060W_C.tif  --FOREST_COVER_LOSSYEAR /projects/my-public-bucket/Data/Harris_et_al_PAPER/_TEST/Hansen_GFC-2022-v1.10_lossyear_10S_060W.tif --BOREAL_AGE /projects/my-public-bucket/Data/Harris_et_al_PAPER/Forest_age_boreal/BOREAL_AGE/10S_060W_BOREAL_AGE.tif --PLANTATIONS /projects/my-public-bucket/Data/Harris_et_al_PAPER/WRI_Planted_Forest_Type_v2/10S_060W_plantation_type_oilpalm_woodfiber_other.tif --IFL_rasters /projects/my-public-bucket/Data/Harris_et_al_PAPER/Intact_Forest_Landscapes/IFL_rasters_2020/10S_060W_IFL_2020.tif --FII_raster /projects/my-public-bucket/Data/Harris_et_al_PAPER/Forest_Integrity_Index/10S_060W_FII.tif --primary_forest_asia /projects/my-public-bucket/Data/Harris_et_al_PAPER/Turubanova_Humid_tropical_primary_forest/10S_060W_Asia.tif --primary_forest_sa /projects/my-public-bucket/Data/Harris_et_al_PAPER/Turubanova_Humid_tropical_primary_forest/10S_060W_SouthAmerica.tif  --primary_forest_africa /projects/my-public-bucket/Data/Harris_et_al_PAPER/Turubanova_Humid_tropical_primary_forest/10S_060W_Africa.tif --EcoCont /projects/my-public-bucket/Data/Harris_et_al_PAPER/EcoCont_v2/10S_060W_EcoCont.tif --tile 10S_060W --output_mask /projects/my-public-bucket/Data/Harris_et_al_PAPER/_TEST/Ages_10S_060W.tif --output_file /projects/my-public-bucket/Data/Harris_et_al_PAPER/_TEST/Classes_10S_060W.tif
######################################################################################

parse = argparse.ArgumentParser(description="Classifies world into GEZ, Continent and Forest Age classes")
parse.add_argument("--FOREST_HEIGHT_2000", help="FOREST_HEIGHT_2000")
parse.add_argument("--FOREST_HEIGHT_2020", help="FOREST_HEIGHT_2020")
parse.add_argument("--FOREST_LOSS", help="FOREST_LOSS")
parse.add_argument("--ESA_FOREST_COVER", help="ESA_FOREST_COVER")
parse.add_argument("--FOREST_COVER_LOSSYEAR", help="FOREST_COVER_LOSSYEAR")
parse.add_argument("--JRC_TMST",help="JRC_TMST")
parse.add_argument("--JRC_DEF",help="JRC_DEF")
parse.add_argument("--JRC_DEG",help="JRC_DEG")
parse.add_argument("--BOREAL_AGE", help="BOREAL_AGE")
parse.add_argument("--PLANTATIONS", help="PLANTATIONS")
parse.add_argument("--GFM", help="GFM")
parse.add_argument("--IFL_rasters", help="IFL_rasters")
parse.add_argument("--FII_raster", help="FII_raster")
parse.add_argument("--primary_forest_asia", help="primary_forest_asia")
parse.add_argument("--primary_forest_sa", help="primary_forest_sa")
parse.add_argument("--primary_forest_africa", help="primary_forest_Africa")
parse.add_argument("--EcoCont", help="EcoCont")
parse.add_argument("--tile", help="Tile")
parse.add_argument("--output_mask", help="output_mask")
parse.add_argument("--output_file", help="Output")
args = parse.parse_args()
    
TILE = args.tile

###################################################################################
###################################################################################
##### PRIMARY AND INTACT FOREST CONDITIONS - GIVEN VALUE 1 ########################
###################################################################################
###################################################################################

@profile

def generate_output():
    PRIMARY = rasterio.open(args.FII_raster).read(1)
    PRIMARY[PRIMARY <= 9600] = 0
    PRIMARY[PRIMARY > 9600] = 1
    
    with rasterio.open(args.IFL_rasters) as src:
        INTACT_FORESTS = src.read(1) # GLAD SUGGESTS TO NOT USE THIS AT ALL 
        PRIMARY[INTACT_FORESTS == 1] = 1
        del(INTACT_FORESTS)
        
    with rasterio.open(args.primary_forest_asia) as src:
        PRIMARY_FOREST_ASIA = src.read(1)
        PRIMARY[PRIMARY_FOREST_ASIA > 0] = 1
        del(PRIMARY_FOREST_ASIA)
        
    with rasterio.open(args.primary_forest_sa) as src:
        PRIMARY_FOREST_SA = src.read(1)
        PRIMARY[PRIMARY_FOREST_SA > 0] = 1
        del(PRIMARY_FOREST_SA)
    
    with rasterio.open(args.primary_forest_africa) as src:
        PRIMARY_FOREST_AFRICA = src.read(1)
        PRIMARY[PRIMARY_FOREST_AFRICA > 0] = 1
        del(PRIMARY_FOREST_AFRICA)
                   
    ##### JRC TRANSITION MAP #####
    JRC_TMST = rasterio.open(args.JRC_TMST).read(1)
    PRIMARY[JRC_TMST == 10] = 1
    del(JRC_TMST)

    ##### FOREST HEIGHT AND HEIGHT LOSS CONDITIONS #####
    
    with rasterio.open(args.FOREST_HEIGHT_2000) as src:
        FOREST_HEIGHT_2000 = src.read(1)
        PRIMARY[FOREST_HEIGHT_2000 < 5] = 0
        del(FOREST_HEIGHT_2000)
        
    with rasterio.open(args.FOREST_HEIGHT_2020) as src:
        FOREST_HEIGHT_2020 = src.read(1)
        PRIMARY[FOREST_HEIGHT_2020 < 5] = 0
        del(FOREST_HEIGHT_2020)

    ##### FOREST TREE COVER LOSS CONDITIONS #####
    FOREST_COVER_LOSSYEAR = rasterio.open(args.FOREST_COVER_LOSSYEAR).read(1)
    PRIMARY[FOREST_COVER_LOSSYEAR > 0] = 0
    del(FOREST_COVER_LOSSYEAR)
    
    ##### JRC DEFORESTATION CONDITIONS #####
    JRC_DEF = rasterio.open(args.JRC_DEF).read(1)
    PRIMARY[JRC_DEF > 1981] = 0
    del(JRC_DEF)
    
    ##### JRC DEGRADATION CONDITIONS #####
    JRC_DEG = rasterio.open(args.JRC_DEG).read(1)
    PRIMARY[JRC_DEG > 1981] = 0
    del(JRC_DEG)

    ##### PLANTATIONS CONDITIONS #####
    PLANTATIONS = rasterio.open(args.PLANTATIONS).read(1)
    PRIMARY[PLANTATIONS > 0] = 0
    del(PLANTATIONS)

    JRC_TMST = rasterio.open(args.JRC_TMST).read(1)
    PRIMARY[(JRC_TMST >= 81) & (JRC_TMST <= 86)] = 0
    del(JRC_TMST)

    GFM = rasterio.open(args.GFM).read(1)
    PRIMARY[(GFM == 32) | (GFM == 40) | (GFM == 53)] = 0
    del(GFM)

    ##### ESA WORLD COVER CONDITION #####
    ESA_FOREST_COVER = rasterio.open(args.ESA_FOREST_COVER).read(1)
    PRIMARY[ESA_FOREST_COVER == 0] = 0
    del(ESA_FOREST_COVER)

    # ###################################################################################
    # ###################################################################################
    # ##### YOUNG FOREST CONDITIONS - GIVEN VALUE 2 #####################################
    # ###################################################################################
    # ###################################################################################
    FOREST_HEIGHT_2000 = rasterio.open(args.FOREST_HEIGHT_2000).read(1)
    FOREST_HEIGHT_2020 = rasterio.open(args.FOREST_HEIGHT_2020).read(1)
    BACKUP_FOREST_HEIGHT_2020 = FOREST_HEIGHT_2020
    FOREST_COVER_LOSSYEAR = rasterio.open(args.FOREST_COVER_LOSSYEAR).read(1)
    
    ##########################################################
    ##### FOREST HEIGHT CONDITIONS (Potapov et al. 2022) #####
    YOUNG = FOREST_HEIGHT_2020
    del(FOREST_HEIGHT_2020)
    YOUNG[BACKUP_FOREST_HEIGHT_2020 < 5] = 0
    YOUNG[(BACKUP_FOREST_HEIGHT_2020 >= 5) & (FOREST_HEIGHT_2000 < 5)] = 2
    YOUNG[(BACKUP_FOREST_HEIGHT_2020 >= 5) & (FOREST_HEIGHT_2000 >= 5) & (FOREST_COVER_LOSSYEAR > 0) & (FOREST_COVER_LOSSYEAR < 18)] = 2
    del(FOREST_COVER_LOSSYEAR)
    JRC_DEF = rasterio.open(args.JRC_DEF).read(1)
    JRC_DEG = rasterio.open(args.JRC_DEG).read(1)
    YOUNG[(BACKUP_FOREST_HEIGHT_2020 >= 5) & (FOREST_HEIGHT_2000 >= 5) & (JRC_DEF > 2000) & (JRC_DEF < 2018)] = 2
    YOUNG[(BACKUP_FOREST_HEIGHT_2020 >= 5) & (FOREST_HEIGHT_2000 >= 5) & (JRC_DEG > 2000) & (JRC_DEG < 2018)] = 2
    del(BACKUP_FOREST_HEIGHT_2020,FOREST_HEIGHT_2000,JRC_DEF,JRC_DEG)
    
    YOUNG[(YOUNG > 2) | (YOUNG < 2)] = 0

    ##### EXCEPTIONS FOR THE BOREAL #####
    BOREAL_AGE = rasterio.open(args.BOREAL_AGE).read(1)
    YOUNG[(BOREAL_AGE <= 20) & (BOREAL_AGE > 0)] = 2 # LAYER ACTIVE ONLY IN BOREAL ZONES
    del(BOREAL_AGE)
                   
    ##### JRC TRANSITION MAP #####
    JRC_TMST = rasterio.open(args.JRC_TMST).read(1)
    YOUNG[JRC_TMST == 32] = 2
    YOUNG[JRC_TMST == 33] = 2
    del(JRC_TMST)

    ##### PLANTATIONS CONDITIONS #####
    PLANTATIONS = rasterio.open(args.PLANTATIONS).read(1)
    YOUNG[PLANTATIONS > 0] = 0
    del(PLANTATIONS)

    JRC_TMST = rasterio.open(args.JRC_TMST).read(1)
    YOUNG[(JRC_TMST >= 81) & (JRC_TMST <= 86)] = 0
    del(JRC_TMST)

    GFM = rasterio.open(args.GFM).read(1)
    YOUNG[(GFM == 32) | (GFM == 40) | (GFM == 53)] = 0
    del(GFM)

    ##### ESA WORLD COVER CONDITION #####
    ESA_FOREST_COVER = rasterio.open(args.ESA_FOREST_COVER).read(1)
    YOUNG[ESA_FOREST_COVER == 0] = 0
    del(ESA_FOREST_COVER)

    ##### PRIMARY FOREST CONDITION #####
    YOUNG[PRIMARY > 0] = 0

    ###################################################################################
    ###################################################################################
    ##### OLD FOREST CONDITIONS - GIVEN VALUE 3 #######################################
    ###################################################################################
    ###################################################################################

    ##### FOREST HEIGHT CONDITIONS #####
    FOREST_HEIGHT_2000 = rasterio.open(args.FOREST_HEIGHT_2000).read(1)
    FOREST_HEIGHT_2020 = rasterio.open(args.FOREST_HEIGHT_2020).read(1)
    
    OLD = FOREST_HEIGHT_2020
    del(FOREST_HEIGHT_2020)
    OLD[OLD < 5] = 0
    OLD[(OLD >= 5) & (FOREST_HEIGHT_2000 < 5)] = 0
    OLD[(OLD >= 5) & (FOREST_HEIGHT_2000 >= 5)] = 3
    del(FOREST_HEIGHT_2000)

    ##### EXCEPTIONS FOR THE BOREAL #####
    BOREAL_AGE = rasterio.open(args.BOREAL_AGE).read(1)
    OLD[(BOREAL_AGE > 20) & (BOREAL_AGE <= 36)] = 3 # LAYER ACTIVE ONLY IN BOREAL ZONES
    del(BOREAL_AGE)
                   
    ##### JRC TRANSITION MAP #####
    JRC_TMST = rasterio.open(args.JRC_TMST).read(1)
    OLD[JRC_TMST == 31] = 3
    del(JRC_TMST)

    ##### FOREST TREE COVER LOSS CONDITIONS #####
    FOREST_COVER_LOSSYEAR = rasterio.open(args.FOREST_COVER_LOSSYEAR).read(1)
    OLD[FOREST_COVER_LOSSYEAR > 0] = 0
    del(FOREST_COVER_LOSSYEAR)
    
    JRC_DEF = rasterio.open(args.JRC_DEF).read(1)
    OLD[JRC_DEF > 2000] = 0
    del(JRC_DEF)
    
    JRC_DEG = rasterio.open(args.JRC_DEG).read(1)
    OLD[JRC_DEG > 2000] = 0
    del(JRC_DEG)

    ##### PLANTATIONS CONDITIONS #####
    PLANTATIONS = rasterio.open(args.PLANTATIONS).read(1)
    OLD[PLANTATIONS > 0] = 0
    del(PLANTATIONS)

    JRC_TMST = rasterio.open(args.JRC_TMST).read(1)
    OLD[(JRC_TMST >= 81) & (JRC_TMST <= 86)] = 0
    del(JRC_TMST)

    GFM = rasterio.open(args.GFM).read(1)
    OLD[(GFM == 32) | (GFM == 40) | (GFM == 53)] = 0
    del(GFM)

    ##### ESA WORLD COVER CONDITION #####
    ESA_FOREST_COVER = rasterio.open(args.ESA_FOREST_COVER).read(1)
    OLD[ESA_FOREST_COVER == 0] = 0
    del(ESA_FOREST_COVER)

    ##### PRIMARY AND YOUNG SECONDARY FOREST CONDITIONS #####
    OLD[PRIMARY > 0] = 0
    OLD[YOUNG > 0] = 0

    ###################################################################################
    ###################################################################################
    ##### AGGREGATING GLOBAL FOREST CLASSES AND OUTPUT TO RASTER ######################
    ###################################################################################
    ###################################################################################

    MASK = PRIMARY + YOUNG + OLD
    del(PRIMARY, YOUNG, OLD)
    
    with rasterio.open(args.EcoCont) as src:
                
        EcoCont = src.read(1)
        kwargs = src.meta
        kwargs.update(dtype="uint16", predictor=2)
        BASE = EcoCont
        EcoCont = EcoCont + (MASK*100)
        EcoCont[EcoCont == BASE] = 0
        EcoCont[(EcoCont == 100) | (EcoCont == 200) | (EcoCont == 300)] = 0
        MASK[EcoCont == 0] = 0
        del(BASE)

        with MemoryFile() as memfile:
            with memfile.open(**kwargs) as mem:
                # Populate the input file with numpy array
                mem.write(EcoCont.reshape(1, EcoCont.shape[0], EcoCont.shape[1]))
                dst_profile = cog_profiles.get("deflate")
                cog_translate(
                    mem,
                    args.output_file,
                    dst_profile,
                    in_memory=True,
                    allow_intermediate_compression=True
                )

        with MemoryFile() as memfile1:
            with memfile1.open(**kwargs) as mem1:
                # Populate the input file with numpy array
                mem1.write(MASK.reshape(1, MASK.shape[0], MASK.shape[1]))
                dst_profile = cog_profiles.get("deflate")
                cog_translate(
                    mem1,
                    args.output_mask,
                    dst_profile,
                    in_memory=True,
                    allow_intermediate_compression=True
                )
            
if __name__ == "__main__":
    generate_output()

# ###################################################################################
# ###################################################################################
# ##### PRIMARY AND INTACT FOREST CONDITIONS - GIVEN VALUE 1 ########################
# ###################################################################################
# ###################################################################################



# FII_raster = rasterio.open(args.FII_raster).read(1)
# INTACT_FORESTS = rasterio.open(args.IFL_rasters).read(1) # GLAD SUGGESTS TO NOT USE THIS AT ALL 
# PRIMARY_FOREST_ASIA = rasterio.open(args.primary_forest_asia).read(1)
# PRIMARY_FOREST_SA = rasterio.open(args.primary_forest_sa).read(1)
# PRIMARY_FOREST_AFRICA = rasterio.open(args.primary_forest_africa).read(1)

# PRIMARY = FII_raster
# del(FII_raster)
# PRIMARY[PRIMARY <= 9600] = 0
# PRIMARY[PRIMARY > 9600] = 1
# PRIMARY[INTACT_FORESTS == 1] = 1
# del(INTACT_FORESTS)
# PRIMARY[PRIMARY_FOREST_ASIA > 0] = 1
# del(PRIMARY_FOREST_ASIA)
# PRIMARY[PRIMARY_FOREST_SA > 0] = 1
# del(PRIMARY_FOREST_SA)
# PRIMARY[PRIMARY_FOREST_AFRICA > 0] = 1
# del(PRIMARY_FOREST_AFRICA)

# ##### FOREST HEIGHT AND HEIGHT LOSS CONDITIONS #####
# FOREST_HEIGHT_2000 = rasterio.open(args.FOREST_HEIGHT_2000).read(1)
# FOREST_HEIGHT_2020 = rasterio.open(args.FOREST_HEIGHT_2020).read(1)
# PRIMARY[FOREST_HEIGHT_2000 < 5] = 0
# PRIMARY[FOREST_HEIGHT_2020 < 5] = 0 
# del(FOREST_HEIGHT_2000,FOREST_HEIGHT_2020)

# ##### FOREST TREE COVER LOSS CONDITIONS #####
# FOREST_COVER_LOSSYEAR = rasterio.open(args.FOREST_COVER_LOSSYEAR).read(1)
# PRIMARY[FOREST_COVER_LOSSYEAR > 0] = 0
# del(FOREST_COVER_LOSSYEAR)

# ##### PLANTATIONS CONDITIONS #####
# PLANTATIONS = rasterio.open(args.PLANTATIONS).read(1)
# PRIMARY[PLANTATIONS > 0] = 0
# del(PLANTATIONS)

# ##### ESA WORLD COVER CONDITION #####
# ESA_FOREST_COVER = rasterio.open(args.ESA_FOREST_COVER).read(1)
# PRIMARY[ESA_FOREST_COVER == 0] = 0
# del(ESA_FOREST_COVER)

# # ###################################################################################
# # ###################################################################################
# # ##### YOUNG FOREST CONDITIONS - GIVEN VALUE 2 #####################################
# # ###################################################################################
# # ###################################################################################

# FOREST_HEIGHT_2000 = rasterio.open(args.FOREST_HEIGHT_2000).read(1)
# FOREST_HEIGHT_2020 = rasterio.open(args.FOREST_HEIGHT_2020).read(1)
# BACKUP_FOREST_HEIGHT_2020 = rasterio.open(args.FOREST_HEIGHT_2020).read(1)
# FOREST_COVER_LOSSYEAR = rasterio.open(args.FOREST_COVER_LOSSYEAR).read(1)

# ####################################
# ##### FOREST HEIGHT CONDITIONS (Potapov et al. 2022) #####
# YOUNG = FOREST_HEIGHT_2020
# YOUNG[BACKUP_FOREST_HEIGHT_2020 < 5] = 0
# YOUNG[(BACKUP_FOREST_HEIGHT_2020 >= 5) & (FOREST_HEIGHT_2000 < 5)] = 2
# YOUNG[(BACKUP_FOREST_HEIGHT_2020 >= 5) & (FOREST_HEIGHT_2000 >= 5) & (FOREST_COVER_LOSSYEAR > 0) & (FOREST_COVER_LOSSYEAR < 18)] = 2
# YOUNG[(YOUNG > 2) | (YOUNG < 2)] = 0

# del(FOREST_COVER_LOSSYEAR,BACKUP_FOREST_HEIGHT_2020,FOREST_HEIGHT_2020,FOREST_HEIGHT_2000)

# ##### EXCEPTIONS FOR THE BOREAL #####
# BOREAL_AGE = rasterio.open(args.BOREAL_AGE).read(1)
# YOUNG[(BOREAL_AGE <= 20) & (BOREAL_AGE > 0)] = 2 # LAYER ACTIVE ONLY IN BOREAL ZONES
# del(BOREAL_AGE)

# ##### PLANTATIONS CONDITIONS #####
# PLANTATIONS = rasterio.open(args.PLANTATIONS).read(1)
# YOUNG[PLANTATIONS > 0] = 0
# del(PLANTATIONS)

# ##### ESA WORLD COVER CONDITION #####
# ESA_FOREST_COVER = rasterio.open(args.ESA_FOREST_COVER).read(1)
# YOUNG[ESA_FOREST_COVER == 0] = 0
# del(ESA_FOREST_COVER)

# ##### PRIMARY AND INTACT FOREST CONDITIONS #####
# YOUNG[PRIMARY > 0] = 0

# ###################################################################################
# ###################################################################################
# ##### OLD FOREST CONDITIONS - GIVEN VALUE 3 #######################################
# ###################################################################################
# ###################################################################################

# ##### FOREST HEIGHT CONDITIONS #####
# FOREST_HEIGHT_2000 = rasterio.open(args.FOREST_HEIGHT_2000).read(1)
# FOREST_HEIGHT_2020 = rasterio.open(args.FOREST_HEIGHT_2020).read(1)
# OLD = FOREST_HEIGHT_2020
# del(FOREST_HEIGHT_2020)
# OLD[OLD < 5] = 0
# OLD[(OLD >= 5) & (FOREST_HEIGHT_2000 < 5)] = 0
# OLD[(OLD >= 5) & (FOREST_HEIGHT_2000 >= 5)] = 3 

# # FOREST_LOSS = rasterio.open(args.FOREST_LOSS).read(1)
# # OLD[FOREST_LOSS > 0] = 0
# # del(FOREST_LOSS)

# ##### EXCEPTIONS FOR THE BOREAL #####
# BOREAL_AGE = rasterio.open(args.BOREAL_AGE).read(1)
# OLD[(BOREAL_AGE > 20) & (BOREAL_AGE <= 36)] = 3 # LAYER ACTIVE ONLY IN BOREAL ZONES
# del(BOREAL_AGE)

# ##### FOREST TREE COVER LOSS CONDITIONS #####
# FOREST_COVER_LOSSYEAR = rasterio.open(args.FOREST_COVER_LOSSYEAR).read(1)
# OLD[FOREST_COVER_LOSSYEAR > 0] = 0
# del(FOREST_COVER_LOSSYEAR)

# ##### PLANTATIONS CONDITIONS #####
# PLANTATIONS = rasterio.open(args.PLANTATIONS).read(1)
# OLD[PLANTATIONS > 0] = 0
# del(PLANTATIONS)

# ##### ESA WORLD COVER CONDITION #####
# ESA_FOREST_COVER = rasterio.open(args.ESA_FOREST_COVER).read(1)
# OLD[ESA_FOREST_COVER == 0] = 0
# del(ESA_FOREST_COVER)

# ##### PRIMARY AND INTACT FOREST CONDITIONS #####
# OLD[PRIMARY > 0] = 0
# OLD[YOUNG > 0] = 0

# ###################################################################################
# ###################################################################################
# ##### AGGREGATING GLOBAL FOREST CLASSES AND OUTPUT TO RASTER ######################
# ###################################################################################
# ###################################################################################

# MASK = PRIMARY + YOUNG + OLD
# EcoCont = rasterio.open(args.EcoCont).read(1)
# kwargs = rasterio.open(args.EcoCont).meta
# kwargs.update(dtype=rasterio.uint16,count=1,compress='lzw') 
# EcoCont = EcoCont + (MASK*100)
# BASE = rasterio.open(args.EcoCont).read(1)
# EcoCont[EcoCont == BASE] = 0
# EcoCont[(EcoCont == 100) | (EcoCont == 200) | (EcoCont == 300)] = 0
# MASK[EcoCont == 0] = 0
# # with rasterio.open(args.output_mask, 'w', **kwargs) as dst: dst.write_band(1, MASK.astype(rasterio.uint16))
# # with rasterio.open(args.output_file, 'w', **kwargs) as dst: dst.write_band(1, EcoCont.astype(rasterio.uint16))

# with MemoryFile() as memfile:
#     with memfile.open(**kwargs) as mem:
#         # Populate the input file with numpy array
#         mem.write(EcoCont.reshape(1, EcoCont.shape[0], EcoCont.shape[1]))

#         dst_profile = cog_profiles.get("deflate")
#         cog_translate(
#             mem,
#             args.output_file,
#             dst_profile,
#             in_memory=False
#         )
        
# with MemoryFile() as memfile1:
#     with memfile1.open(**kwargs) as mem1:
#         # Populate the input file with numpy array
#         mem1.write(MASK.reshape(1, MASK.shape[0], MASK.shape[1]))

#         dst_profile = cog_profiles.get("deflate")
#         cog_translate(
#             mem1,
#             args.output_mask,
#             dst_profile,
#             in_memory=False
#         )