#!/bin/bash

## Use the next line if you need custom env that differs from the container image
# source activate Packages_for_IPCC
basedir=$( cd "$(dirname "$0")" ; pwd -P )   

mkdir -p output

tile=${1}
output_mask=${2}
output_file=${3}
FOREST_HEIGHT_2000="input/2000_${tile}.tif"
FOREST_HEIGHT_2020="input/2020_${tile}.tif"
FOREST_LOSS="input/netloss_${tile}.tif"
ESA_FOREST_COVER="input/WC_${tile}_C.tif"
JRC_TMST="input/${tile}_JRC_Transition_Map.tif"
JRC_DEF="input/${tile}_JRC_Deforestation_Map.tif"
JRC_DEG="input/${tile}_JRC_Degradation_Map.tif"
JRC_FC="input/${tile}_JRC_ForestCover_Map.tif"
FOREST_COVER_LOSSYEAR="input/Hansen_GFC-2023-v1.11_lossyear_${tile}.tif"
BOREAL_AGE="input/${tile}_BOREAL_AGE.tif"
PLANTATIONS="input/C_SDPT_v21_${tile}.tif"
GFM="input/${tile}_GFM.tif"
IFL_rasters="input/${tile}_IFL_2020.tif"
FII_raster="input/${tile}_FII.tif"
primary_forest_asia="input/${tile}_Asia.tif"
primary_forest_sa="input/${tile}_SouthAmerica.tif"
primary_forest_africa="input/${tile}_Africa.tif"
primary_forest_eu="input/${tile}_Africa.tif"
EcoCont="input/PF_Europe__${tile}.tif"

OUTPUTDIR="${PWD}/output"


conda run --live-stream --name python python ${basedir}/IPCC_GEDI_Table4.7.py --tile ${tile} --output_mask output/${output_mask} --output_file output/${output_file} --FOREST_HEIGHT_2000 ${FOREST_HEIGHT_2000} --FOREST_HEIGHT_2020 ${FOREST_HEIGHT_2020} --FOREST_LOSS ${FOREST_LOSS} --ESA_FOREST_COVER ${ESA_FOREST_COVER} --JRC_TMST ${JRC_TMST} --JRC_DEF ${JRC_DEF} --JRC_DEG ${JRC_DEG} --JRC_FC ${JRC_FC} --FOREST_COVER_LOSSYEAR ${FOREST_COVER_LOSSYEAR} --BOREAL_AGE ${BOREAL_AGE} --PLANTATIONS ${PLANTATIONS} --GFM ${GFM} --IFL_rasters ${IFL_rasters} --FII_raster ${FII_raster} --primary_forest_asia ${primary_forest_asia} --primary_forest_sa ${primary_forest_sa} --primary_forest_africa ${primary_forest_africa} --primary_forest_eu ${primary_forest_eu} --EcoCont ${EcoCont}

# python ${basedir}/IPCC_GEDI_Table4.7.py --tile ${tile} --output_mask output/${output_mask} --output_file output/${output_file} --FOREST_HEIGHT_2000 ${FOREST_HEIGHT_2000} --FOREST_HEIGHT_2020 ${FOREST_HEIGHT_2020} --FOREST_LOSS ${FOREST_LOSS} --ESA_FOREST_COVER ${ESA_FOREST_COVER} --JRC_TMST ${JRC_TMST} --JRC_DEF ${JRC_DEF} --JRC_DEG ${JRC_DEG} --FOREST_COVER_LOSSYEAR ${FOREST_COVER_LOSSYEAR} --BOREAL_AGE ${BOREAL_AGE} --PLANTATIONS ${PLANTATIONS} --GFM ${GFM} --IFL_rasters ${IFL_rasters} --FII_raster ${FII_raster} --primary_forest_asia ${primary_forest_asia} --primary_forest_sa ${primary_forest_sa} --primary_forest_africa ${primary_forest_africa} --EcoCont ${EcoCont}