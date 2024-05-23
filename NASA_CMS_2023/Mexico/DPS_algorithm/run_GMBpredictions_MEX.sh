#!/bin/bash

basedir=$( cd "$(dirname "$0")" ; pwd -P )

mkdir -p output

output_csv=${1}
Projects=${2}
mexico="input/ecort08gw_DESECON1_DISS.gpkg"
ccirast="input/cci_mexico_6933.tif"
heirast="input/GLAD_FH_mexico_UINT16_C_6399.tif"
data="input/NFI_CCI_GEDIheights.csv"
INLA_model_fit_v2="input/INLA_model_fit_v2.RData"
MGN2020_INEGI_Urban_Mex_mask_DISS="input/MGN2020_INEGI_Urban_Mex_mask_DISS.gpkg"
All_Products_Comp_Over30_Binary_6933="input/All_Products_Comp_Over30_Binary_6933.tif" 

OUTPUTDIR="${PWD}/output"

conda run --name r --no-capture-output Rscript ${basedir}/GMB_Mexico.R output/${output_csv} ${Projects} ${mexico} ${ccirast} ${heirast} ${data} ${INLA_model_fit_v2} ${MGN2020_INEGI_Urban_Mex_mask_DISS} ${All_Products_Comp_Over30_Binary_6933}