#!/bin/bash

basedir=$( cd "$(dirname "$0")" ; pwd -P )   

mkdir -p output
source activate r

input_file=${1}
split_string=${2}
outfile_string=${3}

OUTPUTDIR="${PWD}/output"

Rscript ${basedir}/LA_Table47.R ${input_file} ${split_string} ${outfile_string}

##### Rscript /projects/ADE_biomass_harmonization/NASA_CMS/NASA_CMS_2023/Hunka_2024_NSD/Layer_Alignment/LA_Table47.R /projects/my-public-bucket/Data/Hunka_et_al_PAPER_NSD/JRC_TFM/JRC_TMST_10S_040E.tif "TMST_" "_JRC_Transition_Map.tif"