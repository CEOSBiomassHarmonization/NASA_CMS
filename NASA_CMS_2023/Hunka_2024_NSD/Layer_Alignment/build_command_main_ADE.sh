#!/bin/bash
set -x
basedir=$( cd "$(dirname "$0")" ; pwd -P )

#install requirements packages
conda env create -f ${basedir}/env_main_ADE.yml

pushd ${HOME}
source activate r

# Do not remove this (PMM Dec 2022)
# source activate icesat2_boreal
# conda install -c conda-forge r-optparse r-ranger -y