#!/bin/bash
set -x
basedir=$( cd "$(dirname "$0")" ; pwd -P )
# Install dependencies
conda env update --name r -f ${basedir}/env_main_ADE.yaml

# Install INLA
conda run --name r --no-capture-output Rscript ${basedir}/install.R
