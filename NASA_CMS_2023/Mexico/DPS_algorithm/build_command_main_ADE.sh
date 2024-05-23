#!/bin/bash
set -x
basedir=$( cd "$(dirname "$0")" ; pwd -P )
# Install dependencies
conda env update --name r -f ${basedir}/env_main_ADE.yaml

# Install INLA
<<<<<<< HEAD
conda run --name r --no-capture-output Rscript ${basedir}/install.R
=======
conda run --name r --no-capture-output Rscript ${basedir}/install.R
>>>>>>> fa24f7168ccaaa40f55441096a436b70c449fb21
