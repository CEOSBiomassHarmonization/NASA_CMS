#!/bin/bash
set -x
basedir=$( cd "$(dirname "$0")" ; pwd -P )
conda env update -n base --solver=libmamba -f ${basedir}/env_main_ADE.yaml

# conda install conda-forge::r-fmesher --yes
# conda install conda-forge::r-exactextractr --yes
# conda install conda-forge::r-sn --yes
# conda install conda-forge::r-MatrixModels --yes
# conda install conda-forge::r-inlabru --yes
# conda install conda-forge::r-terra --yes
# conda install conda-forge::r-dplyr --yes
# conda install conda-forge::r-spdep --yes
# conda install conda-forge::r-sf --yes
# conda install conda-forge::r-fields --yes
# conda install conda-forge::r-Matrix --yes