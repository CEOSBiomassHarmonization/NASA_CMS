#!/bin/bash
set -x

basedir=$( cd "$(dirname "$0")" ; pwd -P )

conda env update -n python --solver=libmamba -f ${basedir}/env_main_ADE.yaml

# conda env update -n python --solver=libmamba -f ${basedir}/env_main_ADE.yaml
# pip3 install pyOpenSSL --upgrade