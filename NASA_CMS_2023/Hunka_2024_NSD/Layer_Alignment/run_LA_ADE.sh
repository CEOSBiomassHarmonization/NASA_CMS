#!/bin/bash

basedir=$( cd "$(dirname "$0")" ; pwd -P )   

mkdir -p output

input_file=${1}
split_string=${2}
outfile_string=${3}

OUTPUTDIR="${PWD}/output"

Rscript ${basedir}/LA_Table4.7.R ${input_file} ${split_string} ${outfile_string}