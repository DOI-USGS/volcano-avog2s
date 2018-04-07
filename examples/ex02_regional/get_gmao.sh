#!/bin/bash

YYYY=$1
MM=$2
DD=$3
HH=$4

wget --password="" ftp://gmao_ops@ftp.nccs.nasa.gov/fp/forecast/Y${YYYY}/M${MM}/D${DD}/H00/GEOS.fp.fcst.inst3_3d_asm_Np.${YYYY}${MM}${DD}_00+${YYYY}${MM}${DD}_${HH}00.V01.nc4
