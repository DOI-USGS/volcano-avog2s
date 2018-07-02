#!/bin/bash

echo "get_gmao.sh:  checking input arguments"
if [ -z $4 ]
then
  echo "Error: Insufficient command-line arguments"
  echo "Usage:  get_gmao.sh YYYY MM DD HH"
  echo "        where HH = 0, 6, 12, or 18"
  exit 1
else
  YYYY=$1
  MM=$2
  DD=$3
  HH=$4
  filename="GEOS.fp.fcst.inst3_3d_asm_Np.${YYYY}${MM}${DD}_00+${YYYY}${MM}${DD}_${HH}00.V01.nc4"
  echo "Downloading file: ${filename}"
fi

WINDROOT="/data/WindFiles"
NASADATAHOME="${WINDROOT}/NASA/GEOS"
#name of directory containing current files
FC_day=${NASADATAHOME}/${yearmonthday}

#go to correct directory
mkdir -p $FC_day
cd $FC_day

wget --password="" ftp://gmao_ops@ftp.nccs.nasa.gov/fp/forecast/Y${YYYY}/M${MM}/D${DD}/H00/GEOS.fp.fcst.inst3_3d_asm_Np.${YYYY}${MM}${DD}_00+${YYYY}${MM}${DD}_${HH}00.V01.nc4
