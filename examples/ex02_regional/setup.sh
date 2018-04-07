#!/bin/bash
ln -s ../../bin/g2s_genSC .
ln -s ../../bin/g2s_Extract_Sonde .
ln -s ../../bin/g2s_Extract_Xsec .
ln -s ../../bin/g2s_Extract_Grid .
ln -s ../../bin/g2s_ResampleAtmos .
ln -s /opt/USGS/AVOG2S/ExternalData/Ap_Forecast/NGDC_NOAA_Archive NGDC
ln -s /opt/USGS/AVOG2S/ExternalData/Windfiles/20161223/GEOS.fp.fcst.inst3_3d_asm_Np.20161223_00+20161223_1200.V01.nc4 .
ln -s /opt/USGS/AVOG2S/ExternalData/Windfiles/20161223/nam.t00z.alaskanest.hiresf12.tm00.grib2.nc .
ln -s /data/TOPO/ETOPO1/ETOPO1_Ice_c_gmt4.nc etopo.nc
