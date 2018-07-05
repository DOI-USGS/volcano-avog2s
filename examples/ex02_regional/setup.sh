#!/bin/bash
ln -s ../../bin/g2s_genSC .
ln -s ../../bin/g2s_Extract_Sonde .
ln -s ../../bin/g2s_Extract_Xsec .
ln -s ../../bin/g2s_Extract_Grid .
ln -s ../../bin/g2s_ResampleAtmos .
ln -s /opt/USGS/AVOG2S/ExternalData/Ap_Forecast/SWPC/Ap_20180411_0005Z.dat Ap.dat
ln -s /opt/USGS/AVOG2S/ExternalData/Ap_Forecast/SWPC/F107_20180411_0005Z.dat F107.dat
ln -s /opt/USGS/AVOG2S/ExternalData/Windfiles/20180410/GEOS.fp.fcst.inst3_3d_asm_Np.20180410_00+20180410_1200.V01.nc4 .
ln -s /opt/USGS/AVOG2S/ExternalData/Windfiles/20180410/nam.t12z.alaskanest.hiresf00.tm00.avo.grib2.nc .
#ln -s /opt/USGS/AVOG2S/ExternalData/Ap_Forecast/SWPC/Ap_20161223_1510Z.dat Ap.dat
#ln -s /opt/USGS/AVOG2S/ExternalData/Ap_Forecast/SWPC/F107_20161223_1510Z.dat F107.dat
#ln -s /opt/USGS/AVOG2S/ExternalData/Windfiles/20161223/GEOS.fp.fcst.inst3_3d_asm_Np.20161223_00+20161223_1200.V01.nc4 .
#ln -s /opt/USGS/AVOG2S/ExternalData/Windfiles/20161223/nam.t00z.alaskanest.hiresf12.tm00.grib2.nc .
