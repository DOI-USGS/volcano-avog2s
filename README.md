AVOG2S
==========

AVOG2S is a collection of programs that can be used to build a Ground-2-Space
model of the atmosphere in a similar style to that outlined by Drob et al.,
(JGR, v108, p4680, 2003, doi:10.1029/2002JD003307) by smoothly merging 1 or 2 
numerical weather prediction (NWP) models (either forecast or reanalysis) with
empirical models of the upper atmosphere.  This Ground-2-Space model can be used
for forward modeling of infrasound propagation.  This model is described in Schwaiger
et al., (submitted) and can be used as both global as well as a regional G2S model.
This model is used by the Alaska Volcano Observatory for infrasound propagation
modeling.

This model relies heavily on externally provided software.  The following 
software packages must be built and installed before the AVOG2S software
suite can be built.
(1) To download and read the NWP data, several USGS libraries must be installed:
   HoursSince https://github.com/hschwaiger-usgs/HoursSince
   projection https://github.com/hschwaiger-usgs/projection
   MetReader  https://github.com/hschwaiger-usgs/MetReader
(2) Upper-atmospheric empirical models are available separately for horizontal winds
    and temperature
   HWM14      http://onlinelibrary.wiley.com/store/10.1002/2014EA000089/asset/supinfo/ess224-sup-0002-supinfo.tgz?v=1&s=2a957ba70b7cf9dd0612d9430076297c3634ea75
   NRLMSIS-00 https://github.com/graziano-giuliani/Meteostuff/tree/master/NRLMSIS_F90
(3) Spectral decomposition of the data utilize two libraries, one for spherical
    harmonic decompositions, and one for Fourier decomposition
   SHTOOLS    https://github.com/SHTOOLS/SHTOOLS/releases
   fftw
(4) Additional libraries needed:
   lapack
   netcdf
   grib2

Data
(1) Ap and F107 values
    Current values available from NOAA Space Weather Prediction Center
      http://services.swpc.noaa.gov/text/wwv.txt
    Archived values available from NOAA National Geophysical Data Center
      ftp://ftp.ngdc.noaa.gov/STP/GEOMAGNETIC_DATA/INDICES/KP_AP/
(2) Numerical Weather Prediction data
    Forecast and reanalysis data available from many sources.  See documentation
    in MetReader.


Preliminary Software Installation
---------------------------------
(1) USGS libraries
  Three libraries from the Ash3d package are needed: projection, HourSince, and
MetReader.  These libraries are currently available at the locations given above.
Installation instructions are given in the repositories for each of these libraries.

(2) Empirical models
  (a) HWM14 is available as supplemental material for Drob et al, 2015
(http://onlinelibrary.wiley.com/doi/10.1002/2014EA000089/abstract).  The build
of AVOG2S expects that this be compiled
as a library with the data files installed in a specified location.  The makefile
src/ExternalDataSoftware/makefile_HWM14.gfortran is provided which can be used to
build the HWM14 library.  The install location can be changed by editing the top
of this file.

To build the library, type:
  cd src/ExternalDataSoftware/
  wget http://onlinelibrary.wiley.com/store/10.1002/2014EA000089/asset/supinfo/ess224-sup-0002-supinfo.tgz?v=1&s=2a957ba70b7cf9dd0612d9430076297c3634ea75
Or open the page http://onlinelibrary.wiley.com/doi/10.1002/2014EA000089/abstract
and download ess224-sup-0002-supinfo.tgz from Supporting Information.

  gunzip ess224-sup-0002-supinfo.tgz
  tar -xvf ess224-sup-0002-supinfo.tar
  cd HWM14
The package makefile builds an example program.  To build as a library, use the makefile
provided by AVOG2S
  make -f ../makefile_HWM14.gfortran libra
  make -f ../makefile_HWM14.gfortran install

To use software that links to this library, you will need to add the following to your
.bash_profile:
HWMPATH=/opt/USGS/AVOG2S/ExternalData/HWM14  # (or whereever you installed the library)
export HWMPATH

If you have difficulty installing HWM14 (e.g. if you are using an older compiler), you
can optionally install HWM07.  The fortran source files and data files are available
from https://spdf.gsfc.nasa.gov/pub/models/old_models_from_nssdc/atmospheric/hwm07/HWM07_all_files/

To build this library, type:
cd src/ExternalDataSoftware/HWM07
wget https://spdf.gsfc.nasa.gov/pub/models/old_models_from_nssdc/atmospheric/hwm07/HWM07_all_files/apexcord.f90
wget https://spdf.gsfc.nasa.gov/pub/models/old_models_from_nssdc/atmospheric/hwm07/HWM07_all_files/apexgrid.dat
wget https://spdf.gsfc.nasa.gov/pub/models/old_models_from_nssdc/atmospheric/hwm07/HWM07_all_files/checkhwm07.f90
wget https://spdf.gsfc.nasa.gov/pub/models/old_models_from_nssdc/atmospheric/hwm07/HWM07_all_files/dwm07b.f90
wget https://spdf.gsfc.nasa.gov/pub/models/old_models_from_nssdc/atmospheric/hwm07/HWM07_all_files/dwm07b_104i.dat
wget https://spdf.gsfc.nasa.gov/pub/models/old_models_from_nssdc/atmospheric/hwm07/HWM07_all_files/hwm07e.f90
wget https://spdf.gsfc.nasa.gov/pub/models/old_models_from_nssdc/atmospheric/hwm07/HWM07_all_files/hwm071308e.dat
wget https://spdf.gsfc.nasa.gov/pub/models/old_models_from_nssdc/atmospheric/hwm07/HWM07_all_files/makefile.ifort
wget https://spdf.gsfc.nasa.gov/pub/models/old_models_from_nssdc/atmospheric/hwm07/HWM07_all_files/readme.txt

  You will have to edit apexcord.f90 (line 76), dwm07b.f90 (line 53), and hwm07e.f90 (line 126)
and change the data file names to include the full path, such as
/opt/USGS/AVOG2S/ExternalData/HWM07/apexgrid.dat

make -f ../makefile_HWM07.gfortran libra
make -f ../makefile_HWM07.gfortran install

This version of HWM does not require an environment variable to be set, but expects
that the data files are at the hardwired location.

  (b) The NRLMSIS-00 code has been ported to fortran 90 and is available from
https://github.com/graziano-giuliani/Meteostuff
After downloading the software bundle into src/ExternalDataSoftware, the code can
be compiled into a library using the makefile provided at:
src/ExternalDataSoftware/makefile_NRLMSIS.gfortran

Note: NRLMSIS from the link above can return erroneous values if a requested altitude
      equals an internal bracketing value.  To correct this, change line 2067 of 
      physics_msis.f90 from 
      if ( dabs(v1-d_1) < nearzero .or. alt > zn2(1) .or. &
        to
      if ( dabs(v1-d_1) < nearzero .or. alt >= zn2(1) .or. &

cd src/ExternalDataSoftware/Meteostuff-master/NRLMSIS_F90
  Add bug fix to source file
vi physics_msis.f90
make -f ../../makefile_NRLMSIS.gfortran all
make -f ../../makefile_NRLMSIS.gfortran install

(3)
 (a) SHTOOLS
  SHTOOLS is available from https://github.com/SHTOOLS/SHTOOLS/releases.  Instructions
for installation are in the source code download
 wget https://github.com/SHTOOLS/SHTOOLS/archive/v4.0.tar.gz
This can be built and installed in the default location.  Take note of where the objects
are installed.  Earlier versions of SHTOOLS installing in /usr/local/SHTOOLS$(version_number).
The latest version simple installs to /usr/local/.

 (b) fftw is available as a distribution package

(4) netcdf4, lapack and grib-api are available as distribution libraries.


AVOG2S Installation
---------------------------------

To compile, edit the makefile to be consistent with the install directory and options
used in the installation of the preliminary software.  Double-check the HWM version
number, the SHTOOL version number if you are using an older version, and whether or not
grib2 support should be built.  Then simply type:

  make all

To install the library, edit the INSTALLDIR variable of the makefile (the
default is /opt/USGS) and type:

  make install

You will need to have write permission in ${INSTALLDIR} or install as root.
This will install the following in ${INSTALLDIR}/bin/:
 g2s_genSC_HWM14   : program to calculate spectral coefficients for G2S model
 g2s_ResampleAtmos : program to reconstruct a gridded atmosphere from the coefficients
 probe_HWT14       : stand-along program to calculate empirical HWM and NRLMSIS values
 g2s_Extract_Sonde : program to resample the gridded atmosphere onto a 1-d profile
 g2s_Extract_Xsec  : program to resample the gridded atmosphere onto a 2-d cross-section
 g2s_Extract_Grid  : program to resample the gridded atmosphere onto a 3-d grid
and the following scripts in ${INSTALLDIR}/ExternalData/Ap_Forecast/
 get_ApFC : script that downloads the current space weather indicies from NOAA
 get_NGDC : script that downloads an archive year of the space weather indicies

Output from the extraction tools is designed to be used with the GeoAc forward
modeling software (https://github.com/LANL-Seismoacoustics/GeoAc).

Please see the user's guide for more information on using this software.

Authors
-------

Hans F. Schwaiger <hschwaiger@usgs.gov>
Alexandra M. Iezzi <amiezzi@alaska.edu>
