##############################################################################
#  Makefile for AVOG2S
#
#    User-specified flags are in this top block
#
###############################################################################

#      This file is a component of the volcanic infrasound monitoring software
#      written at the U.S. Geological Survey by Hans F. Schwaiger (hschwaiger@usgs.gov)
#      and Alexandra M. Iezzi (amiezzi@alaska.edu).  These programs relies on tools
#      developed for the ash transport and dispersion model Ash3d, written at the
#      U.S. Geological Survey by Hans F. Schwaiger (hschwaiger@usgs.gov), Larry G.
#      Mastin (lgmastin@usgs.gov), and Roger P. Denlinger (roger@usgs.gov).

#      The model and its source code are products of the U.S. Federal Government and therefore
#      bear no copyright.  They may be copied, redistributed and freely incorporated 
#      into derivative products.  However as a matter of scientific courtesy we ask that
#      you credit the authors and cite published documentation of this model (below) when
#      publishing or distributing derivative products.

#      Schwaiger, H.F., Alexandra M. Iezzi and David Fee;
#         AVO-G2S:  A modified Ground-to-Space model for volcano monitoring in Alaska,
#         submitted. 

#      We make no guarantees, expressed or implied, as to the usefulness of the software
#      and its documentation for any purpose.  We assume no responsibility to provide
#      technical support to users of this software.

#  SYSTEM specifies which compiler to use
#    Current available options are:
#      gfortran , ifort
#    This variable cannot be left blank
#      
SYSTEM = gfortran
#
#  RUN specifies which collection of compilation flags that should be run
#    Current available options are:
#      DEBUG : includes debugging info and issues warnings
#      PROF  : includes profiling flags with some optimization
#      OPT   : includes optimizations flags for fastest runtime
#    This variable cannot be left blank
RUN = OPT
#
INSTALLDIR=/opt/USGS
#
# DATA FORMATS
#  For each data format you want to include in the library, set the corresponding
#  variable below to 'T'.  Set to 'F' any you do not want compiled or any unavailable
#  Note: netcdf is required.
USEGRIB2 = T
#USEHDF = F   # Note: the hdf reader is not yet functional

# HWM07 from http://nssdcftp.gsfc.nasa.gov/models/atmospheric/hwm07/HWM07_all_files.zip
# NRLMSISE from https://github.com/graziano-giuliani/Meteostuff/tree/master/NRLMSIS_F90
# SHTOOLS  from https://github.com/SHTOOLS/SHTOOLS/releases

# Indicate the version of HWM model (07 or 14)
# Note: 
#       Version 07 can be downloaded from
#   https://spdf.gsfc.nasa.gov/pub/models/old_models_from_nssdc/atmospheric/hwm07/HWM07_all_files/
#       Version 14 can be downloaded from as supplemental software from
#       Drob et al, 2015 (doi:10.1002/2014EA000089)
#   http://onlinelibrary.wiley.com/store/10.1002/2014EA000089/asset/supinfo/ess224-sup-0002-supinfo.tgz?v=1&s=2a957ba70b7cf9dd0612d9430076297c3634ea75
HWMv=14

# Location and version of SHTOOLs installation
#SHTv=2.6
#SHTOOLROOT=/usr/local/SHTOOLS$(SHTv)
#SHTlib=SHTOOLS$(SHTv)
#SHTOOLINC=${SHTOOLROOT}/modules
#  Current version is 4.0 and installs by default to /usr/local
#SHTv=4.0
SHTlib=SHTOOLS
SHTOOLROOT=/usr/local
SHTOOLINC=${SHTOOLROOT}/include

# This is the location of the USGS libraries and include files
# This might be different from the INSTALLDIR path if you want to link to special MetReader
# projection and HoursSince packages.
USGSROOT=/opt/USGS

###############################################################################
#####  END OF USER SPECIFIED FLAGS  ###########################################
###############################################################################

# This is where HWM and NRLMSIS are installed
AVOG2SINSTALLDIR=$(INSTALLDIR)/AVOG2S
HWMPATH=$(AVOG2SINSTALLDIR)/ExternalData/HWM$(HWMv)
HWTLIB=-L$(AVOG2SINSTALLDIR)/lib -lnrlmsis -lHWM$(HWMv) -I$(AVOG2SINSTALLDIR)/include

# Netcdf is required, but MetReader may have been compiled with grib2 support
ifeq ($(USEGRIB2), T)
 # These are the libraries for grib_api
 DATALIBS =  -lnetcdf -lnetcdff -lgrib_api_f90 -lgrib_api
 # These are the libraries for ecCodes
 #DATALIBS =  -lnetcdf -lnetcdff -leccodes -leccodes_f90
 grb2FPPFLAG = -DUSEGRIB2
else
 DATALIBS =  -lnetcdf -lnetcdff
 grb2FPPFLAG=
endif

# location of MetReader, HoursSince and projection
USGSLIBDIR = -L$(USGSROOT)/lib
USGSLIB = -lMetReader -lhourssince -lprojection
USGSINC = -I$(USGSROOT)/include

EXTRALIBS = -L$(SHTOOLROOT)/lib -I$(SHTOOLINC) $(DATALIBS) -l$(SHTlib) -lfftw3 -lfftw3f -lm -llapack

###############################################################################
###########  GNU Fortran Compiler  #############################################
#ifeq ($(SYSTEM), gfortran)
    FCHOME=/usr
    FC = /usr/bin/gfortran

    COMPINC = -I$(FCHOME)/include -I$(FCHOME)/lib64/gfortran/modules
    COMPLIBS = -L$(FCHOME)/lib -L$(FCHOME)/lib64

    #COMPLIBS = -I/usr/lib64/gfortran/modules -I/usr/include -L/usr/lib/x86_64-linux-gnu

    LIBS = $(COMPLIBS) $(USGSLIBDIR) $(USGSINC) $(COMPINC) $(USGSLIB) $(EXTRALIBS)

## Debugging flags
ifeq ($(RUN), DEBUG)
    FFLAGS = -O0 -g3 -Wall -fbounds-check -pedantic -fbacktrace -fimplicit-none -Wunderflow -Wuninitialized -ffpe-trap=invalid,zero,overflow -fdefault-real-8
endif
## Profiling flags
ifeq ($(RUN), PROF)
    FFLAGS = -g -pg -w -fno-math-errno -funsafe-math-optimizations -fno-trapping-math -fno-signaling-nans -fcx-limited-range -fno-rounding-math -fdefault-real-8
endif
## Production run flags
ifeq ($(RUN), OPT)
    FFLAGS = -O3 -w -fno-math-errno -funsafe-math-optimizations -fno-trapping-math -fno-signaling-nans -fcx-limited-range -fno-rounding-math -fdefault-real-8
endif
FPPFLAGS = -x f95-cpp-input -DuseHWM$(HWMv) $(grb2FPPFLAG)

###############################################################################

#$(FC) $(FPPFLAGS) $(FFLAGS) $(HWTLIB) $(COMPLIBS) $(EXTRALIBS) $(USGSINC) $(USGSLIBDIR) $(USGSLIB) -c g2s_genSC.f90
#$(FC) $(FFLAGS) g2s_genSC.o g2s_globvar.o -o g2s_genSC $(HWTLIB) $(COMPLIBS) $(EXTRALIBS) $(USGSINC) $(USGSLIBDIR) $(USGSLIB)


all: g2s_genSC g2s_ResampleAtmos g2s_Extract_Sonde g2s_Extract_Xsec g2s_Extract_Grid probe_HWT

g2s_genSC: g2s_genSC.f90 g2s_globvar.o makefile
	$(FC) $(FPPFLAGS) $(FFLAGS) $(HWTLIB) $(LIBS) -c g2s_genSC.f90
	$(FC) $(FFLAGS) g2s_genSC.o g2s_globvar.o -o g2s_genSC $(HWTLIB) $(LIBS)
	mkdir -p ../bin
	mv g2s_genSC ../bin/
	cp ../bin/g2s_genSC ../bin/g2s_genSC_HWM${HWMv}
probe_HWT: probe_HWT.f90 g2s_globvar.o makefile
	$(FC) $(FPPFLAGS) -DSONDE $(FFLAGS) $(HWTLIB) $(LIBS) -c probe_HWT.f90
	$(FC) $(FFLAGS) probe_HWT.o g2s_globvar.o $(HWTLIB) $(LIBS) -o probe_HWT${HWMv}
	mkdir -p ../bin
	mv probe_HWT${HWMv} ../bin/
g2s_ResampleAtmos: g2s_ResampleAtmos.f90 g2s_globvar.o makefile
	$(FC) $(FFLAGS) $(LIBS) -c g2s_ResampleAtmos.f90
	$(FC) $(FFLAGS) g2s_ResampleAtmos.o g2s_globvar.o $(LIBS) -o g2s_ResampleAtmos
	mkdir -p ../bin
	mv g2s_ResampleAtmos ../bin/
g2s_Extract_Sonde: g2s_Extract.f90 g2s_globvar.o makefile
	$(FC) -x f95-cpp-input -DSONDE $(FFLAGS) $(LIBS) -c g2s_Extract.f90
	$(FC) $(FFLAGS) g2s_Extract.o g2s_globvar.o $(LIBS) -o g2s_Extract_Sonde
	mkdir -p ../bin
	cp g2s_Extract_Sonde $(HOME)/gcc/bin
	mv g2s_Extract_Sonde ../bin/
g2s_Extract_Xsec: g2s_Extract.f90 g2s_globvar.o makefile
	$(FC) -x f95-cpp-input -DXSEC $(FFLAGS) $(LIBS) -c g2s_Extract.f90
	$(FC) $(FFLAGS) g2s_Extract.o g2s_globvar.o $(LIBS) -o g2s_Extract_Xsec
	mkdir -p ../bin
	cp g2s_Extract_Xsec $(HOME)/gcc/bin
	mv g2s_Extract_Xsec ../bin/
g2s_Extract_Grid: g2s_Extract.f90 g2s_globvar.o makefile
	$(FC) -x f95-cpp-input -DGRID $(FFLAGS) $(LIBS) -c g2s_Extract.f90
	$(FC) $(FFLAGS) g2s_Extract.o g2s_globvar.o $(LIBS) -o g2s_Extract_Grid
	mkdir -p ../bin
	cp g2s_Extract_Grid $(HOME)/gcc/bin
	mv g2s_Extract_Grid ../bin/
g2s_globvar.o: g2s_globvar.f90 makefile
	$(FC) $(FFLAGS) g2s_globvar.f90 -c
#checkhwm: checkhwm07.f90
#	$(FC) $(FFLAGS) checkhwm07.f90 $(COMPLIBS) -lHWM${HWMv} -o checkhwm
#testgtd: testgtd7.f90
#	$(FC) $(FFLAGS) testgtd7.f90 $(COMPLIBS) -lnrlmsis -o testgtd
clean:
	rm -f *.mod *.o g2s_genSC g2s_ResampleAtmos g2s_Extract_Sonde g2s_Extract_Xsec

install:
	install -d $(AVOG2SINSTALLDIR)/bin
	install -d $(AVOG2SINSTALLDIR)/ExternalData/Ap_Forecast
	install -d $(AVOG2SINSTALLDIR)/ExternalData/Ap_Forecast/NGDC_NOAA_Archive
	install -d $(AVOG2SINSTALLDIR)/ExternalData/Ap_Forecast/SWPC
	install -m 755 ../bin/g2s_* $(AVOG2SINSTALLDIR)/bin/
	install -m 755 ../bin/probe_HWT* $(AVOG2SINSTALLDIR)/bin/
	install -m 755 ExternalDataSoftware/get_ApFC $(AVOG2SINSTALLDIR)/ExternalData/Ap_Forecast/
	install -m 755 ExternalDataSoftware/get_NGDC $(AVOG2SINSTALLDIR)/ExternalData/Ap_Forecast/
