AVOG2S User Guide
=================

Please see README.md for instructions on installing AVO2G2S.

Once the AVOG2S software is built and installed, the following executables are
copied to $(INSTALLDIR)/bin/ :
 g2s_genSC_HWM14   : program to calculate spectral coefficients for G2S model
 g2s_ResampleAtmos : program to reconstruct a gridded atmosphere from the coefficients
 probe_HWT14       : stand-along program to calculate empirical HWM and NRLMSIS values
 g2s_Extract_Sonde : program to resample the gridded atmosphere onto a 1-d profile
 g2s_Extract_Xsec  : program to resample the gridded atmosphere onto a 2-d cross-section
 g2s_Extract_Grid  : program to resample the gridded atmosphere onto a 3-d grid

The program probe_HWT14 can be used to generate a vertical profile at a point using
just the empirical model.  Running the program without and command-line arguments will
print usage information to stdout before exiting.

The work-flow is as follows.
 (1) Run g2s_genSC_HWM14 to create a spectral coefficient file
 (2) Rebuilt a grid of atmospheric values from the coefficient file
 (3) Extract the needed 1, 2, or 3-d values from the rebuilt grid
 (4) Use the extracted output in an infrasound propagation model such as Art2d or GeoAc.

Step 1: Create the spectral coefficient file
--------------------------------------------

After the default installation, the executable g2s_genSC_HWM14 (also copied to g2s_genSC)
will be copied to /opt/USGS/AVOG2S/bin/.  If you used HWM07 instead of or in addition to
HWM14, then g2s_genSC_HWM07 would be created with the last version built copied to g2s_genSC.

This executable requires a command-line argument.  If none is given, then an informational
usage message is printed to stdout on exit.  The expected usage is to provide the name of
a control file as the command-line argument.  This control file gives all the information
needed to generate the coefficient file.  This file has the following format.

Line 1 : The date and time of the coefficient file
 2016 10 5 4.6        ! year, month, day, hour

Line 2 : The projection information of the needed G2S grid (can be distinct from
the NWP grid).  The details of the projection values are given in the documentation
of the projection library.  The first value is an integer that indicates if Lon/Lat
coordinates are used (1 for 'yes', G2S grid is Lon/Lat, or 0 for a projected grid).
If 0 is given for the first value, then a second integer on the line is expected
identifying the projection to be used: 1 for polar stereographic, 2 for Albers Equal
Area, 3 for UTM, 4 for Lambert conformal conic and 5 for Mercator.  Depending on which
projection is used, additional values are needed to define the full projection.
Examples:
1 4 -107.0 50.0 50.0 50.0 6367.470  Lon/Lat since first value is 1, all other values ignored
0 1 -150.0 90.0 0.933 6371.229      Polar Ster. grid (NAM grid 198 over Alaska)

Line 3: Resolution of G2S grid.  If the G2S grid is Lon/Lat, then the spectral decomposition
is assumed to be global using SHTOOLS.  Two values are expected on this line giving the
degree resolution and the order of the decomposition.  If the G2S grid is projected, then
the resolution term is given in km and is used for both dx and dy.  This is followed by
the order and the positioning information (starting x coordinate, nx, starting y coordinate,
ny).  
Examples:
   For the global case
 0.5 120                                ! ddeg of G2S grid, maxdeg of SH decomp
   For the projected grid over the Aleutians and Bering Sea
 5.950 120 -2172.922 512 -4214.803 340  ! dx,maxdeg,xstart,nx,ystart,ny 

Line 4: Vertical grid information for the empirical model.  Three real values 
are needed: the starting altitude of the HWM/NRLMSIS sampling (in km), the 
maximum altitude and the resolution.  This should overlap with the NWP grid specified
later in the control file.
 20 200 2.5           ! zmin_HWT,zmax_HWT,dz_HWT

Line 5: The number of groups of NWP files.  This refers to the number of meteorological
data sets that will be blended in generating the coefficient file.  Each group of NWP
files may contain one or more files spanning the time given in line 1.  This is useful
if you want to use high-resolution weather forecasting for tropospheric conditions with
a product that gives mesospheric values such as ECMWF or NASA GEOS.
 2                    ! number of windfile groups

Line 6: Format ID and number of files of the first NWP package.  The format ID is the
code expected by the MetReader library for the particular NWP package.  See the MetReader
users manual for the full documentation.  The following products are currently enabled:
    Format_ID = 3  : North American Regional Reanalysis (32 km)
    Format_ID = 12 : 5.9 km forecast over AK (NAM grid 198)
    Format_ID = 13 : 2.95 km forecast over AK (NAM grid 91)
    Format_ID = 20 : 0.5  degree GFS
    Format_ID = 22 : 0.25 degree GFS
    Format_ID = 24 : 1.25 degree NASA Merra Reanalysis
    Format_ID = 25 : 2.5 degree NCEP Reanalysis
    Format_ID = 28 : 0.7 degree ECMWF Interim Reanalysis (ERA)
    Format_ID = 40 : 0.6250 x 0.50 degree NASA GEOS Cp
    Format_ID = 41 : 0.3125 x 0.25 degree NASA GEOS Np
The number of windfiles in this NWP group can be 1, in which case the requested time
must be close to the time associated with the NWP data, or there can be multiple files.
 20 1                 ! format ID of met data1 (20 for new GFS) ,num of windfiles1

Line 7: Vertical grid information needed from the first NWP data.  Similar to line 4
for the empirical model, this line specifies the minimum altitude, maximum altitude, 
and resolution (all in km) to be populated with meteorological data. 
 0.0 25.0 1.0         ! zmin, zmax, dz of Met1

Line 8 - ... : This line starts the list of NWP file in the first group.  The number of
lines read is given in line 6.  Note, if the NCEP reanalysis product is used (format
ID = 25), then only one 'filename' is given, pointing to the root directory of the 
NCEP data.  See example 1 for more information.
 gfs.t00z.pgrb2f00.nc ! windfile name

Subsequent line numbers are not fixed.  Following the first group of NWP files, if
two groups were specified, a second block is read giving information of the 
mesospheric NWP data.
 41 1                 ! format ID of met data2  ,num of windfiles2
 20.0 55.0 2.0        ! zmin, zmax, dz of Met2
 GEOS.fp.fcst.inst3_3d_asm_Np.20161005_00+20161005_0000.V01.nc4 ! windfile name

After the NWP specifications, the two following lines give the location of the files
specifying the Ap and F107 values.  These files are typically soft-linked to files
that are automatically downloaded from NOAA (further details on this are given below),
but they can also be regular files in the current working directory containing the
value of the corresponding index.  If both the Ap and F107 filenames are the same,
the g2s_genSC program interprets this as the name of the directory containing the
NOAA archive of these values.  The archive has a latency of a few months, so this
is only useful for retrospective studies.
 Ap.dat               ! Filename with Ap value
 F107.dat             ! Filename with F107 value

Finally, the output file name for the coefficient file is given on the last line.
 out_GS.nc            ! Output filename of SH coefficients

Here is an example of a control file using a projected grid for a regional Ground-
2-Space model over the Aleutians.  Two forecast products are used, the NAM high-
resolution forecasts over Alaska (0-25km) and the NASA GEOS (20-58km).  Planetary
space weather indicies are from the NOAA Space Weather Prediction Center.

2016 12 23 12.0
0 1 -150.0 90.0 0.933 6371.229
5.950 120 -2172.922 512 -4214.803 340
50 200 2.5
2
12 1 198 3
0.0 25.0 1.0
nam.t00z.alaskanest.hiresf12.tm00.grib2
41 1
20.0 58.0 1.0
GEOS.fp.fcst.inst3_3d_asm_Np.20161223_00+20161223_1200.V01.nc4
Ap.dat
F107.dat
G2S_SC_20161223_12Z_wf12wf41.nc


A global model of the same time using NCEP 2.5 degree NWP data and the archive of
planetary indicies from the NOAA National Geophysical Data Center.

2016 12 23 12.0
1 1 -150.0 90.0 0.933 6371.229
0.5 120
20 200 2.5
1
25 1
0.0 25.0 1.0
NCEP
NGDC
NGDC
G2S_SC_20161223_12Z_wf25.nc


Step 2: Reconstruct a gridded atmosphere from the coefficient file
------------------------------------------------------------------

To reconstruct a gridded atmosphere from the spectral decomposition, run
g2s_ResampleAtmos with the name of the coefficient file as the single command-line
argument.  This will use the default vertical grid from 0-200 km, with a node
spacing of 1 km.

./g2s_ResampleAtmos G2S_SC_20161223_12Z_wf25.nc

Alternatively, you can provide six additional values on the command-line to 
specify three zones to define the grid:

./g2s_ResampleAtmos Coeffic_File.nc nz1 dz1 nz2 dz2 nz3 dz3

For example
  ./g2s_ResampleAtmos Coeffic_File.nc 15 1.0 18 2.0 30 5.0
will resample the atmosphere at a 1.0 km spacing from 0-14 km, then at a 2.0 km
spacing from 16-50 km (18 grid points), and at a 5.0 km spacing from 55-200 km.

This will create three gridded binary files, named with the root of the coefficient
file, but with '_[T,U,V]_res.raw' appended:
G2S_SC_20161223_12Z_wf25_T_res.raw
G2S_SC_20161223_12Z_wf25_U_res.raw
G2S_SC_20161223_12Z_wf25_V_res.raw

Note that these gridded binary files to not contain any information about the
geometry of the grids.  They just are a record of the sequence of values in
the gridded volume.


Step 3: Extract profiles, cross-sections for volumes from the resamples data
----------------------------------------------------------------------------

The three programs for extracting values (in 1,2 and 3-d) from the gridded binary
files are all variants of the source file g2s_Extract.f90, compiled with different
preprocessor directives.  For each of these tools, the program needs to know
the geometry of the raw data, in addition to the information on where to
probe the files.  Each of these programs can be run in one of two ways: 
(1) with the geometry of the raw data assumed and with the probe information
given on the command line.
(2) with the full details given in a control file.

If the extraction tools are run assuming the geometry of the gridded binary files,
then both g2s_ResampleAtmos.f90 and g2s_Extract.f90 must have the same preset values.
These are set at the top of each file, just after the variable declaration.  

For each of these tools, the usage and example control files are written to stdout
if they are run with no command-line arguments.  This usage message will also
include the values assumed if no gridded binary information is given.
The current default values are for a lon/lat grid at 0.5 degree resolution.

1-d: g2s_Extract_Sonde
----------------------
For 1-d values, 5 command-line arguments are required:
x and y coordinate of vertical profile (or lon, lat), the maximum height of the profile
in km, the vertical increment in km, and the root name of the output file

  g2s_Extract_Sonde 190.055 52.8222 180.0 0.2 Clev

Will return a verticle profile above Cleveland Volcano from 0-180 km with a z-spacing
of 0.2 km.  Data will be written to the file Clev0.met with columns for
z [km] : T(z) [K] : u(z) [m/s] : v(z) [m/s] : rho(z) [g/cm3] : p(z) [mbar]

Using a control file provides greater flexibility and is needed to override the
default values.  It has the following format.

190.055 52.8222                          # x (or lon), y (or lat) coordinate of probe
180.0 0.2                                # zmax, dz of output grid
1                                        # LonLat flag (1 for LonLat, 0 for projected)
720 0.5 0.0                              # nx, dx, xstart of gridded binary data
361 0.5 -90.0                            # ny, dy, ystart of gridded binary data
15 1.0                                   # nz1, dz1 \
18 2.0                                   # nz2, dz2 |-> vertical grid description
30 5.0                                   # nz3, dz3 /
G2S_SH_20161223_12Z_wf20wf41_U_res.raw   # filename of gridded binary U values
G2S_SH_20161223_12Z_wf20wf41_V_res.raw   # filename of gridded binary V values
G2S_SH_20161223_12Z_wf20wf41_T_res.raw   # filename of gridded binary T values
Clev                                     # root filename of output
etopo.nc                                 # filename of topography file

If topography is to be used, only the ETOPO1 topography file is currently supported:
https://www.ngdc.noaa.gov/mgg/global/relief/ETOPO1/data/ice_surface/cell_registered/netcdf/ETOPO1_Ice_c_gmt4.grd.gz

2-d: g2s_Extract_Xsec
----------------------
For 2-d values, 7 command-line arguments are required.
x and y coordinate of vertical profile (or lon, lat), the maximum height of the profile
in km, the vertical increment in km, azimuth of the profile, length of the profile in
degrees, and the increment along length in degrees.

When values are passed on the command line for g2s_Extract_Xsec, only parameters for
global gridded binary data is implemented.  For projected grids, a control file is
necessary.

x_in,y_in              #lon,lat (or x,y) of reference point
zmax, dz_prof          #maximum altitude of output, vertical increment
IsLatLon               # 1 for global (lat/lon) grids, 0 for projected
az, len , ds_xsec      #azimuth, length and increment of xsec
             - or - 
x2, y2, len_xsec, dx_xsec # coordinate of in-line point, length and increment of xsec (if projected)
nx, dx, xmin   #number, spacing, and start of gridpoints in x
ny, dy, ymin   #number, spacing, and start of gridpoints in y
nz1, dz1               #number and spacing of low alt points
nz2, dz2               #number and spacing of mid alt points
nz3, dz3               #number and spacing of high alt points
vx_fn                  #file name of resampled vx values
vy_fn                  #file name of resampled vy values
tmp_fn                 #file name of resampled temperature values
out_fn                 #file name of output xsec data for Art2D
out_topo_fn            #file name of output topo data for Art2D
topo_fn                #topography filename

This will create five files designed to be used as input for the infrasound propagation
model, Art2D.

3-d: g2s_Extract_Grid
---------------------
For 3-d values, 9 command-line arguments are required.
   LatStart  : start latitude of grid
   LatEnd    : end latitude of grid
   LatCnt    : number of gridpoints along lat
   LongStart : start longitude of grid
   LongEnd   : end longitude of grid
   LongCnt   : number of gridpoints along lon
   zmax      : maximum altitude of profile
   dz_prof   : vertical increment
   out_root  : root of output filename

This program essentially calls the 1-d profile subroutine over a 2-d grid of lon,lat
coordinates and creates a series of files designed for the the infrasound propagation
model GeoAc.

For example, to extract a grid of values to encompass the region between Cleveland
Volcanao and the Dillingham infrasound array:
     ./g2s_Extract_Grid  50.0 60.0 3 185.0 205.0 5 180.0 0.2 Clev
This produces 15 vertical profile files: Clev[0-14].met;
and two files giving the coordinates of the grid points: Clev.loclat, Clev.loclon

Alternatively, a control file can be used with the following format:

   LatStart,LatEnd,LatCnt    : start, end and count of lat (or y) grid
   LongStart,LongEnd,LongCnt : start, end and count of lon (of x) grid
   zmax, dz_prof             : maximum altitude of output, vertical increment
   IsLatLon          : 1 for global (lat/lon) grids, 0 for projected
   nx, dx, [xmin]    : number and spacing of gridpoints in x, start value for proj grids
   ny, dy, [ymin]    : number and spacing of gridpoints in y, start value for proj grids
   nz1, dz1          : number and spacing of low alt points
   nz2, dz2          : number and spacing of mid alt points
   nz3, dz3          : number and spacing of high alt points
   vx_fn             : file name of resampled vx values
   vy_fn             : file name of resampled vy values
   tmp_fn            : file name of resampled temperature values
   out_fn            : file name of output sonde data for GeoAc
   topo_fn           : (OPTIONAL) topography filename


Step 4: Use output of g2s_Extract_* in infrasound codes
-------------------------------------------------------

The output files can be used by GeoAc by running
GeoAc3D.RngDep -prop Clev Clev.loclon Clev.loclat theta_step=2.0 bounces=2 azimuth=41.0



Examples
--------
(1) Creating a reanalysis G2S model using NCEP 2.5-degree NWP data and archived planetary indicies

The files for this example are located in examples/ex01

Step 1:  First, we need to gather the data that we will need.  For NWP date, this example
will use the NCEP Reanalysis 2.5-degree.  The following script should have been
installed with MetReader and will download the files to the default location at
    /data/Windfiles/NCEP/YYYY
If you ran the example while installing MetReader, you should already have these files.
Note: This will download 6 files, a year of data for each of temperature (air.2016.nc),
      geopotential height (hgt.2016.nc), zonal winds (uwnd.2016.nc), meridonial winds (vwnd.2016.nc),
      vertical velocities (omega.2016.nc) and specific humidity (shum.2016.nc), totalling 2.4 Gb.

  /opt/USGS/bin/autorun_scripts/get_NCEP_50YearReanalysis.sh 2016

Link the data to the current working directory
  ln -s /data/WindFiles/NCEP .

Next, we need the planetary index data.  This example uses the archive from the NOAA National
Geophysical Data Center.

  pushd /opt/USGS/AVOG2S/ExternalData/Ap_Forecast
  ./get_NGDC 2016
  popd
  ln -s /opt/USGS/AVOG2S/ExternalData/Ap_Forecast/NGDC_NOAA_Archive NGDC

Step 2: Generate a coefficient file.  This example directory has a control file 
example1_genSC.ctr:
2016 12 23 12.0                  ! Specifies date and time
1 1 -150.0 90.0 0.933 6371.229   ! Lon/Lat coordinates will be used
0.5 120                          ! 0.5 degree resolution with SphereHarm order of 120
20 200 2.5                       ! Empirical model from 20-200 km
1                                ! One NWP group
25 1                             ! NCEP ID with one filename
0.0 25.0 1.0                     ! NWP data from 0-25km
NCEP                             ! Name of directory with NCEP data
NGDC                             ! Ap name, in this case, the name of directory with NGDC data
NGDC                             ! F107 name, in this case, the name of directory with NGDC data
G2S_SC_20161223_12Z_wf25.nc      ! output spectral coefficient file name

To generate the coefficient file, run

  /opt/USGS/AVOG2S/bin/g2s_genSC example1_genSC.ctr 

Step 3: Create gridded binary files from the coefficient file.

  /opt/USGS/AVOG2S/bin/g2s_ResampleAtmos G2S_SC_20161223_12Z_wf25.nc 15 1.0 18 2.0 30 5.0

Step 4: Extract a 3-d grid from the gridded binary files.

  /opt/USGS/AVOG2S/bin/g2s_Extract_Grid example1_ext3d.ctr

Step 5: Run GeoAc
For a 1-d analysis, run
  ./GeoAc2D -prop Clev0.met theta_min=0.0 theta_max=45.0 theta_step=0.5 azimuth=41 bounces=10 z_src=1.73 CalcAmp=False
This uses the atmospheric profile at the start coordinate of the grid and uses an azimuth of 41, 
corresponding to the azimuth from Cleveland Volcano to the Dillingham infrasound array.

For a 3-d, range-dependent analysis, run
  ./GeoAcGlobal.RngDep -prop Clev Clev.loclat Clev.loclon theta_min=0.0 theta_max=45.0 theta_step=0.5 bounces=10 lat_src=52.8222 lon_src=-169.945 z_src=1.73 CalcAmp=False


