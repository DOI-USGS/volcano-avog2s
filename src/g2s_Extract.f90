!###############################################################################
!
!      This file is a component of the volcanic infrasound monitoring software
!      written at the U.S. Geological Survey by Hans F. Schwaiger (hschwaiger@usgs.gov)
!      and Alexandra M. Iezzi (amiezzi@alaska.edu).  These programs relies on tools
!      developed for the ash transport and dispersion model Ash3d, written at the
!      U.S. Geological Survey by Hans F. Schwaiger (hschwaiger@usgs.gov), Larry G.
!      Mastin (lgmastin@usgs.gov), and Roger P. Denlinger (roger@usgs.gov).
!
!      The model and its source code are products of the U.S. Federal Government and therefore
!      bear no copyright.  They may be copied, redistributed and freely incorporated 
!      into derivative products.  However as a matter of scientific courtesy we ask that
!      you credit the authors and cite published documentation of this model (below) when
!      publishing or distributing derivative products.
!
!      Schwaiger, H.F., Alexandra M. Iezzi and David Fee;
!         AVO-G2S:  A modified Ground-to-Space model for volcano monitoring in Alaska,
!         submitted. 

!      We make no guarantees, expressed or implied, as to the usefulness of the software
!      and its documentation for any purpose.  We assume no responsibility to provide
!      technical support to users of this software.

!##############################################################################
!
!  G2S_Extract
!
!    This source file is used in the makefile to generate three programs:
!     g2s_Extract_Sonde, g2s_Extract_Xsec, and g2s_Extract_Grid
!    for sampling the gridded binary U,V,T files and created 1,2,and 3-d
!    output files, respectively.  These output files can them be used in
!    infrasound propagation software suce as GeoAc or Art2D.
!
!##############################################################################

      PROGRAM G2S_Extract

      use G2S_globvar

      implicit none

      integer :: i,j,k

      real(kind=8) :: x_in,y_in,zmax
      real(kind=8) :: az
      real(kind=8) :: ds_xsec,dz_prof,len_xsec
      integer      :: ns_xsec,nz_prof
      real(kind=8) :: LatStart,LatEnd
      integer      :: LatCnt
      real(kind=8) :: LongStart,LongEnd
      integer      :: LongCnt
      real(kind=8) :: dlon,dlat

      integer      :: ia,ix,iy

      character(len=130) :: lllinebuffer
      integer :: nargs

      integer      :: ilatlonflag
      real(kind=8) :: x2,y2
      real(kind=8) :: tmp1
      integer      :: ioerr

      ! Set all write options to false by default
      logical :: WRITE_SONDE = .false.
      logical :: WRITE_XSEC  = .false.
      logical :: WRITE_GRID  = .false.

      !===================================================
      !  Set the default values here
        ! Vertical node spacing
        ! This one is a bit coarse
      nz1 = 15
      nz2 = 18
      nz3 = 30
      dz1 = 1.0
      dz2 = 2.0
      dz3 = 5.0
        ! Here the three segments all have dz=1.0
      !nz1 = 50
      !nz2 = 50
      !nz3 = 101
      !dz1 = 1.0
      !dz2 = 1.0
      !dz3 = 1.0

       ! Horizontal grid
      IsLatLon  = .true.
      nxmax_g2s = 720
      dx_g2s    = 0.5
      xmin_g2s  = 0.0
      nymax_g2s = 361
      dy_g2s    = 0.5
      ymin_g2s  = -90.0

      !===================================================

      ! One option should be turned on by preprocessor flag
#ifdef SONDE
      WRITE_SONDE = .true.
#endif

#ifdef XSEC
      WRITE_XSEC  = .true.
#endif

#ifdef GRID
      WRITE_GRID  = .true.
#endif

      ! Error-check
      IF((WRITE_SONDE.and.(WRITE_XSEC .or.WRITE_GRID)).or. &
         (WRITE_XSEC .and.(WRITE_SONDE.or.WRITE_GRID)).or. &
         (WRITE_GRID .and.(WRITE_SONDE.or.WRITE_XSEC)))THEN
        write(*,*)"ERROR: Only one preprocessor flag should be given."
        write(*,*)"WRITE_SONDE = ",WRITE_SONDE
        write(*,*)"WRITE_XSEC  = ",WRITE_XSEC
        write(*,*)"WRITE_GRID  = ",WRITE_GRID
        stop
      ENDIF

!     TEST READ COMMAND LINE ARGUMENTS
      nargs = iargc()

      IF(WRITE_SONDE)THEN
          ! For the g2s_Extract_Sonde branch of this program, we need (minimally)
          !   x_in    : longitude or x of sonde point
          !   y_in    : latitude or y of sonde point
          !   zmax    : maximum altitude of profile
          !   dz_prof : vertical increment
          ! With the 4 command-line arguments, we assume
          !  IsLatLon  = .true.
          !  nx      = 720
          !  dx      = 0.5
          !  xmin    = 0.0
          !  ny      = 361
          !  dy      = 0.5
          !  ymin    = -90.0
          !  nz1     = 15
          !  dz1     = 1.0
          !  nz2     = 18
          !  dz2     = 2.0
          !  nz3     = 30
          !  dz3     = 5.0
          !  vx_fn   = 'U_res.raw'
          !  vy_fn   = 'U_res.raw'
          !  tmp_fn  = 'T_res.raw'
          !  topo_fn = 'etopo.nc'
          ! 
          ! Optionally, we can just feed it a command file with these 4 values
          ! along with all values for all these other variables
        if (nargs.eq.5)THEN
          ! parse the four input parameters
          call getarg(1,lllinebuffer)
          read(lllinebuffer,*)x_in
          write(*,*)"x_in",x_in
          call getarg(2,lllinebuffer)
          read(lllinebuffer,*)y_in
          write(*,*)"y_in",y_in
          call getarg(3,lllinebuffer)
          read(lllinebuffer,*)zmax
          write(*,*)"zmax",zmax
          call getarg(4,lllinebuffer)
          read(lllinebuffer,*)dz_prof
          write(*,*)"dz_prof",dz_prof
          call getarg(5,lllinebuffer)
          FILE_OUT_ROOT = adjustl(trim(lllinebuffer))

        elseif (nargs.eq.1)THEN
          ! Assume single argument is a command file with the format
          !   x_in,y_in  : longitude/x, latitude/y of sonde point
          !   zmax, dz_prof : maximum altitude of output, vertical increment
          !   IsLatLon   : 1 for global (lat/lon) grids, 0 for projected
          !   nx, dx, xmin : number, spacing, and start of gridpoints in x
          !   ny, dy, ymin : number, spacing, and start of gridpoints in y
          !   nz1, dz1   : number and spacing of low alt points
          !   nz2, dz2   : number and spacing of mid alt points
          !   nz3, dz3   : number and spacing of high alt points
          !   vx_fn      : file name of resampled vx values
          !   vy_fn      : file name of resampled vy values
          !   tmp_fn     : file name of resampled temperature values
          !   out_fn     : file name of output sonde data for GeoAc
          !   topo_fn    : (OPTIONAL) topography filename
          call getarg(1,lllinebuffer)
          read(lllinebuffer,*)controlfile
          OPEN(UNIT=ct_unit,FILE=controlfile,STATUS='old')
          read(ct_unit,*)x_in,y_in
          write(*,*)"x_in,y_in: ",x_in,y_in
          read(ct_unit,*)zmax, dz_prof
          write(*,*)"zmax, dz_prof: ",zmax, dz_prof
          read(ct_unit,*)ilatlonflag
          IF(ilatlonflag.eq.1)THEN
            IsLatLon = .true.
            ! If grid is global, then assume the start grid is lon=0 and lat=-90
            !  Only read in two values to specify the grid
            read(ct_unit,*)nxmax_g2s, dx_g2s
            xmin_g2s = 0.0
            write(*,*)"nxmax_g2s, dx, xmin: ",nxmax_g2s, dx_g2s, xmin_g2s
            read(ct_unit,*)nymax_g2s, dy_g2s
            ymin_g2s = -90.0
            write(*,*)"nymax_g2s, dy, ymin: ",nymax_g2s, dy_g2s, ymin_g2s
          ELSE
            IsLatLon = .false.
            ! Grid is projected, we need to read in three values (nx,dx,xmin)
            read(ct_unit,*)nxmax_g2s, dx_g2s, xmin_g2s
            write(*,*)"nxmax_g2s, dx_g2s: ",nxmax_g2s, dx_g2s, xmin_g2s
            read(ct_unit,*)nymax_g2s, dy_g2s, ymin_g2s
            write(*,*)"nymax_g2s, dy_g2s: ",nymax_g2s, dy_g2s, ymin_g2s
          ENDIF

          read(ct_unit,*)nz1, dz1
          write(*,*)"nz1, dz1: ",nz1, dz1
          read(ct_unit,*)nz2, dz2
          write(*,*)"nz2, dz2: ",nz2, dz2
          read(ct_unit,*)nz3, dz3
          write(*,*)"nz3, dz3: ",nz3, dz3
          read(ct_unit,*)FILE_U_RES
          write(*,*)"vx_fn = ",FILE_U_RES
          read(ct_unit,*)FILE_V_RES
          write(*,*)"vy_fn = ",FILE_V_RES
          read(ct_unit,*)FILE_T_RES
          write(*,*)"tmp_fn = ",FILE_T_RES
          read(ct_unit,*)FILE_OUT_SONDE
          FILE_OUT_ROOT = adjustl(trim(FILE_OUT_SONDE))
          write(*,*)"out_fn = ",FILE_OUT_ROOT
          read(ct_unit,*,iostat=ioerr)FILE_TOPO
          IF(ioerr.eq.0)THEN
            useTopo = .true.
            write(*,*)"topo_fn = ",FILE_TOPO
          ELSE
            useTopo = .false.
            FILE_TOPO=""
          ENDIF
          close(ct_unit)
        else
          ! Dump useage info to stdout:
          write(*,*)" "
          write(*,*)"No command-line arguments given"
          write(*,*)" "
          write(*,*)"g2s_Extract_Sonde can be run in two modes:"
          write(*,*)"  (1) command-line arguments"
          write(*,*)"    g2s_Extract_Sonde x_in y_in zmax dz_prof"
          write(*,*)"      x_in    : longitude or x of sonde point"
          write(*,*)"      y_in    : latitude or y of sonde point"
          write(*,*)"      zmax    : maximum altitude of profile"
          write(*,*)"      dz_prof : vertical increment"
          write(*,*)"    This mode assumes a global model with 0.5 degree resolution (720x361)"
          write(*,*)"    and 63 irregularly-spaced nodes in z defining the"
          write(*,*)"    resampled files.  The resampled files are assumed to"
          write(*,*)"    be in the current working directory with names:"
          write(*,*)"     [U,V,T]_res.raw"
          write(*,*)"       Note: these binary (raw) files are not portable and"
          write(*,*)"             must have been written on the same platform as"
          write(*,*)"             g2s_Extract_Sonde"
          write(*,*)"    Also assumed is the output filename: InfraAtmos01.met"
          write(*,*)"    This file is in the format suitable for GeoAc."
          write(*,*)"  (2) control file"
          write(*,*)"    g2s_Extract_Sonde input_sonde.ctr"
          write(*,*)"      where input_sonde.ctr has the following format:"
          write(*,*)"      "
          write(*,*)"        x_in,y_in      #lon,lat (or x,y) of sonde point"
          write(*,*)"        zmax, dz_prof  #maximum altitude of output, vertical increment"
          write(*,*)"        IsLatLon       #0 for global (lat,lon) grids, >0 for projected"
          write(*,*)"        nx, dx, xmin   #number, spacing, and start of gridpoints in x"
          write(*,*)"        ny, dy, ymin   #number, spacing, and start of gridpoints in y"
          write(*,*)"        nz1, dz1       #number and spacing of low alt points"
          write(*,*)"        nz2, dz2       #number and spacing of mid alt points"
          write(*,*)"        nz3, dz3       #number and spacing of high alt points"
          write(*,*)"        vx_fn          #file name of resampled vx values"
          write(*,*)"        vy_fn          #file name of resampled vy values"
          write(*,*)"        tmp_fn         #file name of resampled temperature values"
          write(*,*)"        out_fn         #file name of output sonde data for GeoAc"
          write(*,*)"        topo_fn        #(OPTIONAL) topography filename"
          write(*,*)" "
          write(*,*)" Note: the assumed values for the geometry are:"
          write(*,*)"  IsLatLon  = ",IsLatLon
          write(*,*)"         nx = ",nxmax_g2s
          write(*,*)"         dx = ",dx_g2s
          write(*,*)"       xmin = ",xmin_g2s
          write(*,*)"         ny = ",nymax_g2s
          write(*,*)"         dy = ",dy_g2s
          write(*,*)"       ymin = ",ymin_g2s
          write(*,*)"        nz1 = ",nz1
          write(*,*)"        dz1 = ",dz1
          write(*,*)"        nz2 = ",nz2
          write(*,*)"        dz2 = ",dz2
          write(*,*)"        nz3 = ",nz3
          write(*,*)"        dz3 = ",dz3
          write(*,*)" "
          write(*,*)" For example:"
          write(*,*)" "
          write(*,*)"   Cleveland sonde on a command line (global model)"
          write(*,*)"     ./g2s_Extract_Sonde 190.055 52.8222 180.0 0.2 Clev"
          write(*,*)" "
          write(*,*)"   Cleveland sonde using a control file (global model)"
          write(*,*)"        190.055 52.8222"
          write(*,*)"        180.0 0.2"
          write(*,*)"        1"
          write(*,*)"        720 0.5 0.0"
          write(*,*)"        361 0.5 -90.0"
          write(*,*)"        15 1.0"
          write(*,*)"        18 2.0"
          write(*,*)"        30 5.0"
          write(*,*)"        G2S_SH_20161223_12Z_wf20wf41_U_res.raw"
          write(*,*)"        G2S_SH_20161223_12Z_wf20wf41_V_res.raw"
          write(*,*)"        G2S_SH_20161223_12Z_wf20wf41_T_res.raw"
          write(*,*)"        Clev"
          write(*,*)"        etopo.nc"
          write(*,*)" "
          write(*,*)"   Cleveland sonde using a control file (nam198 model)"
          write(*,*)"        -1363.94 -3758.61"
          write(*,*)"        180.0 1.0"
          write(*,*)"        0"
          write(*,*)"        512 5.953 -2172.9221"
          write(*,*)"        340 5.953 -4214.8032"
          write(*,*)"        15 1.0"
          write(*,*)"        18 2.0"
          write(*,*)"        30 5.0"
          write(*,*)"        G2S_FC_20161223_12Z_wf12wf41_U_res.raw"
          write(*,*)"        G2S_FC_20161223_12Z_wf12wf41_V_res.raw"
          write(*,*)"        G2S_FC_20161223_12Z_wf12wf41_T_res.raw"
          write(*,*)"        Clev"
          write(*,*)"        etopo.nc"
          write(*,*)" "
          write(*,*)"Note: to determine the coordinate of interest for the nam198 grid, use:"
          write(*,*)"      proj +proj=stere  +lon_0=210  +lat_0=90 +k_0=0.933 +R=6371.229"
          write(*,*)" "
          stop
        endif
      ELSEIF(WRITE_XSEC)THEN
          ! For the g2s_Extract_XSEC branch of this program, we need
          ! (minimally)
          !   x_in    : longitude or x of center point of cross-section
          !   y_in    : latitude or y of center point of cross-section
          !   zmax    : maximum altitude of cross-section
          !   dz_prof : vertical increment
          !   az      : azimuth of cross-section
          !   len_xsec: profile length in degrees
          !   ds      : length (degree or km if global v.s. projected) increment of cross-section
          ! With the 7 command-line arguments, we assume
          !  IsLatLon= .true.
          !  nx      = 720
          !  ny      = 361
          !  nz1     = 15
          !  dz1     = 1.0
          !  nz2     = 18
          !  dz2     = 2.0
          !  nz3     = 30
          !  dz3     = 5.0
          ! 
          ! Optionally, we can just feed it a command file with these 7 values
          ! along with all values for all these other variables
        if (nargs.eq.7)THEN
          ! parse the seven input parameters
          call getarg(1,lllinebuffer)
          read(lllinebuffer,*)x_in
          write(*,*)"x_in",x_in
          call getarg(2,lllinebuffer)
          read(lllinebuffer,*)y_in
          write(*,*)"y_in",y_in
          call getarg(3,lllinebuffer)
          read(lllinebuffer,*)zmax
          write(*,*)"zmax",zmax
          call getarg(4,lllinebuffer)
          read(lllinebuffer,*)dz_prof
          write(*,*)"dz_prof",dz_prof
          call getarg(5,lllinebuffer)
          read(lllinebuffer,*)az
          write(*,*)"az",az
          call getarg(6,lllinebuffer)
          read(lllinebuffer,*)len_xsec
          write(*,*)"len_xsec",len_xsec
          call getarg(7,lllinebuffer)
          read(lllinebuffer,*)ds_xsec
          write(*,*)"ds_xsec",ds_xsec

        elseif (nargs.eq.1)THEN
          ! Assume single argument is a command file with the format
          !   x_in,y_in              : lon, lat (or x,y) of reference point of cross-section
          !   zmax, dz_prof          : maximum altitude of output, vertical increment
          !   IsLatLon   : 1 for global (lat/lon) grids, 0 for projected
          !   az, len_xsec, ds_xsec  : azimuth, length and increment of xsec (if global)
          !       - or -
          !   x2, y2, len_xsec, dx_xsec : coordinate of in-line point, length and increment of xsec (if projected)
          !   nx, dx, xmin           : number, spacing, and start of gridpoints in x"
          !   ny, dy, ymin           : number, spacing, and start of gridpoints in y"
          !   nz1, dz1               : number and spacing of low alt points
          !   nz2, dz2               : number and spacing of mid alt points
          !   nz3, dz3               : number and spacing of high alt points
          !   vx_fn                  : file name of resampled vx values
          !   vy_fn                  : file name of resampled vy values
          !   tmp_fn                 : file name of resampled temperature values
          !   out_fn                 : file name of output xsec data for Art2D
          !   topo_out_fn            : file name of output topo data for Art2D
          !   topo_fn                : topography filename
          call getarg(1,lllinebuffer)
          read(lllinebuffer,*)controlfile
          OPEN(UNIT=ct_unit,FILE=controlfile,STATUS='old')
          read(ct_unit,*)x_in,y_in
          write(*,*)"x_in,y_in: ",x_in,y_in
          read(ct_unit,*)zmax, dz_prof
          write(*,*)"zmax, dz_prof: ",zmax, dz_prof
          read(ct_unit,*)ilatlonflag
          IF(ilatlonflag.eq.1)THEN
            IsLatLon = .true.
            ! For global grids, xsec are given by azimuth with the input point as the center
            read(ct_unit,*)az, len_xsec, ds_xsec
            write(*,*)"az, len_xsec, ds_xsec =",az, len_xsec, ds_xsec
            ! If grid is global, then assume the start grid is lon=0 and lat=-90
            !  Only read in two values to specify the grid
            read(ct_unit,*)nxmax_g2s, dx_g2s
            xmin_g2s = 0.0
            write(*,*)"nxmax_g2s, dx, xmin: ",nxmax_g2s, dx_g2s, xmin_g2s
            read(ct_unit,*)nymax_g2s, dy_g2s
            ymin_g2s = -90.0
            write(*,*)"nymax_g2s, dy, ymin: ",nymax_g2s, dy_g2s, ymin_g2s
          ELSE
            IsLatLon = .false.
            ! For projected grids, xsec is defined by a second in-line point
            read(ct_unit,*)x2,y2, len_xsec, ds_xsec
            write(*,*)"x2,y2, len_xsec, ds_xsec =",x2,y2, len_xsec, ds_xsec
            ! Grid is projected, we need to read in three values (nx,dx,xmin)
            read(ct_unit,*)nxmax_g2s, dx_g2s, xmin_g2s
            write(*,*)"nxmax_g2s, dx_g2s: ",nxmax_g2s, dx_g2s, xmin_g2s
            read(ct_unit,*)nymax_g2s, dy_g2s, ymin_g2s
            write(*,*)"nymax_g2s, dy_g2s: ",nymax_g2s, dy_g2s, ymin_g2s
          ENDIF
          read(ct_unit,*)nz1, dz1
          write(*,*)"nz1, dz1: ",nz1, dz1
          read(ct_unit,*)nz2, dz2
          write(*,*)"nz2, dz2: ",nz2, dz2
          read(ct_unit,*)nz3, dz3
          write(*,*)"nz3, dz3: ",nz3, dz3
          read(ct_unit,*)FILE_U_RES
          write(*,*)"vx_fn = ",FILE_U_RES
          read(ct_unit,*)FILE_V_RES
          write(*,*)"vy_fn = ",FILE_V_RES
          read(ct_unit,*)FILE_T_RES
          write(*,*)"tmp_fn = ",FILE_T_RES
          read(ct_unit,*)FILE_OUT_XSEC
          FILE_OUT_ROOT = adjustl(trim(FILE_OUT_XSEC))
          write(*,*)"out_fn = ",FILE_OUT_ROOT
          read(ct_unit,*)FILE_LL
          write(*,*)"out_fn = ",FILE_LL
          read(ct_unit,*,iostat=ioerr)FILE_TOPO
          IF(ioerr.eq.0)THEN
            useTopo = .true.
            write(*,*)"topo_fn = ",FILE_TOPO
          ELSE
            useTopo = .false.
            FILE_TOPO=""
          ENDIF
          close(ct_unit)
        else
          ! Dump usage info to stdout:
          write(*,*)" "
          write(*,*)"No command-line arguments given"
          write(*,*)" "
          write(*,*)"g2s_Extract_Xsec can be run in two modes:"
          write(*,*)"  (1) command-line arguments"
          write(*,*)"    g2s_Extract_Xsec x_in y_in zmax dz_prof az len_xsec dx_xsec"
          write(*,*)"      x_in      : longitude or x of reference point of cross-section"
          write(*,*)"      y_in      : latitude or y of reference point of cross-section"
          write(*,*)"      zmax      : maximum altitude of cross-section"
          write(*,*)"      dz_prof   : vertical increment"
          write(*,*)"      az (x2,y2): azimuth of cross-section (for LatLon) or end point (for projected)"
          write(*,*)"      len_xsec  : length (deg or km) of cross-section"
          write(*,*)"      ds_xsec : degree increment of cross-section"
          write(*,*)"    This mode assumes a global 0.5 degree resolution (720x361)"
          write(*,*)"    and 63 irregularly-spaced nodes in z defining the"
          write(*,*)"    resampled files.  The resampled files are assumed to"
          write(*,*)"    be in the current working directory with names:"
          write(*,*)"     [U,V,T]_res.raw"
          write(*,*)"       Note: these binary (raw) files are not portable and"
          write(*,*)"             must have been written on the same platform as"
          write(*,*)"             g2s_Extract_Xsec"
          write(*,*)"    Also assumed is the output filename: InfraAtmos01.env"
          write(*,*)"    This file is in the format suitable for Art2D."
          write(*,*)"  (2) control file"
          write(*,*)"    g2s_Extract_Xsec input_xsec.ctr"
          write(*,*)"      where input_Xsec.ctr has the following format:"
          write(*,*)"      "
          write(*,*)"        x_in,y_in              #lon,lat (or x,y) of reference point"
          write(*,*)"        zmax, dz_prof          #maximum altitude of output, vertical increment"
          write(*,*)"        IsLatLon               # 1 for global (lat/lon) grids, 0 for projected"
          write(*,*)"        az, len , ds_xsec      #azimuth, length and increment of xsec"
          write(*,*)"            - or - "
          write(*,*)"        x2, y2, len_xsec, dx_xsec # coordinate of in-line point, length and increment of xsec (if projected)"
          write(*,*)"        nx, dx, xmin   #number, spacing, and start of gridpoints in x"
          write(*,*)"        ny, dy, ymin   #number, spacing, and start of gridpoints in y"
          write(*,*)"        nz1, dz1               #number and spacing of low alt points"
          write(*,*)"        nz2, dz2               #number and spacing of mid alt points"
          write(*,*)"        nz3, dz3               #number and spacing of high alt points"
          write(*,*)"        vx_fn                  #file name of resampled vx values"
          write(*,*)"        vy_fn                  #file name of resampled vy values"
          write(*,*)"        tmp_fn                 #file name of resampled temperature values"
          write(*,*)"        out_fn                 #file name of output xsec data for Art2D"
          write(*,*)"        out_topo_fn            #file name of output topo data for Art2D"
          write(*,*)"        topo_fn                #topography filename"
          write(*,*)" "
          write(*,*)" Note: the assumed values for the geometry are:"
          write(*,*)"  IsLatLon  = ",IsLatLon
          write(*,*)"         nx = ",nxmax_g2s
          write(*,*)"         dx = ",dx_g2s
          write(*,*)"       xmin = ",xmin_g2s
          write(*,*)"         ny = ",nymax_g2s
          write(*,*)"         dy = ",dy_g2s
          write(*,*)"       ymin = ",ymin_g2s
          write(*,*)"        nz1 = ",nz1
          write(*,*)"        dz1 = ",dz1
          write(*,*)"        nz2 = ",nz2
          write(*,*)"        dz2 = ",dz2
          write(*,*)"        nz3 = ",nz3
          write(*,*)"        dz3 = ",dz3
          write(*,*)" "
          write(*,*)" For example for Cleveland cross-section with az=41 (for Dillingham array)"
          write(*,*)"    ./g2s_Extract_Xsec 190.055 52.8222 180.0 0.2 41.0 47.69 0.01"
          write(*,*)" "
          write(*,*)"   Cleveland cross-section towards Dillingham using a control file (global model)"
          write(*,*)"        190.055 52.8222"
          write(*,*)"        180.0 0.2"
          write(*,*)"        1"
          write(*,*)"        41.0 47.69 0.01"
          write(*,*)"        720 0.5 0.0"
          write(*,*)"        361 0.5 -90.0"
          write(*,*)"        15 1.0"
          write(*,*)"        18 2.0"
          write(*,*)"        30 5.0"
          write(*,*)"        G2S_SH_20161223_12Z_wf20wf41_U_res.raw"
          write(*,*)"        G2S_SH_20161223_12Z_wf20wf41_V_res.raw"
          write(*,*)"        G2S_SH_20161223_12Z_wf20wf41_T_res.raw"
          write(*,*)"        Clev"
          write(*,*)"        ProfLonLat"
          write(*,*)"        etopo.nc"
          write(*,*)" "
          write(*,*)"   Cleveland cross-section towards Dillingham using a control file (nam198 model)"
          write(*,*)"        -1363.94 -3758.61"
          write(*,*)"        180.0 1.0"
          write(*,*)"        0"
          write(*,*)"        -484.26 -3256.79 1000.0 1.0"
          write(*,*)"        512 5.953 -2172.9221"
          write(*,*)"        340 5.953 -4214.8032"
          write(*,*)"        15 1.0"
          write(*,*)"        18 2.0"
          write(*,*)"        30 5.0"
          write(*,*)"        G2S_FC_20161223_12Z_wf12wf41_U_res.raw"
          write(*,*)"        G2S_FC_20161223_12Z_wf12wf41_V_res.raw"
          write(*,*)"        G2S_FC_20161223_12Z_wf12wf41_T_res.raw"
          write(*,*)"        Clev"
          write(*,*)"        ProfLonLat"
          write(*,*)"        etopo.nc"
          write(*,*)" "
          write(*,*)"Note: to determine the coordinate of interest for the nam198 grid, use:"
          write(*,*)"      proj +proj=stere  +lon_0=210  +lat_0=90 +k_0=0.933 +R=6371.229"
          write(*,*)" "

          stop
        endif

      ELSEIF(WRITE_GRID)THEN
          ! For the g2s_Extract_Grid branch of this program, we need
          ! (minimally)
          !   LatStart  : start latitude of grid
          !   LatEnd    : end latitude of grid
          !   LatCnt    : number of gridpoints along lat
          !   LongStart : start longitude of grid
          !   LongEnd   : end longitude of grid
          !   LongCnt   : number of gridpoints along lon
          !   zmax      : maximum altitude of profile
          !   dz_prof   : vertical increment
          !   out_root  : root of output filename
          ! With the 9 command-line arguments, we assume
          !  nx      = 720
          !  ny      = 361
          !  nz1     = 15
          !  dz1     = 1.0
          !  nz2     = 18
          !  dz2     = 2.0
          !  nz3     = 30
          !  dz3     = 5.0
          !  vx_fn   = 'vx_resamp.raw'
          !  vy_fn   = 'vy_resamp.raw'
          !  tmp_fn  = 'temper_resamp.raw'
          !  out_fn  = 'InfraAtmos[].met'
          !  topo_fn = 'etopo.nc'
          ! 
          ! Optionally, we can just feed it a command file with these 9 values
          ! along with all values for all these other variables
        if (nargs.eq.9)THEN
          ! parse the eight input parameters
          call getarg(1,lllinebuffer)
          read(lllinebuffer,*)LatStart
          write(*,*)"LatStart",LatStart
          call getarg(2,lllinebuffer)
          read(lllinebuffer,*)LatEnd
          write(*,*)"LatEnd",LatEnd
          call getarg(3,lllinebuffer)
          read(lllinebuffer,*)LatCnt
          write(*,*)"LatCnt",LatCnt
          call getarg(4,lllinebuffer)
          read(lllinebuffer,*)LongStart
          !IF(LongStart.gt.180.0)THEN
          !  LongStart = LongStart-360.0
          !  write(*,*)"Longitude remapped to domain -180:180"
          !ENDIF
          write(*,*)"LongStart",LongStart
          call getarg(5,lllinebuffer)
          read(lllinebuffer,*)LongEnd
          !IF(LongEnd.gt.180.0)THEN
          !  LongEnd = LongEnd-360.0
          !  write(*,*)"Longitude remapped to domain -180:180"
          !ENDIF
          write(*,*)"LongEnd",LongEnd
          call getarg(6,lllinebuffer)
          read(lllinebuffer,*)LongCnt
          write(*,*)"LongCnt",LongCnt

          call getarg(7,lllinebuffer)
          read(lllinebuffer,*)zmax
          write(*,*)"zmax",zmax
          call getarg(8,lllinebuffer)
          read(lllinebuffer,*)dz_prof
          write(*,*)"dz_prof",dz_prof

          call getarg(9,lllinebuffer)
          FILE_OUT_ROOT = adjustl(trim(lllinebuffer))

        elseif (nargs.eq.1)THEN
          ! Assume single argument is a command file with the format
          !   LatStart,LatEnd,LatCnt    : start, end and count of lat grid
          !   LongStart,LongEnd,LongCnt : start, end and count of lon grid
          !   zmax, dz_prof             : maximum altitude of output, vertical increment
          !   IsLatLon   : 1 for global (lat/lon) grids, 0 for projected
          !   nx, dx     : number and spacing of gridpoints in x
          !   ny, dy     : number and spacing of gridpoints in y
          !   nz1, dz1   : number and spacing of low alt points
          !   nz2, dz2   : number and spacing of mid alt points
          !   nz3, dz3   : number and spacing of high alt points
          !   vx_fn      : file name of resampled vx values
          !   vy_fn      : file name of resampled vy values
          !   tmp_fn     : file name of resampled temperature values
          !   out_fn     : file name of output sonde data for GeoAc
          !   topo_fn    : (OPTIONAL) topography filename
          call getarg(1,lllinebuffer)
          read(lllinebuffer,*)controlfile
          OPEN(UNIT=ct_unit,FILE=controlfile,STATUS='old')
          read(ct_unit,*)LatStart,LatEnd,LatCnt
          write(*,*)"LatStart,LatEnd,LatCnt: ",LatStart,LatEnd,LatCnt
          read(ct_unit,*)LongStart,LongEnd,LongCnt
          write(*,*)"LongStart,LongEnd,LongCnt: ",LongStart,LongEnd,LongCnt
          read(ct_unit,*)zmax, dz_prof
          write(*,*)"zmax, dz_prof: ",zmax, dz_prof

          read(ct_unit,*)ilatlonflag
          IF(ilatlonflag.eq.1)THEN
            IsLatLon = .true.
            ! If grid is global, then assume the start grid is lon=0 and lat=-90
            !  Only read in two values to specify the grid
            read(ct_unit,*)nxmax_g2s, dx_g2s
            xmin_g2s = 0.0
            write(*,*)"nxmax_g2s, dx, xmin: ",nxmax_g2s, dx_g2s, xmin_g2s
            read(ct_unit,*)nymax_g2s, dy_g2s
            ymin_g2s = -90.0
            write(*,*)"nymax_g2s, dy, ymin: ",nymax_g2s, dy_g2s, ymin_g2s
          ELSE
            IsLatLon = .false.
            ! Grid is projected, we need to read in three values (nx,dx,xmin)
            read(ct_unit,*)nxmax_g2s, dx_g2s, xmin_g2s
            write(*,*)"nxmax_g2s, dx_g2s: ",nxmax_g2s, dx_g2s, xmin_g2s
            read(ct_unit,*)nymax_g2s, dy_g2s, ymin_g2s
            write(*,*)"nymax_g2s, dy_g2s: ",nymax_g2s, dy_g2s, ymin_g2s
          ENDIF

          read(ct_unit,*)nz1, dz1
          write(*,*)"nz1, dz1: ",nz1, dz1
          read(ct_unit,*)nz2, dz2
          write(*,*)"nz2, dz2: ",nz2, dz2
          read(ct_unit,*)nz3, dz3
          write(*,*)"nz3, dz3: ",nz3, dz3
          read(ct_unit,*)FILE_U_RES
          write(*,*)"vx_fn = ",FILE_U_RES
          read(ct_unit,*)FILE_V_RES
          write(*,*)"vy_fn = ",FILE_V_RES
          read(ct_unit,*)FILE_T_RES
          write(*,*)"tmp_fn = ",FILE_T_RES
          read(ct_unit,*)FILE_OUT_ROOT
          write(*,*)"out_fn = ",FILE_OUT_ROOT
          read(ct_unit,*,iostat=ioerr)FILE_TOPO
          IF(ioerr.eq.0)THEN
            useTopo = .true.
            write(*,*)"topo_fn = ",FILE_TOPO
          ELSE
            useTopo = .false.
            FILE_TOPO=""
          ENDIF
          close(ct_unit)
        else
          ! Dump useage info to stdout:
          write(*,*)"No command-line arguments given"
          write(*,*)" "
          write(*,*)"g2s_Extract_Grid can be run in two modes:"
          write(*,*)"  (1) command-line arguments"
          write(*,*)"    g2s_Extract_Grid LatStart LatEnd LatCnt LongStart LongEnd LongCnt zmax dz_prof"
          write(*,*)"      LatStart,LatEnd,LatCnt    : start, end and count of lat grid"
          write(*,*)"      LongStart,LongEnd,LongCnt : start, end and count of lon grid"
          write(*,*)"      zmax                      : maximum altitude of profile"
          write(*,*)"      dz_prof                   : vertical increment"
          write(*,*)"    This mode assumes a 0.5 degree resolution (720x361)"
          write(*,*)"    and 63 irregularly-spaced nodes in z defining the"
          write(*,*)"    resampled files.  The resampled files are assumed to"
          write(*,*)"    be in the current working directory with names:"
          write(*,*)"     [U,V,T]_res.raw"
          write(*,*)"       Note: these binary (raw) files are not portable and"
          write(*,*)"             must have been written on the same platform as"
          write(*,*)"             g2s_Extract_Grid"
          write(*,*)"    Also assumed is the output filename: InfraAtmos0.met"
          write(*,*)"    This file is in the format suitable for GeoAc."
          write(*,*)"  (2) control file"
          write(*,*)"    g2s_Extract_Grid input_grid.ctr"
          write(*,*)"      where input_grid.ctr has the following format:"
          write(*,*)"      "
          write(*,*)"        LatStart,LatEnd,LatCnt    #start, end and count of lat grid"
          write(*,*)"        LongStart,LongEnd,LongCnt #start, end and count of lon grid"
          write(*,*)"        zmax, dz_prof  #maximum altitude of output, vertical increment"
          write(*,*)"        IsLatLon       #0 for global (lat/lon) grids, >0 for projected"
          write(*,*)"        nx, dx, xmin   #number, spacing, and start of gridpoints in x"
          write(*,*)"        ny, dy, ymin   #number, spacing, and start of gridpoints in y"
          write(*,*)"        nz1, dz1       #number and spacing of low alt points"
          write(*,*)"        nz2, dz2       #number and spacing of mid alt points"
          write(*,*)"        nz3, dz3       #number and spacing of high alt points"
          write(*,*)"        vx_fn          #file name of resampled vx values"
          write(*,*)"        vy_fn          #file name of resampled vy values"
          write(*,*)"        tmp_fn         #file name of resampled temperature values"
          write(*,*)"        out_fn         #root file name of output sonde data for GeoAc"
          write(*,*)"        topo_fn        #(OPTIONAL) topography filename"
          write(*,*)" "
          write(*,*)" Note: the assumed values for the geometry are:"
          write(*,*)"  IsLatLon  = ",IsLatLon
          write(*,*)"         nx = ",nxmax_g2s
          write(*,*)"         dx = ",dx_g2s
          write(*,*)"         ny = ",nymax_g2s
          write(*,*)"         dy = ",dy_g2s
          write(*,*)"        nz1 = ",nz1
          write(*,*)"        dz1 = ",dz1
          write(*,*)"        nz2 = ",nz2
          write(*,*)"        dz2 = ",dz2
          write(*,*)"        nz3 = ",nz3
          write(*,*)"        dz3 = ",dz3
          write(*,*)" "
          write(*,*)" For example for Cleveland grid extending NE to includ Dillingham:"
          write(*,*)"    ./g2s_Extract_Grid  50.0 60.0 3 185.0 205.0 5 180.0 0.2 Clev"
          write(*,*)" "
          write(*,*)" This will produce the files InfraAtmos[].met and InfraAtmos.lo[]"
          write(*,*)" These files can be used by GeoAc by the commant:"
          write(*,*)" ./GeoAc3D.RngDep -prop InfraAtmos InfraAtmos.locx InfraAtmos.locy theta_step=2.0 bounces=2 azimuth=41.0"

          stop
        endif

      ELSE
        write(*,*)"ERROR: This executable is not flagged as Sonde, Xsec, or Grid."
        write(*,*)"       Recompile using the appropriate pre-processor flag."
        stop
      ENDIF

      allocate(x_g2s_sp(nxmax_g2s))
      allocate(y_g2s_sp(nymax_g2s))

      DO i=1,nxmax_g2s
        x_g2s_sp(i) = xmin_g2s + (i-1)*real(dx_g2s,kind=4)
      ENDDO
      DO j=1,nymax_g2s
        y_g2s_sp(j) = ymin_g2s + (j-1)*real(dx_g2s,kind=4)
      ENDDO

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      data_len = nz1 + nz2 + nz3
      allocate(vx_OUT_dp(nxmax_g2s,nymax_g2s,data_len))
      allocate(vy_OUT_dp(nxmax_g2s,nymax_g2s,data_len))
      allocate(temperature_OUT_dp(nxmax_g2s,nymax_g2s,data_len))
      allocate(new_data_alt(data_len))

      DO k = 1,nz1
        new_data_alt(k) = real((k-1),kind=4)*dz1
      ENDDO
      DO k = 1,nz2
        new_data_alt(k+nz1) = new_data_alt(nz1)+real(k,kind=4)*dz2
      ENDDO
      DO k = 1,nz3
        new_data_alt(k+nz1+nz2) = new_data_alt(nz1+nz2)+real(k,kind=4)*dz3
      ENDDO
      IF(new_data_alt(nz1+nz2+nz3).le.zmax)THEN
        write(*,*)"ERROR:  The requested zmax must be strictly less than the"
        write(*,*)"reconstitued atmosphere"
        write(*,*)"     requested zmax = ",zmax
        write(*,*)" reconstituted zmax = ",new_data_alt(nz1+nz2+nz3)
        stop
      ENDIF

      ! First read in resampled atmosphere generated from G2S_ResampleAtmos
      write(*,*)"Reading output for resampled data files"
      write(*,*)"  Opening ",FILE_U_RES
      open (unit=20,file=FILE_U_RES,access='direct', &
                   recl=2*nxmax_g2s*nymax_g2s*data_len*4)
      read(20,rec=1)(((vx_OUT_dp(i,j,k),i=1,nxmax_g2s),j=1,nymax_g2s),k=1,data_len)
      close (20)
      write(*,*)"  Opening ",FILE_V_RES
      open (unit=20,file=FILE_V_RES,access='direct', &
                   recl=2*nxmax_g2s*nymax_g2s*data_len*4)
      read(20,rec=1)(((vy_OUT_dp(i,j,k),i=1,nxmax_g2s),j=1,nymax_g2s),k=1,data_len)
      close (20)
      write(*,*)"  Opening ",FILE_T_RES
      open (unit=20,file=FILE_T_RES,access='direct', &
                   recl=2*nxmax_g2s*nymax_g2s*data_len*4)
      read(20,rec=1)(((temperature_OUT_dp(i,j,k),i=1,nxmax_g2s),j=1,nymax_g2s),k=1,data_len)
      close (20)

      IF(useTopo)THEN
        write(*,*)"Reading topography file"
        call Read_Topo
      ENDIF

      nz_prof = floor(zmax/dz_prof)+1

      IF(WRITE_SONDE)THEN
        write(*,*)"Extracting sonde data from resampled atmosphere at ",x_in,y_in
        ia = 1
        call Extract_Sonde(x_in,y_in,dz_prof,nz_prof,ia)
      ENDIF

      IF(WRITE_XSEC)THEN
        IF(IsLatLon)THEN
          write(*,*)"Extracting profile from model at azimuth ",az
        ELSE
          write(*,*)"Extracting profile from ",x_in,y_in
          write(*,*)"                     to ",x2,y2
          az = atan2(x2-x_in,y2-y_in)
        ENDIF
        ! Calculate the number of elements along the cross-section
        ns_xsec = floor(len_xsec/ds_xsec)+1

        call Extract_Profile(x_in,y_in,az,       &
                             ds_xsec,ns_xsec,    &
                             dz_prof,nz_prof)
      ENDIF

      IF(WRITE_GRID)THEN
        write(*,*)"Extracting grid of sonde profiles starting at ",LongStart,LatStart
        dlon = (LongEnd-LongStart)/(LongCnt-1)
        dlat = (LatEnd-LatStart)/(LatCnt-1)
        DO iy = 1,LatCnt
          DO ix = 1,LongCnt
            x_in = LongStart + dlon * (ix-1)
            y_in = LatStart  + dlat * (iy-1)
            ia = LongCnt*(iy-1) + ix
            call Extract_Sonde(x_in,y_in,dz_prof,nz_prof,ia)
          ENDDO
        ENDDO

        ! writing coordinate files needed for GeoAc
        FILE_OUT = adjustl(trim(FILE_OUT_ROOT)) // ".loclat"
        OPEN(unit=20,file=FILE_OUT,status='replace')
        DO iy = 1,LatCnt
          write(20,101)LatStart  + dlat * (iy-1)
        ENDDO
 101    FORMAT(F9.3)
        CLOSE(20)

        FILE_OUT = adjustl(trim(FILE_OUT_ROOT)) // ".loclon"
        OPEN(unit=20,file=FILE_OUT,status='replace')
        DO ix = 1,LongCnt
          tmp1 = LongStart  + dlon * (ix-1)
          IF(IsLatLon.and.tmp1.gt.180.0)tmp1=tmp1-360.0
          write(20,102)tmp1
        ENDDO
        CLOSE(20)
 102    FORMAT(F9.3)
      ENDIF

      write(*,*)"Exited normally."

      END PROGRAM G2S_Extract


!##############################################################################
!##############################################################################

      subroutine Read_Topo

      use G2S_globvar
      use netcdf

      implicit none

      integer(kind=2), dimension(:,:) ,allocatable :: dum2d_short
      integer :: nlat,nlon
      real(kind=4), dimension(:)      ,allocatable :: lat_sp
      real(kind=8), dimension(:)      ,allocatable :: lon_dp
      INTEGER :: nSTAT
      INTEGER :: ncid

      INTEGER :: lat_dim_id
      INTEGER :: lon_dim_id

      INTEGER :: lat_var_id
      INTEGER :: lon_var_id
      INTEGER :: topo_var_id
      integer :: i

        ! ETOPO1 (1-minute global topo/batho)
        ! Since this file is about 450Mb, just load the whole thing
        nlat = 10800
        nlon = 21600
        dlon_topo = 1.0/60.0
        dlat_topo = 1.0/60.0

        nSTAT = nf90_open(FILE_TOPO,NF90_NOWRITE,ncid)
        IF(nSTAT.ne.0)write(6,*)'ERROR: nf90_open to read header:', &
                             nf90_strerror(nSTAT)
        IF(nSTAT.ne.0)THEN
          write(6,*)'Could not open ',FILE_TOPO
          write(6,*)'Exiting'
          stop 1
        ENDIF
        ! Get dimensions
        nSTAT = nf90_inq_dimid(ncid,'lon',lon_dim_id)
        nSTAT = nf90_Inquire_Dimension(ncid,lon_dim_id,len=nlon)
        nSTAT = nf90_inq_dimid(ncid,'lat',lat_dim_id)
        nSTAT = nf90_Inquire_Dimension(ncid,lat_dim_id,len=nlat)
        allocate(lon_dp(nlon))
        allocate(lat_sp(nlat))
        allocate(dum2d_short(nlon,nlat))
        allocate(lon_raw(nlon))
        allocate(lat_raw(nlat))
        allocate(topo_raw(nlon,nlat))

        nSTAT = nf90_inq_varid(ncid,'lon',lon_var_id)
        nSTAT = nf90_get_var(ncid,lon_var_id,lon_dp)

        nSTAT = nf90_inq_varid(ncid,'lat',lat_var_id)
        nSTAT = nf90_get_var(ncid,lat_var_id,lat_sp)
        lat_raw = lat_sp
        nSTAT = nf90_inq_varid(ncid,'z',topo_var_id)
        nSTAT = nf90_get_var(ncid,topo_var_id,dum2d_short)

        ! Start topo_raw at prime meridian
        DO i=1,10800
          topo_raw(i,:) = dum2d_short(i+10800,:)
          lon_raw(i)    = lon_dp(i+10800)
        ENDDO
        ! Map western hemisphere onto second half
        DO i=10801,21600
          topo_raw(i,:) = dum2d_short(i-10800,:)
          lon_raw(i)    = lon_dp(i-10800) + 360.0
        ENDDO

        nSTAT = nf90_close(ncid)
        deallocate(lon_dp,lat_sp,dum2d_short)

      end subroutine Read_Topo

!##############################################################################
!##############################################################################

      subroutine Extract_Sonde(x_in,y_in,dz_prof,nz_prof,ia)

      use G2S_globvar

      implicit none

      real(kind=8),INTENT(IN) :: x_in,y_in
      real(kind=8),INTENT(IN) :: dz_prof
      integer     ,INTENT(IN) :: nz_prof
      integer     ,INTENT(IN) :: ia

      real(kind=8) :: alt
      integer      :: ix,iy,ialt
      real(kind=8) :: xfrac,yfrac,zfrac,xc,yc,zc
      real(kind=8) :: a1,a2,a3,a4

      real(kind=8),dimension(nz_prof) :: z_prof
      integer     ,dimension(nz_prof) :: z_prof_idx

      real(kind=8),allocatable,dimension(:) :: temper_sonde
      real(kind=8),allocatable,dimension(:) :: density_sonde
      real(kind=8),allocatable,dimension(:) :: pressure_sonde
      real(kind=8),allocatable,dimension(:) :: Uwind_sonde
      real(kind=8),allocatable,dimension(:) :: Vwind_sonde
      real(kind=8),allocatable,dimension(:) :: Wwind_sonde

      real(kind=8) :: dz_interval
      integer :: k,kk
      real(kind=8) :: vn,ve
      real(kind=8) :: P_Pa, dens_MKS
      character(len=15) :: FILE_EXT

      allocate(temper_sonde(nz_prof))
      allocate(density_sonde(nz_prof))
      allocate(pressure_sonde(nz_prof))
      allocate(Uwind_sonde(nz_prof))
      allocate(Vwind_sonde(nz_prof))
      allocate(Wwind_sonde(nz_prof))

      ! write the file names to the FILE* strings
      IF(ia.le.10)THEN
        write(FILE_EXT,'(i1,a4)') ia-1,'.met'
        FILE_OUT_Sonde = adjustl(trim(FILE_OUT_ROOT)) // adjustl(trim(FILE_EXT))
      ELSEIF(ia.le.100)THEN
        write(FILE_EXT,'(i2,a4)') ia-1,'.met'
        FILE_OUT_Sonde = adjustl(trim(FILE_OUT_ROOT)) // adjustl(trim(FILE_EXT))
      ENDIF

      ! Set up the output z-profile and find which interval of the reconstituted profile
      ! each of the output levels is in
      DO k = 1,nz_prof
        z_prof(k) = (k-1)*dz_prof
        DO kk = 1,data_len-1
          If(new_data_alt(kk).le.z_prof(k).and.new_data_alt(kk+1).gt.z_prof(k))THEN
            z_prof_idx(k) = kk
          ENDIf
        ENDDO
      ENDDO

      ix = floor((x_in-xmin_g2s)/dx_g2s)+1
      iy = floor((y_in-ymin_g2s)/dy_g2s)+1
      IF(ix.gt.nxmax_g2s-1)THEN
        write(*,*)"ix too large",ix,iy
      ENDIF
      IF(iy.gt.nymax_g2s-1)THEN
        write(*,*)"iy too large",ix,iy
      ENDIF

      xfrac=(x_in-x_g2s_sp(ix))/dx_g2s
      yfrac=(y_in-y_g2s_sp(iy))/dy_g2s
      xc = 1.0-xfrac
      yc = 1.0-yfrac
      a1=xc    * yc
      a2=xfrac * yc
      a3=xfrac * yfrac
      a4=xc    * yfrac
      If(a1.lt.0.0.or.&
         a2.lt.0.0.or.&
         a3.lt.0.0.or.&
         a4.lt.0.0)THEN
         write(*,*)"ERROR:"
         write(*,*)"Problem with negative interpolation coefficients"
         write(*,*)"Point is not mapping to expected grid cell."
         write(*,*)"x,y = ",x_in,y_in
         write(*,*)ix,iy
         write(*,*)x_g2s_sp(ix),y_g2s_sp(iy)
         write(*,*)a1
         write(*,*)a2
         write(*,*)a3
         write(*,*)a4
         write(*,*)xfrac,xc
         write(*,*)yfrac,yc
         stop
      ENDIF
      If(a1.gt.1.0.or.&
         a2.gt.1.0.or.&
         a3.gt.1.0.or.&
         a4.gt.1.0)THEN
         write(*,*)"Problem with large a's"
         write(*,*)a1
         write(*,*)a2
         write(*,*)a3
         write(*,*)a4
         write(*,*)xfrac,xc
         write(*,*)yfrac,yc
         stop
      ENDIF

      DO k = 1,nz_prof
        alt = z_prof(k)
        ialt = z_prof_idx(k)
        dz_interval = new_data_alt(ialt+1)-new_data_alt(ialt)
        zfrac = (alt-new_data_alt(ialt))/dz_interval
        zc = 1.0-zfrac
        IF(zfrac.lt.0.0)THEN
          zfrac = 0.0
          zc = 1.0
        ENDIF
        ve =    zc*a1*vx_OUT_dp(ix  ,iy  ,ialt  ) + &
                zc*a2*vx_OUT_dp(ix+1,iy  ,ialt  ) + &
                zc*a3*vx_OUT_dp(ix+1,iy+1,ialt  ) + &
                zc*a4*vx_OUT_dp(ix  ,iy+1,ialt  ) + &
             zfrac*a1*vx_OUT_dp(ix  ,iy  ,ialt+1) + &
             zfrac*a2*vx_OUT_dp(ix+1,iy  ,ialt+1) + &
             zfrac*a3*vx_OUT_dp(ix+1,iy+1,ialt+1) + &
             zfrac*a4*vx_OUT_dp(ix  ,iy+1,ialt+1)

        vn =    zc*a1*vy_OUT_dp(ix  ,iy  ,ialt  ) + &
                zc*a2*vy_OUT_dp(ix+1,iy  ,ialt  ) + &
                zc*a3*vy_OUT_dp(ix+1,iy+1,ialt  ) + &
                zc*a4*vy_OUT_dp(ix  ,iy+1,ialt  ) + &
             zfrac*a1*vy_OUT_dp(ix  ,iy  ,ialt+1) + &
             zfrac*a2*vy_OUT_dp(ix+1,iy  ,ialt+1) + &
             zfrac*a3*vy_OUT_dp(ix+1,iy+1,ialt+1) + &
             zfrac*a4*vy_OUT_dp(ix  ,iy+1,ialt+1)

        Uwind_sonde(k) =  ve
        Vwind_sonde(k) =  vn
        Wwind_sonde(k) = 0.0

        temper_sonde(k) =    zc*a1*temperature_OUT_dp(ix  ,iy  ,ialt  ) + &
                             zc*a2*temperature_OUT_dp(ix+1,iy  ,ialt  ) + &
                             zc*a3*temperature_OUT_dp(ix+1,iy+1,ialt  ) + &
                             zc*a4*temperature_OUT_dp(ix  ,iy+1,ialt  ) + &
                          zfrac*a1*temperature_OUT_dp(ix  ,iy  ,ialt+1) + &
                          zfrac*a2*temperature_OUT_dp(ix+1,iy  ,ialt+1) + &
                          zfrac*a3*temperature_OUT_dp(ix+1,iy+1,ialt+1) + &
                          zfrac*a4*temperature_OUT_dp(ix  ,iy+1,ialt+1)
      ENDDO

      ! Now calculate the dependent variables
      write(*,*)"Calculating pressure"
      DO k = 1,nz_prof
        ! Calculate pressure in mbar or hPa)
        pressure_sonde(k) = 1013.0 * exp(-z_prof(k)/7.0)
      ENDDO
      write(*,*)"Calculating density",nz_prof
      DO k = 1,nz_prof
        ! Calculate density in g/cm3 (Note pressure is converted to Pa)
        P_Pa = 100.0*pressure_sonde(k)
        dens_MKS = P_Pa/(R_GAS_DRYAIR*temper_sonde(k))
        density_sonde(k) = dens_MKS / 1000.0
      ENDDO

      write(*,*)"Now writing file ",FILE_OUT_Sonde

      OPEN(UNIT=12, file=FILE_OUT_Sonde,status='replace')
      ! write met file
      !  Note: The GeoAc manual give the format as 
      !   z [km] : T(z) [K] : u(z) [m/s] : v(z) [m/s] : rho(z) [g/cm3] : p(z) [mbar]
      !  However, the example file gives
      !   z [km] : T(z) [K] : u(z) [m/s] : v(z) [m/s] : rho(z) [g/cm3] : p(z) [mbar/10 = kPa]

      DO k = 1,nz_prof
        write(12,50)z_prof(k),temper_sonde(k),Uwind_sonde(k),Vwind_sonde(k),&
                              density_sonde(k),pressure_sonde(k)
      ENDDO
 50   FORMAT(6E15.8)

      CLOSE(12)

      end subroutine Extract_Sonde


!##############################################################################
!##############################################################################

      subroutine Extract_Profile(x_in,y_in,az, &
                           ds_xsec,ns_xsec,    &
                           dz_prof,nz_prof)

      use G2S_globvar

      implicit none

      real(kind=8) :: x_in,y_in,az
      real(kind=8) :: ds_xsec,dz_prof
      integer      :: ns_xsec,nz_prof

      real(kind=8) :: x,y
      integer      :: ia,ix,iy

      real(kind=8) :: alt
      real(kind=8) :: lon_start,lat_start
      integer      :: ialt
      real(kind=8) :: xfrac,yfrac,zfrac,xc,yc,zc
      real(kind=8) :: a1,a2,a3,a4

      real(kind=8),dimension(ns_xsec) :: x_prof
      real(kind=8),dimension(ns_xsec) :: y_prof
      real(kind=8),dimension(ns_xsec) :: azm_prof
      real(kind=8),dimension(ns_xsec) :: s_prof
      real(kind=8),dimension(ns_xsec) :: topo_prof
      real(kind=8),dimension(nz_prof) :: z_prof
      integer     ,dimension(nz_prof) :: z_prof_idx

      real(kind=8),allocatable,dimension(:,:,:) :: temper_prof
      real(kind=8),allocatable,dimension(:,:,:) :: density_prof
      real(kind=8),allocatable,dimension(:,:,:) :: pressure_prof
      real(kind=8),allocatable,dimension(:,:,:) :: Uwind_prof
      real(kind=8),allocatable,dimension(:,:,:) :: Vwind_prof
      real(kind=8),allocatable,dimension(:,:,:) :: Wwind_prof

      real(kind=8) :: dz_interval,deg_rot,lon_br,lat_br,lon_pole,lat_pole
      real(kind=8) :: az_loc,topo_loc
      integer :: i,j,k,recli,kk
      real(kind=8) :: azr,azrc,saz,caz,sazc,cazc,vn,ve
      real(kind=8) :: P_Pa, dens_MKS
      real(kind=8) :: xs,ys
      integer      :: is,ns_xsec_tmp
      character(len=15) :: FILE_EXT
      character(len=50) :: FILE_OUT_XSEC_U  = "InfraAtmos01.u"
      character(len=50) :: FILE_OUT_XSEC_V  = "InfraAtmos01.v"
      character(len=50) :: FILE_OUT_XSEC_T  = "InfraAtmos01.t"

      ia = 1
      allocate(  temper_prof(2,ns_xsec,nz_prof))
      allocate( density_prof(2,ns_xsec,nz_prof))
      allocate(pressure_prof(2,ns_xsec,nz_prof))
      allocate(   Uwind_prof(2,ns_xsec,nz_prof))
      allocate(   Vwind_prof(2,ns_xsec,nz_prof))
      allocate(   Wwind_prof(2,ns_xsec,nz_prof))

      ! write the file names to the the FILE* strings
      write(FILE_OUT_XSEC,'(a11,i1,a4)') &
                 'InfraAtmos0',ia,'.env'
      write(FILE_LL,'(a11,i1,a4)') &
                 'ProfLonLat0',ia,'.dat'

      ! write the file names to the FILE* strings
      IF(ia.le.10)THEN
        write(FILE_EXT,'(i1,a2)') ia-1,'.u'
        FILE_OUT_XSEC_U = adjustl(trim(FILE_OUT_ROOT)) // adjustl(trim(FILE_EXT))
        write(FILE_EXT,'(i1,a2)') ia-1,'.v'
        FILE_OUT_XSEC_V = adjustl(trim(FILE_OUT_ROOT)) // adjustl(trim(FILE_EXT))
        write(FILE_EXT,'(i1,a2)') ia-1,'.t'
        FILE_OUT_XSEC_T = adjustl(trim(FILE_OUT_ROOT)) // adjustl(trim(FILE_EXT))
      ELSEIF(ia.le.100)THEN
        write(FILE_EXT,'(i2,a2)') ia-1,'.u'
        FILE_OUT_XSEC_U = adjustl(trim(FILE_OUT_ROOT)) // adjustl(trim(FILE_EXT))
        write(FILE_EXT,'(i2,a2)') ia-1,'.v'
        FILE_OUT_XSEC_V = adjustl(trim(FILE_OUT_ROOT)) // adjustl(trim(FILE_EXT))
        write(FILE_EXT,'(i2,a2)') ia-1,'.t'
        FILE_OUT_XSEC_T = adjustl(trim(FILE_OUT_ROOT)) // adjustl(trim(FILE_EXT))
      ENDIF

      ! Calculate x_prof and y_prof
      IF(IsLatLon)THEN
        ! The x_in/y_in will be the center of the profile
        ! To get the lon/lat coordinates of each point in the profile, we rotate
        ! the x_in/y_in point in the direction of az by ds_xsec increments
        ! Get pole of plane defining azimuth vector at this lon/lat
        call get_pole(x_in,y_in,az,lon_pole,lat_pole)
        ! Now get start lon/lat of this profile
        deg_rot = -0.5*(ns_xsec-1)*ds_xsec
        call back_rotate(x_in,y_in,lon_pole,lat_pole,deg_rot,&
                         lon_start,lat_start)
        OPEN(UNIT=12, file=FILE_LL,status='replace')
        DO i = 1,ns_xsec
          deg_rot = i*ds_xsec
          call back_rotate(lon_start,lat_start,lon_pole,lat_pole,deg_rot,lon_br,lat_br)
          If(lon_br.lt.0.0)THEN
            lon_br=lon_br+360.0
          ENDIF
          x_prof(i)=lon_br
          y_prof(i)=lat_br
          call get_local_azim(lon_br,lat_br,lon_pole,lat_pole,az_loc)
          azm_prof(i) = az_loc
          s_prof(i)   = RAD_EARTH*(i-1)*ds_xsec*DEG2RAD
          IF(useTopo)THEN
            call interpolate_topo(lon_br,lat_br,topo_loc)
          ELSE
            topo_loc = 0.0
          ENDIF
          topo_prof(i)= topo_loc/1000.0
          write(12,'(4f25.10)')s_prof(i),deg_rot,x_prof(i),y_prof(i)
        ENDDO
        CLOSE(12)
      ELSE
        ! Here is the projected case
        !  First, use the 'azimuth' from reference point to x2,y2
        !az = atan(x2-x_in,y2-y_in)
        ! User requested ns_xsec values with a spacing of ds_xsec
        ! Check to see if this fits in the domain
        IF(x_in+sin(az)*ns_xsec*ds_xsec.gt.x_g2s_sp(nxmax_g2s).or.&
           x_in+sin(az)*ns_xsec*ds_xsec.lt.x_g2s_sp(        1).or.&
           y_in+cos(az)*ns_xsec*ds_xsec.gt.y_g2s_sp(nymax_g2s).or.&
           y_in+cos(az)*ns_xsec*ds_xsec.gt.y_g2s_sp(        1))THEN
          write(*,*)"Profile length extends outside the domain."
          write(*,*)" Calculating a truncated profile."
        ENDIF
        OPEN(UNIT=12, file=FILE_LL,status='replace')
        DO is = 1,ns_xsec
          xs = x_in+sin(az)*(is-1)*ds_xsec
          ys = y_in+cos(az)*(is-1)*ds_xsec
          IF(xs.gt.x_g2s_sp(nxmax_g2s).or.&
             xs.lt.x_g2s_sp(        1).or.&
             ys.gt.y_g2s_sp(nymax_g2s).or.&
             ys.gt.y_g2s_sp(        1))THEN
            x_prof(is) = xs
            y_prof(is) = ys
            ns_xsec_tmp = is
            azm_prof(is)= az/DEG2RAD
            s_prof(is)  = (is-1)*ds_xsec
            IF(useTopo)THEN
              !call interpolate_topo(lon_br,lat_br,topo_loc)
              topo_loc = 0.0
            ELSE
              topo_loc = 0.0
            ENDIF
            topo_prof(is)= topo_loc/1000.0
            write(12,'(4f25.10)')s_prof(is),s_prof(is),x_prof(is),y_prof(is)
          ELSE
            cycle
          ENDIF
        ENDDO
        CLOSE(12)
        IF(ns_xsec_tmp.lt.ns_xsec)THEN
          ! Redefine the length of the cross section
          ns_xsec = ns_xsec_tmp
        ENDIF
      ENDIF


      DO k = 1,nz_prof
        z_prof(k) = (k-1)*dz_prof
        DO kk = 1,data_len-1
          If(new_data_alt(kk).le.z_prof(k).and.new_data_alt(kk+1).gt.z_prof(k))THEN
            z_prof_idx(k) = kk
          ENDIf
        ENDDO
      ENDDO

      DO i = 1,ns_xsec
        x  = x_prof(i)
        y  = y_prof(i)
        ix = floor((x-x_g2s_sp(1))/dx_g2s)+1
        iy = floor((y-y_g2s_sp(1))/dy_g2s)+1
        IF(ix.gt.nxmax_g2s-1)THEN
          write(*,*)"ix too large",ix,x
        ENDIF
        IF(iy.gt.nymax_g2s)THEN
          write(*,*)"iy too large",iy,y
        ENDIF

        xfrac=(x-x_g2s_sp(ix))/dx_g2s
        yfrac=(y-y_g2s_sp(iy))/dx_g2s
        xc = 1.0-xfrac
        yc = 1.0-yfrac
        a1=xc    * yc
        a2=xfrac * yc
        a3=xfrac * yfrac
        a4=xc    * yfrac
        If(a1.lt.0.0.or.&
           a2.lt.0.0.or.&
           a3.lt.0.0.or.&
           a4.lt.0.0)THEN
           write(*,*)"ERROR:"
           write(*,*)"Problem with negative interpolation coefficients"
           write(*,*)"Point is not mapping to expected grid cell."
           write(*,*)"i = ",i,ns_xsec
           write(*,*)"x,y = ",x,y
           write(*,*)ix,iy
           write(*,*)x_g2s_sp(ix),y_g2s_sp(iy)
           write(*,*)a1
           write(*,*)a2
           write(*,*)a3
           write(*,*)a4
           write(*,*)xfrac,xc
           write(*,*)yfrac,yc
           stop
        ENDIF
        If(a1.gt.1.0.or.&
           a2.gt.1.0.or.&
           a3.gt.1.0.or.&
           a4.gt.1.0)THEN
           write(*,*)"Problem with large a's"
           write(*,*)a1
           write(*,*)a2
           write(*,*)a3
           write(*,*)a4
           write(*,*)xfrac,xc
           write(*,*)yfrac,yc
           stop
        ENDIF

        azr   = azm_prof(i)*DEG2RAD
        azrc  = (90.0-azm_prof(i))*DEG2RAD

        saz = sin(azr)
        caz = cos(azr)
        sazc= sin(azrc)
        cazc= cos(azrc)

        DO k = 1,nz_prof
          alt = z_prof(k)
          ialt = z_prof_idx(k)
          dz_interval = new_data_alt(ialt+1)-new_data_alt(ialt)
          zfrac = (alt-new_data_alt(ialt))/dz_interval
          zc = 1.0-zfrac
          IF(zfrac.lt.0.0)THEN
            zfrac = 0.0
            zc = 1.0
          ENDIF

          ve =    zc*a1*vx_OUT_dp(ix  ,iy  ,ialt  ) + &
                  zc*a2*vx_OUT_dp(ix+1,iy  ,ialt  ) + &
                  zc*a3*vx_OUT_dp(ix+1,iy+1,ialt  ) + &
                  zc*a4*vx_OUT_dp(ix  ,iy+1,ialt  ) + &
               zfrac*a1*vx_OUT_dp(ix  ,iy  ,ialt+1) + &
               zfrac*a2*vx_OUT_dp(ix+1,iy  ,ialt+1) + &
               zfrac*a3*vx_OUT_dp(ix+1,iy+1,ialt+1) + &
               zfrac*a4*vx_OUT_dp(ix  ,iy+1,ialt+1)

          vn =    zc*a1*vy_OUT_dp(ix  ,iy  ,ialt  ) + &
                  zc*a2*vy_OUT_dp(ix+1,iy  ,ialt  ) + &
                  zc*a3*vy_OUT_dp(ix+1,iy+1,ialt  ) + &
                  zc*a4*vy_OUT_dp(ix  ,iy+1,ialt  ) + &
               zfrac*a1*vy_OUT_dp(ix  ,iy  ,ialt+1) + &
               zfrac*a2*vy_OUT_dp(ix+1,iy  ,ialt+1) + &
               zfrac*a3*vy_OUT_dp(ix+1,iy+1,ialt+1) + &
               zfrac*a4*vy_OUT_dp(ix  ,iy+1,ialt+1)

          ! Along Profile winds
          Uwind_prof(1,i,k) =  vn*caz + ve*cazc

          ! Cross Profile winds
          Vwind_prof(1,i,k) = -vn*saz + ve*sazc
          Vwind_prof(1,i,k) = -Vwind_prof(1,i,k)

          Wwind_prof(1,i,k) = 0.0

          temper_prof(1,i,k) =    zc*a1*temperature_OUT_dp(ix  ,iy  ,ialt  ) + &
                                  zc*a2*temperature_OUT_dp(ix+1,iy  ,ialt  ) + &
                                  zc*a3*temperature_OUT_dp(ix+1,iy+1,ialt  ) + &
                                  zc*a4*temperature_OUT_dp(ix  ,iy+1,ialt  ) + &
                               zfrac*a1*temperature_OUT_dp(ix  ,iy  ,ialt+1) + &
                               zfrac*a2*temperature_OUT_dp(ix+1,iy  ,ialt+1) + &
                               zfrac*a3*temperature_OUT_dp(ix+1,iy+1,ialt+1) + &
                               zfrac*a4*temperature_OUT_dp(ix  ,iy+1,ialt+1)
        ENDDO
      ENDDO

      Uwind_prof(2,:,:)  =  Uwind_prof(1,:,:) 
      Vwind_prof(2,:,:)  =  Vwind_prof(1,:,:)
      Wwind_prof(2,:,:)  =  Wwind_prof(1,:,:)
      temper_prof(2,:,:) = temper_prof(1,:,:)

      ! Now calculate the dependent variables
      write(*,*)"Calculating pressure"
      DO k = 1,nz_prof
        ! Calculate pressure in mbar or hPa)
        pressure_prof(:,:,k) = 1013.0 * exp(-z_prof(k)/7.4)
      ENDDO
      write(*,*)"Calculating density",ns_xsec,nz_prof
      DO i = 1,ns_xsec
        DO k = 1,nz_prof
          ! Calculate density in g/cm3 (Note pressure is converted to Pa)
          P_Pa = 100.0*pressure_prof(1,i,k)
          dens_MKS = P_Pa/(R_GAS_DRYAIR*temper_prof(1,i,k))
          density_prof(1,i,k) = dens_MKS / 1000.0

          P_Pa = 100.0*pressure_prof(2,i,k)
          dens_MKS = P_Pa/(R_GAS_DRYAIR*temper_prof(2,i,k))
          density_prof(2,i,k) = dens_MKS / 1000.0
        ENDDO
      ENDDO

      write(*,*)"Now writing files"
    
      ! open file for direct access for 32-bit words
      OPEN(UNIT=12, file=FILE_OUT_XSEC,status='replace', &
           recl=1*4, access='direct',form='unformatted')
      write(12,rec=1) ns_xsec
      write(12,rec=2) nz_prof
      CLOSE(12)

      ! open file for direct access for 64-bit words
      OPEN(UNIT=12, file=FILE_OUT_XSEC, &
           recl=2*4, access='direct',form='unformatted')
      recli = 1
      ! write latitude
      DO i = 1,ns_xsec
        recli = (recli+1)
        write(12,rec=recli)y_prof(i)
      ENDDO
      ! write longitude
      DO i = 1,ns_xsec
        recli = (recli+1)
        write(12,rec=recli)x_prof(i)
      ENDDO
      ! write azimuth
      DO i = 1,ns_xsec
        recli = (recli+1)
        write(12,rec=recli)azm_prof(i)
      ENDDO
      ! write x-distance along profile
      DO i = 1,ns_xsec
        recli = (recli+1)
        write(12,rec=recli)s_prof(i)
      ENDDO
      ! write topography
      DO i = 1,ns_xsec
        recli = (recli+1)
        write(12,rec=recli)topo_prof(i)
      ENDDO
      ! write altitude
      DO k = 1,nz_prof
        recli = (recli+1)
        write(12,rec=recli)z_prof(k)
      ENDDO

      DO j = 1,2
        ! write temperature
        DO i = 1,ns_xsec
          DO k = 1,nz_prof
            recli = (recli+1)
            write(12,rec=recli)temper_prof(j,i,k)
          ENDDO
        ENDDO
        ! write density
        DO i = 1,ns_xsec
          DO k = 1,nz_prof
            recli = (recli+1)
            write(12,rec=recli)density_prof(j,i,k)
          ENDDO
        ENDDO
        ! write pressure
        DO i = 1,ns_xsec
          DO k = 1,nz_prof
            recli = (recli+1)
            write(12,rec=recli)pressure_prof(j,i,k)
          ENDDO
        ENDDO
        ! write Uwind
        DO i = 1,ns_xsec
          DO k = 1,nz_prof
            recli = (recli+1)
            write(12,rec=recli)Uwind_prof(j,i,k)
            !IF(j.eq.1)write(*,*)Uwind_prof(j,i,k)
          ENDDO
        ENDDO
        ! write Vwind
        DO i = 1,ns_xsec
          DO k = 1,nz_prof
            recli = (recli+1)
            write(12,rec=recli)Vwind_prof(j,i,k)
          ENDDO
        ENDDO
        ! write Wwind
        DO i = 1,ns_xsec
          DO k = 1,nz_prof
            recli = (recli+1)
            write(12,rec=recli)Wwind_prof(j,i,k)
          ENDDO
        ENDDO
      ENDDO

      CLOSE(12)

      ! Now dump out the ASCII versions
      write(*,*)"Writing ASCII files of dimensions ",ns_xsec,nz_prof
      OPEN(UNIT=30, file=FILE_OUT_XSEC_U,status='replace')
      OPEN(UNIT=31, file=FILE_OUT_XSEC_V,status='replace')
      OPEN(UNIT=32, file=FILE_OUT_XSEC_T,status='replace')
      DO i = 1,ns_xsec
        DO k = 1,nz_prof
          write(30,'(f10.4)')Uwind_prof(1,i,k)
          write(31,'(f10.4)')Vwind_prof(1,i,k)
          write(32,'(f10.4)')temper_prof(1,i,k)
        ENDDO
      ENDDO

      CLOSE(30)
      CLOSE(31)
      CLOSE(32)

      end subroutine Extract_Profile

!##############################################################################
!##############################################################################

      subroutine interpolate_topo(lon,lat,topo_loc)

       use G2S_globvar

       implicit none

       REAL(KIND=8),INTENT(IN)  :: lon
       REAL(KIND=8),INTENT(IN)  :: lat
       REAL(KIND=8),INTENT(OUT) :: topo_loc

       INTEGER :: ilon,ilat,ilon2,ilat2

       REAL(KIND=8) :: xfrac,yfrac,xc,yc
       REAL(KIND=8) ::  a1,a2,a3,a4,t1,t2,t3,t4

        ilon = floor(lon/dlon_topo)+1
        !IF(Is_y_inverted)THEN
        !  ilat = floor((90.0-lat)/dlat)-1
        !ELSE
          ilat = floor((lat+90.0)/dlat_topo)+1
        !ENDIF

        IF(ilon.eq.21600)THEN
          ilon2 = 1
        ELSEIF(ilon.eq.0)THEN
          ilon = 21600
          ilon2 = 1
        ELSE
          ilon2 = ilon+1
        ENDIF
        IF(ilat.eq.10800)THEN
          ilat2 = 10800
        ELSEIf(ilat.eq.0)THEN
          ilat = 1
          ilat2 = ilat
        ELSE
          ilat2 = ilat+1
        ENDIF

        xfrac=(lon-lon_raw(ilon))/dlon_topo
        yfrac=(lat-lat_raw(ilat))/dlat_topo
        xc = 1.0-xfrac
        yc = 1.0-yfrac
        a1 = xc*yc
        a2 = xfrac*yc
        a3 = xfrac*yfrac
        a4 = yfrac*xc
        t1 = real(topo_raw(ilon ,ilat ),kind=8)
        t2 = real(topo_raw(ilon2,ilat ),kind=8)
        t3 = real(topo_raw(ilon2,ilat2),kind=8)
        t4 = real(topo_raw(ilon ,ilat2),kind=8)

        t1 = max(t1,0.0)
        t2 = max(t2,0.0)
        t3 = max(t3,0.0)
        t4 = max(t4,0.0)

        topo_loc = a1*t1 + a2*t2 + a3*t3 + a4*t4

      end subroutine interpolate_topo

!##############################################################################
!##############################################################################

      function day_of_year(iyear,imonth,iday)

      implicit none

      real(kind=4) :: day_of_year
      integer :: iyear,imonth,iday
      logical :: IsLeapYear
      integer, dimension(12) :: monthdays 

      if  ((mod(iyear,4).eq.0)     .and.                          &
           (mod(iyear,100).ne.0).or.(mod(iyear,400).eq.0))then
        IsLeapYear = .true.
        monthdays = (/31,29,31,30,31,30,31,31,30,31,30,31/)
      else
        IsLeapYear = .false.
        monthdays = (/31,28,31,30,31,30,31,31,30,31,30,31/)
      endif

      IF(imonth.eq.1)THEN
        day_of_year = iday
      ELSE
        day_of_year = sum(monthdays(1:imonth-1)) + iday
      ENDIF

      end function day_of_year

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine get_pole(lon,lat,az,lon_pole,lat_pole)

       use G2S_globvar

       implicit none

       REAL(KIND=8) :: lon,lat
       REAL(KIND=8) :: az
       REAL(KIND=8) :: lon_pole,lat_pole

       REAL(KIND=8),dimension(3,3) :: Rot
       REAL(KIND=8),dimension(3)   :: xvec,bvec,avec,rvec
       REAL(KIND=8) :: phi,lam
       REAL(KIND=8) :: slm,clm,sph,cph
       REAL(KIND=8) :: urx,ury,urz
       REAL(KIND=8) :: omega,somg,comg
       REAL(KIND=8) :: blen,alen,rlen,length,y,deg

       ! Get details of position vector
       phi  = lat*DEG2RAD
       lam  = lon*DEG2RAD
       slm = sin(lam)
       clm = cos(lam)
       sph = sin(phi)
       cph = cos(phi)

       !  % Get unit sphere coordinates of position vector in standard orientation
       xvec(1) = clm*cph
       xvec(2) = slm*cph
       xvec(3) = sph
       ! Here's the vector in the tangent plane pointing along longitude
       ! b = z cross x
       bvec(1) = -xvec(2)
       bvec(2) =  xvec(1)
       bvec(3) = 0.0
       blen = sqrt(dot_product(bvec,bvec))
       bvec = bvec/blen
       ! Now rotate bvec by (az-90) about xvec
       urx = clm*cph
       ury = slm*cph
       urz = sph
       omega = -(az-90.0) * DEG2RAD
       somg = sin(omega)
       comg = cos(omega)
       !% build rotation matrix
       !% (see http://en.wikipedia.org/wiki/Rotation_matrix : Rotation matrix
       !given an axis and an angle)
       Rot(1,1) = urx*urx + (1.0-urx*urx)*comg
       Rot(1,2) = urx*ury*(1.0-comg) - urz*somg
       Rot(1,3) = urx*urz*(1.0-comg) + ury*somg
       Rot(2,1) = urx*ury*(1.0-comg) + urz*somg
       Rot(2,2) = ury*ury + (1.0-ury*ury)*comg
       Rot(2,3) = ury*urz*(1.0-comg) - urx*somg
       Rot(3,1) = urx*urz*(1.0-comg) - ury*somg
       Rot(3,2) = ury*urz*(1.0-comg) + urx*somg
       Rot(3,3) = urz*urz + (1.0-urz*urz)*comg
       ! Now get the azimuth vector
       avec = matmul(Rot,bvec)
       alen = sqrt(dot_product(avec,avec))
       avec = avec/alen

       ! The rotation axis for this lon,lat,az will be given by x cross a   
       rvec(1) = (xvec(2)*avec(3) - xvec(3)*avec(2))
       rvec(2) =-(xvec(1)*avec(3) - xvec(3)*avec(1))
       rvec(3) = (xvec(1)*avec(2) - xvec(2)*avec(1))
       rlen = sqrt(dot_product(rvec,rvec))
       rvec = rvec/rlen

       lat_pole = max(-90.0,min(90.0,asin(rvec(3))/DEG2RAD))
       If((abs(lat_pole-90.0)).gt.1.0e-5)THEN
         length = sqrt(rvec(1)*rvec(1) + rvec(2)*rvec(2) )
         y = abs(rvec(2))
         deg=asin(y/length)/DEG2RAD
         If(rvec(1).ge.0.0.and.rvec(2).ge.0.0)THEN
           ! First quadrant
           !write(*,*)"First quadrant"
           lon_pole = deg
         ELSEIF(rvec(1).lt.0.0.and.rvec(2).ge.0.0)THEN
           ! Second quadrant
           !write(*,*)"Second quadrant"
           lon_pole = 180.0-deg
         ELSEIF(rvec(1).lt.0.0.and.rvec(2).lt.0.0)THEN
           ! Third quadrant
           !write(*,*)"Third quadrant"
           lon_pole = 180.0+deg
         ELSEIF(rvec(1).gt.0.0.and.rvec(2).lt.0.0)THEN
           ! Fourth quadrant
           !write(*,*)"Forth quadrant"
           lon_pole = 360.0-deg
         ELSE
           write(*,*)"Did not find quadrant"
           stop
         ENDIF
       else
         lon_pole = 0.0
       endif

      end subroutine get_pole

!##############################################################################
!##############################################################################

       subroutine get_local_azim(lon_br,lat_br,lon_pole,lat_pole,az_loc)

       use G2S_globvar

       implicit none

       REAL(KIND=8) :: lon_br,lat_br
       REAL(KIND=8) :: lon_pole,lat_pole
       REAL(KIND=8) :: az_loc

       REAL(KIND=8),dimension(3)   :: xvec,zvec,avec,rvec
       REAL(KIND=8) :: phi,lam
       REAL(KIND=8) :: slm,clm,sph,cph
       REAL(KIND=8) :: phir,lamr
       REAL(KIND=8) :: slmr,clmr,sphr,cphr
       REAL(KIND=8) :: zlen

       ! Get details of position vector
       phi  = lat_br*DEG2RAD
       lam  = lon_br*DEG2RAD
       slm = sin(lam)
       clm = cos(lam)
       sph = sin(phi)
       cph = cos(phi)

       !  % Get unit sphere coordinates of position vector in standard
       !  orientation
       xvec(1) = clm*cph
       xvec(2) = slm*cph
       xvec(3) = sph

       ! Get details of rotation vector
       phir  = lat_pole*DEG2RAD
       lamr  = lon_pole*DEG2RAD
       slmr = sin(lamr)
       clmr = cos(lamr)
       sphr = sin(phir)
       cphr = cos(phir)

       !  % Get unit sphere coordinates of position vector in standard
       !  orientation
       rvec(1) = clmr*cphr
       rvec(2) = slmr*cphr
       rvec(3) = sphr

       ! Local azimuth vector is r cross x
       ! The rotation axis for this lon,lat,az will be given by x cross a   
       avec(1) = (rvec(2)*xvec(3) - rvec(3)*xvec(2))
       avec(2) =-(rvec(1)*xvec(3) - rvec(3)*xvec(1))
       avec(3) = (rvec(1)*xvec(2) - rvec(2)*xvec(1))

       zvec(1) = 0.0
       zvec(2) = 0.0
       zvec(3) = 1.0
       ! Project zvec into local tangent plane
       zvec = zvec - dot_product(zvec,xvec)*xvec
       zlen = sqrt(dot_product(zvec,zvec))
       zvec = zvec/zlen
       az_loc = acos(min(max(dot_product(zvec,avec),-1.0),1.0))/DEG2RAD

       end subroutine get_local_azim

!##############################################################################
!##############################################################################

       subroutine back_rotate(lon,lat,    &
                              lon_pole,lat_pole,deg_rot, &
                              lon_br,lat_br)

       use G2S_globvar

       implicit none

       REAL(KIND=8) :: lon,lat
       REAL(KIND=8) :: lon_pole,lat_pole
       REAL(KIND=8) :: deg_rot
       REAL(KIND=8) :: lon_br,lat_br

       REAL(KIND=8),dimension(3,3) :: Rot
       REAL(KIND=8),dimension(3)   :: xvec,xnew

       REAL(KIND=8) :: clm,clmr,comg,cph,cphr
       REAL(KIND=8) :: lam,lamr,omega
       REAL(KIND=8) :: phi,phir
       REAL(KIND=8) :: slm,slmr,somg,sph,sphr
       REAL(KIND=8) :: urx,ury,urz
       REAL(KIND=8) :: x_usph,y_usph,z_usph
       REAL(KIND=8) :: xcolat,xlm,xph

       ! back_rotate rotates the point lon,lat by an angle deg_rot about
       ! the axis lon_pole,lat_pole
       ! returns the lon,lat coordinates of the rotated point

       ! Get details of position vector
       phi  = lat*DEG2RAD
       lam  = lon*DEG2RAD
       slm = sin(lam)
       clm = cos(lam)
       sph = sin(phi)
       cph = cos(phi)

       ! Get details of rotation axis
       phir = lat_pole*DEG2RAD
       lamr = lon_pole*DEG2RAD
       slmr = sin(lamr)
       clmr = cos(lamr)
       sphr = sin(phir)
       cphr = cos(phir)
       urx = clmr*cphr
       ury = slmr*cphr
       urz = sphr

       omega = deg_rot * DEG2RAD
       somg = sin(omega)
       comg = cos(omega)

       !% build rotation matrix
       !% (see http://en.wikipedia.org/wiki/Rotation_matrix : Rotation matrix
       !given an axis and an angle)
       Rot(1,1) = urx*urx + (1.0-urx*urx)*comg
       Rot(1,2) = urx*ury*(1.0-comg) - urz*somg
       Rot(1,3) = urx*urz*(1.0-comg) + ury*somg
       Rot(2,1) = urx*ury*(1.0-comg) + urz*somg
       Rot(2,2) = ury*ury + (1.0-ury*ury)*comg
       Rot(2,3) = ury*urz*(1.0-comg) - urx*somg
       Rot(3,1) = urx*urz*(1.0-comg) - ury*somg
       Rot(3,2) = ury*urz*(1.0-comg) + urx*somg
       Rot(3,3) = urz*urz + (1.0-urz*urz)*comg

       !  % Get unit sphere coordinates in standard orientation
       x_usph = clm*cph
       y_usph = slm*cph
       z_usph = sph
       !  % Back-rotate coordinates
       xvec(1) = x_usph
       xvec(2) = y_usph
       xvec(3) = z_usph
       xnew = matmul(Rot,xvec)
       !  % Get lon/lat for this new position
       xcolat = acos(min(max(xnew(3),-1.0),1.0))
       xph    = (90.0*DEG2RAD)-xcolat
       xlm    = atan2(xnew(2),xnew(1))

       lon_br = xlm/DEG2RAD
       lat_br = xph/DEG2RAD

       end subroutine back_rotate

