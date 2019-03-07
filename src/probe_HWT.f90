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
!         AVO-G2S:  A modified, open-source Ground-to-Space atmospheric specifications
!           for infrasound model; submitted.

!      We make no guarantees, expressed or implied, as to the usefulness of the software
!      and its documentation for any purpose.  We assume no responsibility to provide
!      technical support to users of this software.

!##############################################################################
!
!  Probe_HWT
!
!    This is a stand-alone program that exports a vertical profile of 
!    Z, T, U, V, RelHum, and P using the empirical models HWM[07/14] and
!    NRLMSISE-00.  This program takes 10 command-line arguments:
!     
!       year : integer
!       month: integer
!       day  : integer
!       hour : real
!       lon  : longitude of sonde point
!       lat  : latitude of sonde point
!       ap   : Ap planetary index
!       f107 : F107 value
!       zmax : maximum altitude of profile
!       dz   : vertical increment
!
!    probe_HWT 2013 5 13 19.6 -161.887 55.42 4.0 153.5 200.0 2.5
!
!    Output is written to the file sonde_ztuvrp.dat
!
!##############################################################################

      program probe_HWT

      use G2S_globvar

      implicit none

      integer :: k

      real(kind=8) :: lon_in,lat_in,zmax
      real(kind=8) :: dz_prof
      integer      :: nz_prof
      integer :: kp1,kp2,kp3,kp4,kp5,kp6,kp7,kp8,kp9
      integer :: ap1,ap2,ap3,ap4,ap5,ap6,ap7,ap8,ap9
      integer :: C9,ssnum,f107q
      real(kind=4) :: Cp
      integer :: day
      integer :: ihour,iminute,isecond
      real(kind=8) :: alt
      real(kind=4) :: u_g2s,v_g2s,temperature,density,pressure
      integer :: day_of_year

      character(len=130) :: lllinebuffer
      integer :: nargs
      integer :: iy,im,id,ibar,inum

      real(kind=8),dimension(:),allocatable   :: z_prof

!     TEST READ COMMAND LINE ARGUMENTS
      nargs = iargc()

        ! For the probe_HWT branch of this program, we need (minimally)
        !   inyear  : year
        !   inmonth : month
        !   inday   : day
        !   inhour  : hour (as float)
        !   lon_in  : longitude of sonde point
        !   lat_in  : latitude of sonde point
        !   ap      : Ap value
        !   f107    : F107 value
        !   zmax    : maximum altitude of profile
        !   dz_prof : vertical increment
        ! Optionally, we can just feed it a command file with these 4 values
        ! along with all values for all these other variables
      if (nargs.eq.10)then
        ! parse the eight input parameters
        write(G2S_global_info,*)"Parsing argument list: "
        call getarg(1,lllinebuffer)
        read(lllinebuffer,*)inyear
        write(G2S_global_info,*)"inyear",inyear
        call getarg(2,lllinebuffer)
        read(lllinebuffer,*)inmonth
        write(G2S_global_info,*)"inmonth",inmonth
        call getarg(3,lllinebuffer)
        read(lllinebuffer,*)inday
        write(G2S_global_info,*)"inday",inday
        call getarg(4,lllinebuffer)
        read(lllinebuffer,*)inhour
        write(G2S_global_info,*)"inhour",inhour
        call getarg(5,lllinebuffer)
        read(lllinebuffer,*)lon_in
        write(G2S_global_info,*)"lon_in",lon_in
        call getarg(6,lllinebuffer)
        read(lllinebuffer,*)lat_in
        write(G2S_global_info,*)"lat_in",lat_in
        call getarg(7,lllinebuffer)
        read(lllinebuffer,*)ap
        write(G2S_global_info,*)"ap",ap
        call getarg(8,lllinebuffer)
        read(lllinebuffer,*)f107
        write(G2S_global_info,*)"f107",f107
        call getarg(9,lllinebuffer)
        read(lllinebuffer,*)zmax
        write(G2S_global_info,*)"zmax",zmax
        call getarg(10,lllinebuffer)
        read(lllinebuffer,*)dz_prof
        write(G2S_global_info,*)"dz_prof",dz_prof

      elseif (nargs.eq.1)then
        ! Assume single argument is a command file with the format
        !   inyear,inmonth,inday,inhour :  year, month, day, hour
        !   lon_in,lat_in : longitude, latitude of sonde point
        !   apfile        : Ap filename
        !   f107file      : F107 filename
        !   zmax, dz_prof : maximum altitude of output, vertical increment
        call getarg(1,lllinebuffer)
        read(lllinebuffer,*)controlfile
        open(unit=ct_unit,file=controlfile,status='old')
        read(ct_unit,*)inyear,inmonth,inday,inhour
        write(G2S_global_info,*)"inyear,inmonth,inday,inhour",inyear,inmonth,inday,inhour
        read(ct_unit,*)lon_in,lat_in
        write(G2S_global_info,*)"lon_in,lat_in: ",lon_in,lat_in
        read(ct_unit,'(a130)')apfile
        read(ct_unit,'(a130)')f107file
        read(ct_unit,*)zmax, dz_prof
        write(G2S_global_info,*)"zmax, dz_prof: ",zmax, dz_prof

        if(apfile.ne.f107file)then
          ! Open and read individual Ap and F107files
          open(unit=ap_unit,file=apfile,status='old')
          read(ap_unit,*)ap
          close(ap_unit)
          open(unit=f107_unit,file=f107file,status='old')
          read(f107_unit,*)f107
          close(f107_unit)
        else
          write(G2S_global_info,*)"Ap and F107 are equivalent, assume it is a directory to archive."
          ! Open archive file and get find the line for the requested day
          write(apfile,115)trim(adjustl(apfile)),'/',inyear
 115      format(a4,a1,i4)
          open(unit=ap_unit,file=apfile,status='old')
          read(ap_unit,101)iy,im,id,ibar,inum,&
                           kp1,kp2,kp3,kp4,kp5,kp6,kp7,kp8,kp9,&
                           ap1,ap2,ap3,ap4,ap5,ap6,ap7,ap8,ap9,&
                           Cp,C9,ssnum,f107,f107q
          do while(im.ne.inmonth.or.id.ne.inday)
            read(ap_unit,101)iy,im,id,ibar,inum,&
                             kp1,kp2,kp3,kp4,kp5,kp6,kp7,kp8,kp9,&
                             ap1,ap2,ap3,ap4,ap5,ap6,ap7,ap8,ap9,&
                             Cp,C9,ssnum,f107,f107q
          enddo
          if(inhour.le.3.0)then
            ap=ap1
          elseif(inhour.le.6.0)then
            ap=ap2
          elseif(inhour.le.9.0)then
            ap=ap3
          elseif(inhour.le.12.0)then
            ap=ap4
          elseif(inhour.le.15.0)then
            ap=ap5
          elseif(inhour.le.18.0)then
            ap=ap6
          elseif(inhour.le.21.0)then
            ap=ap7
          elseif(inhour.lt.24.0)then
            ap=ap8
          endif

 101      format(3i2,i4,9i2,10i3,f3.1,i1,i3,f5.1,i1)
        endif
        close(f107_unit)
        write(G2S_global_info,*)"Using an ap value of  : ",ap
        write(G2S_global_info,*)"Using a solar flux of : ",f107

      else
        ! Dump useage info to stdout:
        write(G2S_global_info,*)"No command-line arguments given"
        write(G2S_global_info,*)" "
        write(G2S_global_info,*)"probe_HWT can be run in two modes:"
        write(G2S_global_info,*)"  (1) command-line arguments"
        write(G2S_global_info,*)"    probe_HWT followed by 10 command-line arguments"
        write(G2S_global_info,*)"      year : integer"
        write(G2S_global_info,*)"      month: integer"
        write(G2S_global_info,*)"      day  : integer"
        write(G2S_global_info,*)"      hour : real"
        write(G2S_global_info,*)"      lon  : longitude of sonde point"
        write(G2S_global_info,*)"      lat  : latitude of sonde point"
        write(G2S_global_info,*)"      ap   : Ap planetary index"
        write(G2S_global_info,*)"      f107 : F107 value"
        write(G2S_global_info,*)"      zmax : maximum altitude of profile"
        write(G2S_global_info,*)"      dz   : vertical increment"
        write(G2S_global_info,*)"  (2) control file"
        write(G2S_global_info,*)"    probe_HWT input_HWTprobe.ctr"
        write(G2S_global_info,*)"      where input_HWTprobe.ctr has the following format:"
        write(G2S_global_info,*)"      "
        write(G2S_global_info,*)"        inyear,inmonth,inday,inhour  # time of sonde"
        write(G2S_global_info,*)"        lon_in,lat_in  # longitude, latitude of sonde point"
        write(G2S_global_info,*)"        ap_filename    # name of Ap file"
        write(G2S_global_info,*)"        f107_filename  # name of F107 file"
        write(G2S_global_info,*)"        zmax, dz_prof  # maximum altitude of output, vertical increment"
        write(G2S_global_info,*)" "
        write(G2S_global_info,*)" For example for Pavlof sonde"
        write(G2S_global_info,*)"    ./probe_HWT 2013 5 13 19.6 -161.887 55.42 4.0 153.5 200.0 2.5"

        stop 1
      endif

      start_year = inyear
      start_month= inmonth
      start_day  = inday
      day        = day_of_year(start_year,start_month,start_day)
      ihour      = floor(inhour)
      iminute    = floor((inhour-ihour)*60.0)
      isecond    = floor((inhour - ihour - iminute/60.0)*3600.0)

      nz_prof = floor(zmax/dz_prof)+1
      allocate(z_prof(nz_prof))
      do k = 1,nz_prof
        z_prof(k) = (k-1)*dz_prof
      enddo

      write(G2S_global_info,*)"Calculating sonde data from HWT at ",lon_in,lat_in
      open(55,file='sonde_ztuvrp.dat',status='replace')
      do k = 1,nz_prof
        alt = z_prof(k)
        call Get_WindTempRhoP_Empir(lon_in,lat_in,alt,              &
                    start_year,day,ihour,iminute,isecond, &
                    ap,f107,                              &
                    u_g2s,v_g2s,temperature,density,pressure)
        pressure = 1013.0_4 * exp(-real(alt,kind=4)/7.4_4) ! Eq 1.5 of Seinfeld/Pandis
        density  = (pressure*100.0_4) / (temperature*287.058_4)/1000.0_4
        !write(*,*)lon_in,lat_in,alt,start_year,day,ihour,iminute,isecond,ap,f107,u_g2s,v_g2s,temperature
        write(55,"(6E15.5)")alt,temperature,u_g2s,v_g2s,density,pressure

      enddo
      close(55)

      write(G2S_global_info,*)"probe_HWT exited normally."

      end program probe_HWT


!##############################################################################
!##############################################################################

!##############################################################################
!##############################################################################
!
!     Get_WindTempRho_Empir
!
!     This subroutine is the interface to the empirical models HWM07 and
!     NRLMSISE.  This takes as arguments, the full coordinates (x,y,z,t) as well
!     as the space-weather indices Ap and F107, then returns the Vx, Vy, and T
!     values.
!
!     No global variables are filled.  Vx,Vy, and T are returned through the
!     argument list.
!
!##############################################################################

      subroutine Get_WindTempRhoP_Empir(lon,lat,alt,          &
                              year,day,ihour,iminute,isecond, &
                              ap,f107,                        &
                              u_g2s,v_g2s,temperature,density,pressure)

      ! Here are the NRLMSISE-00 modules
      use utils_constants
      use physics_msis

#ifdef useHWM07
      use NEWmodel
#endif

      implicit none

      integer,intent(in)       :: day,year
      integer,intent(in)       :: ihour,iminute,isecond
      real(kind=8),intent(in)  :: lat,lon,alt
      real(kind=4),intent(in)  :: ap,f107
      real(kind=4),intent(out) :: u_g2s,v_g2s,temperature,density,pressure

      ! Variables for HWM07 use single-precision reals
      real(4) :: sec_sp            ! UT(SEC)
      real(4) :: alt_sp            ! ALTITUDE(KM)
      real(4) :: glat_sp,glon_sp   ! GEODETIC LATITUDE(DEG) LONGITUDE(DEG)
      real(4) :: stl_sp            ! Solar time local: Not used
      real(4) :: F107A_sp          ! Not used
      real(4) :: F107_sp           ! Not used
      real(4) :: ap_sp(2)          ! 1:Not used; 2: CURRENT 3HR ap INDEX
      real(4) :: w_sp(2)           ! 1: MERIDIONAL WIND (m/sec + Northward)
                                   ! 2: ZONAL WIND (m/sec + Eastward)
      !real(4) :: qw_sp(2)          ! 
      !real(4) :: dw_sp(2)          ! 1: MERIDIONAL DISTURBANCE WIND (m/sec +
      !Geo. Northward)
                                   ! 2: ZONAL DISTURBANCE WIND (m/sec + Geo.
                                   ! Eastward)
      real(4) :: ut_sp             ! UT(HOURS)
      real(4) :: apqt_sp(2)        ! 1:Not used; 2: CURRENT 3HR ap INDEX
      real(4) :: pershift          ! 
      integer :: IYD               ! YEAR AND DAY AS YYDDD
      !integer :: ialt,istl

      ! Variables for NRLMSISE-00 use double-precision reals
      real(dp) :: ut_dp            ! UT(SEC)
      real(dp) :: alt_dp           ! ALTITUDE(KM)
      real(dp) :: xlat_dp          ! GEODETIC LATITUDE(DEG)
      real(dp) :: xlon_dp          ! GEODETIC LONGITUDE(DEG)
      real(dp) :: xlst_dp          ! LOCAL APPARENT SOLAR TIME(HRS)
      real(dp) :: f107a_dp         ! 81 day AVERAGE OF F10.7 FLUX (centered onday DDD)
      real(dp) :: f107_dp          ! DAILY F10.7 FLUX FOR PREVIOUS DAY
      !real(dp) :: ap_dp            ! MAGNETIC INDEX(DAILY)

      real(dp) , dimension(7) :: aph_dp
      real(dp) , dimension(9) :: d_dp  ! D(1) - HE NUMBER DENSITY(CM-3)
                                       ! D(2) - O NUMBER DENSITY(CM-3)
                                       ! D(3) - N2 NUMBER DENSITY(CM-3)
                                       ! D(4) - O2 NUMBER DENSITY(CM-3)
                                       ! D(5) - AR NUMBER DENSITY(CM-3)
                                       ! D(6) - TOTAL MASS DENSITY(GM/CM3)
                                       ! [includes anomalous oxygen]
                                       ! D(7) - H NUMBER DENSITY(CM-3)
                                       ! D(8) - N NUMBER DENSITY(CM-3)
                                       ! D(9) - Anomalous oxygen NUMBER
                                       ! DENSITY(CM-3)
      real(dp) , dimension(2) :: t_dp  ! T(1) - EXOSPHERIC TEMPERATURE
                                       ! T(2) - TEMPERATURE AT ALT

      real(4) :: vx_new, vy_new

      ! get winds
      iyd = 1000*mod(year,100) + day
      ut_sp = real(ihour,kind=4) + real(iminute,kind=4)/60.0_4 + &
              real(isecond,kind=4)/3600.0_4
      sec_sp = ut_sp * 3600.0_4
      glat_sp = real(lat,kind=4)
      glon_sp = real(lon,kind=4)
      alt_sp = real(alt,kind=4)
      stl_sp = pershift(ut_sp + glon_sp/15.0_4, (/0.0_4, 24.0_4/) )
      ap_sp(1) = real(0.0,kind=4)
      ap_sp(2) = real(ap,kind=4)
      apqt_sp(1) = real(0.0,kind=4)
      apqt_sp(2) = real(-1.0,kind=4)

      f107_sp  = real(f107,kind=4)
      f107a_sp = f107_sp

#ifdef useHWM07
      call HWM07(iyd,sec_sp,alt_sp,glat_sp,glon_sp,stl_sp,f107a_sp,f107_sp,ap_sp,w_sp)
#endif
#ifdef useHWM14
      call hwm14(iyd,sec_sp,alt_sp,glat_sp,glon_sp,stl_sp,f107a_sp,f107_sp,ap_sp,w_sp)
#endif
      vx_new = w_sp(2)
      vy_new = w_sp(1)

      ! Get temperature
      ut_dp    = real(ihour,kind=8)*3600.0 + real(iminute,kind=8)*60.0 + &
                 real(isecond,kind=8)
      alt_dp   = real(alt,kind=8)
      xlat_dp  = real(lat,kind=8)
      xlon_dp  = real(lon,kind=8)
      xlst_dp  = real(stl_sp,kind=8)
      f107_dp  = real(f107,kind=8)
      f107a_dp = f107_dp
      aph_dp   = real(ap,kind=8)

      call gtd7(iyd,ut_dp,alt_dp,xlat_dp,xlon_dp,xlst_dp,f107a_dp,f107_dp,aph_dp,48, &
                d_dp,t_dp)
      u_g2s       = vx_new
      v_g2s       = vy_new
      temperature = real(t_dp(2),kind=4)
      density     = real(d_dp(6),kind=4)
      pressure    = real((d_dp(6)*1000.0_8) & ! convert to kg/m3
                    * t_dp(2)               & ! in K
                    * 287.05_8              & ! R = Specific gas const dry air J/kgK
                    * 0.01_8,kind=4)          ! convert Pa to mbar

      end subroutine Get_WindTempRhoP_Empir

!##############################################################################
!##############################################################################


      function day_of_year(iyear,imonth,iday)

      implicit none

      integer :: day_of_year
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

      if(imonth.eq.1)then
        day_of_year = iday
      else
        day_of_year = sum(monthdays(1:imonth-1)) + iday
      endif

      end function day_of_year

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

