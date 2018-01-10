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
!  g2s_genSC
!
!  This program reads one or more groups of NWP files and smoothly merges them
!  with empirical values taken from HWM14 and NRLMSIS-00.  The spectral (either
!  spherical harmonic or Fourier) coefficients of this smoothed atmospheric
!  conditions is written out to a netcdf file.
!
!##############################################################################

      program G2S_genSC

      use G2S_globvar
      use MetReader
      use projection

      implicit none

      integer :: day
      integer :: ihour,iminute,isecond
      real(kind=4) :: lat,lon,alt
      real(kind=4) :: u_g2s,v_g2s,temperature
      integer :: day_of_year

      integer :: i,j,k
      integer :: ios

      character    :: testchar
      character(len=130) :: lllinebuffer
      integer :: nargs
      integer :: iy,im,id
      integer :: ap1,ap2,ap3,ap4,ap5,ap6,ap7,ap8 !,ap9

      integer :: idf,igrid,iw
      integer :: ivalue3,ivalue4
      real(kind=8) :: HS_hours_since_baseyear
      real(kind=8) :: Met_needed_StartHour
      real(kind=8) :: Simtime_in_hours
      integer :: BaseYear = 1900
      logical :: useLeap  = .true.

      logical :: IsPeriodic
      integer :: ivar

      character(len=71)  :: linebuffer
      integer :: ioerr

      integer      :: ilatlonflag
      integer      :: iprojflag
      real(kind=8) :: k0_scale
      real(kind=8) :: phi0    != 90.0_ip        ! latitude of projection point
      real(kind=8) :: lam0 != -135.0_ip   ! longitude of projection point
      real(kind=8) :: lam1,phi1
      real(kind=8) :: lam2,phi2
      real(kind=8) :: radius_earth

      real(kind=8) :: HoursIntoInterval
      real(kind=8) :: Interval_Frac

      real(kind=8) :: ptlon,ptlat
      real(kind=4) :: rotang
      real(kind=8) :: x_in,x_out,de_x
      real(kind=8) :: y_in,y_out,de_y

      ! Make user MetReader is using the same calendar
      MR_BaseYear = BaseYear
      MR_useLeap  = useLeap

!     TEST READ COMMAND LINE ARGUMENTS
      nargs = iargc()
      if (nargs.ge.1) then
        call getarg(1,lllinebuffer)
        ! If the input argument is a number, assume it is the forecast hour and
        ! iwf1 = 20 (the GFS 0.5-degree forecast files), otherwise,
        ! assume the argument is an input control file
        read(lllinebuffer,'(a1)')testchar
        ! Assuming a standard ASCII character set
        if(ichar(testchar).ge.48.and.ichar(testchar).le.57)then
          write(G2S_global_info,*)"Testing first character of command-line argument..."
          write(G2S_global_info,*)"Input arg #1 is a number : ",testchar
          write(G2S_global_info,*)"Assume this is a forecast hour."
          write(G2S_global_info,*)"This section needs to be customized to your system: exiting"
          stop 1

         ! This section can be set up with system-specific paths so that only
         ! a forecast hour is needed for command-line arguments

          !read(lllinebuffer,*)fc_hour
          !iwf1 = 20
          !if(fc_hour.ne. 0.and.fc_hour.ne. 3.and.fc_hour.ne. 6.and.fc_hour.ne. 9.and.&
          !   fc_hour.ne.12.and.fc_hour.ne.15.and.fc_hour.ne.18.and.fc_hour.ne.21.and.&
          !   fc_hour.ne.24.and.fc_hour.ne.27.and.fc_hour.ne.30)then
          !  write(windfile,121)fc_hour
          !  ! Edit this format line to point to your windfile
          !  !  An example is below
! 121        format('/data/Windfiles/gfs/latest/latest.f',i2,'.nc')
          !else
          !  write(G2S_global_info,*)"fc must be 0,3,6,9,12,15,18,21,24,27, or 30."
          !  stop 1
          !endif
          !! Edit this format line to point to your Ap and F107 files
          !apfile   = "/opt/USGS/AVOG2S/ExternalData/Ap_Forecast/Ap.dat"
          !f107file = "/opt/USGS/AVOG2S/ExternalData/Ap_Forecast/F107.dat"
          !! Edit this format line to point to your coefficient files
          !file_SH = "/data/Ground2Space/Model_SH/G2S_SH.nc"

          !! Read the latest 3-hourly Ap megnetic index which are downloaded and
          !! extracted from http://www.swpc.noaa.gov/ftpdir/latest/wwv.txt
          !! www.swpc.noaa.gov/products/geophysical-alert-wwv-text
          !open(unit=ap_unit,file=apfile,status='old')
          !read(ap_unit,*)ap
          !close(ap_unit)
          !write(G2S_global_info,*)"Using an ap value of  : ",ap
      
          !! Read the latest 10.7cm solar flux index, also from www.txt
          !open(unit=f107_unit,file=f107file,status='old')
          !read(f107_unit,*)f107
          !close(f107_unit)
          !write(G2S_global_info,*)"Using a solar flux of : ",f107

        else
          ! Loading an input control file
          !  Control file format example
          !     2008 8 8 4.6    ! year, month, day, hour
          !     1 4 -107.0 50.0 50.0 50.0 6367.470    ! IsLatLon, iprojflag, ...
          !     0.5  120        ! G2S_ddeg Spec_MaxOrder
          !     20 200 2.5      ! zmin_HWT,zmax_HWT,dz_HWT
          !     1               ! number of groups of winddfiles
          !     20 1            ! iwf_Met1, nwindfiles
          !     0.0 25.0 1.0    ! zmin_Met1 zmax_Met1, dz_Met1
          !     /data/windfiles/iw5_iwf25_NCEP_Reanalysis1
          !     /data/Ground2Space/Ap_Forecast/NGDC_NOAA_Archive/2008
          !     /data/Ground2Space/Ap_Forecast/NGDC_NOAA_Archive/2008
          !     out_GS.nc

          read(lllinebuffer,*)controlfile

          write(G2S_global_info,*)"Input arg #1 is NOT a number : ",controlfile
          write(G2S_global_info,*)"Assume this is a control file."

          open(unit=ct_unit,file=controlfile,status='old')
          read(ct_unit,*)inyear,inmonth,inday,inhour
          !READ PROJECTION PARAMETERS
          read(ct_unit,'(a80)') Comp_projection_line
          read(Comp_projection_line,*)ilatlonflag,iprojflag
          if (ilatlonflag.eq.0) then
            ! expecting input variables to be in the same projection as
            ! specified by iprojflag and parameters
            IsLatLon          = .false.
          else
            ! expecting input variables to be in lat/lon
           IsLatLon          = .true.
          endif
          !SET PROJECTION PARAMETERS
          if (IsLatLon.eqv..false.)then
            call PJ_Set_Proj_Params(Comp_projection_line)
            ilatlonflag  = PJ_ilatlonflag
            iprojflag    = PJ_iprojflag
            k0_scale     = PJ_k0
            radius_earth = PJ_radius_earth
            lam0         = PJ_lam0
            lam1         = PJ_lam1
            lam2         = PJ_lam2
            phi0         = PJ_phi0
            phi1         = PJ_phi1
            phi2         = PJ_phi2
          endif

          if (IsLatLon.eqv..true.)then
              ! We are reading a Lon/Lat file, use global Spherical Harmonics
              ! and assume lon goes from 0-360 with nxmax_g2s set by dx_g2s
              !            lat goes from -90-90 with nymax_g2s set by dx_g2s
            read(ct_unit,*)dx_g2s,Spec_MaxOrder
            xmin_g2s =   0.0
            ymin_g2s = -90.0
            dy_g2s   = dx_g2s
            nxmax_g2s = nint(360.0/dx_g2s)
            nymax_g2s = nint(180.0/dy_g2s) + 1
          else
              ! This is a projected grid, use Fourier decomposition
              ! We need the extra info to specify the computational grid
            read(ct_unit,*)dx_g2s,Spec_MaxOrder,xmin_g2s,nxmax_g2s,ymin_g2s,nymax_g2s
            dy_g2s = dx_g2s
          endif

          read(ct_unit,*)zmin_HWT, zmax_HWT, dz_HWT
          read(ct_unit,*)nwindfile_groups
          idf = 2
          if(nwindfile_groups.eq.1.or.nwindfile_groups.eq.2)then
            read(ct_unit,'(a71)')linebuffer
              ! Assume we can at least read the windfile format number
            read(linebuffer,*)iwf1
              ! Try to read the number of windfiles as well
            read(linebuffer,*,iostat=ioerr)iwf1,nwindfiles1
            if (ioerr.ne.0)then
              ! number of windfiles not present, assume 1
              nwindfiles1 = 1
            else
              ! Succeeded in reading the two required values, try for three
              read(linebuffer,*,iostat=ioerr) iwf1,nwindfiles1, ivalue3
              if (ioerr.eq.0)then
                ! Success reading three values, try for four
                igrid = ivalue3
                read(linebuffer,*,iostat=ioerr) iwf1,nwindfiles1, igrid, ivalue4
                if (ioerr.eq.0)then
                  ! Success!, set data format (ascii, netcdf, hdf, grib)
                  idf = ivalue4
                endif
              else
                igrid = 0
              endif
            endif
            if(iwf1.eq.0)then
              ! If iwindformat = 0, then the input file is a not a known format
              ! Read an extra line given the name of a template file.
              if(idf.ne.2)then
                write(G2S_global_info,*)" Currently only netcdf reader implemented, resetting idf to 2"
                idf = 2
              endif
              read(10,'(a80)')linebuffer
              read(linebuffer,'(a80)') MR_iwf_template
              call MR_Read_Met_Template
            endif
            read(ct_unit,*) zmin_Met1, zmax_Met1, dz_Met1
            allocate(windfile1(nwindfiles1))
            do iw = 1,nwindfiles1
              read(ct_unit,'(a130)')windfile1(iw)
            enddo
          else
            write(G2S_global_error,*)"ERROR: number of windfile groups in neither 1 nor 2"
            write(G2S_global_error,*)"       nwindfile_groups = ",nwindfile_groups
            stop 1
          endif

          if(nwindfile_groups.eq.2)then
            read(ct_unit,'(a71)')linebuffer
              ! Assume we can at least read the windfile format number
            read(linebuffer,*)iwf2
              ! Try to read the number of windfiles as well
            read(linebuffer,*,iostat=ioerr)iwf2,nwindfiles2
            if (ioerr.ne.0)then
              ! number of windfiles not present, assume 1
              nwindfiles2 = 1
            endif
            if(iwf2.eq.0)then
              ! If iwindformat = 0, then the input file is a not a known format
              ! Read an extra line given the name of a template file.
              if(idf.ne.2)then
                write(G2S_global_info,*)" Currently only netcdf reader implemented, resetting idf to 2"
                idf = 2
              endif
              read(10,'(a80)')linebuffer
              read(linebuffer,'(a80)') MR_iwf_template
              call MR_Read_Met_Template
            endif
            read(ct_unit,*) zmin_Met2, zmax_Met2, dz_Met2
            allocate(windfile2(nwindfiles2))
            do iw = 1,nwindfiles2
              read(ct_unit,'(a130)')windfile2(iw)
            enddo
          endif
          read(ct_unit,'(a130)')apfile
          read(ct_unit,'(a130)')f107file
          read(ct_unit,'(a130)')file_SH
          write(G2S_global_info,*)"Year  = ",inyear
          write(G2S_global_info,*)"Month = ",inmonth
          write(G2S_global_info,*)"Day   = ",inday
          write(G2S_global_info,*)"Hour  = ",inhour
          write(G2S_global_info,*)"dx_g2s",dx_g2s
          write(G2S_global_info,*)"Spec_MaxOrder",Spec_MaxOrder
          write(G2S_global_info,*)"xmin_g2s",xmin_g2s
          write(G2S_global_info,*)"nxmax_g2s",nxmax_g2s
          write(G2S_global_info,*)"ymin_g2s",ymin_g2s
          write(G2S_global_info,*)"nymax_g2s",nymax_g2s
          write(G2S_global_info,*)"zmin_HWT  = ",zmin_HWT
          write(G2S_global_info,*)"zmax_HWT  = ",zmax_HWT
          write(G2S_global_info,*)"dz_HWT    = ",dz_HWT
          write(G2S_global_info,*)"nwindgroups= ",nwindfile_groups
          write(G2S_global_info,*)"windformat1= ",iwf1
          write(G2S_global_info,*)"nwindfiles1= ",nwindfiles1
          write(G2S_global_info,*)"zmin_Met1  = ", zmin_Met1
          write(G2S_global_info,*)"zmax_Met1  = ", zmax_Met1
          write(G2S_global_info,*)"dz_Met1    = ", dz_Met1
          write(G2S_global_info,*)"Windfile   = ",trim(adjustl(windfile1(1)))
          if(nwindfile_groups.eq.2)then
            write(G2S_global_info,*)"windformat2= ",iwf2
            write(G2S_global_info,*)"nwindfiles2= ",nwindfiles2
            write(G2S_global_info,*)"zmin_Met1  = ", zmin_Met2
            write(G2S_global_info,*)"zmax_Met1  = ", zmax_Met2
            write(G2S_global_info,*)"dz_Met1    = ", dz_Met2
            write(G2S_global_info,*)'-',trim(adjustl(windfile2(1))),'-'
          endif
          write(G2S_global_info,*)"Ap file   = ",trim(adjustl(apfile))
          write(G2S_global_info,*)"F107 file = ",trim(adjustl(f107file))
          write(G2S_global_info,*)"outfile   = ",trim(adjustl(file_SH))
          !write(G2S_global_info,*)apfile.eq.f107file
          close(ct_unit)

          if(apfile.ne.f107file)then
            ! Open and read individual Ap and F107files
            open(unit=ap_unit,file=apfile,status='old')
            read(ap_unit,*)ap
            close(ap_unit)
            open(unit=f107_unit,file=f107file,status='old')
            read(f107_unit,*)f107
            close(f107_unit)
          else
            write(G2S_global_info,*)"Ap and F107 are equivalant, assume it is a directory to archive."
            ! Open archive file and get find the line for the requested day
            write(apfile,115)trim(adjustl(apfile)),'/',inyear
 115        format(a4,a1,i4)
            open(unit=ap_unit,file=apfile,status='old')
            read(ap_unit,'(a71)')linebuffer
            read(linebuffer,102)iy,im,id, &
                                ap1,ap2,ap3,ap4,ap5,ap6,ap7,ap8, &
                                f107
 102        format(3i2,25x,8i3,10x,f5.1,1x)

            do while(im.ne.inmonth.or.id.ne.inday)
              read(ap_unit,'(a71)')linebuffer
              read(linebuffer,102)iy,im,id, &
                                  ap1,ap2,ap3,ap4,ap5,ap6,ap7,ap8, &
                                  f107
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

! 101        format(3i2,i4,9i2,10i3,f3.1,i1,i3,f5.1,i1)
          endif
          close(f107_unit)
          write(G2S_global_info,*)"Using an ap value of  : ",ap
          write(G2S_global_info,*)"Using a solar flux of : ",f107

        endif
      else
        ! No command-line arguments given, print usage and exit

        write(G2S_global_info,*)"No command-line arguments given"
        write(G2S_global_info,*)"Either run g2s with a forecast hour or a control file."
        write(G2S_global_info,*)" "
        write(G2S_global_info,*)" The control file is of the following format:"
        write(G2S_global_info,*)"----------------------------------------------------------------------"
        write(G2S_global_info,*)"2016 10 5 4.6        ! year, month, day, hour"
        write(G2S_global_info,*)"1 4 -107.0 50.0 50.0 50.0 6367.470    ! IsLatLon, iprojflag, ..."
        write(G2S_global_info,*)"0.5 120              ! ddeg of G2S grid, maxdeg of SH decomp"
        write(G2S_global_info,*)"20 200 2.5           ! zmin_HWT,zmax_HWT,dz_HWT"
        write(G2S_global_info,*)"2                    ! number of windfile groups"
        write(G2S_global_info,*)"20 1                 ! format ID of met data1 (20 for new GFS) ,num of windfiles1"
        write(G2S_global_info,*)"0.0 25.0 1.0         ! zmin, zmax, dz of Met1"
        write(G2S_global_info,*)"gfs.t00z.pgrb2f00.nc ! windfile name"
        write(G2S_global_info,*)"41 1                 ! format ID of met data2  ,num of windfiles2"
        write(G2S_global_info,*)"20.0 55.0 2.0        ! zmin, zmax, dz of Met2"
        write(G2S_global_info,*)"GEOS.fp.fcst.inst3_3d_asm_Np.20161005_00+20161005_0000.V01.nc4 ! windfile name"
        write(G2S_global_info,*)"Ap.dat               ! Filename with Ap value"
        write(G2S_global_info,*)"F107.dat             ! Filename with F107 value"
        write(G2S_global_info,*)"out_GS.nc            ! Output filename of SH coefficients"
        write(G2S_global_info,*)"                              "
        write(G2S_global_info,*)"  Note: If the G2S grid is projected, line 3 of the control file should be"
        write(G2S_global_info,*)"5.950 120 -2172.922 512 -4214.803 340  ! dx,maxdeg,xstart,nx,ystart,ny"
        write(G2S_global_info,*)"----------------------------------------------------------------------"

        stop 1

      endif
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Finished parsing command line / command file parameters
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! First working on Met1
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      write(G2S_global_info,*)"------------------------------"
      write(G2S_global_info,*)"  Getting data from met file 1"
      write(G2S_global_info,*)"------------------------------"
      write(G2S_global_info,*)"   **First set up grid.**"

      start_year = inyear
      start_month= inmonth
      start_day  = inday
      day     = day_of_year(start_year,start_month,start_day)

      ! Use MetReader library
      if(iwf1.eq.20.or.iwf1.eq.21.or.iwf1.eq.26)then
        iw      = 4   ! read from multiple files
        igrid   = 4   ! 0.5 degree global
        !idf     = 2   ! netcdf
      elseif(iwf1.eq.22)then
        iw      = 4   ! read from single file
        igrid   = 193 ! 0.25 degree global
        !idf     = 2   ! netcdf
      elseif(iwf1.eq.24)then
        iw      = 3   ! read from single file
        igrid   = 1024 ! MERRA
        !idf     = 2   ! netcdf
      elseif(iwf1.eq.25)then
        iw      = 5   ! read from multiple files
        igrid   = 2   ! 0.5 degree global
        !idf     = 2   ! netcdf
      elseif(iwf1.eq.28)then
        iw      = 4   ! read from multiple files
        igrid   = 170   ! 0.5 degree global
        idf     = 2   ! netcdf
      elseif(iwf1.eq.40)then
        iw      = 4   ! read from single file
        igrid   = 1040 ! NASA GMAO Cp
        !idf     = 2   ! netcdf
      elseif(iwf1.eq.41)then
        iw      = 4   ! read from single file
        igrid   = 1041 ! NASA GMAO Np
        !idf     = 2   ! netcdf
      elseif(iwf1.eq.3.or.iwf1.eq.4)then
        iw      = 3   ! read from single file
        igrid   = 221 ! NARR 32-km
        idf     = 2   ! netcdf
      elseif(iwf1.eq.12)then
        iw      = 4   ! read from single file
        igrid   = 198 ! nam 5.9km AK
        !idf     = 2   ! netcdf
      elseif(iwf1.eq.13)then
        iw      = 4   ! read from single file
        igrid   = 91 ! nam 2.95km AK
        !idf     = 2   ! netcdf
      else
        write(G2S_global_error,*)"ERROR: Only iwf2 = 3,4,12,13,20,21,22,24,25,26,28,40,41 implemented."
        stop 1
      endif
      Met_needed_StartHour = HS_hours_since_baseyear(inyear,inmonth,inday,inhour,MR_BaseYear,MR_useLeap)
      Simtime_in_hours     = 0.0
      call MR_Allocate_FullMetFileList(iw,iwf1,igrid,idf,nwindfiles1)
      do iw = 1,nwindfiles1
        MR_windfiles(iw) = trim(ADJUSTL(windfile1(iw)))
      enddo
        ! Now check for existance and compatibility with simulation time requirements
      call MR_Read_Met_DimVars(inyear)

        ! Note: the Met files will always be on the g2s grid, possibly resampled
        ! from the native resolution within MetReader

      if (IsLatLon.eqv..true.)then
        IsPeriodic = .true.
      else
        IsPeriodic = .false.
      endif

      allocate(x_g2s_sp(nxmax_g2s))
      allocate(y_g2s_sp(nymax_g2s))
      do i = 1,nxmax_g2s
        x_g2s_sp(i) = real(xmin_g2s + (i-1)*dx_g2s,kind=4)
      enddo
      do i = 1,nymax_g2s
        y_g2s_sp(i) = real(ymin_g2s + (i-1)*dy_g2s,kind=4)
      enddo

      nzmax_Met1 = nint((zmax_Met1-zmin_Met1)/dz_Met1) + 1
      allocate(z_Met1_sp(nzmax_Met1))
        ! Fill z array
      do i = 1,nzmax_Met1
        z_Met1_sp(i) = real(zmin_Met1 + (i-1)*dz_Met1,kind=4)
      enddo

      call MR_Set_CompProjection(IsLatLon,iprojflag,lam0,phi0,phi1,phi2,k0_scale,radius_earth)

      call MR_Initialize_Met_Grids(nxmax_g2s,nymax_g2s,nzmax_Met1, &
                                   x_g2s_sp(1:nxmax_g2s),        &
                                   y_g2s_sp(1:nymax_g2s),        &
                                   z_Met1_sp(1:nzmax_Met1),        &
                                   IsPeriodic)
      call MR_Set_Met_Times(Met_needed_StartHour, Simtime_in_hours)

      if (IsLatLon_CompGrid.eqv..false.)then
        ! For projected grids, we will need the proj to LL mapping in order to
        ! use the HWT routines.  Also, we'll need to rotate the returned wind
        ! vectors back to the projected grid
        allocate(Lon_of_proj_node_dp(nxmax_g2s,nymax_g2s))
        allocate(Lat_of_proj_node_dp(nxmax_g2s,nymax_g2s))
        allocate(Theta_of_proj_node_dp(nxmax_g2s,nymax_g2s))
        do i=1,nxmax_g2s
          x_in = x_g2s_sp(i)
          do j=1,nymax_g2s
            y_in = y_g2s_sp(j)
            call PJ_proj_inv(x_in, y_in, iprojflag, &
                           lam0,phi0,phi1,phi2,k0_scale,radius_earth, &
                           x_out,y_out)
            Lon_of_proj_node_dp(i,j) = x_out
            Lat_of_proj_node_dp(i,j) = y_out
          enddo
        enddo
        ! Now find the rotation angle
        do i=1,nxmax_g2s
          x_in = x_g2s_sp(i)
          do j=1,nymax_g2s
            y_in = y_g2s_sp(j)
              ! Get lon/lat of point in question
            ptlon = Lon_of_proj_node_dp(i,j)
            ptlat = Lat_of_proj_node_dp(i,j)
              ! Get projected coordinate of de at the currrent point
            call PJ_proj_for(ptlon+1.0_8/60.0_8,ptlat, iprojflag, &
                       lam0,phi0,phi1,phi2,k0_scale,radius_earth, &
                       x_out,y_out)
            de_x = x_out-x_in
            de_y = y_out-y_in
              ! Now recover the angle between de and x (of map grid)
            Theta_of_proj_node_dp(i,j) = atan(de_y/de_x)
          enddo
        enddo
      endif

      ! These variables hold the Met values interpolated onto the g2s grid
      allocate(vx_Met_loc_sp(nxmax_g2s,nymax_g2s,nzmax_Met1))
      allocate(vy_Met_loc_sp(nxmax_g2s,nymax_g2s,nzmax_Met1))
      allocate(temperature_Met_loc_sp(nxmax_g2s,nymax_g2s,nzmax_Met1))

      Spec_MaxOrder = min(Spec_MaxOrder,nymax_g2s/2-1)
      ihour     = start_hour
      iminute   = start_min
      isecond   = start_sec

      !Load the 3d arrays for uwind, vwind, temperature, geopotential
      write(G2S_global_info,*)"   **Now read in the state variables from met file.**"
      write(G2S_global_info,*)"     Find MR_iMetStep_Now for ",Met_needed_StartHour
      MR_iMetStep_Now = 1
      do i = 1,MR_MetSteps_Total-1
        write(G2S_global_info,*)Met_needed_StartHour,MR_MetStep_Hour_since_baseyear(i:i+1)
        if(Met_needed_StartHour.ge.MR_MetStep_Hour_since_baseyear(i).and.&
           Met_needed_StartHour.lt.MR_MetStep_Hour_since_baseyear(i+1))then
          MR_iMetStep_Now = i
          write(G2S_global_info,*)i,MR_MetStep_Hour_since_baseyear(i)
          !write(G2S_global_info,*)MR_iMetStep_Now,TimeNow_fromRefTime,&
          !    MR_MetStep_Hour_since_baseyear(MR_iMetStep_Now)
          cycle
        endif
      enddo

      HoursIntoInterval = Met_needed_StartHour - MR_MetStep_Hour_since_baseyear(MR_iMetStep_Now)
      MR_ForecastInterval = MR_MetStep_Interval(MR_iMetStep_Now)
      Interval_Frac = HoursIntoInterval / MR_ForecastInterval

      ! This subroutine sets both last and next geoH arrays so call with
      MR_iMetStep_Now = 1
      call MR_Read_HGT_arrays(MR_iMetStep_Now)

      allocate(tmp3d_1_sp(nxmax_g2s,nymax_g2s,nzmax_Met1))
      allocate(tmp3d_2_sp(nxmax_g2s,nymax_g2s,nzmax_Met1))

      if(Map_case.eq.1.or.Map_case.eq.2)then
        ! Either both the comp and met1 grids are LL (Map_case = 1)
        ! or they are both the same projection (Map_case = 2) so
        ! we can read the velocity components individually and interpolate onto
        ! the computational (g2s) grid
        ivar = 2 ! Vx
        call MR_Read_3d_Met_Variable_to_CompGrid(ivar,MR_iMetStep_Now)
        tmp3d_1_sp(1:nxmax_g2s,:,:) = MR_dum3d_compH(1:nxmax_g2s,:,:)
        call MR_Read_3d_Met_Variable_to_CompGrid(ivar,MR_iMetStep_Now+1)
        tmp3d_2_sp(1:nxmax_g2s,:,:) = MR_dum3d_compH(1:nxmax_g2s,:,:)
        vx_Met_loc_sp(1:nxmax_g2s,1:nymax_g2s,1:nzmax_Met1) = &
                tmp3d_1_sp(1:nxmax_g2s,1:nymax_g2s,1:nzmax_Met1) + &
               (tmp3d_2_sp(1:nxmax_g2s,1:nymax_g2s,1:nzmax_Met1) - &
                 tmp3d_1_sp(1:nxmax_g2s,1:nymax_g2s,1:nzmax_Met1)) * &
                                         real(Interval_Frac,kind=4)
  
        ivar = 3 ! Vy
        call MR_Read_3d_Met_Variable_to_CompGrid(ivar,MR_iMetStep_Now)
        tmp3d_1_sp(1:nxmax_g2s,:,:) = MR_dum3d_compH(1:nxmax_g2s,:,:)
        call MR_Read_3d_Met_Variable_to_CompGrid(ivar,MR_iMetStep_Now+1)
        tmp3d_2_sp(1:nxmax_g2s,:,:) = MR_dum3d_compH(1:nxmax_g2s,:,:)
        vy_Met_loc_sp(1:nxmax_g2s,1:nymax_g2s,1:nzmax_Met1) = &
                tmp3d_1_sp(1:nxmax_g2s,1:nymax_g2s,1:nzmax_Met1) + &
               (tmp3d_2_sp(1:nxmax_g2s,1:nymax_g2s,1:nzmax_Met1) - &
                 tmp3d_1_sp(1:nxmax_g2s,1:nymax_g2s,1:nzmax_Met1)) * &
                                         real(Interval_Frac,kind=4)
      else
        allocate(tmp3d_1_2_sp(nxmax_g2s,nymax_g2s,nzmax_Met1))
        allocate(tmp3d_2_2_sp(nxmax_g2s,nymax_g2s,nzmax_Met1))

        ! If this is not a global case (Met_Case=1) then we can assume the comp
        ! grid is projected.  If the met grid is also projected, and the same
        ! projection, then the above branch also applied.  So either we have met
        ! data that is projected differently that comp data (Map_case=5) or we
        ! have LL Met data and projected comp grid.
        if(Map_Case.eq.5)then
          ! In this case, we first need to rotate the projected Met data to
          ! Earth-relative on the MetP grid
          call MR_Rotate_UV_GR2ER_Met(MR_iMetStep_Now)
        endif
        call MR_Rotate_UV_ER2GR_Comp(MR_iMetStep_Now)
        ! U for comp grid is stored in MR_dum3d_compH
        tmp3d_1_sp(1:nxmax_g2s,:,:) = MR_dum3d_compH(1:nxmax_g2s,:,:)
        ! V for comp grid is stored in MR_dum3d_compH_2
        tmp3d_1_2_sp(1:nxmax_g2s,:,:) = MR_dum3d_compH_2(1:nxmax_g2s,:,:)

        if(Map_Case.eq.5)then
          call MR_Rotate_UV_GR2ER_Met(MR_iMetStep_Now+1)
        endif
        call MR_Rotate_UV_ER2GR_Comp(MR_iMetStep_Now+1)
        ! U for comp grid is stored in MR_dum3d_compH
        tmp3d_2_sp(1:nxmax_g2s,:,:) = MR_dum3d_compH(1:nxmax_g2s,:,:)
        ! V for comp grid is stored in MR_dum3d_compH_2
        tmp3d_2_2_sp(1:nxmax_g2s,:,:) = MR_dum3d_compH_2(1:nxmax_g2s,:,:)

        vx_Met_loc_sp(1:nxmax_g2s,1:nymax_g2s,1:nzmax_Met1) = &
                tmp3d_1_sp(1:nxmax_g2s,1:nymax_g2s,1:nzmax_Met1) + &
               (tmp3d_2_sp(1:nxmax_g2s,1:nymax_g2s,1:nzmax_Met1) - &
                 tmp3d_1_sp(1:nxmax_g2s,1:nymax_g2s,1:nzmax_Met1)) * &
                                         real(Interval_Frac,kind=4)

        vy_Met_loc_sp(1:nxmax_g2s,1:nymax_g2s,1:nzmax_Met1) = &
                tmp3d_1_2_sp(1:nxmax_g2s,1:nymax_g2s,1:nzmax_Met1) + &
               (tmp3d_2_2_sp(1:nxmax_g2s,1:nymax_g2s,1:nzmax_Met1) - &
                 tmp3d_1_2_sp(1:nxmax_g2s,1:nymax_g2s,1:nzmax_Met1)) * &
                                         real(Interval_Frac,kind=4)
        deallocate(tmp3d_1_2_sp)
        deallocate(tmp3d_2_2_sp)

      endif
      ivar = 5 ! Temperature
      call MR_Read_3d_Met_Variable_to_CompGrid(ivar,MR_iMetStep_Now)
      tmp3d_1_sp(1:nxmax_g2s,:,:) = MR_dum3d_compH(1:nxmax_g2s,:,:)
      call MR_Read_3d_Met_Variable_to_CompGrid(ivar,MR_iMetStep_Now+1)
      tmp3d_2_sp(1:nxmax_g2s,:,:) = MR_dum3d_compH(1:nxmax_g2s,:,:)
      temperature_Met_loc_sp(1:nxmax_g2s,1:nymax_g2s,1:nzmax_Met1) = &
              tmp3d_1_sp(1:nxmax_g2s,1:nymax_g2s,1:nzmax_Met1) + &
             (tmp3d_2_sp(1:nxmax_g2s,1:nymax_g2s,1:nzmax_Met1) - &
               tmp3d_1_sp(1:nxmax_g2s,1:nymax_g2s,1:nzmax_Met1)) * &
                                         real(Interval_Frac,kind=4)
      deallocate(tmp3d_1_sp,tmp3d_2_sp)


      if(write_test_output)then
        ! Dump out map at bottom of Met grid
        open(unit=30,file='plan_Met.dat',status='replace')
        do i=1,nxmax_g2s
          do j=1,nymax_g2s
            write(30,*)x_g2s_sp(i),y_g2s_sp(j),temperature_Met_loc_sp(i,j,1),&
                      vx_Met_loc_sp(i,j,1),vy_Met_loc_sp(i,j,1)
          enddo
        enddo
        close(30)
        ! Dump out equatorial cross-section of Met grid
        open(unit=40,file='xsec_Met.dat',status='replace')
        do i=1,nxmax_g2s
          do k=1,nzmax_Met1
            write(40,*)x_g2s_sp(i),y_g2s_sp(jout),z_Met1_sp(k),&
                       temperature_Met_loc_sp(i,jout,k),&
                       vx_Met_loc_sp(i,jout,k),vy_Met_loc_sp(i,jout,k)
          enddo
        enddo
        close(40)
        ! Dump out meridonal cross-section of Met grid
        open(unit=40,file='ysec_Met.dat',status='replace')
        do j=1,nymax_g2s
          do k=1,nzmax_Met1
            write(40,*)x_g2s_sp(iout),y_g2s_sp(j),z_Met1_sp(k),& 
                       temperature_Met_loc_sp(iout,j,k),&
                       vx_Met_loc_sp(iout,j,k),vy_Met_loc_sp(iout,j,k)
          enddo
        enddo
        close(40)
      endif

      if (IsLatLon_CompGrid.eqv..true.)then
        write(G2S_global_info,*)"   **Calculating spherical harmonics of met data.**"
        call Get_Met_SH(1)
      else
        write(G2S_global_info,*)"   **Calculating Fourier decomposition of met data.**"
        call Get_Met_FS(1)
      endif

      deallocate(vx_Met_loc_sp)
      deallocate(vy_Met_loc_sp)
      deallocate(temperature_Met_loc_sp)

      if(nwindfile_groups.eq.2)then
        call MR_Reset_Memory()
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Now working on Met2
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        write(G2S_global_info,*)"------------------------------"
        write(G2S_global_info,*)"  Getting data from met file 2"
        write(G2S_global_info,*)"------------------------------"
        write(G2S_global_info,*)"   **First set up grid.**"
  
        start_year = inyear
        start_month= inmonth
        start_day  = inday
        day     = day_of_year(start_year,start_month,start_day)
  
        ! Use MetReader library
        if(iwf2.eq.20.or.iwf2.eq.21.or.iwf2.eq.26)then
          iw      = 4   ! read from multiple files
          igrid   = 4   ! 0.5 degree global
          idf     = 2   ! netcdf
        elseif(iwf2.eq.22)then
          iw      = 4   ! read from single file
          igrid   = 193 ! 0.25 degree global
          idf     = 2   ! netcdf
        elseif(iwf2.eq.24)then
          iw      = 3   ! read from single file
          igrid   = 1024 ! MERRA
          idf     = 2   ! netcdf
        elseif(iwf2.eq.25)then
          iw      = 5   ! read from multiple files
          igrid   = 2   ! 0.5 degree global
          idf     = 2   ! netcdf
        elseif(iwf2.eq.28)then
          iw      = 4   ! read from multiple files
          igrid   = 170   ! 0.5 degree global
          idf     = 2   ! netcdf
        elseif(iwf2.eq.40)then
          iw      = 4   ! read from single file
          igrid   = 1040 ! MERRA
          idf     = 2   ! netcdf
        elseif(iwf2.eq.41)then
          iw      = 4   ! read from single file
          igrid   = 1041 ! MERRA
          idf     = 2   ! netcdf
        elseif(iwf2.eq.3)then
          iw      = 3   ! read from single file
          igrid   = 221 ! NARR 32-km
          idf     = 2   ! netcdf
        elseif(iwf2.eq.12)then
          iw      = 4   ! read from single file
          igrid   = 198 ! nam 5.9km AK
          idf     = 2   ! netcdf
        elseif(iwf2.eq.13)then
          iw      = 4   ! read from single file
          igrid   = 91 ! nam 5.9km AK
          idf     = 2   ! netcdf
        else
          write(G2S_global_error,*)"ERROR: Only iwf2 = 3,12,13,20,21,22,24,25,26,28,40,41 implemented."
          stop 1
        endif
        Met_needed_StartHour = HS_hours_since_baseyear(inyear,inmonth,inday,inhour,MR_BaseYear,MR_useLeap)
        Simtime_in_hours     = 0.0
        call MR_Allocate_FullMetFileList(iw,iwf2,igrid,idf,nwindfiles2)
        do iw = 1,nwindfiles2
          MR_windfiles(iw) = trim(ADJUSTL(windfile2(iw)))
        enddo
          ! Now check for existance and compatibility with simulation time requirements
        call MR_Read_Met_DimVars(inyear)

        nzmax_Met2 = nint((zmax_Met2-zmin_Met2)/dz_Met2) + 1
          ! Note: the Met files will always be on the g2s grid, possibly resampled
          ! from the native resolution within MetReader
        allocate(z_Met2_sp(nzmax_Met2))
  
          ! Fill z array
        do i = 1,nzmax_Met2
          z_Met2_sp(i) = real(zmin_Met2 + (i-1)*dz_Met2,kind=4)
        enddo
  
        call MR_Set_CompProjection(IsLatLon,iprojflag,lam0,phi0,phi1,phi2,k0_scale,radius_earth)

        call MR_Initialize_Met_Grids(nxmax_g2s,nymax_g2s,nzmax_Met2,&
                                x_g2s_sp(1:nxmax_g2s), &
                                y_g2s_sp(1:nymax_g2s), &
                                z_Met2_sp(1:nzmax_Met2)    , &
                                IsPeriodic)
        call MR_Set_Met_Times(Met_needed_StartHour, Simtime_in_hours)
 
        ! These variables hold the Met values interpolated onto the g2s grid
        allocate(vx_Met_loc_sp(nxmax_g2s,nymax_g2s,nzmax_Met2))
        allocate(vy_Met_loc_sp(nxmax_g2s,nymax_g2s,nzmax_Met2))
        allocate(temperature_Met_loc_sp(nxmax_g2s,nymax_g2s,nzmax_Met2))
  
        ihour     = start_hour
        iminute   = start_min
        isecond   = start_sec
  
        !Load the 3d arrays for uwind, vwind, temperature, geopotential
        write(G2S_global_info,*)"   **Now read in the state variables from met file.**"
        MR_iMetStep_Now = 1
        do i = 1,MR_MetSteps_Total-1
          if(Met_needed_StartHour.ge.MR_MetStep_Hour_since_baseyear(i).and.&
             Met_needed_StartHour.lt.MR_MetStep_Hour_since_baseyear(i+1))then
            MR_iMetStep_Now = i
            !write(G2S_global_info,*)MR_iMetStep_Now,TimeNow_fromRefTime,MetStep_Hour_since_1900(MR_iMetStep_Now)
            cycle
          endif
        enddo
        HoursIntoInterval = Met_needed_StartHour - MR_MetStep_Hour_since_baseyear(MR_iMetStep_Now)
        MR_ForecastInterval = MR_MetStep_Interval(MR_iMetStep_Now)
        Interval_Frac = HoursIntoInterval / MR_ForecastInterval
  
        ! This subroutine sets both last and next geoH arrays so call with MR_iMetStep_Now
        MR_iMetStep_Now = 1
        call MR_Read_HGT_arrays(MR_iMetStep_Now,.true.)
 
        allocate(tmp3d_1_sp(nxmax_g2s,nymax_g2s,nzmax_Met2))
        allocate(tmp3d_2_sp(nxmax_g2s,nymax_g2s,nzmax_Met2))
 
        if(Map_case.eq.1.or.Map_case.eq.2)then
          ! Either both the comp and met2 grids are LL (Map_case = 1)
          ! or they are both the same projection (Map_case = 2) so
          ! we can read the velocity components individually and interpolate onto
          ! the computational (g2s) grid

          ivar = 2 ! Vx
          call MR_Read_3d_Met_Variable_to_CompGrid(ivar,MR_iMetStep_Now)
          tmp3d_1_sp(1:nxmax_g2s,:,:) = MR_dum3d_compH(1:nxmax_g2s,:,:)
          call MR_Read_3d_Met_Variable_to_CompGrid(ivar,MR_iMetStep_Now+1)
          tmp3d_2_sp(1:nxmax_g2s,:,:) = MR_dum3d_compH(1:nxmax_g2s,:,:)
          vx_Met_loc_sp(1:nxmax_g2s,1:nymax_g2s,1:nzmax_Met2) = &
                  tmp3d_1_sp(1:nxmax_g2s,1:nymax_g2s,1:nzmax_Met2) + &
                 (tmp3d_2_sp(1:nxmax_g2s,1:nymax_g2s,1:nzmax_Met2) - &
                   tmp3d_1_sp(1:nxmax_g2s,1:nymax_g2s,1:nzmax_Met2)) * &
                                         real(Interval_Frac,kind=4)
          ivar = 3 ! Vy
          call MR_Read_3d_Met_Variable_to_CompGrid(ivar,MR_iMetStep_Now)
          tmp3d_1_sp(1:nxmax_g2s,:,:) = MR_dum3d_compH(1:nxmax_g2s,:,:)
          call MR_Read_3d_Met_Variable_to_CompGrid(ivar,MR_iMetStep_Now+1)
          tmp3d_2_sp(1:nxmax_g2s,:,:) = MR_dum3d_compH(1:nxmax_g2s,:,:)
          vy_Met_loc_sp(1:nxmax_g2s,1:nymax_g2s,1:nzmax_Met2) = &
                  tmp3d_1_sp(1:nxmax_g2s,1:nymax_g2s,1:nzmax_Met2) + &
                 (tmp3d_2_sp(1:nxmax_g2s,1:nymax_g2s,1:nzmax_Met2) - &
                   tmp3d_1_sp(1:nxmax_g2s,1:nymax_g2s,1:nzmax_Met2)) * &
                                         real(Interval_Frac,kind=4)
        else
          allocate(tmp3d_1_2_sp(nxmax_g2s,nymax_g2s,nzmax_Met2))
          allocate(tmp3d_2_2_sp(nxmax_g2s,nymax_g2s,nzmax_Met2))
  
          if(Map_Case.eq.5)then
            ! In this case, we first need to rotate the projected Met data to
            ! Earth-relative on the MetP grid
            call MR_Rotate_UV_GR2ER_Met(MR_iMetStep_Now)
          endif
          ! In this case, the comp and met2 grids differ so we will need to rotate
          ! the velocity vectors
          call MR_Rotate_UV_ER2GR_Comp(MR_iMetStep_Now)
          ! U for comp grid is stored in MR_dum3d_compH
          tmp3d_1_sp(1:nxmax_g2s,:,:) = MR_dum3d_compH(1:nxmax_g2s,:,:)
          ! V for comp grid is stored in MR_dum3d_compH_2
          tmp3d_1_2_sp(1:nxmax_g2s,:,:) = MR_dum3d_compH_2(1:nxmax_g2s,:,:)


          if(Map_Case.eq.5)then
            ! In this case, we first need to rotate the projected Met data to
            ! Earth-relative on the MetP grid
            call MR_Rotate_UV_GR2ER_Met(MR_iMetStep_Now+1)
          endif
          call MR_Rotate_UV_ER2GR_Comp(MR_iMetStep_Now+1)
          ! U for comp grid is stored in MR_dum3d_compH
          tmp3d_2_sp(1:nxmax_g2s,:,:) = MR_dum3d_compH(1:nxmax_g2s,:,:)
          ! V for comp grid is stored in MR_dum3d_compH_2
          tmp3d_2_2_sp(1:nxmax_g2s,:,:) = MR_dum3d_compH_2(1:nxmax_g2s,:,:)
  
          vx_Met_loc_sp(1:nxmax_g2s,1:nymax_g2s,1:nzmax_Met2) = &
                  tmp3d_1_sp(1:nxmax_g2s,1:nymax_g2s,1:nzmax_Met2) + &
                 (tmp3d_2_sp(1:nxmax_g2s,1:nymax_g2s,1:nzmax_Met2) - &
                   tmp3d_1_sp(1:nxmax_g2s,1:nymax_g2s,1:nzmax_Met2)) * &
                                         real(Interval_Frac,kind=4)
  
          vy_Met_loc_sp(1:nxmax_g2s,1:nymax_g2s,1:nzmax_Met2) = &
                  tmp3d_1_2_sp(1:nxmax_g2s,1:nymax_g2s,1:nzmax_Met2) + &
                 (tmp3d_2_2_sp(1:nxmax_g2s,1:nymax_g2s,1:nzmax_Met2) - &
                   tmp3d_1_2_sp(1:nxmax_g2s,1:nymax_g2s,1:nzmax_Met2)) * &
                                         real(Interval_Frac,kind=4)
          deallocate(tmp3d_1_2_sp)
          deallocate(tmp3d_2_2_sp)
  
        endif

        ivar = 5 ! Temperature
        call MR_Read_3d_Met_Variable_to_CompGrid(ivar,MR_iMetStep_Now)
        tmp3d_1_sp(1:nxmax_g2s,:,:) = MR_dum3d_compH(1:nxmax_g2s,:,:)
        call MR_Read_3d_Met_Variable_to_CompGrid(ivar,MR_iMetStep_Now+1)
        tmp3d_2_sp(1:nxmax_g2s,:,:) = MR_dum3d_compH(1:nxmax_g2s,:,:)
        temperature_Met_loc_sp(1:nxmax_g2s,1:nymax_g2s,1:nzmax_Met2) = &
                tmp3d_1_sp(1:nxmax_g2s,1:nymax_g2s,1:nzmax_Met2) + &
               (tmp3d_2_sp(1:nxmax_g2s,1:nymax_g2s,1:nzmax_Met2) - &
                 tmp3d_1_sp(1:nxmax_g2s,1:nymax_g2s,1:nzmax_Met2)) * &
                                         real(Interval_Frac,kind=4)
        deallocate(tmp3d_1_sp,tmp3d_2_sp)

        if(write_test_output)then
          ! Dump out map at bottom of Met grid
          open(unit=30,file='plan_Met2.dat',status='replace')
          do i=1,nxmax_g2s
            do j=1,nymax_g2s
              write(30,*)x_g2s_sp(i),y_g2s_sp(j),temperature_Met_loc_sp(i,j,1),&
                        vx_Met_loc_sp(i,j,1),vy_Met_loc_sp(i,j,1)
            enddo
          enddo
          close(30)
          ! Dump out equatorial cross-section of Met grid
          open(unit=40,file='xsec_Met2.dat',status='replace')
          do i=1,nxmax_g2s
            do k=1,nzmax_Met2
              write(40,*)x_g2s_sp(i),y_g2s_sp(j),z_Met2_sp(k),&
                         temperature_Met_loc_sp(i,jout,k),&
                         vx_Met_loc_sp(i,jout,k),vy_Met_loc_sp(i,jout,k)
            enddo
          enddo
          close(40)
          ! Dump out meridonal cross-section of Met grid
          open(unit=40,file='ysec_Met2.dat',status='replace')
          do j=1,nymax_g2s
            do k=1,nzmax_Met2
              write(40,*)x_g2s_sp(iout),y_g2s_sp(j),z_Met2_sp(k),&
                         temperature_Met_loc_sp(iout,j,k),&
                         vx_Met_loc_sp(iout,j,k),vy_Met_loc_sp(iout,j,k)
            enddo
          enddo
          close(40)

        endif
 
        if (IsLatLon_CompGrid.eqv..true.)then
          write(G2S_global_info,*)"   **Calculating spherical harmonics of met data.**"
          call Get_Met_SH(2)
        else
          write(G2S_global_info,*)"   **Calculating Fourier decomposition of met data.**"
          call Get_Met_FS(2)
        endif
        deallocate(vx_Met_loc_sp)
        deallocate(vy_Met_loc_sp)
        deallocate(temperature_Met_loc_sp)
      else
        nzmax_Met2 = 0
      endif ! nwindfiles.eq.2

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Now working on HWT
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      nzmax_HWT = int((zmax_HWT-zmin_HWT)/dz_HWT) + 1
      write(G2S_global_info,*)" Allocating HWT grids of size: ",nxmax_g2s,nymax_g2s,nzmax_HWT
      allocate(z_HWT_g2s_sp(nzmax_HWT))
      allocate(vx_HWT_sp(nxmax_g2s,nymax_g2s,nzmax_HWT))
      allocate(vy_HWT_sp(nxmax_g2s,nymax_g2s,nzmax_HWT))
      allocate(temperature_HWT_sp(nxmax_g2s,nymax_g2s,nzmax_HWT))
      full_data_len = nzmax_HWT+nzmax_Met1+nzmax_Met2
      allocate(data_alt(full_data_len))

      do k = 1,nzmax_HWT
        z_HWT_g2s_sp(k) = zmin_HWT + (k-1) * dz_HWT
      enddo

      if(LOAD_HWT)then
        write(G2S_global_info,*)"Reading precomputed empirical vx,vy,temperature."
        open (unit=20,file=file_U_HWT,access='direct', &
                     recl=nxmax_g2s*nymax_g2s*nzmax_HWT,iostat=ios,status='old')
        if(ios.ne.0)then
          write(G2S_global_error,*)"ERROR: Could not open ",file_U_HWT
          stop 1
        endif
        read(20,rec=1) (((vx_HWT_sp(i,j,k),i=1,nxmax_g2s),j=1,nymax_g2s),k=1,nzmax_HWT)
        close (20)

        open (unit=20,file=file_V_HWT,access='direct', &
                     recl=nxmax_g2s*nymax_g2s*nzmax_HWT,iostat=ios,status='old')
        if(ios.ne.0)then
          write(G2S_global_error,*)"ERROR: Could not open ",file_U_HWT
          stop 1
        endif
        read(20,rec=1)(((vy_HWT_sp(i,j,k),i=1,nxmax_g2s),j=1,nymax_g2s),k=1,nzmax_HWT)
        close (20)

        open (unit=20,file=file_T_HWT,access='direct', &
                     recl=nxmax_g2s*nymax_g2s*nzmax_HWT,iostat=ios,status='old')
        if(ios.ne.0)then
          write(G2S_global_error,*)"ERROR: Could not open ",file_T_HWT
          stop 1
        endif
        read(20,rec=1)(((temperature_HWT_sp(i,j,k),i=1,nxmax_g2s),j=1,nymax_g2s),k=1,nzmax_HWT)
        close (20)
      else

#ifdef useHWM07
        write(G2S_global_info,*)"Computing vx,vy,temperature from HWM07 and MSIS"
#endif
#ifdef useHWM14
        write(G2S_global_info,*)"Computing vx,vy,temperature from HWM14 and MSIS"
#endif


        do i = 1,nxmax_g2s
          write(G2S_global_info,*)"Calculating x=",i," of ",nxmax_g2s
          do j = 1,nymax_g2s
            if (IsLatLon_CompGrid.eqv..true.)then
              lon = x_g2s_sp(i)
              lat = y_g2s_sp(j)
            else
              ! Need to submit lon/lat of this projected point to the HWT
              ! calculator
              lon    = Lon_of_proj_node_dp(i,j)
              lat    = Lat_of_proj_node_dp(i,j)
              rotang = Theta_of_proj_node_dp(i,j)
            endif

            do k = 1,nzmax_HWT
              alt = z_HWT_g2s_sp(k)
              call Get_WindTemp_Empir(lon,lat,alt,              &
                          start_year,day,ihour,iminute,isecond, &
                          ap,f107,                              &
                          u_g2s,v_g2s,temperature)
              if (IsLatLon_CompGrid.eqv..true.)then
                vx_HWT_sp(i,j,k)          = u_g2s
                vy_HWT_sp(i,j,k)          = v_g2s
              else
                ! Now we need to rotate the earth relative components to grid
                ! relative
                vx_HWT_sp(i,j,k) = u_g2s * cos(rotang) - v_g2s * sin(rotang)
                vy_HWT_sp(i,j,k) = u_g2s * sin(rotang) + v_g2s * cos(rotang)
              endif
              temperature_HWT_sp(i,j,k) = temperature
            enddo
          enddo
        enddo
      endif

      if(write_test_output)then
        ! Dump out map at bottom of HWT grid (~20km)
        open(unit=31,file='plan_HWT.dat',status='replace')
        do i=1,nxmax_g2s
          do j=1,nymax_g2s
            write(31,*)x_g2s_sp(i),y_g2s_sp(j),temperature_HWT_sp(i,j,1),&
                      vx_HWT_sp(i,j,1),vy_HWT_sp(i,j,1)
          enddo
        enddo
        close(31)
        ! Dump out equatorial cross-section of HWT grid
        open(unit=41,file='xsec_HWT.dat',status='replace')
        do i=1,nxmax_g2s
          do k=1,nzmax_HWT
            write(41,*)x_g2s_sp(i),y_g2s_sp(jout),z_HWT_g2s_sp(k),&
                       temperature_HWT_sp(i,jout,k),&
                       vx_HWT_sp(i,jout,k),vy_HWT_sp(i,jout,k)
          enddo
        enddo
        close(41)
        ! Dump out meridonal cross-section of HWT grid
        open(unit=41,file='ysec_HWT.dat',status='replace')
        do j=1,nymax_g2s
          do k=1,nzmax_HWT
            write(41,*)x_g2s_sp(iout),y_g2s_sp(j),z_HWT_g2s_sp(k),&
                       temperature_HWT_sp(iout,j,k),&
                       vx_HWT_sp(iout,j,k),vy_HWT_sp(iout,j,k)
          enddo
        enddo
        close(41)

      endif

      if(write_test_output)then
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !  Write HWT data to raw files for later 
        write(G2S_global_info,*)"Writing output for HWT files"
        open (unit=20,file=file_U_HWT,access='direct', &
                     recl=nxmax_g2s*nymax_g2s*nzmax_HWT*4)
        write(20,rec=1)(((vx_HWT_sp(i,j,k),i=1,nxmax_g2s),j=1,nymax_g2s),k=1,nzmax_HWT)
        close (20)
        open (unit=20,file=file_V_HWT,access='direct', &
                     recl=nxmax_g2s*nymax_g2s*nzmax_HWT*4)
        write(20,rec=1)(((vy_HWT_sp(i,j,k),i=1,nxmax_g2s),j=1,nymax_g2s),k=1,nzmax_HWT)
        close (20)
        open (unit=20,file=file_T_HWT,access='direct', &
                     recl=nxmax_g2s*nymax_g2s*nzmax_HWT*4)
        write(20,rec=1)(((temperature_HWT_sp(i,j,k),i=1,nxmax_g2s),j=1,nymax_g2s),k=1,nzmax_HWT)
        close (20)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      endif

      if (IsLatLon_CompGrid.eqv..true.)then
        write(G2S_global_info,*)"   Calculating spherical harmonics of HWT data."
        call Get_HWT_SH
      else
        write(G2S_global_info,*)"   Calculating Fourier decomposition of HWT data."
        call Get_HWT_FS
      endif
      write(G2S_global_info,*)"Fit B-Spline along vertical coordinate"
      call Merge_HWT_Met
      write(G2S_global_info,*)"Writing spectra to netcdf file."
      call Write_Spec_nc

      write(G2S_global_info,*)"Exited normally."
      end program G2S_genSC


!##############################################################################
!##############################################################################
!
!  SUBROUTINES
!
!##############################################################################
!
!     Get_Met_SH
!
!     This subroutine takes the data from the met file and calculates
!     the coefficients of the spherical harmonic decomposition.
!     NOTE:  a scalar spherical harmonic decomposition is performed using the 
!            SHTOOLS library.  This is as opposed to the vector spherical
!            harmonic decomposition of the horizontal winds that Drob
!            recommends.  The scalar decomposition does not ensure solenoidal
!            winds remain divergence-free, but it seems to be good enough.
!     Global variables filled are:
!        vx1_SH_dp          
!        vy1_SH_dp          
!        temperature1_SH_dp 
!
!##############################################################################

      subroutine Get_Met_SH(imet)

      use G2S_globvar
      use shtools

      implicit none

      integer :: imet

      real(kind=8),dimension(:,:,:),allocatable ::  cilm
      real(kind=8),dimension(:,:)  ,allocatable ::  grid
      real(kind=8),dimension(:,:)  ,allocatable ::  griddh
      real(kind=8) :: interval
      integer :: n, lmax2, j, k
      integer :: nzmax_Met

      write(G2S_global_info,*)"Inside Get_Met_SH",imet

      if(imet.eq.1)then
        nzmax_Met = nzmax_Met1
        allocate(         vx1_SH_dp(2, Spec_MaxOrder+1, Spec_MaxOrder+1,nzmax_Met))
        allocate(         vy1_SH_dp(2, Spec_MaxOrder+1, Spec_MaxOrder+1,nzmax_Met))
        allocate(temperature1_SH_dp(2, Spec_MaxOrder+1, Spec_MaxOrder+1,nzmax_Met))
        vx1_SH_dp          = 0.0_8
        vy1_SH_dp          = 0.0_8
        temperature1_SH_dp = 0.0_8
      else
        nzmax_Met = nzmax_Met2
        allocate(         vx1_2_SH_dp(2, Spec_MaxOrder+1, Spec_MaxOrder+1,nzmax_Met))
        allocate(         vy1_2_SH_dp(2, Spec_MaxOrder+1, Spec_MaxOrder+1,nzmax_Met))
        allocate(temperature1_2_SH_dp(2, Spec_MaxOrder+1, Spec_MaxOrder+1,nzmax_Met))
        vx1_2_SH_dp          = 0.0_8
        vy1_2_SH_dp          = 0.0_8
        temperature1_2_SH_dp = 0.0_8
      endif

      n = int(nxmax_g2s*0.5)
      
      interval = real(dx_g2s,kind=8)

      allocate(cilm(2, Spec_MaxOrder+1, Spec_MaxOrder+1))
      allocate(grid(  nymax_g2s+1, nxmax_g2s+1))
      allocate(griddh(nxmax_g2s, nxmax_g2s))

      cilm   = 0.0_8
      grid   = 0.0_8
      griddh = 0.0_8

      do k = 1,nzmax_Met
        write(G2S_global_info,*)"Decomposing Spherical Harmonics for vx: ",k
        griddh = 0.0
        do j = 1,n
          griddh(j,1:nxmax_g2s) = vx_Met_loc_sp(1:nxmax_g2s,j,k)
        enddo
        ! Calculate spherical harmonic coefficients out to degree maxdeg
        call SHExpandDH(griddh, &  ! i contains the gridded data (n,2*n)
                        n,      &  ! i number of samples in lat of gridded data
                        cilm,   &  ! o real SH coefficients (2, LMAX_CALC+1,LMAX_CALC+1)
                        lmax2,  &  ! o max SH bandwidth
                        sampling=2, &        ! i 2=equal spacing
                        LMAX_CALC=Spec_MaxOrder-1 & ! i max SH degree
                        )
        if(imet.eq.1)then
          vx1_SH_dp(1:2, 1:Spec_MaxOrder+1, 1:Spec_MaxOrder+1, k) = &
               cilm(1:2, 1:Spec_MaxOrder+1, 1:Spec_MaxOrder+1)
        else
          vx1_2_SH_dp(1:2, 1:Spec_MaxOrder+1, 1:Spec_MaxOrder+1, k) = &
                 cilm(1:2, 1:Spec_MaxOrder+1, 1:Spec_MaxOrder+1)
        endif
      enddo

      do k = 1,nzmax_Met
        write(G2S_global_info,*)"Decomposing Spherical Harmonics for vy: ",k
        griddh = 0.0
        do j = 1,n
          griddh(j,1:nxmax_g2s) = vy_Met_loc_sp(1:nxmax_g2s,j,k)
        enddo
        ! Calculate spherical harmonic coefficients out to degree maxdeg
        call SHExpandDH(griddh, &  ! i contains the gridded data (n,2*n)
                        n,      &  ! i number of samples in lat of gridded data
                        cilm,   &  ! o real SH coefficients (2,LMAX_CALC+1,LMAX_CALC+1)
                        lmax2,  &  ! o max SH bandwidth
                        sampling=2, &        ! i 2=equal spacing
                        LMAX_CALC=Spec_MaxOrder-1 & ! i max SH degree
                        )
        if(imet.eq.1)then
          vy1_SH_dp(1:2, 1:Spec_MaxOrder+1, 1:Spec_MaxOrder+1, k) = &
               cilm(1:2, 1:Spec_MaxOrder+1, 1:Spec_MaxOrder+1)
        else
          vy1_2_SH_dp(1:2, 1:Spec_MaxOrder+1, 1:Spec_MaxOrder+1, k) = &
                 cilm(1:2, 1:Spec_MaxOrder+1, 1:Spec_MaxOrder+1)
        endif
      enddo

      do k = 1,nzmax_Met
        write(G2S_global_info,*)"Decomposing Spherical Harmonics for temperature: ",k
        griddh = 0.0
        do j = 1,n
          griddh(j,1:nxmax_g2s) = temperature_Met_loc_sp(1:nxmax_g2s,j,k)
        enddo
        ! Calculate spherical harmonic coefficients out to degree maxdeg
        call SHExpandDH(griddh, &  ! i contains the gridded data (n,2*n)
                        n,      &  ! i number of samples in lat of gridded data
                        cilm,   &  ! o real SH coefficients (2,LMAX_CALC+1,LMAX_CALC+1)
                        lmax2,  &  ! o max SH bandwidth
                        sampling=2, &        ! i 2=equal spacing
                        LMAX_CALC=Spec_MaxOrder-1 & ! i max SH degree
                        )
        if(imet.eq.1)then
          temperature1_SH_dp(1:2, 1:Spec_MaxOrder+1, 1:Spec_MaxOrder+1, k) = &
                        cilm(1:2, 1:Spec_MaxOrder+1, 1:Spec_MaxOrder+1)
        else
          temperature1_2_SH_dp(1:2, 1:Spec_MaxOrder+1, 1:Spec_MaxOrder+1, k) = &
                          cilm(1:2,1:Spec_MaxOrder+1, 1:Spec_MaxOrder+1)
        endif
      enddo

      deallocate(cilm)
      deallocate(grid)
      deallocate(griddh)

      write(G2S_global_info,*)"    Leaving Get_Met_SH"

      end subroutine Get_Met_SH

!##############################################################################
!##############################################################################
!
!     Get_HWT_SH
!
!     This subroutine takes the data from the HWT models and calculates
!     the coefficients of the spherical harmonic decomposition.
!     NOTE:  a scalar spherical harmonic decomposition is performed using the 
!            SHTOOLS library.  This is as opposed to the vector spherical
!            harmonic decomposition of the horizontal winds that Drob
!            recommends.  The scalar decomposition does not ensure solenoidal
!            winds remain divergence-free, but it seems to be good enough.
!     Global variables filled are:
!        vx2_SH_dp          
!        vy2_SH_dp          
!        temperature2_SH_dp 
!
!##############################################################################

      subroutine Get_HWT_SH

      use G2S_globvar
      use SHTOOLS

      implicit none

      real(kind=8),dimension(:,:,:),allocatable ::  cilm
      real(kind=8),dimension(:,:)  ,allocatable ::  grid
      real(kind=8),dimension(:,:)  ,allocatable ::  griddh
      real(kind=8) :: interval
      integer :: n, lmax2, j, k

      allocate(         vx2_SH_dp(2, Spec_MaxOrder+1, Spec_MaxOrder+1,nzmax_HWT))
      allocate(         vy2_SH_dp(2, Spec_MaxOrder+1, Spec_MaxOrder+1,nzmax_HWT))
      allocate(temperature2_SH_dp(2, Spec_MaxOrder+1, Spec_MaxOrder+1,nzmax_HWT))

      interval = real(dx_g2s,kind=8)

      allocate(cilm(2, Spec_MaxOrder+1, Spec_MaxOrder+1))
      allocate(grid(  nymax_g2s+1, nxmax_g2s+1))
      allocate(griddh(nxmax_g2s, nxmax_g2s))

      n = int(nxmax_g2s*0.5)

      do k = 1,nzmax_HWT
        !write(G2S_global_info,*)"Decomposing Spherical Harmonics for vx: ",k
        griddh = 0.0
        do j = 1,n
          griddh(j,1:nxmax_g2s) = vx_HWT_sp(1:nxmax_g2s,j,k)
        enddo
        ! Calculate spherical harmonic coefficients out to degree maxdeg
        call SHExpandDH(griddh, &  ! i contains the gridded data (n,2*n)
                        n,      &  ! i number of samples in lat of gridded data
                        cilm,   &  ! o real SH coefficients (2,LMAX_CALC+1,LMAX_CALC+1)
                        lmax2,  &  ! o max SH bandwidth
                        sampling=2, &        ! i 2=equal spacing
                        LMAX_CALC=Spec_MaxOrder-1 & ! i max SH degree
                        )
        vx2_SH_dp(1:2, 1:Spec_MaxOrder+1, 1:Spec_MaxOrder+1, k) = &
             cilm(1:2, 1:Spec_MaxOrder+1, 1:Spec_MaxOrder+1)
      enddo

      do k = 1,nzmax_HWT
        !write(G2S_global_info,*)"Decomposing Spherical Harmonics for vy: ",k
        griddh = 0.0
        do j = 1,n
          griddh(j,1:nxmax_g2s) = vy_HWT_sp(1:nxmax_g2s,j,k)
        enddo
        ! Calculate spherical harmonic coefficients out to degree maxdeg
        call SHExpandDH(griddh, &  ! i contains the gridded data (n,2*n)
                        n,      &  ! i number of samples in lat of gridded data
                        cilm,   &  ! o real SH coefficients(2,LMAX_CALC+1,LMAX_CALC+1)
                        lmax2,  &  ! o max SH bandwidth
                        sampling=2, &        ! i 2=equal spacing
                        LMAX_CALC=Spec_MaxOrder-1 & ! i max SH degree
                        )
        vy2_SH_dp(1:2, 1:Spec_MaxOrder+1, 1:Spec_MaxOrder+1, k) = &
             cilm(1:2, 1:Spec_MaxOrder+1, 1:Spec_MaxOrder+1)
      enddo

      do k = 1,nzmax_HWT
        !write(G2S_global_info,*)"Decomposing Spherical Harmonics for vx: ",k
        griddh = 0.0
        do j = 1,n
          griddh(j,1:nxmax_g2s) = temperature_HWT_sp(1:nxmax_g2s,j,k)
        enddo
        ! Calculate spherical harmonic coefficients out to degree maxdeg
        call SHExpandDH(griddh, &  ! i contains the gridded data (n,2*n)
                        n,      &  ! i number of samples in lat of gridded data
                        cilm,   &  ! o real SH coefficients (2,LMAX_CALC+1,LMAX_CALC+1)
                        lmax2,  &  ! o max SH bandwidth
                        sampling=2, &        ! i 2=equal spacing
                        LMAX_CALC=Spec_MaxOrder-1 & ! i max SH degree
                        )
        temperature2_SH_dp(1:2, 1:Spec_MaxOrder+1, 1:Spec_MaxOrder+1, k) = &
                      cilm(1:2, 1:Spec_MaxOrder+1, 1:Spec_MaxOrder+1)
      enddo

      deallocate(cilm)
      deallocate(grid)
      deallocate(griddh)

      end subroutine Get_HWT_SH

!##############################################################################
!
!     Get_Met_FS
!
!     This subroutine takes the data from the met file and calculates
!     the coefficients of the Fourier decomposition.
!     Global variables filled are:
!        vx1_FS_C          
!        vy1_FS_C          
!        temperature1_FS_C 
!
!##############################################################################

      subroutine Get_Met_FS(imet)

      use G2S_globvar

      implicit none

      integer :: imet

      integer :: k
      integer :: nzmax_Met

      real(kind=8)   ,dimension(:,:),allocatable :: gridU_in
      real(kind=8)   ,dimension(:,:),allocatable :: gridV_in
      real(kind=8)   ,dimension(:,:),allocatable :: gridT_in
      complex(kind=8),dimension(:,:),allocatable :: gridU_out
      complex(kind=8),dimension(:,:),allocatable :: gridV_out
      complex(kind=8),dimension(:,:),allocatable :: gridT_out
      integer(kind=8)                            :: plan_for2Dfft_U
      integer(kind=8)                            :: plan_for2Dfft_V
      integer(kind=8)                            :: plan_for2Dfft_T


      !write(G2S_global_info,*)"Inside Get_Met_FS",imet

      allocate(gridU_in( nxmax_g2s      , nymax_g2s))
      allocate(gridV_in( nxmax_g2s      , nymax_g2s))
      allocate(gridT_in( nxmax_g2s      , nymax_g2s))
      allocate(gridU_out(nxmax_g2s/2 + 1, nymax_g2s))
      allocate(gridV_out(nxmax_g2s/2 + 1, nymax_g2s))
      allocate(gridT_out(nxmax_g2s/2 + 1, nymax_g2s))
      call dfftw_plan_dft_r2c_2d(plan_for2Dfft_U, nxmax_g2s, nymax_g2s, &
                                 gridU_in , gridU_out,                  &
                                 FFTW_MEASURE)
      call dfftw_plan_dft_r2c_2d(plan_for2Dfft_V, nxmax_g2s, nymax_g2s, &
                                 gridV_in , gridV_out,                  &
                                 FFTW_MEASURE)
      call dfftw_plan_dft_r2c_2d(plan_for2Dfft_T, nxmax_g2s, nymax_g2s, &
                                 gridT_in , gridT_out,                  &
                                 FFTW_MEASURE)
      if(imet.eq.1)then
        nzmax_Met = nzmax_Met1
        write(G2S_global_info,*)"Allocating arrays of size ",nxmax_g2s/2+1, nymax_g2s,nzmax_Met
        allocate(         vx1_FS_C(nxmax_g2s/2+1, nymax_g2s,nzmax_Met))
        allocate(         vy1_FS_C(nxmax_g2s/2+1, nymax_g2s,nzmax_Met))
        allocate(temperature1_FS_C(nxmax_g2s/2+1, nymax_g2s,nzmax_Met))
        vx1_FS_C          = 0.0_8
        vy1_FS_C          = 0.0_8
        temperature1_FS_C = 0.0_8
      else
        nzmax_Met = nzmax_Met2
        allocate(         vx1_2_FS_C(nxmax_g2s/2+1, nymax_g2s,nzmax_Met))
        allocate(         vy1_2_FS_C(nxmax_g2s/2+1, nymax_g2s,nzmax_Met))
        allocate(temperature1_2_FS_C(nxmax_g2s/2+1, nymax_g2s,nzmax_Met))
        vx1_2_FS_C          = 0.0_8
        vy1_2_FS_C          = 0.0_8
        temperature1_2_FS_C = 0.0_8
      endif

      do k = 1,nzmax_Met
        write(G2S_global_info,*)"Decomposing Fourier coefficients for U,V,T: ",k
        gridU_in = 0.0
        gridV_in = 0.0
        gridT_in = 0.0
        gridU_in(1:nxmax_g2s,1:nymax_g2s) = real(vx_Met_loc_sp(1:nxmax_g2s,1:nymax_g2s,k),kind=8)
        gridV_in(1:nxmax_g2s,1:nymax_g2s) = real(vy_Met_loc_sp(1:nxmax_g2s,1:nymax_g2s,k),kind=8)
        gridT_in(1:nxmax_g2s,1:nymax_g2s) = real(temperature_Met_loc_sp(1:nxmax_g2s,1:nymax_g2s,k),kind=8)

        ! Calculate Fourier coefficients
        call dfftw_execute_dft_r2c(plan_for2Dfft_U,gridU_in,gridU_out)
        call dfftw_execute_dft_r2c(plan_for2Dfft_V,gridV_in,gridV_out)
        call dfftw_execute_dft_r2c(plan_for2Dfft_T,gridT_in,gridT_out)

        if(imet.eq.1)then
          vx1_FS_C(:,:,k)          = gridU_out(:,:)
          vy1_FS_C(:,:,k)          = gridV_out(:,:)
          temperature1_FS_C(:,:,k) = gridT_out(:,:)
        else
          vx1_2_FS_C(:,:,k)          = gridU_out(:,:)
          vy1_2_FS_C(:,:,k)          = gridV_out(:,:)
          temperature1_2_FS_C(:,:,k) = gridT_out(:,:)
        endif
      enddo

      deallocate(gridU_in)
      deallocate(gridV_in)
      deallocate(gridT_in)
      deallocate(gridU_out)
      deallocate(gridV_out)
      deallocate(gridT_out)

      call dfftw_destroy_plan(plan_for2Dfft_U)
      call dfftw_destroy_plan(plan_for2Dfft_V)
      call dfftw_destroy_plan(plan_for2Dfft_T)

      end subroutine Get_Met_FS

!##############################################################################
!
!     Get_HWT_FS
!
!     This subroutine takes the data from the HWT model and calculates
!     the coefficients of the Fourier decomposition.
!     Global variables filled are:
!        vx2_FS_C          
!        vy2_FS_C          
!        temperature2_FS_C 
!
!##############################################################################

      subroutine Get_HWT_FS

      use G2S_globvar

      implicit none

      integer :: k

      real(kind=8)   ,dimension(:,:),allocatable :: gridU_in
      real(kind=8)   ,dimension(:,:),allocatable :: gridV_in
      real(kind=8)   ,dimension(:,:),allocatable :: gridT_in
      complex(kind=8),dimension(:,:),allocatable :: gridU_out
      complex(kind=8),dimension(:,:),allocatable :: gridV_out
      complex(kind=8),dimension(:,:),allocatable :: gridT_out
      integer(kind=8)                            :: plan_for2Dfft_U
      integer(kind=8)                            :: plan_for2Dfft_V
      integer(kind=8)                            :: plan_for2Dfft_T

      !write(G2S_global_info,*)"Inside Get_HWT_FS"

      allocate(gridU_in( nxmax_g2s      , nymax_g2s))
      allocate(gridV_in( nxmax_g2s      , nymax_g2s))
      allocate(gridT_in( nxmax_g2s      , nymax_g2s))
      allocate(gridU_out(nxmax_g2s/2 + 1, nymax_g2s))
      allocate(gridV_out(nxmax_g2s/2 + 1, nymax_g2s))
      allocate(gridT_out(nxmax_g2s/2 + 1, nymax_g2s))
      call dfftw_plan_dft_r2c_2d(plan_for2Dfft_U, nxmax_g2s, nymax_g2s, &
                                 gridU_in , gridU_out,                  &
                                 FFTW_MEASURE)
      call dfftw_plan_dft_r2c_2d(plan_for2Dfft_V, nxmax_g2s, nymax_g2s, &
                                 gridV_in , gridV_out,                  &
                                 FFTW_MEASURE)
      call dfftw_plan_dft_r2c_2d(plan_for2Dfft_T, nxmax_g2s, nymax_g2s, &
                                 gridT_in , gridT_out,                  &
                                 FFTW_MEASURE)
      write(G2S_global_info,*)"Allocating arrays of size ",nxmax_g2s/2+1, nymax_g2s,nzmax_HWT
      allocate(         vx2_FS_C(nxmax_g2s/2+1, nymax_g2s,nzmax_HWT))
      allocate(         vy2_FS_C(nxmax_g2s/2+1, nymax_g2s,nzmax_HWT))
      allocate(temperature2_FS_C(nxmax_g2s/2+1, nymax_g2s,nzmax_HWT))
      vx2_FS_C          = 0.0_8
      vy2_FS_C          = 0.0_8
      temperature2_FS_C = 0.0_8

      do k = 1,nzmax_HWT
        write(G2S_global_info,*)"Decomposing Fourier coefficients for U,V,T: ",k
        gridU_in = 0.0
        gridV_in = 0.0
        gridT_in = 0.0
        gridU_in(1:nxmax_g2s,1:nymax_g2s) = real(vx_HWT_sp(1:nxmax_g2s,1:nymax_g2s,k),kind=8)
        gridV_in(1:nxmax_g2s,1:nymax_g2s) = real(vy_HWT_sp(1:nxmax_g2s,1:nymax_g2s,k),kind=8)
        gridT_in(1:nxmax_g2s,1:nymax_g2s) = real(temperature_HWT_sp(1:nxmax_g2s,1:nymax_g2s,k),kind=8)

        ! Calculate Fourier coefficients
        call dfftw_execute_dft_r2c(plan_for2Dfft_U,gridU_in,gridU_out)
        call dfftw_execute_dft_r2c(plan_for2Dfft_V,gridV_in,gridV_out)
        call dfftw_execute_dft_r2c(plan_for2Dfft_T,gridT_in,gridT_out)

        vx2_FS_C(:,:,k)          = gridU_out(:,:)
        vy2_FS_C(:,:,k)          = gridV_out(:,:)
        temperature2_FS_C(:,:,k) = gridT_out(:,:)
      enddo

      deallocate(gridU_in)
      deallocate(gridV_in)
      deallocate(gridT_in)
      deallocate(gridU_out)
      deallocate(gridV_out)
      deallocate(gridT_out)

      call dfftw_destroy_plan(plan_for2Dfft_U)
      call dfftw_destroy_plan(plan_for2Dfft_V)
      call dfftw_destroy_plan(plan_for2Dfft_T)

      end subroutine Get_HWT_FS


!##############################################################################

!##############################################################################
!##############################################################################
!
!     Merge_HWT_Met
!
!     This subroutine takes the SH coefficients for both Met and HWT, and for
!     each x,y node, solves for the best-fit B-spline along the z-direction.
!     One point to note is that the (var)[1,2,3]_SH_dp all use different
!     vertical ranges:
!       (var)1_SH_dp : met file  ( 0 ->  30 km or so)
!       (var)2_SH_dp : HWT file  (20 -> 200 km or so)
!       (var)3_SH_dp : Bspline file  (around 73 knot points which must be
!                                     consistent in the profiler!)
!
!     Global variables filled are:
!       vx3_SH_dp(2, maxdeg+1, maxdeg+1,model_len))
!       vy3_SH_dp(2, maxdeg+1, maxdeg+1,model_len))
!       temperature3_SH_dp(2, maxdeg+1, maxdeg+1,model_len))
!
!##############################################################################

      subroutine Merge_HWT_Met

      use G2S_globvar
      use SHTOOLS
      use MetReader

      implicit none

      real(kind=8) :: data_r(full_data_len)
      real(kind=8) :: data_i(full_data_len)

      integer :: i,j,k,ii,iii,jj,ioff,Mfull
      integer :: ll
      integer :: model_len
      real(kind=8) :: alt

      logical :: First_time

      integer          M, N, NRHS
      integer          LDA,LDB
      integer          LWMAX
      integer          INFO, LWORK, RANK
      real             RCOND
      integer,dimension(:) ,allocatable :: IWORK    ! ( 3*M*0+11*M )
      real(kind=8),dimension(:,:),allocatable :: A        ! ( LDA, N )
      real(kind=8),dimension(:,:),allocatable :: B        ! ( LDB, NRHS )
      real(kind=8),dimension(:)  ,allocatable :: S        ! ( M )
      real(kind=8),dimension(:)  ,allocatable :: JPVT     ! ( N )
      real(kind=8),dimension(:)  ,allocatable :: WORK     ! ( LWMAX )

      real(kind=8),dimension(nzmax_Met1) :: z_Met1_dp
      real(kind=8),dimension(nzmax_Met2) :: z_Met2_dp

      real(kind=8),dimension(:,:,:),allocatable :: mm
      real(kind=8),dimension(:,:,:),allocatable :: mmm

      complex(kind=8) :: dum_C

      z_Met1_dp = z_Met1_sp
      if(nwindfile_groups.eq.2)then
        z_Met2_dp = z_Met2_sp
      endif

      data_len  = full_data_len
      Mfull     = Nknots+P_ord+P_ord
      model_len = Nknots+P_ord

      First_time = .true.

      if (IsLatLon_CompGrid.eqv..true.)then
        x_order = Spec_MaxOrder+1
        y_order = Spec_MaxOrder+1
        allocate(         vx3_SH_dp(2, x_order,y_order,model_len))
        allocate(         vy3_SH_dp(2, x_order,y_order,model_len))
        allocate(temperature3_SH_dp(2, x_order,y_order,model_len))
        !allocate(              cilm(2, x_order,y_order))
      else
        x_order = nxmax_g2s/2+1
        y_order = nymax_g2s
        allocate(         vx3_FS_C(x_order,y_order,model_len))
        allocate(         vy3_FS_C(x_order,y_order,model_len))
        allocate(temperature3_FS_C(x_order,y_order,model_len))
      endif
      !allocate(grid(  nymax_g2s+1, nxmax_g2s+1))
      !allocate(griddh(nxmax_g2s, nxmax_g2s))

      ! set lapack values
      M     = data_len
      N     = model_len
      NRHS  = 2
      LDA   = max(1,M)
      LDB   = max(1,max(M,N))
      LWMAX = 20000
      RCOND = -1.0
      allocate(IWORK( 3*M*0+11*M ))
      allocate(A( LDA, N ))
      allocate(B( LDB, NRHS ))
      allocate(JPVT( N ))
      allocate(S( M ))
      allocate(WORK ( LWMAX ))

      allocate(mm(Mfull,P_ord+1,data_len))

      allocate(mmm(Mfull,P_ord+1,63))
      allocate(new_data_alt(63))
      do k = 1,15
        new_data_alt(k) = (k-1)*1.0_4
      enddo
      do k = 1,18
        new_data_alt(k+15) = 14.0_4+(k)*2.0_4
      enddo
      do k = 1,30
        new_data_alt(k+15+18) = 50.0_4+(k)*5.0_4
      enddo

      ! Define knot vector
      do i = 1,Mfull
        if(i.le.P_ord)then
          knots(i) = iknots(1)
        elseif(i.ge.Nknots+P_ord)then
          knots(i) = iknots(Nknots)
        else
          knots(i) = iknots(i-P_ord)
        endif
      enddo

      ! Set the altitudes for the HWT data
      !   Start at around 15km or 20 km
      data_alt(1:nzmax_HWT) = z_HWT_g2s_sp(1:nzmax_HWT)
      data_alt(nzmax_HWT+1:nzmax_HWT+nzmax_Met1)   = real(z_Met1_dp(1:nzmax_Met1),kind=4)
      if(nwindfile_groups.eq.2)then
        data_alt(nzmax_HWT+nzmax_Met1+1:nzmax_HWT+nzmax_Met1+nzmax_Met2) = &
                  real(z_Met2_dp(1:nzmax_Met2),kind=4)
      endif
      ! Now get B-Spline basis for the HWT data
      !% Generate coefficient matrix of B-Splines
      ioff = 1
      mm = 0.0
      do k = 1,data_len
        alt = data_alt(k)
        do j=1,Mfull-1
          if (alt.eq.knots(j+1).and.alt.gt.knots(j)) then
            mm (j+1,1,k) = 1.0
          elseif (alt.lt.knots(j+1).and.alt.ge.knots(j)) then
            mm(j,1,k) = 1.0
          endif
        enddo
      enddo
      !% higher order basis functions
      do ii = 2,P_ord+1
        iii = ii-1
        do k=1,data_len
          alt = data_alt(k)
          do j=1,Mfull-iii
            if (knots(j+iii).gt.knots(j)) then
              mm(j,ii,k) = (alt - knots(j) ) / ( knots(j+iii)-knots(j)) * mm(j,ii-1,k)
            endif
            if (j+iii+1.le.Mfull) then
              if (knots(j+iii+1)>knots(j+1))then
                mm(j,ii,k) = mm(j,ii,k) + ( knots(j+iii+1) - alt) /       &
                                          ( knots(j+iii+1) - knots(j+1)) * &
                                             mm(j+1,ii-1,k)
              endif
            endif
          enddo
        enddo
      enddo

      allocate(coeff(data_len,model_len))
      do k = 1,data_len
        do j = 1,Mfull-1
          if (j.ge.model_len) then
            coeff(k,model_len) = 0.0;
          else
            coeff(k,j)         = real(mm(j,P_ord+1,k),kind=4)
          endif
        enddo
      enddo

      ! Start preparing for the least-squares best fit
      ! Fill the data vectors
      do i=1,x_order
        do j=1,y_order
          do ll=1,3
            select case (ll)
            case(1) 
              do k = 1,nzmax_HWT
                if(IsLatLon_CompGrid.eqv..true.)then
                  data_r(k) = vx2_SH_dp(1,i,j,k) ! real HWT
                  data_i(k) = vx2_SH_dp(2,i,j,k) ! imag HWT
                else
                  data_r(k) =  real(vx2_FS_C(i,j,k)) ! real HWT
                  data_i(k) = aimag(vx2_FS_C(i,j,k)) ! imag HWT
                endif
              enddo
              do k = 1,nzmax_Met1
                if(IsLatLon_CompGrid.eqv..true.)then
                  data_r(nzmax_HWT+k) = vx1_SH_dp(1,i,j,k) ! real Met
                  data_i(nzmax_HWT+k) = vx1_SH_dp(2,i,j,k) ! imag Met
                else
                  data_r(nzmax_HWT+k) =  real(vx1_FS_C(i,j,k)) ! real Met
                  data_i(nzmax_HWT+k) = aimag(vx1_FS_C(i,j,k)) ! imag Met
                endif
              enddo
              if(nwindfile_groups.eq.2)then
                do k = 1,nzmax_Met2
                  if(IsLatLon_CompGrid.eqv..true.)then
                    data_r(nzmax_HWT+nzmax_Met1+k) = vx1_2_SH_dp(1,i,j,k) ! real Met
                    data_i(nzmax_HWT+nzmax_Met1+k) = vx1_2_SH_dp(2,i,j,k) ! imag Met
                  else
                    data_r(nzmax_HWT+nzmax_Met1+k) =  real(vx1_2_FS_C(i,j,k)) ! real Met
                    data_i(nzmax_HWT+nzmax_Met1+k) = aimag(vx1_2_FS_C(i,j,k)) ! imag Met
                  endif
                enddo
              endif
            case(2)
              do k = 1,nzmax_HWT
                if(IsLatLon_CompGrid.eqv..true.)then
                  data_r(k) = vy2_SH_dp(1,i,j,k) ! real HWT
                  data_i(k) = vy2_SH_dp(2,i,j,k) ! imag HWT
                else
                  data_r(k) =  real(vy2_FS_C(i,j,k)) ! real HWT
                  data_i(k) = aimag(vy2_FS_C(i,j,k)) ! imag HWT
                endif
              enddo
              do k = 1,nzmax_Met1
                if(IsLatLon_CompGrid.eqv..true.)then
                  data_r(nzmax_HWT+k) = vy1_SH_dp(1,i,j,k) ! real Met
                  data_i(nzmax_HWT+k) = vy1_SH_dp(2,i,j,k) ! imag Met
                else
                  data_r(nzmax_HWT+k) =  real(vy1_FS_C(i,j,k)) ! real Met
                  data_i(nzmax_HWT+k) = aimag(vy1_FS_C(i,j,k)) ! imag Met
                endif
              enddo
              if(nwindfile_groups.eq.2)then
                do k = 1,nzmax_Met2
                  if(IsLatLon_CompGrid.eqv..true.)then
                    data_r(nzmax_HWT+nzmax_Met1+k) = vy1_2_SH_dp(1,i,j,k) ! real Met
                    data_i(nzmax_HWT+nzmax_Met1+k) = vy1_2_SH_dp(2,i,j,k) ! imag Met
                  else
                    data_r(nzmax_HWT+nzmax_Met1+k) =  real(vy1_2_FS_C(i,j,k)) ! real Met
                    data_i(nzmax_HWT+nzmax_Met1+k) = aimag(vy1_2_FS_C(i,j,k)) ! imag Met
                  endif
                enddo
              endif
            case(3)
              do k = 1,nzmax_HWT
                if(IsLatLon_CompGrid.eqv..true.)then
                  data_r(k) = temperature2_SH_dp(1,i,j,k) ! real HWT
                  data_i(k) = temperature2_SH_dp(2,i,j,k) ! imag HWT
                else
                  data_r(k) =  real(temperature2_FS_C(i,j,k)) ! real HWT
                  data_i(k) = aimag(temperature2_FS_C(i,j,k)) ! imag HWT
                endif
              enddo
              do k = 1,nzmax_Met1
                if(IsLatLon_CompGrid.eqv..true.)then
                  data_r(nzmax_HWT+k) = temperature1_SH_dp(1,i,j,k) ! real Met
                  data_i(nzmax_HWT+k) = temperature1_SH_dp(2,i,j,k) ! imag Met
                else
                  data_r(nzmax_HWT+k) =  real(temperature1_FS_C(i,j,k)) ! real Met
                  data_i(nzmax_HWT+k) = aimag(temperature1_FS_C(i,j,k)) ! imag Met
                endif
              enddo
              if(nwindfile_groups.eq.2)then
                do k = 1,nzmax_Met2
                  if(IsLatLon_CompGrid.eqv..true.)then
                    data_r(nzmax_HWT+nzmax_Met1+k) = temperature1_2_SH_dp(1,i,j,k) ! real Met
                    data_i(nzmax_HWT+nzmax_Met1+k) = temperature1_2_SH_dp(2,i,j,k) ! imag Met
                  else
                    data_r(nzmax_HWT+nzmax_Met1+k) =  real(temperature1_2_FS_C(i,j,k)) ! real Met
                    data_i(nzmax_HWT+nzmax_Met1+k) = aimag(temperature1_2_FS_C(i,j,k)) ! imag Met
                  endif
                enddo
              endif
            end select

            do k = 1,data_len
              ! set up to operate on both real and imag vectors
              B(k,1) = data_r(k)
              B(k,2) = data_i(k)
              do jj = 1,model_len
                A(k,jj) = coeff(k,jj)
              enddo
            enddo

            if(First_time)then
              LWORK = -1
              CALL DGELSD( M, N, NRHS, A, LDA, B, LDB, S, RCOND, RANK, WORK, &
                           LWORK, IWORK, INFO )
              LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
              First_time = .false.

              ! Reset A and B
              do k = 1,data_len
                ! set up to operate on both real and imag vectors
                B(k,1) = data_r(k)
                B(k,2) = data_i(k)
                do jj = 1,model_len
                  A(k,jj) = coeff(k,jj)
                enddo
              enddo

            endif

            ! Now solve the least squared problem
            !  dGELSY computes the minimum-norm solution to a real linear least
            !    squares problem:  minimize || A * X - B ||
            !    using a complete orthogonal factorization of A.  A is an M-by-N
            !    matrix which may be rank-deficient.
            !call dgelsy(m, n, nrhs, a, lda, b, ldb, jpvt, rcond, rank, work, lwork, info)

            !  dGELSD computes the minimum-norm solution to a real linear least
            !    squares problem:  minimize 2-norm(| b - A*x |)
            !    using the singular value decomposition (SVD) of A. A is an M-by-N
            !    matrix which may be rank-deficient.
            !call dgelss(m, n, nrhs, a, lda, b, ldb, s, rcond, rank, work, lwork, iwork, info)
            !    Here's using the SVD via divide and conquere
            call dgelsd(m, n, nrhs, a, lda, b, ldb, s, rcond, rank, work, lwork, iwork, info)

            if( INFO.GT.0 ) then
              write(G2S_global_info,*)'The algorithm computing SVD failed to converge;'
              write(G2S_global_info,*)'the least squares solution could not be computed.'
              stop 1
            endif

            select case (ll)
            case(1)
              if(IsLatLon_CompGrid.eqv..true.)then
                vx3_SH_dp(1,i,j,1:model_len) = B(1:model_len,1)
                vx3_SH_dp(2,i,j,1:model_len) = B(1:model_len,2)
              else
                do k = 1,model_len
                  dum_C = complex(B(k,1),B(k,2))
                  vx3_FS_C(i,j,k) = dum_C
                enddo
              endif
            case(2)
              if(IsLatLon_CompGrid.eqv..true.)then
                vy3_SH_dp(1,i,j,1:model_len) = B(1:model_len,1)
                vy3_SH_dp(2,i,j,1:model_len) = B(1:model_len,2)
              else
                do k = 1,model_len
                  dum_C = complex(B(k,1),B(k,2))
                  vy3_FS_C(i,j,k) = dum_C
                enddo
              endif
            case(3)
              if(IsLatLon_CompGrid.eqv..true.)then
                temperature3_SH_dp(1,i,j,1:model_len) = B(1:model_len,1)
                temperature3_SH_dp(2,i,j,1:model_len) = B(1:model_len,2)
              else
                do k = 1,model_len
                  dum_C = complex(B(k,1),B(k,2))
                  temperature3_FS_C(i,j,k) = dum_C
                enddo
              endif
            end select

          enddo
        enddo
      enddo

      deallocate(IWORK,A,B,S,JPVT,WORK)
      deallocate(mm)

      end subroutine Merge_HWT_Met

!##############################################################################
!##############################################################################
!
!     Write_Spec_nc
!
!     This subroutine takes the merged SH coefficients at the knot points
!     and writes out a netcdf file.
!
!##############################################################################

      subroutine Write_Spec_nc()

      use G2S_globvar
      use MetReader
      use G2S_NCvars
      use netcdf

      implicit none

      integer :: nSTAT
      integer :: ncid

      real(kind=8),dimension(:)      ,allocatable :: dum1d_out
      real(kind=8),dimension(:,:,:,:),allocatable :: dum4d_out
      character(len=3),dimension(5) :: odim_names
      character(len=30) :: cdf_title

      odim_names(1) = "r"  ! real/imag
      odim_names(2) = "z"  ! alt
      odim_names(3) = "y"  ! order in y
      odim_names(4) = "x"  ! order in x
      odim_names(5) = "k"  ! knot

      ! Create and open netcdf file
      cdf_title = "G2S atmosphere model"
      nSTAT = nf90_create(file_SH,nf90_clobber, ncid)
      if(nSTAT.ne.NF90_NOERR) &
          write(G2S_global_log,*)'ERROR: create file_SH: ', &
                              nf90_strerror(nSTAT)
      nSTAT = nf90_put_att(ncid,nf90_global,"Title",cdf_title)
      if(nSTAT.ne.NF90_NOERR) &
          write(G2S_global_log,*)'ERROR: put_att title: ',nf90_strerror(nSTAT)

      nSTAT = nf90_put_att(ncid,nf90_global,"Forecast_Hour",fc_hour)
      if(nSTAT.ne.NF90_NOERR) &
          write(G2S_global_log,*)'ERROR: put_att fc_hour: ',nf90_strerror(nSTAT)

      nSTAT = nf90_put_att(ncid,nf90_global,"year",inyear)
      if(nSTAT.ne.NF90_NOERR) &
          write(G2S_global_log,*)'ERROR: put_att year: ',nf90_strerror(nSTAT)
      nSTAT = nf90_put_att(ncid,nf90_global,"month",inmonth)
      if(nSTAT.ne.NF90_NOERR) &
          write(G2S_global_log,*)'ERROR: put_att month: ',nf90_strerror(nSTAT)
      nSTAT = nf90_put_att(ncid,nf90_global,"day",inday)
      if(nSTAT.ne.NF90_NOERR) &
          write(G2S_global_log,*)'ERROR: put_att day: ',nf90_strerror(nSTAT)
      nSTAT = nf90_put_att(ncid,nf90_global,"hour",inhour)
      if(nSTAT.ne.NF90_NOERR) &
          write(G2S_global_log,*)'ERROR: put_att hour: ',nf90_strerror(nSTAT)

      nSTAT = nf90_put_att(ncid,nf90_global,"Num_windfile_groups",nwindfile_groups)
      if(nSTAT.ne.NF90_NOERR) &
          write(G2S_global_log,*)'ERROR: put_att nwindfiles: ',nf90_strerror(nSTAT)
      nSTAT = nf90_put_att(ncid,nf90_global,"iwf1",iwf1)
      if(nSTAT.ne.NF90_NOERR) &
          write(G2S_global_log,*)'ERROR: put_att iwf1: ',nf90_strerror(nSTAT)
      nSTAT = nf90_put_att(ncid,nf90_global,"Num_windfiles1",nwindfiles1)
      if(nSTAT.ne.NF90_NOERR) &
          write(G2S_global_log,*)'ERROR: put_att nwindfiles1: ',nf90_strerror(nSTAT)

      nSTAT = nf90_put_att(ncid,nf90_global,"windfile1",windfile1(1))
      if(nSTAT.ne.NF90_NOERR) &
          write(G2S_global_log,*)'ERROR: put_att windfile1: ',nf90_strerror(nSTAT)
      if(nwindfile_groups.eq.2)then
        nSTAT = nf90_put_att(ncid,nf90_global,"iwf2",iwf2)
        if(nSTAT.ne.NF90_NOERR) &
            write(G2S_global_log,*)'ERROR: put_att iwf2: ',nf90_strerror(nSTAT)
        nSTAT = nf90_put_att(ncid,nf90_global,"Num_windfiles2",nwindfiles2)
        if(nSTAT.ne.NF90_NOERR) &
            write(G2S_global_log,*)'ERROR: put_att nwindfiles2: ',nf90_strerror(nSTAT)
        nSTAT = nf90_put_att(ncid,nf90_global,"windfile2",windfile2(1))
        if(nSTAT.ne.NF90_NOERR) &
            write(G2S_global_log,*)'ERROR: put_att windfile2: ',nf90_strerror(nSTAT)
      endif

      nSTAT = nf90_put_att(ncid,nf90_global,"Ap",ap)
      if(nSTAT.ne.NF90_NOERR) &
          write(G2S_global_log,*)'ERROR: put_att Ap: ',nf90_strerror(nSTAT)
      nSTAT = nf90_put_att(ncid,nf90_global,"F107",f107)
      if(nSTAT.ne.NF90_NOERR) &
          write(G2S_global_log,*)'ERROR: put_att F107: ',nf90_strerror(nSTAT)

        ! Add projection line from control file
      nSTAT = nf90_put_att(ncid,nf90_global,"proj",Comp_projection_line)
      if(nSTAT.ne.NF90_NOERR) &
          write(G2S_global_log,*)'ERROR: put_att Comment: ',nf90_strerror(nSTAT)
      if(IsLatLon_CompGrid.eqv..false.)then
        nSTAT = nf90_put_att(ncid,nf90_global,"xmin_g2s",xmin_g2s)
        if(nSTAT.ne.NF90_NOERR) &
            write(G2S_global_log,*)'ERROR: put_att xmin_g2s: ',nf90_strerror(nSTAT)
        nSTAT = nf90_put_att(ncid,nf90_global,"ymin_g2s",ymin_g2s)
        if(nSTAT.ne.NF90_NOERR) &
            write(G2S_global_log,*)'ERROR: put_att ymin_g2s: ',nf90_strerror(nSTAT)
      endif

      nSTAT = nf90_put_att(ncid,nf90_global,"dx_g2s",dx_g2s)
      if(nSTAT.ne.NF90_NOERR) &
          write(G2S_global_log,*)'ERROR: put_att dx_g2s: ',nf90_strerror(nSTAT)
      nSTAT = nf90_put_att(ncid,nf90_global,"dy_g2s",dy_g2s)
      if(nSTAT.ne.NF90_NOERR) &
          write(G2S_global_log,*)'ERROR: put_att dy_g2s: ',nf90_strerror(nSTAT)
      nSTAT = nf90_put_att(ncid,nf90_global,"Spec_MaxOrder",Spec_MaxOrder)
      if(nSTAT.ne.NF90_NOERR) &
          write(G2S_global_log,*)'ERROR: put_att Spec_MaxOrder: ',nf90_strerror(nSTAT)


      nSTAT = nf90_put_att(ncid,nf90_global,"Spline_Order",P_ord)
      if(nSTAT.ne.NF90_NOERR) &
          write(G2S_global_log,*)'ERROR: put_att P_ord: ',nf90_strerror(nSTAT)
      nSTAT = nf90_put_att(ncid,nf90_global,"NKnots",Nknots)
      if(nSTAT.ne.NF90_NOERR) &
          write(G2S_global_log,*)'ERROR: put_att NKnots: ',nf90_strerror(nSTAT)

      nSTAT = nf90_put_att(ncid,nf90_global,"nx_g2s",nxmax_g2s)
      if(nSTAT.ne.NF90_NOERR) &
          write(G2S_global_log,*)'ERROR: put_att nx_g2s: ',nf90_strerror(nSTAT)
      nSTAT = nf90_put_att(ncid,nf90_global,"ny_g2s",nymax_g2s)
      if(nSTAT.ne.NF90_NOERR) &
          write(G2S_global_log,*)'ERROR: put_att ny_g2s: ',nf90_strerror(nSTAT)

      nSTAT = nf90_put_att(ncid,nf90_global,"nx_Met1",nxmax_g2s)
      if(nSTAT.ne.NF90_NOERR) &
          write(G2S_global_log,*)'ERROR: put_att nx_Met1: ',nf90_strerror(nSTAT)
      nSTAT = nf90_put_att(ncid,nf90_global,"ny_Met1",nymax_g2s)
      if(nSTAT.ne.NF90_NOERR) &
          write(G2S_global_log,*)'ERROR: put_att ny_Met1: ',nf90_strerror(nSTAT)
      nSTAT = nf90_put_att(ncid,nf90_global,"nz_Met1",nzmax_Met1)
      if(nSTAT.ne.NF90_NOERR) &
          write(G2S_global_log,*)'ERROR: put_att nz_Met1: ',nf90_strerror(nSTAT)
      nSTAT = nf90_put_att(ncid,nf90_global,"zmin_Met1",zmin_Met1)
      if(nSTAT.ne.NF90_NOERR) &
          write(G2S_global_log,*)'ERROR: put_att zmin_Met1: ',nf90_strerror(nSTAT)
      nSTAT = nf90_put_att(ncid,nf90_global,"zmax_Met1",zmax_Met1)
      if(nSTAT.ne.NF90_NOERR) &
          write(G2S_global_log,*)'ERROR: put_att zmax_Met1: ',nf90_strerror(nSTAT)
      nSTAT = nf90_put_att(ncid,nf90_global,"dz_Met1",dz_Met1)
      if(nSTAT.ne.NF90_NOERR) &
          write(G2S_global_log,*)'ERROR: put_att dz_Met1: ',nf90_strerror(nSTAT)

      if(nwindfile_groups.eq.2)then
        nSTAT = nf90_put_att(ncid,nf90_global,"nx_Met2",nxmax_g2s)
        if(nSTAT.ne.NF90_NOERR) &
            write(G2S_global_log,*)'ERROR: put_att nx_Met2: ',nf90_strerror(nSTAT)
        nSTAT = nf90_put_att(ncid,nf90_global,"ny_Met2",nymax_g2s)
        if(nSTAT.ne.NF90_NOERR) &
            write(G2S_global_log,*)'ERROR: put_att ny_Met2: ',nf90_strerror(nSTAT)
        nSTAT = nf90_put_att(ncid,nf90_global,"nz_Met2",nzmax_Met2)
        if(nSTAT.ne.NF90_NOERR) &
            write(G2S_global_log,*)'ERROR: put_att nz_Met2: ',nf90_strerror(nSTAT)
        nSTAT = nf90_put_att(ncid,nf90_global,"zmin_Met2",zmin_Met2)
        if(nSTAT.ne.NF90_NOERR) &
            write(G2S_global_log,*)'ERROR: put_att zmin_Met2: ',nf90_strerror(nSTAT)
        nSTAT = nf90_put_att(ncid,nf90_global,"zmax_Met2",zmax_Met2)
        if(nSTAT.ne.NF90_NOERR) &
            write(G2S_global_log,*)'ERROR: put_att zmax_Met2: ',nf90_strerror(nSTAT)
        nSTAT = nf90_put_att(ncid,nf90_global,"dz_Met2",dz_Met2)
        if(nSTAT.ne.NF90_NOERR) &
            write(G2S_global_log,*)'ERROR: put_att dz_Met2: ',nf90_strerror(nSTAT)
      endif

      nSTAT = nf90_put_att(ncid,nf90_global,"nx_HWT",nxmax_g2s)
      if(nSTAT.ne.NF90_NOERR) &
          write(G2S_global_log,*)'ERROR: put_att nx_HWT: ',nf90_strerror(nSTAT)
      nSTAT = nf90_put_att(ncid,nf90_global,"ny_HWT",nymax_g2s)
      if(nSTAT.ne.NF90_NOERR) &
          write(G2S_global_log,*)'ERROR: put_att ny_HWT: ',nf90_strerror(nSTAT)
      nSTAT = nf90_put_att(ncid,nf90_global,"nz_HWT",nzmax_HWT)
      if(nSTAT.ne.NF90_NOERR) &
          write(G2S_global_log,*)'ERROR: put_att nz_HWT: ',nf90_strerror(nSTAT)
      nSTAT = nf90_put_att(ncid,nf90_global,"zmin_HWT",zmin_HWT)
      if(nSTAT.ne.NF90_NOERR) &
          write(G2S_global_log,*)'ERROR: put_att zmin: ',nf90_strerror(nSTAT)
      nSTAT = nf90_put_att(ncid,nf90_global,"zmax_HWT",zmax_HWT)
      if(nSTAT.ne.NF90_NOERR) &
          write(G2S_global_log,*)'ERROR: put_att zmax: ',nf90_strerror(nSTAT)
      nSTAT = nf90_put_att(ncid,nf90_global,"dz_HWT",dz_HWT)
      if(nSTAT.ne.NF90_NOERR) &
          write(G2S_global_log,*)'ERROR: put_att dz_HWT: ',nf90_strerror(nSTAT)


      nSTAT = nf90_def_dim(ncid,odim_names(1),2,r_dim_id)
      if(nSTAT.ne.NF90_NOERR) &
          write(G2S_global_log,*)'ERROR: ',nf90_strerror(nSTAT)
      nSTAT = nf90_def_dim(ncid,odim_names(2),Nknots+P_ord,z_dim_id)
      if(nSTAT.ne.NF90_NOERR) &
          write(G2S_global_log,*)'ERROR: ',nf90_strerror(nSTAT)
      !nSTAT = nf90_def_dim(ncid,odim_names(3),Spec_MaxOrder+1,y_dim_id)
      nSTAT = nf90_def_dim(ncid,odim_names(3),y_order,y_dim_id)
      if(nSTAT.ne.NF90_NOERR) &
          write(G2S_global_log,*)'ERROR: ',nf90_strerror(nSTAT)
      !nSTAT = nf90_def_dim(ncid,odim_names(4),Spec_MaxOrder+1,x_dim_id)
      nSTAT = nf90_def_dim(ncid,odim_names(4),x_order,x_dim_id)
      if(nSTAT.ne.NF90_NOERR) &
          write(G2S_global_log,*)'ERROR: ',nf90_strerror(nSTAT)
      
      nSTAT = nf90_def_var(ncid,"knot",nf90_double,  &
            (/z_dim_id/), &
              kn_var_id)
      nSTAT = nf90_put_att(ncid,kn_var_id,"long_name",&
                           "Height of B-spline knot point")
      if(nSTAT.ne.NF90_NOERR)write(G2S_global_log,*)'ERROR: put_att: ',nf90_strerror(nSTAT)
      nSTAT = nf90_put_att(ncid,kn_var_id,"units","km")
      if(nSTAT.ne.NF90_NOERR)write(G2S_global_log,*)'ERROR: put_att: ',nf90_strerror(nSTAT)

      nSTAT = nf90_def_var(ncid,"vx_sh",nf90_double,  &
            (/x_dim_id,y_dim_id,z_dim_id,r_dim_id/), &
              vx_var_id)
      nSTAT = nf90_put_att(ncid,vx_var_id,"long_name",&
                           "SH Coefficient for zonal velocity (Vx)")
      if(nSTAT.ne.NF90_NOERR)write(G2S_global_log,*)'ERROR: put_att: ',nf90_strerror(nSTAT)

      nSTAT = nf90_def_var(ncid,"vy_sh",nf90_double,  &
            (/x_dim_id,y_dim_id,z_dim_id,r_dim_id/), &
              vy_var_id)
      nSTAT = nf90_put_att(ncid,vy_var_id,"long_name",&
                           "SH Coefficient for meridional velocity (Vy)")
      if(nSTAT.ne.NF90_NOERR)write(G2S_global_log,*)'ERROR: put_att: ',nf90_strerror(nSTAT)

      nSTAT = nf90_def_var(ncid,"tp_sh",nf90_double,  &
            (/x_dim_id,y_dim_id,z_dim_id,r_dim_id/), &
              temp_var_id)
      nSTAT = nf90_put_att(ncid,temp_var_id,"long_name",&
                           "SH Coefficient for temperature")
      if(nSTAT.ne.NF90_NOERR)write(G2S_global_log,*)'ERROR: put_att: ',nf90_strerror(nSTAT)

        ! Leaving define mode.
      nSTAT = nf90_enddef(ncid)

      allocate(dum1d_out(Nknots+P_ord))
      !allocate(dum4d_out(Spec_MaxOrder+1, Spec_MaxOrder+1,Nknots+P_ord,2))
      allocate(dum4d_out(x_order, y_order,Nknots+P_ord,2))

      dum1d_out(1:Nknots+P_ord) = real(knots(1:Nknots+P_ord),kind=8)
      nSTAT=nf90_put_var(ncid,kn_var_id,dum1d_out,(/1/))
      if(nSTAT.ne.NF90_NOERR) &
        write(G2S_global_log,*)'ERROR: put_var knot: ',nf90_strerror(nSTAT)

        ! Vx
      if (IsLatLon_CompGrid.eqv..true.)then
        dum4d_out(:,:,:,1) = vx3_SH_dp(1,1:x_order,1:y_order,1:Nknots+P_ord)
        dum4d_out(:,:,:,2) = vx3_SH_dp(2,1:x_order,1:y_order,1:Nknots+P_ord)
      else
        dum4d_out(:,:,:,1) =  real(vx3_FS_C(1:x_order,1:y_order,1:Nknots+P_ord))
        dum4d_out(:,:,:,2) = aimag(vx3_FS_C(1:x_order,1:y_order,1:Nknots+P_ord))
      endif
      nSTAT=nf90_put_var(ncid,vx_var_id,dum4d_out,(/1,1,1,1/))
      if(nSTAT.ne.NF90_NOERR) &
        write(G2S_global_log,*)'ERROR: put_var Vx: ',nf90_strerror(nSTAT)

        ! Vy
      if (IsLatLon_CompGrid.eqv..true.)then
        dum4d_out(:,:,:,1) = vy3_SH_dp(1,1:x_order,1:y_order,1:Nknots+P_ord)
        dum4d_out(:,:,:,2) = vy3_SH_dp(2,1:x_order,1:y_order,1:Nknots+P_ord)
      else
        dum4d_out(:,:,:,1) =  real(vy3_FS_C(1:x_order,1:y_order,1:Nknots+P_ord))
        dum4d_out(:,:,:,2) = aimag(vy3_FS_C(1:x_order,1:y_order,1:Nknots+P_ord))
      endif
      !dum4d_out(:,:,:,1) = vy3_SH_dp(1,1:Spec_MaxOrder+1,1:Spec_MaxOrder+1,1:Nknots+P_ord)
      !dum4d_out(:,:,:,2) = vy3_SH_dp(2,1:Spec_MaxOrder+1,1:Spec_MaxOrder+1,1:Nknots+P_ord)
      nSTAT=nf90_put_var(ncid,vy_var_id,dum4d_out,(/1,1,1,1/))
      if(nSTAT.ne.NF90_NOERR) &
        write(G2S_global_log,*)'ERROR: put_var Vy: ',nf90_strerror(nSTAT)

        ! Temp
      if (IsLatLon_CompGrid.eqv..true.)then
        dum4d_out(:,:,:,1) = temperature3_SH_dp(1,1:x_order,1:y_order,1:Nknots+P_ord)
        dum4d_out(:,:,:,2) = temperature3_SH_dp(2,1:x_order,1:y_order,1:Nknots+P_ord)
      else
        dum4d_out(:,:,:,1) =  real(temperature3_FS_C(1:x_order,1:y_order,1:Nknots+P_ord))
        dum4d_out(:,:,:,2) = aimag(temperature3_FS_C(1:x_order,1:y_order,1:Nknots+P_ord))
      endif
      !dum4d_out(:,:,:,1) = temperature3_SH_dp(1,1:Spec_MaxOrder+1,1:Spec_MaxOrder+1,1:Nknots+P_ord)
      !dum4d_out(:,:,:,2) = temperature3_SH_dp(2,1:Spec_MaxOrder+1,1:Spec_MaxOrder+1,1:Nknots+P_ord)
      nSTAT=nf90_put_var(ncid,temp_var_id,dum4d_out,(/1,1,1,1/))
      if(nSTAT.ne.NF90_NOERR) &
        write(G2S_global_log,*)'ERROR: put_var Temp: ',nf90_strerror(nSTAT)

      ! Close file
      nSTAT = nf90_close(ncid)

      end subroutine Write_Spec_nc

!##############################################################################
!##############################################################################
!
!     Get_WindTemp_Empir
!
!     This subroutine is the interface to the empirical models HWM07 and
!     NRLMSIS.  This takes as arguments, the full coordinates (x,y,z,t) as well
!     as the space-weather indicies Ap and F107, then returns the Vx, Vy, and T
!     values.
!
!     No global variables are filled.  Vx,Vy, and T are returned through the
!     argument list.
!
!##############################################################################

      subroutine Get_WindTemp_Empir(lon,lat,alt,              &
                              year,day,ihour,iminute,isecond, &
                              ap,f107,                        &
                              u_g2s,v_g2s,temperature)

      ! Here are the NRLMSISE-00 modules
      use utils_constants
      use physics_msis

#ifdef useHWM07
      use NEWmodel
#endif
      implicit none

      integer,intent(in)       :: day,year
      integer,intent(in)       :: ihour,iminute,isecond
      real(kind=4),intent(in)  :: lat,lon,alt,ap,f107
      real(kind=4),intent(out) :: u_g2s,v_g2s,temperature

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
      !real(4) :: dw_sp(2)          ! 1: MERIDIONAL DISTURBANCE WIND (m/sec + Geo. Northward)
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
      real(dp) :: f107a_dp         ! 81 day AVERAGE OF F10.7 FLUX (centered on day DDD)
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

      call gtd7(day,ut_dp,alt_dp,xlat_dp,xlon_dp,xlst_dp,f107a_dp,f107_dp,aph_dp,48,d_dp,t_dp)

      u_g2s       = vx_new
      v_g2s       = vy_new
      temperature = real(t_dp(2),kind=4)

      end subroutine Get_WindTemp_Empir

!##############################################################################
!##############################################################################
!
!     day_of_year
!
!     This function takes the year, month and day as integers and returns the
!     integer day-of-year.
!
!##############################################################################

      function day_of_year(iyear,imonth,iday)

      implicit none

      !real(kind=4) :: day_of_year
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

