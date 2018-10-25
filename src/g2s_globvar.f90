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
!  G2S_globvar
!
!  Module containing the global variables used by AVOG2S
!
!##############################################################################

      module G2S_globvar

      implicit none

      real(kind=4), parameter :: DEG2RAD      = 1.7453292519943295e-2_4
      real(kind=4), parameter :: RAD_EARTH    = 6371.229_4
      real(kind=4), parameter :: R_GAS_DRYAIR = 286.98_4    ! Specific gas constant of R=286.98 J /(kg K)

      integer, parameter ::     P_ord   = 3    ! This is the order of the B-Spline
!      integer, parameter ::     Nknots  = 28   ! Number of knots in the B-Splines
      integer, parameter ::     Nknots  = 37   ! Number of knots in the B-Splines

        ! These are useful for writing out cross-sections of the Met1/2, HWT, and reconctructed data
      logical, parameter ::     write_test_output = .false.
      integer, parameter ::     jout              = 180   ! This is for a equatorial cross-section
      integer, parameter ::     iout              = 366   ! This is x=0 for the subgrid of the nam198 used in Fig.4
                                                          !    x ~= 0 in the projection (lon=210)
      !integer, parameter ::     iout              = 421   ! This is lat=210 in the GFS 0.5-degree


      integer        :: G2S_global_essential    = 6
      integer        :: G2S_global_production   = 6
      integer        :: G2S_global_debug        = 6
      integer        :: G2S_global_info         = 6
      integer        :: G2S_global_log          = 9
      integer        :: G2S_global_error        = 0

      ! These are the variables specifying the restructured atmosphere
      integer      :: nz1
      integer      :: nz2
      integer      :: nz3
      real(kind=4) :: dz1
      real(kind=4) :: dz2
      real(kind=4) :: dz3

      real(kind=4) :: dx_g2s,dy_g2s        ! Resolution of the g2s grid (deg or km)
      integer      :: Spec_MaxOrder        ! Order of the spectral decomposition (for SH or Fourier)
      integer      :: x_order
      integer      :: y_order
      character(len=80)  :: Comp_projection_line

      integer :: nwindfile_groups = 1      ! number of MET datasets to use (1 or 2)
      integer :: nwindfiles1  = 1          ! number of windfiles in group 1
      integer :: nwindfiles2  = 1          ! number of windfiles in group 2
      character(len=130)  :: controlfile   ! name of control file to read.
      character(len=130),allocatable,dimension(:)  :: windfile1     ! name of 3D wind file to read (low-altitude)
      character(len=130),allocatable,dimension(:)  :: windfile2     ! name of 3D wind file to read (mid-altitude)
      character(len=130)  :: apfile        ! name of ap file to read.
      character(len=130)  :: f107file      ! name of ap file to read.

      integer :: ct_unit   = 14
      integer :: ap_unit   = 15
      integer :: f107_unit = 16

      integer :: iwf1,iwf2

      real(kind=4),dimension(:,:,:)     ,allocatable :: vx_Met_loc_sp
      real(kind=4),dimension(:,:,:)     ,allocatable :: vy_Met_loc_sp
      real(kind=4),dimension(:,:,:)     ,allocatable :: temperature_Met_loc_sp
      real(kind=4),dimension(:,:,:)     ,allocatable :: tmp3d_1_sp
      real(kind=4),dimension(:,:,:)     ,allocatable :: tmp3d_2_sp
      real(kind=4),dimension(:,:,:)     ,allocatable :: tmp3d_1_2_sp
      real(kind=4),dimension(:,:,:)     ,allocatable :: tmp3d_2_2_sp

      real(kind=4),dimension(:,:,:)   ,allocatable :: vx_HWT_sp
      real(kind=4),dimension(:,:,:)   ,allocatable :: vy_HWT_sp
      real(kind=4),dimension(:,:,:)   ,allocatable :: temperature_HWT_sp

      ! Holds the Met Spherical Harmonic decomposition
      real(kind=8),dimension(:,:,:,:)   ,allocatable :: vx1_SH_dp
      real(kind=8),dimension(:,:,:,:)   ,allocatable :: vy1_SH_dp
      real(kind=8),dimension(:,:,:,:)   ,allocatable :: temperature1_SH_dp
      ! Holds the Upper-atmos Met Spherical Harmonic decomposition
      real(kind=8),dimension(:,:,:,:)   ,allocatable :: vx1_2_SH_dp
      real(kind=8),dimension(:,:,:,:)   ,allocatable :: vy1_2_SH_dp
      real(kind=8),dimension(:,:,:,:)   ,allocatable :: temperature1_2_SH_dp
      ! Holds the HWT Spherical Harmonic decomposition
      real(kind=8),dimension(:,:,:,:)   ,allocatable :: vx2_SH_dp
      real(kind=8),dimension(:,:,:,:)   ,allocatable :: vy2_SH_dp
      real(kind=8),dimension(:,:,:,:)   ,allocatable :: temperature2_SH_dp
      ! Holds the Interpolated Spherical Harmonic decomposition
      real(kind=8),dimension(:,:,:,:)   ,allocatable :: vx3_SH_dp
      real(kind=8),dimension(:,:,:,:)   ,allocatable :: vy3_SH_dp
      real(kind=8),dimension(:,:,:,:)   ,allocatable :: temperature3_SH_dp

      ! This is used to flag projected (Fourier) vs global (spherical harmonic)
      logical :: IsLatLon
      ! Holds the Met Fourier decomposition
      complex(kind=8),dimension(:,:,:)   ,allocatable :: vx1_FS_C
      complex(kind=8),dimension(:,:,:)   ,allocatable :: vy1_FS_C
      complex(kind=8),dimension(:,:,:)   ,allocatable :: temperature1_FS_C
      ! Holds the Upper-atmos Met Fourier decomposition
      complex(kind=8),dimension(:,:,:)   ,allocatable :: vx1_2_FS_C
      complex(kind=8),dimension(:,:,:)   ,allocatable :: vy1_2_FS_C
      complex(kind=8),dimension(:,:,:)   ,allocatable :: temperature1_2_FS_C
      ! Holds the HWT Fourier decomposition
      complex(kind=8),dimension(:,:,:)   ,allocatable :: vx2_FS_C
      complex(kind=8),dimension(:,:,:)   ,allocatable :: vy2_FS_C
      complex(kind=8),dimension(:,:,:)   ,allocatable :: temperature2_FS_C
      ! Holds the Interpolated Fourier decomposition
      complex(kind=8),dimension(:,:,:)   ,allocatable :: vx3_FS_C
      complex(kind=8),dimension(:,:,:)   ,allocatable :: vy3_FS_C
      complex(kind=8),dimension(:,:,:)   ,allocatable :: temperature3_FS_C

      real(kind=8),dimension(:,:),allocatable :: Lon_of_proj_node_dp
      real(kind=8),dimension(:,:),allocatable :: Lat_of_proj_node_dp
      real(kind=8),dimension(:,:),allocatable :: Theta_of_proj_node_dp

      real(kind=8),dimension(:,:,:),allocatable :: vx_OUT_dp
      real(kind=8),dimension(:,:,:),allocatable :: vy_OUt_dp
      real(kind=8),dimension(:,:,:),allocatable :: temperature_OUT_dp

      integer :: start_year,start_month,start_day
      integer :: start_hour,start_min,start_sec

      integer :: inyear      = 0
      integer :: inmonth     = 0
      integer :: inday       = 0
      real(kind=8) :: inhour = 0.0

      real(kind=4),dimension(Nknots)         :: iknots
      real(kind=4),dimension(Nknots+2*P_ord) :: knots

      real(kind=4),dimension(:)     ,allocatable :: data_alt
      real(kind=4),dimension(:)     ,allocatable :: new_data_alt
      real(kind=4),dimension(:,:)   ,allocatable :: coeff

      integer :: full_data_len

      ! B-spline knot (28) points suitable for thermosphere (similare to NRL-G2S)
!      data iknots/0.0,  5.0, 10.0, 15.0, 20.0, &
!                 25.0, 30.0, 35.0, 40.0, 45.0, &
!                 50.0, 55.0, 60.0, 65.0, 70.0, &
!                 75.0, 80.0, 85.0, 90.0, 95.0, &
!                100.0, 105.0, 110.0, 117.5, 125.0, &
!                135.0, 150.0, 225.0/

      ! B-spline knot points (37) suitable for thermosphere with higher resolution in
      ! the troposphere
      data iknots/0.0,  2.5, 5.0, 7.5, 10.0, 12.5, 15.0, 17.5, 20.0, &
                 22.5, 25.0, 27.5, 30.0, 32.5, 35.0, 37.5, 40.0, 42.5, 45.0, &
                 50.0, 55.0, 60.0, 65.0, 70.0, &
                 75.0, 80.0, 85.0, 90.0, 95.0, &
                100.0, 105.0, 110.0, 117.5, 125.0, &
                135.0, 150.0, 225.0/

      ! B-spline knot points (37) suitable for exosphere
!      data iknots/0.0,  5.0, 10.0, 15.0, 20.0, &
!                 25.0, 30.0, 35.0, 40.0, 45.0, &
!                 50.0, 55.0, 60.0, 65.0, 70.0, &
!                 75.0, 80.0, 85.0, 90.0, 95.0, &
!                100.0, 105.0, 110.0, 117.5, 125.0, &
!                135.0, 150.0, 175.0, 225.0, 275.0, &
!                300.0, 375.0, 450.0, 525.0, 600.0, &
!                675.0, 750.0/


      integer :: fc_hour = 0
      real(kind=4) :: ap,f107

      logical :: LOAD_HWT    = .false.

      character(len=50) :: file_SH    = "out20_GS.nc"
      character(len=50) :: file_T_HWT = "temper_HWT.raw"
      character(len=50) :: file_U_HWT = "vx_HWT.raw"
      character(len=50) :: file_V_HWT = "vy_HWT.raw"
      character(len=50) :: file_T_RES = "T_res.raw"
      character(len=50) :: file_U_RES = "U_res.raw"
      character(len=50) :: file_V_RES = "V_res.raw"
      character(len=50) :: file_TOPO  = "ETOPO1_Ice_c_gmt4.nc"
      character(len=50) :: file_OUT_ROOT  = "InfraAtmos"
      character(len=50) :: file_OUT_SONDE = "InfraAtmos01.met"
      character(len=50) :: file_OUT_XSEC  = "InfraAtmos01.env"
      character(len=50) :: file_OUT_GRID  = "Volc"
      character(len=50) :: file_OUT
      character(len=50) :: file_LL        = "ProfLonLat01.dat"

      integer :: data_len

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Variables only used in topography for extract
      logical      :: useTopo
      real(kind=8) :: dlat_topo
      real(kind=8) :: dlon_topo
      real(kind=8)   , dimension(:)   ,allocatable :: lat_raw
      real(kind=8)   , dimension(:)   ,allocatable :: lon_raw
      integer(kind=2), dimension(:,:) ,allocatable :: topo_raw
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      logical :: IsThere     ! variable used for checking existance of files.

      integer      :: nxmax_g2s ! number of x points in the g2s grid
      integer      :: nymax_g2s ! number of y points in the g2s grid
      real(kind=4) :: xmin_g2s,ymin_g2s
      integer      :: nzmax_Met1
      real(kind=4) :: zmin_Met1, zmax_Met1, dz_Met1
      integer      :: nzmax_Met2
      real(kind=4) :: zmin_Met2, zmax_Met2, dz_Met2
      integer      :: nzmax_HWT
      real(kind=4) :: zmin_HWT,  zmax_HWT,  dz_HWT
      real(kind=4),dimension(:)      ,allocatable :: x_g2s_sp    ! G2S grid
      real(kind=4),dimension(:)      ,allocatable :: y_g2s_sp    ! G2S grid
      real(kind=4),dimension(:)      ,allocatable :: z_Met1_sp ! Met1 is resampled onto this
      real(kind=4),dimension(:)      ,allocatable :: z_Met2_sp ! Met2 ''
      real(kind=4),dimension(:)      ,allocatable :: z_HWT_g2s_sp  ! z of HWT grid
      real(kind=8),dimension(:)      ,allocatable :: x_g2s_dp    ! G2S grid
      real(kind=8),dimension(:)      ,allocatable :: y_g2s_dp    ! G2S grid

      integer FFTW_MEASURE
      parameter (FFTW_MEASURE=0)
      real(kind=8) :: ForecastInterval

      end module G2S_globvar

!##############################################################################
!##############################################################################

      module G2S_NCvars

      implicit none

      integer :: z_dim_id
      integer :: y_dim_id
      integer :: x_dim_id
      integer :: r_dim_id

      integer :: vx_var_id      = 0
      integer :: vy_var_id      = 0
      integer :: temp_var_id    = 0

      integer :: kn_var_id      = 0

      end module G2S_NCvars


