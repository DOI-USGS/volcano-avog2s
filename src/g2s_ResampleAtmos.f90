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
!  G2S_ResampleAtmos
!
!    This program reads a coeffiecient file and generates three gridded binary
!    files for U, V, and T that can be used the the g2s_Extract_[Sonde,Xsec,Grid]
!    programs.  Vertical spacing of grid can either be hardwired below (must be
!    consistent with G2S_Extract.f90), or given on the command-line.
!
!##############################################################################

      PROGRAM G2S_ResampleAtmos

      use G2S_globvar

      implicit none

      integer :: nargs
      character(len=50) :: llinebuffer
      integer :: root_len 


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
       !===================================================


!     TEST READ COMMAND LINE ARGUMENTS
      nargs = iargc()
      if (nargs.ge.1) then
        call getarg(1,llinebuffer)
        FILE_SH = trim(adjustl(llinebuffer))
        ! If a filename is provided, strip off the ".nc" and append
        ! the labels for TUV
        root_len = len(trim(FILE_SH)) - 3
        FILE_T_RES = FILE_SH(1:root_len) // "_T_res.raw"
        FILE_U_RES = FILE_SH(1:root_len) // "_U_res.raw"
        FILE_V_RES = FILE_SH(1:root_len) // "_V_res.raw"
      else
        write(*,*)"No command-line arguments given."
        write(*,*)" Normal usage is to provide the coefficient file as an argument."
        write(*,*)"Assuming default SH file name and"
        write(*,*)"using default Res file names."
      endif
      if (nargs.ge.2) then
        ! Additional command-line arguments are provided.  Assume these are
        ! nz1 dz1 nz2 dz2 nz3 dz3
        if (nargs.lt.7)then
          write(*,*)"Multiple command-line argument provided, but not enough to"
          write(*,*)"specify the grid."
          write(*,*)"Expected format:"
          write(*,*)"g2s_ResampleAtmos Coeffic_File.nc nz1 dz1 nz2 dz2 nz3 dz3"
          stop 1
        else
          ! Read the next six arguments
          call getarg(2,llinebuffer)
          read(llinebuffer,*)nz1
          call getarg(3,llinebuffer)
          read(llinebuffer,*)dz1
          call getarg(4,llinebuffer)
          read(llinebuffer,*)nz2
          call getarg(5,llinebuffer)
          read(llinebuffer,*)dz2
          call getarg(6,llinebuffer)
          read(llinebuffer,*)nz3
          call getarg(7,llinebuffer)
          read(llinebuffer,*)dz3
        endif
      else
        ! Use the default vertical grid
        write(*,*)"Using the default values for vertical grid."
        write(*,*)"nz1 = ",nz1
        write(*,*)"dz1 = ",dz1
        write(*,*)"nz2 = ",nz2
        write(*,*)"dz2 = ",dz2
        write(*,*)"nz3 = ",nz3
        write(*,*)"dz3 = ",dz3
      endif
      write(*,*)"Opening SH coefficient file: ",FILE_SH
      call Read_SH_nc()
 
      write(*,*)&
        "Recreating Atmos based on spectra at native lat/lon grid."
      write(*,*)"Number of levels = ",nz1+nz2+nz3
      call Resample_Atmos

      end program G2S_ResampleAtmos

!##############################################################################
!##############################################################################
!
!     Read_SH_nc
!
!     This subroutine reads the merged SH coefficients file generated by g2s
!
!     Global variables filled are:
!       vx3_SH_dp(2, Spec_MaxOrder+1, Spec_MaxOrder+1,model_len))
!       vy3_SH_dp(2, Spec_MaxOrder+1, Spec_MaxOrder+1,model_len))
!       temperature3_SH_dp(2, Spec_MaxOrder+1, Spec_MaxOrder+1,model_len))
!
!##############################################################################

      subroutine Read_SH_nc()

      use G2S_globvar
      use G2S_NCvars
      use netcdf

      implicit none

      integer :: nSTAT
      integer :: ncid

      real(kind=8),dimension(:)      ,allocatable :: dum1d_out
      real(kind=8),dimension(:,:,:,:),allocatable :: dum4d_out
      character(len=3),dimension(5) :: odim_names
      integer :: model_len
      integer :: i,j,k
      integer :: ilatlonflag
      integer :: iprojflag

      odim_names(1) = "r"  ! real/imag
      odim_names(2) = "z"  ! alt
      odim_names(3) = "y"  ! order in y
      odim_names(4) = "x"  ! order in x
      odim_names(5) = "k"  ! knot

      INQUIRE( FILE=adjustl(trim(FILE_SH)), EXIST=IsThere )
      IF(.not.IsThere)THEN
        write(*,*)"ERROR: Could not file input file."
        write(*,*)adjustl(trim(FILE_SH))
        stop 1
      ENDIF

      ! Open existing netcdf file
      nSTAT = nf90_open(FILE_SH,NF90_NOWRITE,ncid)
      IF(nSTAT.ne.0)write(6,*)'ERROR: nf90_open to read header:', &
                           nf90_strerror(nSTAT)
      IF(nSTAT.ne.0)THEN
        write(6,*)'Could not open ',FILE_SH
        write(6,*)'Exiting'
        stop 1
      ENDIF

      nSTAT = nf90_get_att(ncid,nf90_global,"nx_g2s",nxmax_g2s)
      nSTAT = nf90_get_att(ncid,nf90_global,"ny_g2s",nymax_g2s)
      nSTAT = nf90_get_att(ncid,nf90_global,"dx_g2s",dx_g2s)
      IF(nSTAT.ne.0)THEN
        write(*,*)"WARNING: Attribute dx_g2s not found in Spectral Coefficient file."
        write(*,*)"         Assuming 0.5"
        dx_g2s = 0.5
      ENDIF
      nSTAT = nf90_get_att(ncid,nf90_global,"dy_g2s",dy_g2s)
      IF(nSTAT.ne.0)THEN
        write(*,*)"WARNING: Attribute dy_g2s not found in Spectral Coefficient file."
        write(*,*)"         Assuming 0.5"
        dy_g2s = 0.5
      ENDIF
      nSTAT = nf90_get_att(ncid,nf90_global,"Spec_MaxOrder",Spec_MaxOrder)
      IF(nSTAT.ne.0)THEN
        write(*,*)&
          "WARNING: Attribute Spec_MaxOrder not found in Spectral Coefficient file."
        write(*,*)"         Assuming 120"
        Spec_MaxOrder = 120
      ENDIF
      nSTAT = nf90_get_att(ncid,nf90_global,"proj",Comp_projection_line)
      IF(nSTAT.ne.0)THEN
        write(*,*)&
          "WARNING: Attribute proj not found in Spectral Coefficient file."
        write(*,*)"         Assuming Lat/Lon grid"
        Comp_projection_line = "1 4 -107.0 50.0 50.0 50.0 6367.470"
      ENDIF
      read(Comp_projection_line,*)ilatlonflag,iprojflag
      if (ilatlonflag.eq.0) then
        ! expecting input variables to be in the same projection as
        ! specified by iprojflag and parameters
        IsLatLon          = .false.
        x_order = nxmax_g2s/2+1
        y_order = nymax_g2s
        nSTAT = nf90_get_att(ncid,nf90_global,"xmin_g2s",xmin_g2s)
        IF(nSTAT.ne.0)THEN
          write(*,*)"ERROR: grid is projected, but no xmin_g2s in SH file"
          stop
        ENDIF
        nSTAT = nf90_get_att(ncid,nf90_global,"ymin_g2s",ymin_g2s)
        IF(nSTAT.ne.0)THEN
          write(*,*)"ERROR: grid is projected, but no ymin_g2s in SH file"
          stop
        ENDIF
      else
        ! expecting input variables to be in lat/lon
        x_order = Spec_MaxOrder+1
        y_order = Spec_MaxOrder+1
        IsLatLon          = .true.
      endif

      allocate(x_g2s_sp(nxmax_g2s))
      allocate(y_g2s_sp(nymax_g2s))

      IF(IsLatLon)THEN
        DO i=1,nxmax_g2s
          x_g2s_sp(i) = (i-1)*real(dx_g2s,kind=4)
        ENDDO
        DO j=1,nymax_g2s
          y_g2s_sp(j) = -90.0_4+(j-1)*real(dx_g2s,kind=4)
        ENDDO
      ELSE
        DO i=1,nxmax_g2s
          x_g2s_sp(i) = xmin_g2s + (i-1)*real(dx_g2s,kind=4)
        ENDDO
        DO j=1,nymax_g2s
          y_g2s_sp(j) = ymin_g2s + (j-1)*real(dy_g2s,kind=4)
        ENDDO
      ENDIF
      !nSTAT = nf90_put_att(ncid,nf90_global,"Spline_Order",P_ord)
      !P_ord = 3
      
      nSTAT = nf90_inq_dimid(ncid,odim_names(1),r_dim_id)
       ! r_dim length = 2
      nSTAT = nf90_inq_dimid(ncid,odim_names(2),z_dim_id)
!      nSTAT = nf90_Inquire_Dimension(ncid,z_dim_id,len=Nknots)
!      Nknots = Nknots - P_ord
      nSTAT = nf90_inq_dimid(ncid,odim_names(3),y_dim_id)
      !nSTAT = nf90_Inquire_Dimension(ncid,y_dim_id,len=Spec_MaxOrder)
      !Spec_MaxOrder = Spec_MaxOrder - 1
      nSTAT = nf90_inq_dimid(ncid,odim_names(4),x_dim_id)
        ! Assume x length is the same as y length

      !nSTAT = nf90_def_var(ncid,"knot",nf90_double,  &
      !      (/z_dim_id/), &
      !        kn_var_id)

      model_len = Nknots+P_ord
      allocate(dum1d_out(model_len))
      allocate(dum4d_out(x_order, y_order,Nknots+P_ord,2))
      IF(IsLatLon)THEN
        allocate(         vx3_SH_dp(2, x_order, y_order,model_len))
        allocate(         vy3_SH_dp(2, x_order, y_order,model_len))
        allocate(temperature3_SH_dp(2, x_order, y_order,model_len))
      ELSE
        allocate(         vx3_FS_C(x_order,y_order,model_len))
        allocate(         vy3_FS_C(x_order,y_order,model_len))
        allocate(temperature3_FS_C(x_order,y_order,model_len))
      ENDIF

      nSTAT = nf90_inq_varid(ncid,"knot",kn_var_id)
      nSTAT = nf90_get_var(ncid,kn_var_id,dum1d_out, &
             start = (/1/),       &
             count = (/model_len/))
      knots(1:model_len) = real(dum1d_out(1:model_len),kind=4)
      knots(model_len:Nknots+2*P_ord) = knots(model_len)
      iknots(1:Nknots) = real(dum1d_out(P_ord+1:model_len),kind=4)

      nSTAT = nf90_inq_varid(ncid,"vx_sh",vx_var_id)
      IF(nSTAT.ne.0) &
         write(6,*)'ERROR: inq_varid:Vx ',nf90_strerror(nSTAT)
      dum4d_out = 0.0
      nSTAT = nf90_get_var(ncid,vx_var_id,dum4d_out, &
             start = (/1,1,1,1/),       &
             count = (/x_order, y_order,Nknots+P_ord,2/))
      IF(nSTAT.ne.0) &
         write(6,*)'ERROR: get_var:Vx ',nf90_strerror(nSTAT)
      IF(IsLatLon)THEN
        vx3_SH_dp(1,:,:,:) = dum4d_out(:,:,:,1)
        vx3_SH_dp(2,:,:,:) = dum4d_out(:,:,:,2)
      ELSE
        DO i=1,x_order
          DO j=1,y_order
            DO k=1,model_len
              vx3_FS_C(i,j,k) = complex(dum4d_out(i,j,k,1),dum4d_out(i,j,k,2))
            ENDDO
          ENDDO
        ENDDO
      ENDIF

      nSTAT = nf90_inq_varid(ncid,"vy_sh",vy_var_id)
      IF(nSTAT.ne.0) &
         write(6,*)'ERROR: inq_varid:Vy ',nf90_strerror(nSTAT)
      dum4d_out = 0.0
      nSTAT = nf90_get_var(ncid,vy_var_id,dum4d_out, &
             start = (/1,1,1,1/),       &
             count = (/x_order, y_order,Nknots+P_ord,2/))
      IF(nSTAT.ne.0) &
         write(6,*)'ERROR: get_var:Vy ',nf90_strerror(nSTAT)
      IF(IsLatLon)THEN
        vy3_SH_dp(1,:,:,:) = dum4d_out(:,:,:,1)
        vy3_SH_dp(2,:,:,:) = dum4d_out(:,:,:,2)
      ELSE
        DO i=1,x_order
          DO j=1,y_order
            DO k=1,model_len
              vy3_FS_C(i,j,k) = complex(dum4d_out(i,j,k,1),dum4d_out(i,j,k,2))
            ENDDO
          ENDDO
        ENDDO
      ENDIF

      nSTAT = nf90_inq_varid(ncid,"tp_sh",temp_var_id)
      IF(nSTAT.ne.0) &
         write(6,*)'ERROR: inq_varid:Temperature ',nf90_strerror(nSTAT)
      dum4d_out = 0.0
      nSTAT = nf90_get_var(ncid,temp_var_id,dum4d_out, &
             start = (/1,1,1,1/),       &
             count = (/x_order, y_order,Nknots+P_ord,2/))
      IF(nSTAT.ne.0) &
         write(6,*)'ERROR: get_var:Temperature ',nf90_strerror(nSTAT)
      IF(IsLatLon)THEN
        temperature3_SH_dp(1,:,:,:) = dum4d_out(:,:,:,1)
        temperature3_SH_dp(2,:,:,:) = dum4d_out(:,:,:,2)
      ELSE
        DO i=1,x_order
          DO j=1,y_order
            DO k=1,model_len
              temperature3_FS_C(i,j,k) = complex(dum4d_out(i,j,k,1),dum4d_out(i,j,k,2))
            ENDDO
          ENDDO
        ENDDO
      ENDIF
      ! Close file
      nSTAT = nf90_close(ncid)

      end subroutine Read_SH_nc

!##############################################################################
!##############################################################################
!
!     Resample_Atmos
!
!     This subroutine reads the merged SH coefficients file generated by g2s
!
!     Global variables filled are:
!       vx3_SH_dp(2, Spec_MaxOrder+1, Spec_MaxOrder+1,model_len))
!       vy3_SH_dp(2, Spec_MaxOrder+1, Spec_MaxOrder+1,model_len))
!       temperature3_SH_dp(2, Spec_MaxOrder+1, Spec_MaxOrder+1,model_len))
!
!##############################################################################

      subroutine Resample_Atmos

      use G2S_globvar
      use SHTOOLS

      implicit none

      real(kind=8) :: cilm( 2, x_order, y_order)
      real(kind=8) :: grid(  nxmax_g2s+1, nxmax_g2s+1)
      real(kind=8) :: interval
      integer :: i, j, k, ioff, ii,iii,nlat,nlong
      integer :: Mfull,model_len,lmax
      real(kind=8) :: alt

      real(kind=8),dimension(:,:,:),allocatable :: mm
      real(kind=8),dimension(:),allocatable :: model_r,model_i

      real(kind=8)   ,dimension(:,:),allocatable :: gridU_in
      real(kind=8)   ,dimension(:,:),allocatable :: gridV_in
      real(kind=8)   ,dimension(:,:),allocatable :: gridT_in
      complex(kind=8),dimension(:,:),allocatable :: gridU_out
      complex(kind=8),dimension(:,:),allocatable :: gridV_out
      complex(kind=8),dimension(:,:),allocatable :: gridT_out
      integer(kind=8)                            :: plan_bak2Dfft_U
      integer(kind=8)                            :: plan_bak2Dfft_V
      integer(kind=8)                            :: plan_bak2Dfft_T

      interval = real(dx_g2s,kind=8)
      data_len  = nz1 + nz2 + nz3
      Mfull     = Nknots+P_ord+P_ord
      model_len = Nknots+P_ord

      IF(.not.IsLatLon)THEN
        allocate(gridU_in( nxmax_g2s      , nymax_g2s))
        allocate(gridV_in( nxmax_g2s      , nymax_g2s))
        allocate(gridT_in( nxmax_g2s      , nymax_g2s))
        allocate(gridU_out(nxmax_g2s/2 + 1, nymax_g2s))
        allocate(gridV_out(nxmax_g2s/2 + 1, nymax_g2s))
        allocate(gridT_out(nxmax_g2s/2 + 1, nymax_g2s))
        call dfftw_plan_dft_c2r_2d(plan_bak2Dfft_U, nxmax_g2s, nymax_g2s, &
                                   gridU_out , gridU_in,                  &
                                   FFTW_MEASURE)
        call dfftw_plan_dft_c2r_2d(plan_bak2Dfft_V, nxmax_g2s, nymax_g2s, &
                                   gridV_out , gridV_in,                  &
                                   FFTW_MEASURE)
        call dfftw_plan_dft_c2r_2d(plan_bak2Dfft_T, nxmax_g2s, nymax_g2s, &
                                   gridT_out , gridT_in,                  &
                                   FFTW_MEASURE)
      ENDIF

      allocate(vx_OUT_dp(nxmax_g2s,nymax_g2s,data_len))
      allocate(vy_OUT_dp(nxmax_g2s,nymax_g2s,data_len))
      allocate(temperature_OUT_dp(nxmax_g2s,nymax_g2s,data_len))

      ! Set up for recreating B-Splines
      allocate(mm(Mfull,P_ord+1,data_len))
      allocate(new_data_alt(data_len))

      DO k = 1,nz1
        new_data_alt(k) = (k-1)*dz1
      ENDDO
      DO k = 1,nz2
        new_data_alt(k+nz1) = new_data_alt(nz1) + (k)*dz2
      ENDDO
      DO k = 1,nz3
        new_data_alt(k+nz1+nz2) = new_data_alt(nz1+nz2) + (k)*dz3
      ENDDO

      ! Define knot vector
      DO i = 1,Mfull
        If(i.le.P_ord)THEN
          knots(i) = iknots(1)
        ELSEIF(i.ge.Nknots+P_ord)THEN
          knots(i) = iknots(Nknots)
        ELSE
          knots(i) = iknots(i-P_ord)
        ENDIF
      ENDDO

      ! Now get B-Spline basis for the HWT data
      !% Generate coefficient matrix of B-Splines
      ioff = 1
      mm = 0.0
      DO k = 1,data_len
        alt = new_data_alt(k)
        DO j=1,Mfull-1
          IF (alt.eq.knots(j+1).and.alt.gt.knots(j)) THEN
            mm (j+1,1,k) = 1.0
          ELSEIf (alt.lt.knots(j+1).and.alt.ge.knots(j)) THEN
            mm(j,1,k) = 1.0
          ENDIF
        ENDDO
      ENDDO
      !% higher order basis functions
      DO ii = 2,P_ord+1
        iii = ii-1
        DO k=1,data_len
          alt = new_data_alt(k)
          DO j=1,Mfull-iii
            IF (knots(j+iii).gt.knots(j)) THEN
              mm(j,ii,k) = (alt - knots(j) ) / ( knots(j+iii)-knots(j)) * mm(j,ii-1,k)
            ENDIF
            IF (j+iii+1.le.Mfull) THEN
              IF (knots(j+iii+1)>knots(j+1))THEN
                mm(j,ii,k) = mm(j,ii,k) + ( knots(j+iii+1) - alt) /       &
                                          ( knots(j+iii+1) - knots(j+1)) * &
                                             mm(j+1,ii-1,k)
              ENDIF
            ENDIF
          ENDDO
        ENDDO
      ENDDO

      allocate(model_r(model_len))
      allocate(model_i(model_len))

      DO k = 1,data_len
        write(*,*)k,data_len,new_data_alt(k)

        ! Vx
        DO i = 1,x_order
          DO j = 1,y_order
            IF(IsLatLon)THEN
              model_r = vx3_SH_dp(1,i,j, :)
              model_i = vx3_SH_dp(2,i,j, :)
            ELSE
              model_r =  real(vx3_FS_C(i,j, :))
              model_i = aimag(vx3_FS_C(i,j, :))
            ENDIF
            cilm(1,i,j) = dot_product(mm(1:model_len,P_ord+1,k),model_r(:))
            cilm(2,i,j) = dot_product(mm(1:model_len,P_ord+1,k),model_i(:))
            IF(.not.IsLatLon)THEN
              gridU_out(i,j) = complex(cilm(1,i,j),cilm(2,i,j))
            ENDIF
          ENDDO
        ENDDO
        ! Now resample the field from the Spectral expansion and write out
        IF(IsLatLon)THEN
          ! Use Spherical harmonic coefficients
          lmax = Spec_MaxOrder-1
          call MakeGrid2D(grid,  & ! o resampled map ()
                          cilm,  & ! i real SH coefficients (2,LMAX,LMAX)
                          lmax,  & ! i max SH degree of cilm
                          interval, & ! i lat/lon spacing of grid
                          nlat,  & ! number of output lat samples
                          nlong  & ! number of output lon samples
                          )
          DO i=1,nlong-1
            vx_OUT_dp(i,1:nymax_g2s,k) = grid(1:nymax_g2s,i)
          ENDDO
        ELSE
          ! Use Fourier coefficients
          call dfftw_execute_dft_c2r(plan_bak2Dfft_U,gridU_out,gridU_in)
          gridU_in = gridU_in/(nxmax_g2s*nymax_g2s)
          vx_OUT_dp(1:nxmax_g2s,1:nymax_g2s,k) = gridU_in(1:nxmax_g2s,1:nymax_g2s)
        ENDIF

        ! Vy
        DO i = 1,x_order
          DO j = 1,y_order
            IF(IsLatLon)THEN
              model_r = vy3_SH_dp(1,i,j, :)
              model_i = vy3_SH_dp(2,i,j, :)
            ELSE
              model_r =  real(vy3_FS_C(i,j, :))
              model_i = aimag(vy3_FS_C(i,j, :))
            ENDIF
            cilm(1,i,j) = dot_product(mm(1:model_len,P_ord+1,k),model_r(:))
            cilm(2,i,j) = dot_product(mm(1:model_len,P_ord+1,k),model_i(:))
            IF(.not.IsLatLon)THEN
              gridV_out(i,j) = complex(cilm(1,i,j),cilm(2,i,j))
            ENDIF
          ENDDO
        ENDDO
        ! Now resample the field from the SH expansion and write out.
        IF(IsLatLon)THEN
          ! Use Spherical harmonic coefficients

          lmax = Spec_MaxOrder-1
          call MakeGrid2D(grid,  & ! o resampled map ()
                          cilm,  & ! i real SH coefficients (2,LMAX,LMAX)
                          lmax,  & ! i max SH degree of cilm
                          interval, & ! i lat/lon spacing of grid
                          nlat,  & ! number of output lat samples
                          nlong  & ! number of output lon samples
                          )
          DO i=1,nlong-1
            vy_OUT_dp(i,1:nymax_g2s,k) = grid(1:nymax_g2s,i)
          ENDDO
        ELSE
          ! Use Fourier coefficients
          call dfftw_execute_dft_c2r(plan_bak2Dfft_V,gridV_out,gridV_in)
          gridV_in = gridV_in/(nxmax_g2s*nymax_g2s)
          vy_OUT_dp(1:nxmax_g2s,1:nymax_g2s,k) = gridV_in(1:nxmax_g2s,1:nymax_g2s)
        ENDIF

        ! Temperature
        DO i = 1,x_order
          DO j = 1,y_order
            IF(IsLatLon)THEN
              model_r = temperature3_SH_dp(1,i,j, :)
              model_i = temperature3_SH_dp(2,i,j, :)
            ELSE
              model_r =  real(temperature3_FS_C(i,j, :))
              model_i = aimag(temperature3_FS_C(i,j, :))
            ENDIF
            cilm(1,i,j) = dot_product(mm(1:model_len,P_ord+1,k),model_r(:))
            cilm(2,i,j) = dot_product(mm(1:model_len,P_ord+1,k),model_i(:))
            IF(.not.IsLatLon)THEN
              gridT_out(i,j) = complex(cilm(1,i,j),cilm(2,i,j))
            ENDIF
          ENDDO
        ENDDO
        ! Now resample the field from the SH expansion and write out
        IF(IsLatLon)THEN
          ! Use Spherical harmonic coefficients
          lmax = Spec_MaxOrder-1
          call MakeGrid2D(grid,  & ! o resampled map ()
                          cilm,  & ! i real SH coefficients (2,LMAX,LMAX)
                          lmax,  & ! i max SH degree of cilm
                          interval, & ! i lat/lon spacing of grid
                          nlat,  & ! number of output lat samples
                          nlong  & ! number of output lon samples
                          )
          DO i=1,nlong-1
            temperature_OUT_dp(i,1:nymax_g2s,k) = grid(1:nymax_g2s,i)
          ENDDO
        ELSE
          ! Use Fourier coefficients
          call dfftw_execute_dft_c2r(plan_bak2Dfft_T,gridT_out,gridT_in)
          gridT_in = gridT_in/(nxmax_g2s*nymax_g2s)
          temperature_OUT_dp(1:nxmax_g2s,1:nymax_g2s,k) = gridT_in(1:nxmax_g2s,1:nymax_g2s)
        ENDIF
      ENDDO

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !  Write Resampled data to raw files for later 
      write(*,*)"Writing output for resampled data files"
      open (unit=20,file=FILE_U_RES,access='direct', &
                   recl=2*nxmax_g2s*nymax_g2s*data_len*4)
      write(20,rec=1)(((vx_OUT_dp(i,j,k),i=1,nxmax_g2s),j=1,nymax_g2s),k=1,data_len)
      close (20)
      open (unit=20,file=FILE_V_RES,access='direct', &
                   recl=2*nxmax_g2s*nymax_g2s*data_len*4)
      write(20,rec=1)(((vy_OUT_dp(i,j,k),i=1,nxmax_g2s),j=1,nymax_g2s),k=1,data_len)
      close (20)
      open (unit=20,file=FILE_T_RES,access='direct', &
                   recl=2*nxmax_g2s*nymax_g2s*data_len*4)
      write(20,rec=1)(((temperature_OUT_dp(i,j,k),i=1,nxmax_g2s),j=1,nymax_g2s),k=1,data_len)
      close (20)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if(write_test_output)then
        ! Dump out map at bottom of Met grid
        OPEN(UNIT=32,FILE='plan_res.dat',STATUS='replace')
        DO i=1,nxmax_g2s
          DO j=1,nymax_g2s
            write(32,*)x_g2s_sp(i),y_g2s_sp(j),         &
                       real(temperature_OUT_dp(i,j,1),kind=4),  &
                       real(         vx_OUT_dp(i,j,1),kind=4),           &
                       real(         vy_OUT_dp(i,j,1),kind=4)
          ENDDO
        ENDDO
        CLOSE(32)
        ! Dump out equitorial cross-section of Met grid
        OPEN(UNIT=42,FILE='xsec_res.dat',STATUS='replace')
        DO i=1,nxmax_g2s
          DO k=1,data_len
            write(42,*)x_g2s_sp(i),real(new_data_alt(k),kind=4), &
                       real(temperature_OUT_dp(i,jout,k),kind=4),     &
                       real(         vx_OUT_dp(i,jout,k),kind=4),              &
                       real(         vy_OUT_dp(i,jout,k),kind=4)
          ENDDO
        ENDDO
        CLOSE(42)
        ! Dump out meridonal cross-section of Met grid
        OPEN(UNIT=42,FILE='ysec_res.dat',STATUS='replace')
        DO j=1,nymax_g2s
          DO k=1,data_len
            write(42,*)y_g2s_sp(j),real(new_data_alt(k),kind=4), &
                       real(temperature_OUT_dp(iout,j,k),kind=4),     &
                       real(         vx_OUT_dp(iout,j,k),kind=4),              &
                       real(         vy_OUT_dp(iout,j,k),kind=4)
          ENDDO
        ENDDO
        CLOSE(42)

      endif

      end subroutine Resample_Atmos

