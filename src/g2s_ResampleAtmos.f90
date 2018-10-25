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
!  G2S_ResampleAtmos
!
!    This program reads a coeffiecient file and generates three gridded binary
!    files for U, V, and T that can be used the the g2s_Extract_[Sonde,Xsec,Grid]
!    programs.  Vertical spacing of grid can either be hardwired below (must be
!    consistent with G2S_Extract.f90), or given on the command-line.
!
!##############################################################################

      program G2S_ResampleAtmos

      use G2S_globvar

      implicit none

      integer :: nargs
      character(len=50) :: llinebuffer
      integer :: root_len 


       !===================================================
       !  Set the default values here
        ! Vertical node spacing
          ! This one is a bit coarse
        !nz1 = 15
        !nz2 = 18
        !nz3 = 30
        !dz1 = 1.0
        !dz2 = 2.0
        !dz3 = 5.0
          ! This one is good for operations
        !nz1 = 26
        !nz2 = 50
        !nz3 = 50
        !dz1 = 1.0
        !dz2 = 1.5
        !dz3 = 2.0
          ! Here the three segments all have dz=1.0
        nz1 = 50
        nz2 = 50
        nz3 = 101
        dz1 = 1.0
        dz2 = 1.0
        dz3 = 1.0

       !===================================================


!     TEST READ COMMAND LINE ARGUMENTS
      nargs = iargc()
      if (nargs.ge.1) then
        call getarg(1,llinebuffer)
        file_SH = trim(adjustl(llinebuffer))
        ! If a filename is provided, strip off the ".nc" and append
        ! the labels for TUV
        root_len = len(trim(file_SH)) - 3
        file_T_RES = file_SH(1:root_len) // "_T_res.raw"
        file_U_RES = file_SH(1:root_len) // "_U_res.raw"
        file_V_RES = file_SH(1:root_len) // "_V_res.raw"
      else
        write(G2S_global_info,*)"No command-line arguments given."
        write(G2S_global_info,*)" Normal usage is to provide the coefficient file as an argument."
        write(G2S_global_info,*)"Assuming default SH file name and"
        write(G2S_global_info,*)"using default Res file names."
      endif
      if (nargs.ge.2) then
        ! Additional command-line arguments are provided.  Assume these are
        ! nz1 dz1 nz2 dz2 nz3 dz3
        if (nargs.lt.7)then
          write(G2S_global_info,*)"Multiple command-line argument provided, but not enough to"
          write(G2S_global_info,*)"specify the grid."
          write(G2S_global_info,*)"Expected format:"
          write(G2S_global_info,*)"g2s_ResampleAtmos Coeffic_File.nc nz1 dz1 nz2 dz2 nz3 dz3"
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
        write(G2S_global_info,*)"Using the default values for vertical grid."
        write(G2S_global_info,*)"nz1 = ",nz1
        write(G2S_global_info,*)"dz1 = ",dz1
        write(G2S_global_info,*)"nz2 = ",nz2
        write(G2S_global_info,*)"dz2 = ",dz2
        write(G2S_global_info,*)"nz3 = ",nz3
        write(G2S_global_info,*)"dz3 = ",dz3
      endif
      write(G2S_global_info,*)"Opening SH coefficient file: ",file_SH
      call Read_SH_nc()
 
      write(G2S_global_info,*)&
        "Recreating Atmos based on spectra at native lat/lon grid."
      write(G2S_global_info,*)"Number of levels = ",nz1+nz2+nz3
      call Resample_Atmos

      write(G2S_global_info,*)"G2S_ResampleAtmos exited normally."

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

      INQUIRE( file=adjustl(trim(file_SH)), EXIST=IsThere )
      if(.not.IsThere)then
        write(G2S_global_error,*)"ERROR: Could not find input file."
        write(G2S_global_error,*)adjustl(trim(file_SH))
        stop 1
      endif

      ! Open existing netcdf file
      inquire(file=file_SH, exist=IsThere )
      if(.not.IsThere)then
        ! Coefficient file does not exits, issue error message and exit
        write(G2S_global_error,*)"ERROR: The coefficient file does not exist"
        stop 1
      endif
      nSTAT = nf90_open(file_SH,NF90_NOWRITE,ncid)
      if(nSTAT.ne.NF90_NOERR)then
        write(G2S_global_error,*)'ERROR: nf90_open:',nf90_strerror(nSTAT)
        write(G2S_global_info,*)'Could not open ',file_SH
        write(G2S_global_info,*)'Exiting'
        stop 1
      endif

      nSTAT = nf90_get_att(ncid,nf90_global,"nx_g2s",nxmax_g2s)
      if(nSTAT.ne.NF90_NOERR)then
        write(G2S_global_error,*)'ERROR: nf90_get_att:',nf90_strerror(nSTAT)
        stop 1 
      endif
      nSTAT = nf90_get_att(ncid,nf90_global,"ny_g2s",nymax_g2s)
      if(nSTAT.ne.NF90_NOERR)then
        write(G2S_global_error,*)'ERROR: nf90_get_att:',nf90_strerror(nSTAT)
        stop 1 
      endif
      nSTAT = nf90_get_att(ncid,nf90_global,"dx_g2s",dx_g2s)
      if(nSTAT.ne.NF90_NOERR)then
        write(G2S_global_error,*)'ERROR: nf90_get_att:',nf90_strerror(nSTAT)
        write(G2S_global_info,*)"WARNING: Attribute dx_g2s not found in Spectral Coefficient file."
        write(G2S_global_info,*)"         Assuming 0.5"
        dx_g2s = 0.5
      endif
      nSTAT = nf90_get_att(ncid,nf90_global,"dy_g2s",dy_g2s)
      if(nSTAT.ne.NF90_NOERR)then
        write(G2S_global_error,*)'ERROR: nf90_get_att:',nf90_strerror(nSTAT)
        write(G2S_global_info,*)"WARNING: Attribute dy_g2s not found in Spectral Coefficient file."
        write(G2S_global_info,*)"         Assuming 0.5"
        dy_g2s = 0.5
      endif
      nSTAT = nf90_get_att(ncid,nf90_global,"Spec_MaxOrder",Spec_MaxOrder)
      if(nSTAT.ne.NF90_NOERR)then
        write(G2S_global_error,*)'ERROR: nf90_get_att:',nf90_strerror(nSTAT)
        write(G2S_global_info,*)&
          "WARNING: Attribute Spec_MaxOrder not found in Spectral Coefficient file."
        write(G2S_global_info,*)"         Assuming 120"
        Spec_MaxOrder = 120
      endif
      nSTAT = nf90_get_att(ncid,nf90_global,"proj",Comp_projection_line)
      if(nSTAT.ne.NF90_NOERR)then
        write(G2S_global_error,*)'ERROR: nf90_get_att:',nf90_strerror(nSTAT)
        write(G2S_global_info,*)&
          "WARNING: Attribute proj not found in Spectral Coefficient file."
        write(G2S_global_info,*)"         Assuming Lat/Lon grid"
        Comp_projection_line = "1 4 -107.0 50.0 50.0 50.0 6367.470"
      endif
      read(Comp_projection_line,*)ilatlonflag,iprojflag
      if (ilatlonflag.eq.0) then
        ! expecting input variables to be in the same projection as
        ! specified by iprojflag and parameters
        IsLatLon          = .false.
        x_order = nxmax_g2s/2+1
        y_order = nymax_g2s
        nSTAT = nf90_get_att(ncid,nf90_global,"xmin_g2s",xmin_g2s)
        if(nSTAT.ne.NF90_NOERR)then
          write(G2S_global_error,*)'ERROR: nf90_get_att:',nf90_strerror(nSTAT)
          write(G2S_global_error,*)"ERROR: grid is projected, but no xmin_g2s in SH file"
          stop 1
        endif
        nSTAT = nf90_get_att(ncid,nf90_global,"ymin_g2s",ymin_g2s)
        if(nSTAT.ne.NF90_NOERR)then
          write(G2S_global_error,*)'ERROR: nf90_get_att:',nf90_strerror(nSTAT)
          write(G2S_global_error,*)"ERROR: grid is projected, but no ymin_g2s in SH file"
          stop 1
        endif
      else
        ! expecting input variables to be in lat/lon
        x_order = Spec_MaxOrder+1
        y_order = Spec_MaxOrder+1
        IsLatLon          = .true.
      endif

      allocate(x_g2s_sp(nxmax_g2s))
      allocate(y_g2s_sp(nymax_g2s))

      if(IsLatLon)then
        do i=1,nxmax_g2s
          x_g2s_sp(i) = (i-1)*real(dx_g2s,kind=4)
        enddo
        do j=1,nymax_g2s
          y_g2s_sp(j) = -90.0_4+(j-1)*real(dx_g2s,kind=4)
        enddo
      else
        do i=1,nxmax_g2s
          x_g2s_sp(i) = xmin_g2s + (i-1)*real(dx_g2s,kind=4)
        enddo
        do j=1,nymax_g2s
          y_g2s_sp(j) = ymin_g2s + (j-1)*real(dy_g2s,kind=4)
        enddo
      endif
      !nSTAT = nf90_put_att(ncid,nf90_global,"Spline_Order",P_ord)
      !P_ord = 3
      
      nSTAT = nf90_inq_dimid(ncid,odim_names(1),r_dim_id)
      if(nSTAT.ne.NF90_NOERR)then
        write(G2S_global_error,*)'ERROR: nf90_inq_dimid:',nf90_strerror(nSTAT)
        stop 1
      endif
       ! r_dim length = 2
      nSTAT = nf90_inq_dimid(ncid,odim_names(2),z_dim_id)
      if(nSTAT.ne.NF90_NOERR)then
        write(G2S_global_error,*)'ERROR: nf90_inq_dimid:',nf90_strerror(nSTAT)
        stop 1
      endif
!      nSTAT = nf90_Inquire_Dimension(ncid,z_dim_id,len=Nknots)
!      Nknots = Nknots - P_ord
      nSTAT = nf90_inq_dimid(ncid,odim_names(3),y_dim_id)
      if(nSTAT.ne.NF90_NOERR)then
        write(G2S_global_error,*)'ERROR: nf90_inq_dimid:',nf90_strerror(nSTAT)
        stop 1
      endif
      !nSTAT = nf90_Inquire_Dimension(ncid,y_dim_id,len=Spec_MaxOrder)
      !Spec_MaxOrder = Spec_MaxOrder - 1
      nSTAT = nf90_inq_dimid(ncid,odim_names(4),x_dim_id)
      if(nSTAT.ne.NF90_NOERR)then
        write(G2S_global_error,*)'ERROR: nf90_inq_dimid:',nf90_strerror(nSTAT)
        stop 1
      endif
        ! Assume x length is the same as y length

      !nSTAT = nf90_def_var(ncid,"knot",nf90_double,  &
      !      (/z_dim_id/), &
      !        kn_var_id)

      model_len = Nknots+P_ord
      allocate(dum1d_out(model_len))
      allocate(dum4d_out(x_order, y_order,Nknots+P_ord,2))
      if(IsLatLon)then
        allocate(         vx3_SH_dp(2, x_order, y_order,model_len))
        allocate(         vy3_SH_dp(2, x_order, y_order,model_len))
        allocate(temperature3_SH_dp(2, x_order, y_order,model_len))
      else
        allocate(         vx3_FS_C(x_order,y_order,model_len))
        allocate(         vy3_FS_C(x_order,y_order,model_len))
        allocate(temperature3_FS_C(x_order,y_order,model_len))
      endif

      nSTAT = nf90_inq_varid(ncid,"knot",kn_var_id)
      if(nSTAT.ne.NF90_NOERR)then
        write(G2S_global_error,*)'ERROR: nf90_inq_varid:',nf90_strerror(nSTAT)
        stop 1
      endif
      nSTAT = nf90_get_var(ncid,kn_var_id,dum1d_out, &
             start = (/1/),       &
             count = (/model_len/))
      if(nSTAT.ne.NF90_NOERR)then
        write(G2S_global_error,*)'ERROR: nf90_get_var:',nf90_strerror(nSTAT)
        stop 1
      endif
      knots(1:model_len) = real(dum1d_out(1:model_len),kind=4)
      knots(model_len:Nknots+2*P_ord) = knots(model_len)
      iknots(1:Nknots) = real(dum1d_out(P_ord+1:model_len),kind=4)

      nSTAT = nf90_inq_varid(ncid,"vx_sh",vx_var_id)
      if(nSTAT.ne.NF90_NOERR) &
         write(G2S_global_error,*)'ERROR: inq_varid:Vx ',nf90_strerror(nSTAT)
      dum4d_out = 0.0
      nSTAT = nf90_get_var(ncid,vx_var_id,dum4d_out, &
             start = (/1,1,1,1/),       &
             count = (/x_order, y_order,Nknots+P_ord,2/))
      if(nSTAT.ne.NF90_NOERR) &
         write(G2S_global_error,*)'ERROR: get_var:Vx ',nf90_strerror(nSTAT)
      if(IsLatLon)then
        vx3_SH_dp(1,:,:,:) = dum4d_out(:,:,:,1)
        vx3_SH_dp(2,:,:,:) = dum4d_out(:,:,:,2)
      else
        do i=1,x_order
          do j=1,y_order
            do k=1,model_len
              vx3_FS_C(i,j,k) = complex(dum4d_out(i,j,k,1),dum4d_out(i,j,k,2))
            enddo
          enddo
        enddo
      endif

      nSTAT = nf90_inq_varid(ncid,"vy_sh",vy_var_id)
      if(nSTAT.ne.NF90_NOERR) &
         write(G2S_global_error,*)'ERROR: inq_varid:Vy ',nf90_strerror(nSTAT)
      dum4d_out = 0.0
      nSTAT = nf90_get_var(ncid,vy_var_id,dum4d_out, &
             start = (/1,1,1,1/),       &
             count = (/x_order, y_order,Nknots+P_ord,2/))
      if(nSTAT.ne.NF90_NOERR) &
         write(G2S_global_error,*)'ERROR: get_var:Vy ',nf90_strerror(nSTAT)
      if(IsLatLon)then
        vy3_SH_dp(1,:,:,:) = dum4d_out(:,:,:,1)
        vy3_SH_dp(2,:,:,:) = dum4d_out(:,:,:,2)
      else
        do i=1,x_order
          do j=1,y_order
            do k=1,model_len
              vy3_FS_C(i,j,k) = complex(dum4d_out(i,j,k,1),dum4d_out(i,j,k,2))
            enddo
          enddo
        enddo
      endif

      nSTAT = nf90_inq_varid(ncid,"tp_sh",temp_var_id)
      if(nSTAT.ne.NF90_NOERR) &
         write(G2S_global_error,*)'ERROR: inq_varid:Temperature ',nf90_strerror(nSTAT)
      dum4d_out = 0.0
      nSTAT = nf90_get_var(ncid,temp_var_id,dum4d_out, &
             start = (/1,1,1,1/),       &
             count = (/x_order, y_order,Nknots+P_ord,2/))
      if(nSTAT.ne.NF90_NOERR) &
         write(G2S_global_error,*)'ERROR: get_var:Temperature ',nf90_strerror(nSTAT)
      if(IsLatLon)then
        temperature3_SH_dp(1,:,:,:) = dum4d_out(:,:,:,1)
        temperature3_SH_dp(2,:,:,:) = dum4d_out(:,:,:,2)
      else
        do i=1,x_order
          do j=1,y_order
            do k=1,model_len
              temperature3_FS_C(i,j,k) = complex(dum4d_out(i,j,k,1),dum4d_out(i,j,k,2))
            enddo
          enddo
        enddo
      endif
      ! Close file
      nSTAT = nf90_close(ncid)
      if(nSTAT.ne.NF90_NOERR)then
        write(G2S_global_error,*)'ERROR: nf90_close:',nf90_strerror(nSTAT)
        stop 1
      endif

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

      if(.not.IsLatLon)then
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
      endif

      allocate(vx_OUT_dp(nxmax_g2s,nymax_g2s,data_len))
      allocate(vy_OUT_dp(nxmax_g2s,nymax_g2s,data_len))
      allocate(temperature_OUT_dp(nxmax_g2s,nymax_g2s,data_len))

      ! Set up for recreating B-Splines
      allocate(mm(Mfull,P_ord+1,data_len))
      allocate(new_data_alt(data_len))

      do k = 1,nz1
        new_data_alt(k) = (k-1)*dz1
      enddo
      do k = 1,nz2
        new_data_alt(k+nz1) = new_data_alt(nz1) + (k)*dz2
      enddo
      do k = 1,nz3
        new_data_alt(k+nz1+nz2) = new_data_alt(nz1+nz2) + (k)*dz3
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

      ! Now get B-Spline basis for the HWT data
      !% Generate coefficient matrix of B-Splines
      ioff = 1
      mm = 0.0
      do k = 1,data_len
        alt = new_data_alt(k)
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
          alt = new_data_alt(k)
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

      allocate(model_r(model_len))
      allocate(model_i(model_len))

      do k = 1,data_len
        write(G2S_global_info,*)k,data_len,new_data_alt(k)

        ! Vx
        do i = 1,x_order
          do j = 1,y_order
            if(IsLatLon)then
              model_r = vx3_SH_dp(1,i,j, :)
              model_i = vx3_SH_dp(2,i,j, :)
            else
              model_r =  real(vx3_FS_C(i,j, :))
              model_i = aimag(vx3_FS_C(i,j, :))
            endif
            cilm(1,i,j) = dot_product(mm(1:model_len,P_ord+1,k),model_r(:))
            cilm(2,i,j) = dot_product(mm(1:model_len,P_ord+1,k),model_i(:))
            if(.not.IsLatLon)then
              gridU_out(i,j) = complex(cilm(1,i,j),cilm(2,i,j))
            endif
          enddo
        enddo
        ! Now resample the field from the Spectral expansion and write out
        if(IsLatLon)then
          ! Use Spherical harmonic coefficients
          lmax = Spec_MaxOrder-1
          call MakeGrid2D(grid,  & ! o resampled map ()
                          cilm,  & ! i real SH coefficients (2,LMAX,LMAX)
                          lmax,  & ! i max SH degree of cilm
                          interval, & ! i lat/lon spacing of grid
                          nlat,  & ! number of output lat samples
                          nlong  & ! number of output lon samples
                          )
          do i=1,nlong-1
            vx_OUT_dp(i,1:nymax_g2s,k) = grid(1:nymax_g2s,i)
          enddo
        else
          ! Use Fourier coefficients
          call dfftw_execute_dft_c2r(plan_bak2Dfft_U,gridU_out,gridU_in)
          gridU_in = gridU_in/(nxmax_g2s*nymax_g2s)
          vx_OUT_dp(1:nxmax_g2s,1:nymax_g2s,k) = gridU_in(1:nxmax_g2s,1:nymax_g2s)
        endif

        ! Vy
        do i = 1,x_order
          do j = 1,y_order
            if(IsLatLon)then
              model_r = vy3_SH_dp(1,i,j, :)
              model_i = vy3_SH_dp(2,i,j, :)
            else
              model_r =  real(vy3_FS_C(i,j, :))
              model_i = aimag(vy3_FS_C(i,j, :))
            endif
            cilm(1,i,j) = dot_product(mm(1:model_len,P_ord+1,k),model_r(:))
            cilm(2,i,j) = dot_product(mm(1:model_len,P_ord+1,k),model_i(:))
            if(.not.IsLatLon)then
              gridV_out(i,j) = complex(cilm(1,i,j),cilm(2,i,j))
            endif
          enddo
        enddo
        ! Now resample the field from the SH expansion and write out.
        if(IsLatLon)then
          ! Use Spherical harmonic coefficients

          lmax = Spec_MaxOrder-1
          call MakeGrid2D(grid,  & ! o resampled map ()
                          cilm,  & ! i real SH coefficients (2,LMAX,LMAX)
                          lmax,  & ! i max SH degree of cilm
                          interval, & ! i lat/lon spacing of grid
                          nlat,  & ! number of output lat samples
                          nlong  & ! number of output lon samples
                          )
          do i=1,nlong-1
            vy_OUT_dp(i,1:nymax_g2s,k) = grid(1:nymax_g2s,i)
          enddo
        else
          ! Use Fourier coefficients
          call dfftw_execute_dft_c2r(plan_bak2Dfft_V,gridV_out,gridV_in)
          gridV_in = gridV_in/(nxmax_g2s*nymax_g2s)
          vy_OUT_dp(1:nxmax_g2s,1:nymax_g2s,k) = gridV_in(1:nxmax_g2s,1:nymax_g2s)
        endif

        ! Temperature
        do i = 1,x_order
          do j = 1,y_order
            if(IsLatLon)then
              model_r = temperature3_SH_dp(1,i,j, :)
              model_i = temperature3_SH_dp(2,i,j, :)
            else
              model_r =  real(temperature3_FS_C(i,j, :))
              model_i = aimag(temperature3_FS_C(i,j, :))
            endif
            cilm(1,i,j) = dot_product(mm(1:model_len,P_ord+1,k),model_r(:))
            cilm(2,i,j) = dot_product(mm(1:model_len,P_ord+1,k),model_i(:))
            if(.not.IsLatLon)then
              gridT_out(i,j) = complex(cilm(1,i,j),cilm(2,i,j))
            endif
          enddo
        enddo
        ! Now resample the field from the SH expansion and write out
        if(IsLatLon)then
          ! Use Spherical harmonic coefficients
          lmax = Spec_MaxOrder-1
          call MakeGrid2D(grid,  & ! o resampled map ()
                          cilm,  & ! i real SH coefficients (2,LMAX,LMAX)
                          lmax,  & ! i max SH degree of cilm
                          interval, & ! i lat/lon spacing of grid
                          nlat,  & ! number of output lat samples
                          nlong  & ! number of output lon samples
                          )
          do i=1,nlong-1
            temperature_OUT_dp(i,1:nymax_g2s,k) = grid(1:nymax_g2s,i)
          enddo
        else
          ! Use Fourier coefficients
          call dfftw_execute_dft_c2r(plan_bak2Dfft_T,gridT_out,gridT_in)
          gridT_in = gridT_in/(nxmax_g2s*nymax_g2s)
          temperature_OUT_dp(1:nxmax_g2s,1:nymax_g2s,k) = gridT_in(1:nxmax_g2s,1:nymax_g2s)
        endif
      enddo

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !  Write Resampled data to raw files for later 
      write(G2S_global_info,*)"Writing output for resampled data files"
      open (unit=20,file=file_U_RES,access='direct', &
                   recl=2*nxmax_g2s*nymax_g2s*data_len*4)
      write(20,rec=1)(((vx_OUT_dp(i,j,k),i=1,nxmax_g2s),j=1,nymax_g2s),k=1,data_len)
      close (20)
      open (unit=20,file=file_V_RES,access='direct', &
                   recl=2*nxmax_g2s*nymax_g2s*data_len*4)
      write(20,rec=1)(((vy_OUT_dp(i,j,k),i=1,nxmax_g2s),j=1,nymax_g2s),k=1,data_len)
      close (20)
      open (unit=20,file=file_T_RES,access='direct', &
                   recl=2*nxmax_g2s*nymax_g2s*data_len*4)
      write(20,rec=1)(((temperature_OUT_dp(i,j,k),i=1,nxmax_g2s),j=1,nymax_g2s),k=1,data_len)
      close (20)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if(write_test_output)then
        ! Dump out map at bottom of Met grid
        open(unit=32,file='plan_res.dat',status='replace')
        do i=1,nxmax_g2s
          do j=1,nymax_g2s
            write(32,*)x_g2s_sp(i),                             &
                       y_g2s_sp(j),                             &
                       real(temperature_OUT_dp(i,j,1),kind=4),  &
                       real(         vx_OUT_dp(i,j,1),kind=4),  &
                       real(         vy_OUT_dp(i,j,1),kind=4)
          enddo
        enddo
        close(32)
        ! Dump out equitorial cross-section of Met grid
        open(unit=42,file='xsec_res.dat',status='replace')
        do i=1,nxmax_g2s
          do k=1,data_len
            write(42,*)x_g2s_sp(i),                               &
                       y_g2s_sp(jout),                            &
                       real(new_data_alt(k),kind=4),              &
                       real(temperature_OUT_dp(i,jout,k),kind=4), &
                       real(         vx_OUT_dp(i,jout,k),kind=4), &
                       real(         vy_OUT_dp(i,jout,k),kind=4)
          enddo
        enddo
        close(42)
        ! Dump out meridonal cross-section of Met grid
        open(unit=42,file='ysec_res.dat',status='replace')
        do j=1,nymax_g2s
          do k=1,data_len
            write(42,*)x_g2s_sp(iout),                            &
                       y_g2s_sp(j),                               &
                       real(new_data_alt(k),kind=4),              &
                       real(temperature_OUT_dp(iout,j,k),kind=4), &
                       real(         vx_OUT_dp(iout,j,k),kind=4), &
                       real(         vy_OUT_dp(iout,j,k),kind=4)
          enddo
        enddo
        close(42)

      endif

      end subroutine Resample_Atmos

