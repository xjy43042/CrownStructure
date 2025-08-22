program aggr_500m_to_0p1deg
!====================================================================
! Aggregate regional 5x5° crown structure data (CW/CD/AR) 
! to global 0.1° resolution.
!
! Example: Take CH90-based prediction as input.
!
! Usage:
!   gfortran -g -mcmodel=large -fbounds-check -o aggr aggr_CWCD_globe.F90 \
!            -I/usr/include -lnetcdff -lnetcdf
!
! Input:
!   Regional 5x5° NetCDF tiles (Crown structure + PCT_PFT)
!
! Output:
!   Global 0.1° NetCDF composite
!
! Authors:
!   Jiayi Xiang; Wenzong Dong
!====================================================================
  use netcdf
  implicit none

  !=============================
  ! Parameters
  !=============================
  integer, parameter :: r8 = selected_real_kind(12)   ! double precision
  integer, parameter :: r4 = 4                        ! single precision

  ! Coarse output grid: 0.1° resolution
  real(r8), parameter :: coarse_delta = 0.1_r8
  integer, parameter :: nPFT = 8                      ! number of plant functional types

  ! Fine input grid: 500m ≈ 15s (~0.0041667°)
  integer, parameter :: nx = 1200
  integer, parameter :: ny = 1200
  real(r8), parameter :: fine_delta = 1.0_r8 / 240.0_r8

  real(r8), parameter :: missing_value = -999.0_r8
  real(r8), parameter :: lat_max =  90.0_r8
  real(r8), parameter :: lat_min = -90.0_r8
  real(r8), parameter :: lon_min = -180.0_r8
  real(r8), parameter :: lon_max =  180.0_r8

  ! Mathematical constants
  real(r8) :: pi, deg2rad, re

  ! Derived global coarse grid size
  integer :: nlon, nlat, XY3D(3)

  !=============================
  ! Arrays
  !=============================
  ! Accumulators (lon, lat, PFT)
  real(r8), allocatable :: CW_accum(:,:,:), CD_accum(:,:,:), AR_accum(:,:,:)
  real(r8), allocatable :: weight_sum(:,:,:), pct_f(:,:,:)

  ! Fine grid buffers
  real(r8), allocatable :: lat_fine(:), lon_fine(:), area_fine(:,:)
  real(r8), allocatable :: CW(:,:,:), CD(:,:,:), AR(:,:,:)

  ! Coarse grid coordinates
  real(r8), allocatable :: lon_c(:), lat_c(:)

  !=============================
  ! Loop & I/O
  !=============================
  integer :: i, j, i_fine_lon, j_fine_lat, pft_idx, pft_pct_idx
  integer :: ncid, lat_id, lon_id, CW_id, CD_id, AR_id, pct_id, ncid1
  integer :: ncout, lon_dim_out, lat_dim_out, pft_dim_out
  integer :: lon_varid, lat_varid, cr_varid, cd_varid, ar_varid
  integer :: status
  character(len=256) :: filename, pctname, outfile
  logical :: fileExists, pctExists

  ! Indices & temporaries
  integer :: i_coarse_lon, j_coarse_lat
  real(r8) :: lonval, latval
  real(r8) :: lon_west, lon_east, lat_south, lat_north

  ! Regional definition file
  integer :: reg(4), iostatus
  character(len=256) :: REGFILE
  character(len=20) :: reg1, reg2, reg3, reg4

  !=============================
  ! Initialize constants
  !=============================
  pi      = 4.0_r8 * atan(1.0_r8)
  deg2rad = pi / 180.0_r8
  re      = 6.37122e6_r8   ! Earth radius (m)

  !=============================
  ! Define coarse grid
  !=============================
  nlon = int((lon_max - lon_min) / coarse_delta + 0.5_r8)
  nlat = int((lat_max - lat_min) / coarse_delta + 0.5_r8)

  allocate(lon_c(nlon), lat_c(nlat))
  do i = 1, nlon
     lon_c(i) = lon_min + (real(i,r8) - 0.5_r8) * coarse_delta
  end do
  do j = 1, nlat
     lat_c(j) = lat_max - (real(j,r8) - 0.5_r8) * coarse_delta
  end do

  print *, "Coarse grid lon range=", lon_c(1), lon_c(nlon)
  print *, "Coarse grid lat range=", lat_c(1), lat_c(nlat)

  !=============================
  ! Allocate accumulators
  !=============================
  allocate(CW_accum(nlon,nlat,nPFT))
  allocate(CD_accum(nlon,nlat,nPFT))
  allocate(AR_accum(nlon,nlat,nPFT))
  allocate(weight_sum(nlon,nlat,nPFT))
  CW_accum = 0.0_r8; CD_accum = 0.0_r8; AR_accum = 0.0_r8; weight_sum = 0.0_r8

  ! Fine grid buffers
  allocate(lat_fine(ny), lon_fine(nx))
  allocate(area_fine(nx,ny))
  allocate(CW(nx,ny,nPFT), CD(nx,ny,nPFT), AR(nx,ny,nPFT), pct_f(nx,ny,16))

  !=============================
  ! Open region definition file
  !=============================
  REGFILE = "/tera10/yuanhua/xiangjy/factors_new/make/aggr/final/reg_5x5"
  open(unit=11, file=REGFILE, form="formatted", status="old", action="read")

  !=============================
  ! Loop through regional tiles
  !=============================
  do
     read(11, *, iostat=iostatus) reg
     if (iostatus /= 0) exit  ! EOF

     ! Convert integers to strings
     write(reg1,"(i4)") reg(1)
     write(reg2,"(i4)") reg(2)
     write(reg3,"(i4)") reg(3)
     write(reg4,"(i4)") reg(4)

     !-----------------------------------
     ! Load PCT_PFT file
     !-----------------------------------
     write(pctname,'(A,A,"_",A,"_",A,"_",A,".MOD2020.nc")') &
          '/tera10/yuanhua/xiangjy/factors_new/rawdata/PFTLAI_5x5/RG_', &
          trim(adjustL(reg1)), trim(adjustL(reg2)), trim(adjustL(reg3)), trim(adjustL(reg4))

     inquire(file=trim(pctname), exist=pctExists)
     if (.not. pctExists) then
        print *, 'PCT_PFT file does not exist: ', trim(pctname)
        cycle
     else
        print *, 'Loading PCT_PFT file: ', pctname
     end if

     call nccheck(nf90_open(pctname, NF90_NOWRITE, ncid1))
     call nccheck(nf90_inq_varid(ncid1, 'PCT_PFT', pct_id))
     call nccheck(nf90_get_var  (ncid1, pct_id, pct_f))
     call nccheck(nf90_close(ncid1))

     !-----------------------------------
     ! Load CW/CD/AR canopy structure file
     !-----------------------------------
     write(filename,'(A,A,"_",A,"_",A,"_",A,".CanopyStructure_15s_CH90.nc")') &
          '/tera10/yuanhua/xiangjy/factors_new/make/PreRegion/output/RG_', &
          trim(adjustL(reg1)), trim(adjustL(reg2)), trim(adjustL(reg3)), trim(adjustL(reg4))

     inquire(file=trim(filename), exist=fileExists)
     if (fileExists) print *, 'Loading CWCD file: ', filename

     call nccheck(nf90_open(trim(filename), NF90_NOWRITE, ncid))
     call nccheck(nf90_inq_varid(ncid,'lat',lat_id)); call nccheck(nf90_get_var(ncid, lat_id, lat_fine))
     call nccheck(nf90_inq_varid(ncid,'lon',lon_id)); call nccheck(nf90_get_var(ncid, lon_id, lon_fine))

     call nccheck(nf90_inq_varid(ncid,'CROWN_WIDTH' , CW_id)); call nccheck(nf90_get_var(ncid, CW_id, CW))
     call nccheck(nf90_inq_varid(ncid,'CROWN_DEPTH' , CD_id)); call nccheck(nf90_get_var(ncid, CD_id, CD))
     call nccheck(nf90_inq_varid(ncid,'ASPECT_RATIO', AR_id)); call nccheck(nf90_get_var(ncid, AR_id, AR))

     !-----------------------------------
     ! Compute fine-cell area (spherical)
     !-----------------------------------
     do i_fine_lon = 1, nx
        lon_west = lon_fine(i_fine_lon) - 0.5_r8 * fine_delta
        lon_east = lon_fine(i_fine_lon) + 0.5_r8 * fine_delta
        do j_fine_lat = 1, ny
           lat_south = lat_fine(j_fine_lat) - 0.5_r8 * fine_delta
           lat_north = lat_fine(j_fine_lat) + 0.5_r8 * fine_delta
           area_fine(i_fine_lon,j_fine_lat) = (lon_east-lon_west)*deg2rad * &
                                              (sin(lat_north*deg2rad)-sin(lat_south*deg2rad)) * re*re
        end do
     end do

     !-----------------------------------
     ! Accumulate to regional coarse grid
     !-----------------------------------
     do i_fine_lon = 1, nx
        lonval = lon_fine(i_fine_lon)
        i_coarse_lon = int( floor( (lonval - lon_min)/coarse_delta ) ) + 1
        if (i_coarse_lon < 1 .or. i_coarse_lon > nlon) cycle
        do j_fine_lat = 1, ny
           latval = lat_fine(j_fine_lat)
           j_coarse_lat = int( floor( (lat_max - latval)/coarse_delta ) ) + 1
           if (j_coarse_lat < 1 .or. j_coarse_lat > nlat) cycle

           do pft_idx = 1, nPFT
              pft_pct_idx = pft_idx + 1  ! pct_f indexing (skip first column)
              if (pct_f(i_fine_lon,j_fine_lat,pft_pct_idx) < 0.01_r8) cycle
              if (.not. is_valid(CW(i_fine_lon,j_fine_lat,pft_idx))) cycle
              if (.not. is_valid(CD(i_fine_lon,j_fine_lat,pft_idx))) cycle
              if (.not. is_valid(AR(i_fine_lon,j_fine_lat,pft_idx))) cycle

              CW_accum(i_coarse_lon,j_coarse_lat,pft_idx) = CW_accum(i_coarse_lon,j_coarse_lat,pft_idx) + &
                     CW(i_fine_lon,j_fine_lat,pft_idx) * area_fine(i_fine_lon,j_fine_lat) * pct_f(i_fine_lon,j_fine_lat,pft_pct_idx)
              CD_accum(i_coarse_lon,j_coarse_lat,pft_idx) = CD_accum(i_coarse_lon,j_coarse_lat,pft_idx) + &
                     CD(i_fine_lon,j_fine_lat,pft_idx) * area_fine(i_fine_lon,j_fine_lat) * pct_f(i_fine_lon,j_fine_lat,pft_pct_idx)
              AR_accum(i_coarse_lon,j_coarse_lat,pft_idx) = AR_accum(i_coarse_lon,j_coarse_lat,pft_idx) + &
                     AR(i_fine_lon,j_fine_lat,pft_idx) * area_fine(i_fine_lon,j_fine_lat) * pct_f(i_fine_lon,j_fine_lat,pft_pct_idx)
              weight_sum(i_coarse_lon,j_coarse_lat,pft_idx) = weight_sum(i_coarse_lon,j_coarse_lat,pft_idx) + &
                     area_fine(i_fine_lon,j_fine_lat) * pct_f(i_fine_lon,j_fine_lat,pft_pct_idx)
           end do
        end do
     end do

     call nccheck(nf90_close(ncid))
  end do  ! end of region loop

  !=============================
  ! Normalize accumulators (pct-area-weighted mean)
  !=============================
  do i = 1, nlon
     do j = 1, nlat
        do pft_idx = 1, nPFT
           if (weight_sum(i,j,pft_idx) > 0.0_r8) then
              CW_accum(i,j,pft_idx) = CW_accum(i,j,pft_idx) / weight_sum(i,j,pft_idx)
              CD_accum(i,j,pft_idx) = CD_accum(i,j,pft_idx) / weight_sum(i,j,pft_idx)
              AR_accum(i,j,pft_idx) = AR_accum(i,j,pft_idx) / weight_sum(i,j,pft_idx)
           else
              CW_accum(i,j,pft_idx) = missing_value
              CD_accum(i,j,pft_idx) = missing_value
              AR_accum(i,j,pft_idx) = missing_value
           end if
        end do
     end do
  end do

  !=============================
  ! Write output global NetCDF
  !=============================
  outfile = '/tera10/yuanhua/xiangjy/factors_new/make/aggr/final/output/Global_CrownStructure_0.1deg_CH90.nc'
  print*, 'Creating global output: ', trim(outfile)

  call nccheck(nf90_create(trim(outfile), NF90_NETCDF4, ncout))
  call nccheck(nf90_def_dim(ncout, 'lon', nlon, lon_dim_out))
  call nccheck(nf90_def_dim(ncout, 'lat', nlat, lat_dim_out))
  call nccheck(nf90_def_dim(ncout, 'PFT', nPFT, pft_dim_out))

  call nccheck(nf90_def_var(ncout, 'lon', NF90_DOUBLE, (/lon_dim_out/), lon_varid, deflate_level=9))
  call nccheck(nf90_put_att(ncout, lon_varid, 'long_name', 'Longitude'))
  call nccheck(nf90_put_att(ncout, lon_varid, 'units', 'degrees_east'))

  call nccheck(nf90_def_var(ncout, 'lat', NF90_DOUBLE, (/lat_dim_out/), lat_varid, deflate_level=9))
  call nccheck(nf90_put_att(ncout, lat_varid, 'long_name', 'Latitude'))
  call nccheck(nf90_put_att(ncout, lat_varid, 'units', 'degrees_north'))

  XY3D = (/lon_dim_out,lat_dim_out,pft_dim_out/)
  call nccheck(nf90_def_var(ncout, 'CROWN_WIDTH', NF90_FLOAT, XY3D, cr_varid, deflate_level=9))
  call nccheck(nf90_put_att(ncout, cr_varid, '_FillValue', real(missing_value,r4)))
  call nccheck(nf90_put_att(ncout, cr_varid, 'long_name', 'crown width predicted by & 
  the 90th percentile canopy-top height'))
  call nccheck(nf90_put_att(ncout, cr_varid, 'units','m'))

  call nccheck(nf90_def_var(ncout, 'CROWN_DEPTH', NF90_FLOAT, XY3D, cd_varid, deflate_level=9))
  call nccheck(nf90_put_att(ncout, cd_varid, '_FillValue', real(missing_value,r4)))
  call nccheck(nf90_put_att(ncout, cd_varid, 'long_name', 'crown depth predicted by & 
  the 90th percentile canopy-top height'))
  call nccheck(nf90_put_att(ncout, cd_varid, 'units','m'))

  call nccheck(nf90_def_var(ncout, 'ASPECT_RATIO', NF90_FLOAT, XY3D, ar_varid, deflate_level=9))
  call nccheck(nf90_put_att(ncout, ar_varid, '_FillValue', real(missing_value,r4)))
  call nccheck(nf90_put_att(ncout, ar_varid, 'long_name', 'crown depth to width ratio & 
  predicted by the 90th percentile canopy-top height'))

  call nccheck(nf90_enddef(ncout))

  call nccheck(nf90_put_var(ncout, lon_varid, lon_c))
  call nccheck(nf90_put_var(ncout, lat_varid, lat_c))
  call nccheck(nf90_put_var(ncout, cr_varid, real(CW_accum,r4)))
  call nccheck(nf90_put_var(ncout, cd_varid, real(CD_accum,r4)))
  call nccheck(nf90_put_var(ncout, ar_varid, real(AR_accum,r4)))

  ! Global attributes
  call nccheck(nf90_put_att(ncout, NF90_GLOBAL, 'Title',      'Land surface model input crown morphological structure data'))
  call nccheck(nf90_put_att(ncout, NF90_GLOBAL, 'resolution', '0.1 degree, 3600x1800 (lon x lat) global'))
  call nccheck(nf90_put_att(ncout, NF90_GLOBAL, 'coordinate', 'Geographic, degrees longitude and latitude'))
  call nccheck(nf90_put_att(ncout, NF90_GLOBAL, 'source',     'CanopyStructure_15s_CH90 tiles (500m)'))
  call nccheck(nf90_put_att(ncout, NF90_GLOBAL, 'Authors',    'Yongjiu Dai land group at Sun Yat-sen University'))
  call nccheck(nf90_put_att(ncout, NF90_GLOBAL, 'Address',    'School of Atmospheric Sciences, & 
  Sun Yat-sen University, Zhuhai, China'))
  call nccheck(nf90_put_att(ncout, NF90_GLOBAL, 'Contact',    'Jiayi Xiang (xiangjy8@mail2.sysu.edu.cn), & 
  Hua Yuan (yuanh25@mail.sysu.edu.cn)'))

  call nccheck(nf90_close(ncout))
  print*, 'Successfully created '//trim(outfile)//'!'

  !=============================
  ! Cleanup
  !=============================
  if (allocated(CW_accum))  deallocate(CW_accum)
  if (allocated(CD_accum))  deallocate(CD_accum)
  if (allocated(AR_accum))  deallocate(AR_accum)
  if (allocated(weight_sum)) deallocate(weight_sum)
  if (allocated(lon_c))     deallocate(lon_c)
  if (allocated(lat_c))     deallocate(lat_c)
  if (allocated(lat_fine))  deallocate(lat_fine)
  if (allocated(lon_fine))  deallocate(lon_fine)
  if (allocated(area_fine)) deallocate(area_fine)
  if (allocated(CW))        deallocate(CW)
  if (allocated(CD))        deallocate(CD)
  if (allocated(AR))        deallocate(AR)

contains

  !-----------------------------
  subroutine nccheck(status)
    integer, intent(in) :: status
    if (status /= nf90_noerr) then
       write(*,*) 'NetCDF error: ', trim(nf90_strerror(status))
       stop 2
    end if
  end subroutine nccheck

  !-----------------------------
  logical function is_valid(x)
    use, intrinsic :: ieee_arithmetic
    real(r8), intent(in) :: x
    is_valid = (.not. ieee_is_nan(x)) .and. (abs(x - missing_value) > 1.0e-6_r8)
  end function is_valid

end program aggr_500m_to_0p1deg
