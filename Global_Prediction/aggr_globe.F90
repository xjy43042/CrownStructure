PROGRAM merge_var
 
!====================================================================
! Aggregate regional 5x5° crown structure data (CW/CD/AR) to global 0.5° resolution
!   Take CHavg-based prediction as example 
! Usage: 
!   gfortran -g -mcmodel=large -fbounds-check -o aggr aggr_CWCD_globe.F90 -I/usr/include -lnetcdff -lnetcdf
! Input:  
!   5x5° regional NetCDF files
! Output: 
!   0.5° global composite
!====================================================================
USE netcdf

IMPLICIT NONE

! Grid dimensions
integer , parameter :: nx = 1200 , ny = 1200    ! High-res input grid
integer , parameter :: nxo= 720, nyo= 360      ! Low-res output grid (0.5 degree)
integer , parameter :: r8 = selected_real_kind(12)
real(r8), parameter :: delta = 1./2.           ! Grid spacing
real(r8), parameter :: reso  = 1./2.           ! Output resolution (0.5 degree)

! Data arrays
real(r8), dimension(1200) :: lat, lon, lats, latn, lonw, lone  ! Grid coordinates
real(r8), dimension(1200, 1200) :: area, lc                   ! Area and land cover
real(r8), dimension(1200, 1200, 8) :: CR, CD, AR             ! Crown properties (high-res)
real(r8), dimension(nxo , nyo, 8) :: CR_, CD_, AR_, wgt       ! Aggregated properties (low-res)
real(r8), dimension(nyo) :: slat                              ! Output latitude
real(r8), dimension(nxo) :: slon                               ! Output longitude

! NetCDF variables
integer :: ncid, lat_id, lon_id, CR_id, CD_id, AR_id, lc_id
integer :: lat_dimid, lon_dimid, lc_dimid

! Loop and index variables
integer :: i, j, ix, iy, io, jo, ilc
integer :: ei, ej, si, sj, imon, argn
integer :: reglat,reglon,reglon_,sreglat,sreglon,ereglat,ereglon
integer :: XY3D(3), reg(4)

! File handling variables
character(len=256) :: reg1, reg2, reg3, reg4
character(len=256) :: iyear, cmon, lndname
REAL(r8) :: dx, dy, deg2rad, pi, re
logical  :: fileExists

! Region boundaries (global)
sreglat =   90
ereglat =  -90
sreglon = -180
ereglon =  180

! Get year from command line if provided
argn = IARGC()
IF (argn > 0) THEN
    CALL getarg(1, iyear)
ENDIF

! Constants for coordinate calculations
pi = 4._r8*atan(1.)
deg2rad = pi/180._r8
re = 6.37122e6 * 0.001 ! Earth radius (km -> m)

! Main processing loop over 5x5 degree regions
DO reglat = sreglat, ereglat, -5
    DO reglon_ = sreglon, ereglon, 5

       ! Handle longitude wrapping
       reglon = reglon_
       IF (reglon_ > 180) reglon = reglon_ - 360
       
       ! Define current region boundaries
       reg(1) = reglat      ! North
       reg(2) = reglon      ! West 
       reg(3) = reglat - 5  ! South
       reg(4) = reglon + 5  ! East

       ! Create region identifier strings
       write(reg1, "(i4)") reg(1)
       write(reg2, "(i4)") reg(2)
       write(reg3, "(i4)") reg(3)
       write(reg4, "(i4)") reg(4)

       ! Initialize arrays
       CR(:,:,:) = 0.
       CD(:,:,:) = 0.
       AR(:,:,:) = 0.
       
       ! Build input filename
       lndname = '/tera10/yuanhua/xiangjy/factors_new/make/output/CH90/' // &
                  'RG_' // trim(adjustL(reg1)) // '_' // trim(adjustL(reg2)) // '_' // &
                  trim(adjustL(reg3)) // '_' //trim(adjustL(reg4)) // '.CanopyStructure_15s_CH90.nc'
                  
       ! Check if file exists
       inquire (file=lndname, exist=fileExists)

       IF (fileExists) THEN
          print*, 'Processing '//trim(lndname)

          ! Open and read NetCDF file
          CALL nccheck( nf90_open(trim(lndname), nf90_nowrite, ncid) )
          CALL nccheck( nf90_inq_varid(ncid, "CROWN_WIDTH" , CR_id) )
          CALL nccheck( nf90_get_var  (ncid, CR_id, CR   ) )
          CALL nccheck( nf90_inq_varid(ncid, "CROWN_DEPTH" , CD_id) )
          CALL nccheck( nf90_get_var  (ncid, CD_id, CD   ) )
          CALL nccheck( nf90_inq_varid(ncid, "ASPECT_RATIO" , AR_id) )
          CALL nccheck( nf90_get_var  (ncid, AR_id, AR   ) )
          CALL nccheck( nf90_inq_varid(ncid, "lat" , lat_id) )
          CALL nccheck( nf90_get_var  (ncid, lat_id, lat   ) )
          CALL nccheck( nf90_inq_varid(ncid, "lon" , lon_id) )
          CALL nccheck( nf90_get_var  (ncid, lon_id, lon   ) )
          CALL nccheck( nf90_close(ncid) )

          ! Calculate grid cell boundaries
          latn(:) = lat(:) + delta/2
          lats(:) = lat(:) - delta/2
          lonw(:) = lon(:) - delta/2
          lone(:) = lon(:) + delta/2

          ! Calculate grid cell areas (m²)
          DO i = 1, ny
             dx = (lone(1)-lonw(1))*deg2rad
             dy = sin(latn(i)*deg2rad) - sin(lats(i)*deg2rad)
             area(:,i) = dx*dy*re*re
          ENDDO

          ! Aggregate high-res data to low-res grid
          DO ix = 1, nx
             DO iy = 1, ny
                ! Find output grid indices
                jo = NINT((90.-lat(iy))/reso+0.5)
                io = NINT((lon(ix)+180.)/reso+0.5)
                
                ! Store aggregated values
                CR_(io,jo,:) = CR(ix,iy,:)
                CD_(io,jo,:) = CD(ix,iy,:)
                AR_(io,jo,:) = AR(ix,iy,:)
             ENDDO
          ENDDO
       ENDIF
    ENDDO
ENDDO

! Create output latitude/longitude coordinates
DO ix = 1, nyo
    slat(ix) = 90 - reso*(ix-1) - reso/2
ENDDO
DO iy = 1, nxo
    slon(iy) = -180 + reso*(iy-1) + reso/2
ENDDO

! Create output NetCDF file
lndname = '/tera10/yuanhua/xiangjy/factors_new/make/aggr/Global_CanopyStructure_0.5deg_CH90.nc'
print*, 'Making Global Data'

CALL nccheck( nf90_create (trim(lndname), NF90_NETCDF4, ncid) )

! Define dimensions
CALL nccheck( nf90_def_dim(ncid, "PFT", 8, lc_dimid) )
CALL nccheck( nf90_def_dim(ncid, "lat", nyo, lat_dimid) )
CALL nccheck( nf90_def_dim(ncid, "lon", nxo, lon_dimid) )

! Define coordinate variables
CALL nccheck( nf90_def_var(ncid, "lat", NF90_DOUBLE, lat_dimid, lat_id, deflate_level=6) )
CALL nccheck( nf90_put_att(ncid, lat_id, "long_name", "Latitude") )
CALL nccheck( nf90_put_att(ncid, lat_id, "units", "degrees_north") )

CALL nccheck( nf90_def_var(ncid, "lon", NF90_DOUBLE, lon_dimid, lon_id, deflate_level=6) )
CALL nccheck( nf90_put_att(ncid, lon_id, "long_name", "Longitude") )
CALL nccheck( nf90_put_att(ncid, lon_id, "units", "degrees_east") )

! Define 3D data variables
XY3D = (/lon_dimid, lat_dimid, lc_dimid/)

! Crown width
CALL nccheck( nf90_def_var(ncid, "CROWN_WIDTH", NF90_DOUBLE, XY3D, CR_id, deflate_level=6) )
CALL nccheck( nf90_put_att(ncid, CR_id, "long_name", "crown width predicted by the average canopy-top height") )
CALL nccheck( nf90_put_att(ncid, CR_id, "units", "m") )
CALL nccheck( nf90_put_att(ncid, CR_id, "_FillValue", -999._r8) )

! Crown depth
CALL nccheck( nf90_def_var(ncid, "CROWN_DEPTH", NF90_DOUBLE, XY3D, CD_id, deflate_level=6) )
CALL nccheck( nf90_put_att(ncid, CD_id, "long_name", "crown depth predicted by the average canopy-top height") )
CALL nccheck( nf90_put_att(ncid, CD_id, "units", "m") )
CALL nccheck( nf90_put_att(ncid, CD_id, "_FillValue", -999._r8) )

! Aspect ratio
CALL nccheck( nf90_def_var(ncid, "ASPECT_RATIO", NF90_DOUBLE, XY3D, AR_id, deflate_level=6) )
CALL nccheck( nf90_put_att(ncid, AR_id, "long_name", "crown depth to width ratio predicted by the average canopy-top height") )
CALL nccheck( nf90_put_att(ncid, AR_id, "units", "None") )
CALL nccheck( nf90_put_att(ncid, AR_id, "_FillValue", -999._r8) )

! Write data to file
CALL nccheck( nf90_inq_varid(ncid, "lat", lat_id) )
CALL nccheck( nf90_put_var(ncid, lat_id, slat) )
CALL nccheck( nf90_inq_varid(ncid, "lon", lon_id) )
CALL nccheck( nf90_put_var(ncid, lon_id, slon) )
CALL nccheck( nf90_inq_varid(ncid, "CROWN_WIDTH", CR_id) )
CALL nccheck( nf90_put_var(ncid, CR_id, CR_) )
CALL nccheck( nf90_inq_varid(ncid, "CROWN_DEPTH", CD_id) )
CALL nccheck( nf90_put_var(ncid, CD_id, CD_) )
CALL nccheck( nf90_inq_varid(ncid, "ASPECT_RATIO", AR_id) )
CALL nccheck( nf90_put_var(ncid, AR_id, AR_) )

! Add global attributes
CALL nccheck( nf90_put_att(ncid, NF90_GLOBAL, 'Title', 'Land surface model input canopy morphological structure data') )
CALL nccheck( nf90_put_att(ncid, NF90_GLOBAL, 'resolution', '0.5 degree, 720x360 (lonxlat) global') )
CALL nccheck( nf90_put_att(ncid, NF90_GLOBAL, 'coordinate', 'Geographic, degrees longitude and latitude') )
CALL nccheck( nf90_put_att(ncid, NF90_GLOBAL, 'Authors', 'Yongjiu Dai land group at Sun Yat-sen University') )
CALL nccheck( nf90_put_att(ncid, NF90_GLOBAL, 'Address', 'School of Atmospheric Sciences, Sun Yat-sen University, Zhuhai, China') )
CALL nccheck( nf90_put_att(ncid, NF90_GLOBAL, 'Contact', 'Jiayi Xiang (xiangjy8@mail2.sysu.edu.cn), Hua Yuan (yuanh25@mail.sysu.edu.cn)') )

CALL nccheck( nf90_close(ncid) )
print*, 'Successfully created '//trim(lndname)//'!'

CONTAINS
  ! NetCDF error handling subroutine
  SUBROUTINE nccheck(status)
      INTEGER, INTENT(IN) :: status
      IF (status /= nf90_noerr) THEN
         print *, trim(nf90_strerror(status))
         stop 2
      END IF
   END SUBROUTINE nccheck
END PROGRAM merge_var