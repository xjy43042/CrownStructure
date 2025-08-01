Program Aggre

 USE omp_lib
 USE netcdf

!-------------------------------------
! Program to aggregate high-resolution Canopy Height data to 500m resolution
! Usage: 
!   gfortran -mcmodel=large -fbounds-check -g -o aggr aggr_CH_500m.F90 -I/usr/include -lnetcdf -lnetcdff
!   ./aggr 40 -5 35 0
! Input: 
!   Region coordinates (lat/lon bounds)
! Output:
!   Aggregated canopy height statistics at 500m resolution
!-------------------------------------
! Precision and array size parameters
INTEGER , parameter :: r8   = selected_real_kind(12)  ! Double precision
INTEGER , parameter :: nlat = 1200 , nlon = 1200      ! Output grid dimensions (500m)
INTEGER , parameter :: glat = 18554 , glon = 18554    ! Input GFCC grid dimensions (30m)

! Resolution parameters (in degrees)
REAL(r8), parameter :: gdelta = 0.0002694837         ! GFCC 30m resolution
REAL(r8), parameter :: hdelta = 0.00025              ! Target 30m resolution
REAL(r8), parameter :: delta = 1._r8/240._r8         ! Output 500m resolution
REAL    , parameter :: spval = -999.                 ! Missing value flag

! Data arrays
REAL(r8), allocatable, dimension(:)   :: lat, lon, nclat, nclon
INTEGER , allocatable, dimension(:,:) :: cntdata      ! Count of valid pixels
INTEGER , allocatable, dimension(:,:) :: sumdata      ! Sum of values for averaging
REAL(r8), allocatable, dimension(:,:) :: outdata      ! Output mean values
INTEGER , allocatable, dimension(:,:) :: cntdata1     ! Secondary counter
INTEGER , allocatable, dimension(:,:) :: sumdata1     ! Secondary sum
REAL(r8), allocatable, dimension(:,:) :: outdata1      ! Secondary output
REAL(r8), allocatable, dimension(:,:) :: outdata2     ! 98th percentile
REAL(r8), allocatable, dimension(:,:) :: outdata3     ! 95th percentile 
REAL(r8), allocatable, dimension(:,:) :: outdata4     ! 90th percentile
REAL(r8), allocatable, dimension(:,:) :: outdata5     ! Maximum value
REAL(r8) :: gfcc(glon,glat), lon1(glon), lat1(glat)  ! GFCC input data

! Index and loop variables
INTEGER  :: argn                                     ! Argument counter
INTEGER(kind=4), allocatable, dimension(:,:) :: band1 ! Input data band
INTEGER(kind=4), allocatable, dimension(:,:) :: max_value ! Max values
INTEGER(kind=4), allocatable, dimension(:,:,:) :: validata ! Valid data storage
INTEGER(kind=4), allocatable, dimension(:) :: temp_data   ! Temporary array

! Geographic calculation variables
REAL(r8) :: dlat, dlon, clat, clon
INTEGER  :: tlat, tlon  ! Target grid indices
INTEGER  :: nyo , nxo   ! Input dimensions
INTEGER  :: xid , yid   ! GFCC grid indices
INTEGER  :: i   , j     ! Loop counters

! NetCDF file handles
INTEGER  :: ncid, tree_id, lat_id, lon_id, nid, nclat_id, nclon_id, band_id
INTEGER  :: ncid1, lat_id1, lon_id1, gfcc_id ! GFCC file handles
INTEGER  :: xx  , yy      ! Dimension IDs
INTEGER  :: index90, index95, index98  ! Percentile indices

! Region parameters from command line
CHARACTER (len=6)  :: creg1, creg2, creg3, creg4
CHARACTER(len=6) :: str_reg1, str_reg2, str_reg3, str_reg4
INTEGER :: reg1, reg2, reg3, reg4  ! Region bounds
LOGICAL :: ncExists, gfccExists    ! File existence flags

! File paths
CHARACTER(len=255) :: ncfile, filename, infile
CHARACTER(len=255) :: outdir= '/tera10/yuanhua/xiangjy/CHM_Tolan/mkGlobe/new/aggr/RG_5x5/'

 ! Get command line arguments
 argn = IARGC()
 IF (argn > 1) THEN
    CALL getarg(1, creg1)
    CALL getarg(2, creg2)
    CALL getarg(3, creg3)
    CALL getarg(4, creg4)
    read(creg1,*) reg1
    read(creg2,*) reg2
    read(creg3,*) reg3
    read(creg4,*) reg4
 ENDIF
 print*,'reg1, lat_max = ', reg1
 print*,'reg2, lon_min = ', reg2
 print*,'reg3, lat_min = ', reg3
 print*,'reg4, lon_max = ', reg4

 ! Allocate and initialize arrays
 ALLOCATE( nclat(nlat) )
 ALLOCATE( nclon(nlon) )

 ALLOCATE( sumdata(nlon,nlat) )
 ALLOCATE( cntdata(nlon,nlat) )
 ALLOCATE( outdata(nlon,nlat) )
 ALLOCATE( sumdata1(nlon,nlat) )
 ALLOCATE( cntdata1(nlon,nlat) )
 ALLOCATE( outdata1(nlon,nlat) )
 ALLOCATE( outdata2(nlon,nlat) )
 ALLOCATE( outdata3(nlon,nlat) )
 ALLOCATE( outdata4(nlon,nlat) )
 ALLOCATE(max_value(nlon, nlat))
 ALLOCATE(validata(nlon, nlat, 6000))

 ! Initialize arrays
 validata(:,:,:) = -999.
 sumdata(:,:) = 0
 cntdata(:,:) = 0
 outdata(:,:) = -999.
 max_value(:,:) = -999.

 sumdata1(:,:) = 0
 cntdata1(:,:) = 0
 outdata1(:,:) = -999.
 outdata2(:,:) = -999.
 outdata3(:,:) = -999.
 outdata4(:,:) = -999.

 ! Calculate output grid coordinates
 nclat(1) = delta/2 + reg3 ! First latitude center point: 0.5delta+latmin
 DO i=2,nlat,1
    nclat(i) = (i-1)*delta + nclat(1)
 ENDDO

 nclon(1) = delta/2 + reg2
 DO i=2,nlon,1
    nclon(i) = (i-1)*delta + nclon(1)
 ENDDO

 ! Load GFCC 30m reference data
 WRITE(str_reg1, '(I0)') reg1
 WRITE(str_reg2, '(I0)') reg2
 WRITE(str_reg3, '(I0)') reg3
 WRITE(str_reg4, '(I0)') reg4
 
 filename = '/stu01/linwy20/MKSRF30/GFCC/GFCC30mFill/RG_'&
         //TRIM(adjustL(str_reg1))//'_' &
         //TRIM(adjustL(str_reg2))//'_' &
         //TRIM(adjustL(str_reg3))//'_' &
         //TRIM(adjustL(str_reg4))//'.' &
         //'GFCC30mFill_2015.nc'
 print*, 'Loading GFCC file: ', filename
 inquire (file=filename, exist=gfccExists)
 IF (gfccExists) THEN
    print*, 'gfcc file exist!'
 ENDIF
 CALL check( nf90_open(filename, NF90_NOWRITE, ncid1  ) )
 ! Read GFCC data
 CALL check( nf90_inq_varid(ncid1, 'lat' , lat_id1) )
 CALL check( nf90_get_var  (ncid1, lat_id1, lat1   ) )
 CALL check( nf90_inq_varid(ncid1, 'lon' , lon_id1) )
 CALL check( nf90_get_var  (ncid1, lon_id1, lon1   ) )
 CALL check( nf90_inq_varid(ncid1, 'PCT_Tree', gfcc_id) )
 CALL check( nf90_get_var  (ncid1, gfcc_id   , gfcc   ) )
 CALL check( nf90_close(ncid1) )

 ! Process input canopy height files
 infile = '/tera10/yuanhua/xiangjy/CHM_Tolan/mkGlobe/new/aggr/CHtxt/RG_'&
 //TRIM(adjustL(str_reg1))//'_' &
 //TRIM(adjustL(str_reg2))//'_' &
 //TRIM(adjustL(str_reg3))//'_' &
 //TRIM(adjustL(str_reg4))//'_' &
 //'CHnc.txt'

 OPEN(11, file=trim(infile))

 DO WHILE(.true.)
    READ(11,*,end=100) ncfile                     
    !  PRINT*, ncfile
    ! Open canopy height file
    CALL check( nf90_open('/tera10/yuanhua/xiangjy/CHM_Tolan/mkGlobe/new/nc_30mmax/'& 
    //trim(ncfile), NF90_NOWRITE, nid  ) )

    ! Get dimensions
    CALL check( nf90_inq_dimid        (nid, 'lat' , lat_id ) )
    CALL check( nf90_inquire_dimension(nid, lat_id, len=nyo) )
    CALL check( nf90_inq_dimid        (nid, 'lon' , lon_id ) )
    CALL check( nf90_inquire_dimension(nid, lon_id, len=nxo) )

    ALLOCATE( lat(nyo) )
    ALLOCATE( lon(nxo) )
    ALLOCATE( band1(nxo,nyo) )

    CALL check( nf90_inq_varid(nid, 'lat' , lat_id) )
    CALL check( nf90_get_var  (nid, lat_id, lat   ) )
    CALL check( nf90_inq_varid(nid, 'lon' , lon_id) )
    CALL check( nf90_get_var  (nid, lon_id, lon   ) )
    CALL check( nf90_inq_varid(nid, 'Band1', band_id) )
    CALL check( nf90_get_var  (nid, band_id, band1  ) )
    CALL check( nf90_close(nid) )

    ! Process each pixel
    DO i=1,nxo
       DO j=1,nyo
            ! Calculate pixel center coordinates
            clat = lat(j) + hdelta/2.
            clon = lon(i) + hdelta/2.

            ! Find target 500m grid cell
            tlat = NINT((clat-reg3)/delta+0.5)
            tlon = NINT((clon-reg2)/delta+0.5)
            ! Find corresponding CH/GFCC 30m grid cell
            yid = NINT((reg1-clat)/gdelta+0.5)
            xid = NINT((clon-reg2)/gdelta+0.5)

            ! Skip pixels outside region bounds
            IF (clat<reg3 .or. clat>reg1 .or. clon<reg2 .or. clon>reg4) THEN
               CYCLE
            ELSE
               ! Only process valid pixels (GFCC tree cover > 0 and height between 2-255m)
               IF (gfcc(xid,yid)>0 .and. band1(i,j)>=2 .and. band1(i,j)<=255) THEN ! =0是空值，不需要求平均
                  ! Weighted sum using GFCC tree cover as weights
                  sumdata(tlon,tlat) = sumdata(tlon,tlat) + band1(i,j)*gfcc(xid,yid)
                  cntdata(tlon,tlat) = cntdata(tlon,tlat) + 1*gfcc(xid,yid)
 
                  ! Unweighted sum for secondary statistics
                  sumdata1(tlon,tlat) = sumdata1(tlon,tlat) + band1(i,j)
                  cntdata1(tlon,tlat) = cntdata1(tlon,tlat) + 1
 
                  ! Track maximum value
                  max_value(tlon,tlat) = MAX(max_value(tlon,tlat), band1(i,j))
                  
                  ! Store values for percentile calculation
                  validata(tlon,tlat,cntdata1(tlon,tlat)) = band1(i,j)
               ENDIF
            ENDIF
       ENDDO
    ENDDO
   !print*, "end"
    
   ! Clean up
   DEALLOCATE( lat   )
   DEALLOCATE( lon   )
   DEALLOCATE( band1 )
ENDDO

100 close(11)

! Calculate final statistics
DO i=1,nlon,1
   DO j=1,nlat,1
      ! Calculate weighted mean (GFCC tree cover weighted)
      IF (cntdata(i,j) .gt. 0) THEN
         outdata(i,j) = sumdata(i,j)*1. / cntdata(i,j)
      ENDIF

      ! Calculate unweighted statistics
      IF (cntdata1(i,j) .gt. 0) THEN 
         ! Unweighted mean
         outdata1(i,j) = sumdata1(i,j)*1. / cntdata1(i,j)   

         ! Calculate percentiles
         ALLOCATE(temp_data(cntdata1(i,j)))
         temp_data = validata(i,j,1:cntdata1(i,j))
         CALL quicksort_int32(temp_data)    

         index98 = NINT(0.98 * REAL(cntdata1(i,j))) ! 98th percentile
         index95 = NINT(0.95 * REAL(cntdata1(i,j))) ! 95th percentile
         index90 = NINT(0.90 * REAL(cntdata1(i,j))) ! 90th percentile
         outdata2(i,j) = temp_data(index98)
         outdata3(i,j) = temp_data(index95)
         outdata4(i,j) = temp_data(index90)
         DEALLOCATE(temp_data)
      ENDIF
   ENDDO
ENDDO 

 ! Write output NetCDF file
 CALL check( nf90_create(trim(outdir)//'RG_'&
 //TRIM(adjustL(str_reg1))//'_' &
 //TRIM(adjustL(str_reg2))//'_' &
 //TRIM(adjustL(str_reg3))//'_' &
 //TRIM(adjustL(str_reg4))//'.' &
 //'CanopyHeight_15s.nc', NF90_NETCDF4, ncid) )


 CALL check( nf90_def_dim(ncid, 'lat', nlat, xx) )
 CALL check( nf90_def_dim(ncid, 'lon', nlon, yy) )

 ! Define variables and attributes
 CALL check( nf90_def_var(ncid, 'lat'   , NF90_DOUBLE    , (/xx/), nclat_id, deflate_level=6) )
 CALL check( nf90_put_att(ncid, nclat_id, 'standard_name', 'latitude'     ) )
 CALL check( nf90_put_att(ncid, nclat_id, 'long_name'    , 'latitude'     ) )
 CALL check( nf90_put_att(ncid, nclat_id, 'units'        , 'degrees_north') )

 CALL check( nf90_def_var(ncid, 'lon'   , NF90_DOUBLE    , (/yy/), nclon_id, deflate_level=6) )
 CALL check( nf90_put_att(ncid, nclon_id, 'standard_name', 'longitude'   ) )
 CALL check( nf90_put_att(ncid, nclon_id, 'long_name'    , 'longitude'   ) )
 CALL check( nf90_put_att(ncid, nclon_id, 'units'        , 'degrees_east') )

 CALL check( nf90_def_var(ncid, 'CH_gfccmean', NF90_FLOAT  , (/yy,xx/), tree_id, deflate_level=6) )
 CALL check( nf90_put_att(ncid, tree_id   , 'long_name' , &
            'canopy height more than 2m, weighted with PCT_Tree'    ) )
 CALL check( nf90_put_att(ncid, tree_id   , 'units'     , 'm'               ) )
 CALL check( nf90_put_att(ncid, tree_id   , '_FillValue', spval             ) )

 CALL check( nf90_def_var(ncid, 'CH_ge2mean', NF90_FLOAT  , (/yy,xx/), tree_id, deflate_level=6) )
 CALL check( nf90_put_att(ncid, tree_id   , 'long_name' , &
            'canopy height more than 2m, with PCT_Tree greater than 0'    ) )
 CALL check( nf90_put_att(ncid, tree_id   , 'units'     , 'm'               ) )
 CALL check( nf90_put_att(ncid, tree_id   , '_FillValue', spval             ) )

 CALL check( nf90_def_var(ncid, 'CH_ge2max', NF90_FLOAT  , (/yy,xx/), tree_id, deflate_level=6) )
 CALL check( nf90_put_att(ncid, tree_id   , 'long_name' , &
            'the maximum canopy height of the 500m grid, with PCT_Tree greater than 0'    ) )
 CALL check( nf90_put_att(ncid, tree_id   , 'units'     , 'm'               ) )
 CALL check( nf90_put_att(ncid, tree_id   , '_FillValue', spval             ) )

 CALL check( nf90_def_var(ncid, 'CH_98', NF90_FLOAT  , (/yy,xx/), tree_id, deflate_level=6) )
 CALL check( nf90_put_att(ncid, tree_id   , 'long_name' , &
            'the 98% maximum canopy height of the 500m grid, with PCT_Tree greater than 0'    ) )
 CALL check( nf90_put_att(ncid, tree_id   , 'units'     , 'm'               ) )
 CALL check( nf90_put_att(ncid, tree_id   , '_FillValue', spval             ) )

 CALL check( nf90_def_var(ncid, 'CH_95', NF90_FLOAT  , (/yy,xx/), tree_id, deflate_level=6) )
 CALL check( nf90_put_att(ncid, tree_id   , 'long_name' , &
             'the 95% maximum canopy height of the 500m grid, with PCT_Tree greater than 0'    ) )
 CALL check( nf90_put_att(ncid, tree_id   , 'units'     , 'm'               ) )
 CALL check( nf90_put_att(ncid, tree_id   , '_FillValue', spval             ) )

 CALL check( nf90_def_var(ncid, 'CH_90', NF90_FLOAT  , (/yy,xx/), tree_id, deflate_level=6) )
 CALL check( nf90_put_att(ncid, tree_id   , 'long_name' , &
            'the 90% maximum canopy height of the 500m grid, with PCT_Tree greater than 0'    ) )
 CALL check( nf90_put_att(ncid, tree_id   , 'units'     , 'm'               ) )
 CALL check( nf90_put_att(ncid, tree_id   , '_FillValue', spval             ) )

 CALL check( nf90_enddef(ncid) )
 
 ! Write data to file
 CALL check( nf90_inq_varid(ncid, 'lat' , nclat_id) )
 CALL check( nf90_put_var  (ncid, nclat_id, nclat ) )

 CALL check( nf90_inq_varid(ncid, 'lon' , nclon_id) )
 CALL check( nf90_put_var  (ncid, nclon_id, nclon ) )

 CALL check( nf90_inq_varid(ncid, 'CH_gfccmean', tree_id) )
 CALL check( nf90_put_var  (ncid, tree_id   , outdata) )

 CALL check( nf90_inq_varid(ncid, 'CH_ge2mean', tree_id) )
 CALL check( nf90_put_var  (ncid, tree_id   , outdata1) )

 CALL check( nf90_inq_varid(ncid, 'CH_ge2max', tree_id) )
 CALL check( nf90_put_var  (ncid, tree_id   , max_value) )

 CALL check( nf90_inq_varid(ncid, 'CH_98', tree_id) )
 CALL check( nf90_put_var  (ncid, tree_id   , outdata2) )

 CALL check( nf90_inq_varid(ncid, 'CH_95', tree_id) )
 CALL check( nf90_put_var  (ncid, tree_id   , outdata3) )

 CALL check( nf90_inq_varid(ncid, 'CH_90', tree_id) )
 CALL check( nf90_put_var  (ncid, tree_id   , outdata4) )

 CALL check( nf90_put_att(ncid, NF90_GLOBAL, 'Title'  , &
                     'Canopy structure prediction model input average canopy-top height data') )
 CALL check( nf90_put_att(ncid, NF90_GLOBAL, 'resolution'  , '500m, 1200x1200 (lonxlat) regional grid') )
 CALL check( nf90_put_att(ncid, NF90_GLOBAL, 'coordinate'  , 'Geographic, degrees longitude and latitude') )
 CALL check( nf90_put_att(ncid, NF90_GLOBAL, 'source: ', 'Tolan et al., 2024 (10.1016/j.rse.2023.113888)') )
 CALL check( nf90_put_att(ncid, NF90_GLOBAL, 'Authors', 'Yongjiu Dai land group at Sun Yat-sen University' ) )
 CALL check( nf90_put_att(ncid, NF90_GLOBAL, 'Address', 'School of Atmospheric Sciences, &
                    Sun Yat-sen University, Zhuhai, China') )
 CALL check( nf90_put_att(ncid, NF90_GLOBAL, 'Contact', 'Jiayi Xiang (xiangjy8@mail2.sysu.edu.cn), &
 Hua Yuan (yuanh25@mail.sysu.edu.cn)') )
 CALL check( nf90_close(ncid) )

 ! Clean up memory
 DEALLOCATE(nclat  )
 DEALLOCATE(nclon  )
 DEALLOCATE(sumdata)
 DEALLOCATE(cntdata)
 DEALLOCATE(outdata)
 DEALLOCATE(sumdata1)
 DEALLOCATE(cntdata1)
 DEALLOCATE(outdata2)
 DEALLOCATE(outdata3)
 DEALLOCATE(outdata4)
 DEALLOCATE(max_value)
 DEALLOCATE(validata)

 contains
 ! NetCDF error checking subroutine
 SUBROUTINE check(status)
    INTEGER, intent (in) :: status

    IF(status /= nf90_noerr) THEN
      print *, trim(nf90_strerror(status))
      stop "Stopped"
    END IF
 END SUBROUTINE check

 ! Quick sort implementation for percentile calculation
 recursive SUBROUTINE quicksort_int32 (A)
     IMPLICIT NONE
 
     INTEGER, intent(inout) :: A(:)
 
     ! Local variables
     INTEGER:: nA      
     INTEGER :: left, right
     INTEGER :: pivot
     INTEGER :: marker
     INTEGER :: itemp
     
     nA = SIZE(A)
     IF (nA > 1) THEN
 
     pivot = A(nA/2)
     left  = 0
     right = nA + 1
 
     DO WHILE (left < right)
         right = right - 1
         DO WHILE (A(right) > pivot)
             right = right - 1
         ENDDO
 
         left = left + 1
         DO WHILE (A(left) < pivot)
             left = left + 1
         ENDDO
 
         IF (left < right) THEN
             itemp    = A(left)
             A(left)  = A(right)
             A(right) = itemp
         ENDIF
     ENDDO
 
     marker = right
 
     CALL quicksort_int32 (A(1:marker))
     CALL quicksort_int32 (A(marker+1:nA))
 
     ENDIF
 
 END SUBROUTINE quicksort_int32
 
END PROGRAM Aggre
