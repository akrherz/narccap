PROGRAM interpb

   USE all_io
   USE date_pack
   USE diags
   USE header_data
   USE interp
   USE util
   USE xyll

   IMPLICIT NONE

   !  NAMELIST declarations

   !  RECORD0

   CHARACTER (LEN=132) , DIMENSION(100) :: input_file

   !  RECORD1

   INTEGER :: start_year , start_month , start_day , start_hour ,  start_minute , start_second , start_frac
   INTEGER ::   end_year ,   end_month ,   end_day ,   end_hour ,    end_minute ,   end_second , end_frac
   INTEGER :: interval

   !  RECORD2

   REAL , DIMENSION(1000) :: pressure_bu_no_sfc_Pa

   !  RECORD3

   LOGICAL :: print_info

   !  RECORD4

   CHARACTER (LEN=20) :: uv_interp_method , t_interp_method , moist_interp_method , &
                         height_interp_method , p_interp_method , else_interp_method ,  &
                         uv_extrap_up , t_extrap_up , moist_extrap_up , height_extrap_up , p_extrap_up , else_extrap_up ,  &
                         uv_extrap_low , t_extrap_low , moist_extrap_low , height_extrap_low , p_extrap_low , else_extrap_low
  
   !  Local declarations.

   !  Generic loop counters.

   INTEGER :: i , j , k

   !  Time loop.

   INTEGER :: loop_count , iprocess , itimes , loop

   !  Compute pressure levels, etc.

   INTEGER :: kxp , k_loop
   REAL , ALLOCATABLE , DIMENSION(:) :: newcoord

   !  I/O local variables.

   INTEGER :: input_count 
   INTEGER , PARAMETER :: input_sigma = 50

   !  Domain sizes, specifications from big header.

   INTEGER :: imx , jmx , kxs
   INTEGER :: index , domain
   REAL :: p00 , ts0 , tlp , tiso , ds , ptop

   !  Vertcal interplation variables

   INTEGER :: crsdot

   INTEGER :: int_method
   INTEGER , PARAMETER :: linear=1 , log=2 , kappa=3

   INTEGER :: extrap_up , extrap_low
   INTEGER , PARAMETER :: constant=1 , extrapolate=2

   !  Variable locations in the large arrays

   INTEGER :: index_pp , index_t , index_qv , index_pd , index_pc , index_z , index_rh , index_p
   INTEGER :: index_psc , index_ter , index_psfc , index_psfc2 , index_slp , index_slpc1 , index_slpd1 , index_fat1001
   INTEGER :: index_latd , index_lond , index_pat1001 , index_slpd
   INTEGER :: index_sigma

   !  Output unit info.

   INTEGER , PARAMETER :: output_mmoutp       = 20
   INTEGER , PARAMETER :: output_regridder    = 21
   INTEGER , PARAMETER :: output_intermediate = 22

   !  Timing info.

   REAL :: start_time_loop , end_time_loop , &
           start_time      , end_time

   !  Map info, cross point lat/lon lower left corner point, for use in computing the
   !  dot point lat/lon arrays.

   REAL :: lat1c , lon1c

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   !  Initialize timing for whole routine.

#ifdef F95
   CALL cpu_time ( start_time ) 
#endif

   !  Read in the NAMELIST info.

   CALL do_namelist ( input_file , &
                      start_year , start_month , start_day , start_hour , &
                      start_minute , start_second , start_frac , &
                      end_year ,   end_month ,   end_day ,   end_hour , &
                      end_minute ,   end_second ,   end_frac , &
                      interval , &
                      pressure_bu_no_sfc_Pa , &
                      print_info , &
                      uv_interp_method , t_interp_method , moist_interp_method , &
                      height_interp_method , p_interp_method , else_interp_method ,  &
                      uv_extrap_up , t_extrap_up , moist_extrap_up , height_extrap_up , p_extrap_up , else_extrap_up ,  &
                      uv_extrap_low , t_extrap_low , moist_extrap_low , height_extrap_low , p_extrap_low , else_extrap_low )


   !  The date info was contained in the namelist.  When to start and
   !  when to end, interval.

   WRITE ( start_date ,  &
           '(I4.4,"-",I2.2,"-",I2.2,"_",I2.2,":",I2.2,":",I2.2)') &
           start_year , start_month , start_day , start_hour , &
           start_minute , start_second

   current_date = start_date

   WRITE ( end_date ,  &
           '(I4.4,"-",I2.2,"-",I2.2,"_",I2.2,":",I2.2,":",I2.2)') &
           end_year , end_month , end_day , end_hour , &
           end_minute , end_second

   PRINT '(A)','Time periods to process'
   loop_count = 0
   which_times : DO

      loop_count = loop_count + 1
      IF ( loop_count .GT. 1000 ) THEN
         PRINT '(A)','Seems like a lot of loops.'
         STOP 'lots_of_loops_check_dates'
      END IF

      IF ( current_date .EQ. end_date ) THEN
         PRINT '(A,A,A)','Computed valid time ',current_date,'.  This is the ending time.'
         iprocess = loop_count
         EXIT which_times
      ELSE IF ( current_date .LT. end_date ) THEN
            PRINT '(A,A,A)','Computed valid time ',current_date,'.'
      ELSE
         iprocess = loop_count -1
         EXIT which_times
      END IF

      CALL geth_newdate ( new_date , current_date , interval )
      current_date = new_date

   END DO which_times

   !  Compute the number of pressure layers.

   kxp = 0
   find_pressure : DO k_loop = 1 , 1000
      IF ( pressure_bu_no_sfc_Pa(k_loop) .LT. 0.5 ) THEN
         kxp = k_loop - 1
         EXIT find_pressure
      END IF
   END DO find_pressure
   IF ( kxp .LT. 1 ) THEN
      PRINT '(A)','Error in finding how many pressure levels there are.'
      PRINT *,pressure_bu_no_sfc_Pa
      STOP 'how_many_pressure_levels'
   END IF

   !  ALLOCATE space for the pressure array.  Note that kxp now includes the surface level.

   ALLOCATE ( newcoord(kxp+1) )
   newcoord(1)     = 100100
   newcoord(2:kxp+1) = pressure_bu_no_sfc_Pa(1:kxp)
   kxp = kxp + 1

   !  Check to see if most things are OK with these selected pressure values.

   CALL check ( kxp , newcoord )

   !  This is some more I/O set up.  This OPENS the sigma level
   !  data set with logical unit 50 (input_sigma), and with the name specified
   !  in the NAMELIST file.

   input_count = 1
   CALL get_fg_file ( input_file(input_count) , input_sigma )

   !  Once we have the sigma level data ready to read, we
   !  can procees.  We process the data in the big header first.
   !  The big header data comes out of this routine, and is available
   !  via the USE association for header_data.

   CALL read_bh ( input_sigma )

   !  Some of this big header information is important.  Particularly, we
   !  pick up the horizontal size of the domain.  This allows us to go ahead
   !  and ALLOCATE some space a bit later.

   imx   = bhi(16,1)
   jmx   = bhi(17,1)
   kxs   = bhi(12,bhi(1,1))
   index = bhi(1,1)
   domain= bhi(13,1)

   ds  = bhr( 9,1)
   ptop = bhr(2,2)
   p00 = bhr(2,5)
   ts0 = bhr(3,5)
   tlp = bhr(4,5)
   tiso = bhr(5,5)

   !  Test requested pressure levels against ptop from header.

   IF      ( newcoord(kxp) .LT. ptop ) THEN
      PRINT '(A,F8.2,A)','PTOP from header = ',ptop,' Pa.'
      PRINT '(A,F8.2,A)','Minimum requested pressure = ',newcoord(kxp),' Pa.'
      STOP 'Namelist_press_ptop_mismatch'
   END IF

   !  Some more I/O things.  This sets up the names of the output files
   !  for the INTERPB program.  UNIT=20 is the MMOUTP_DOMAINx file.
   !  The bhi(13,1) is the domain ID, so that we can use a number for
   !  the "x" in DOMAINx.

   CALL open_out_file ( output_mmoutp , output_regridder , domain ) 

   !  Loop over all of the time periods in the pressure level input
   !  data.  We only pick the times that match the requested values
   !  from the NAMELIST file.  The number of vertical pressure levels
   !  is permitted to change between times, so the ALLOCATION of 3D
   !  pressure data occurs when we know how big the vertical component
   !  will be.

   ALLOCATE ( alls_3d(max3d) )
   ALLOCATE ( alls_2d(max2d) )
   ALLOCATE ( alls_1d(max1d) )

   ALLOCATE ( allp_3d(max3d) )
   ALLOCATE ( allp_2d(max2d) )
   ALLOCATE ( allp_1d(max1d) )

   itimes = 0
   current_date = start_date
   time_loop : DO

#ifdef F95
      CALL cpu_time ( start_time_loop ) 
#endif

      !  READ the data for this time.  This is a loop over an unknown
      !  number of variables.  We READ flags, headers and data until we
      !  find the flag that signifies the end of the time.  The data for
      !  the all_3d, all_2d and all_1d is just the input data from the
      !  sigma level data set.  This is available from USE association
      !  of the all_io module.

      DO loop = 1 , max3d
         NULLIFY(alls_3d(loop)%array)
         NULLIFY(allp_3d(loop)%array)
      END DO
      DO loop = 1 , max2d
         NULLIFY(alls_2d(loop)%array)
         NULLIFY(allp_2d(loop)%array)
      END DO
      DO loop = 1 , max1d
         NULLIFY(alls_1d(loop)%array)
         NULLIFY(allp_1d(loop)%array)
      END DO

      bhi(1,1) = index

      CALL read_data ( input_sigma , input_file , input_count , kxp , newcoord , interval )

      !  If this is MM5 input, and the forecast was "dry" (no Qv output), then we need to 
      !  create an artificial Qv field of zero.  This allows us to avoid any additional 
      !  logical tests as to the existence of Qv. 

      IF ( ( index .EQ. 11 ) .AND. ( bhi(3,13) .EQ. 1 ) ) THEN
         CALL add_qv
      END IF

      !  If this is MM5 input, and IFSNOW=0 with slab model, then we need to 
      !  create an artificial snowcovr field of zero.  This allows us to avoid any additional 
      !  logical tests as to the existence of snowcovr. 

![DEH] No thank you, we already have SNOWCOVR
!      IF ( ( index .EQ. 11 ) .AND. ( bhi(16,13) .EQ. 0 ) ) THEN
!         print*,"NO SNOWCOVR field found - add a fake field of 0"
!         CALL add_snowcovr
!      END IF

      !  Where are the important arrays located?  We assume that these locations will stay the same
      !  throughout the time periods.

      IF ( itimes .EQ. 0 ) THEN
         CALL assign_index ( index_pp , index_t , index_qv , &
                             index_psc , index_ter , & 
                             index_sigma ) 
      END IF

      !  There are not too many diagnostics, but there are a few that we must 
      !  compute.

      !  Compute 3d total pressure, and fill in the small header.

      nums_3d = nums_3d + 1
      nump_3d = nump_3d + 1
      index_pc = nums_3d
      ALLOCATE ( alls_3d(nums_3d)%array(imx,jmx,kxs) )
      ALLOCATE ( allp_3d(nump_3d)%array(imx,jmx,kxp) )
      alls_3d(nums_3d)%small_header = alls_3d(nums_3d-1)%small_header
      alls_3d(nums_3d)%small_header%staggering(1:1) = 'C'
      alls_3d(nums_3d)%small_header%name(1:8) = 'P3DCRS  '
      alls_3d(nums_3d)%small_header%units(1:25) = 'Pa                       '
      alls_3d(nums_3d)%small_header%description(1:46) = 'Pressure on cross points                      '
      allp_3d(nump_3d)%small_header = alls_3d(nums_3d)%small_header
      allp_3d(nump_3d)%small_header%ordering(1:4) = 'YXP '
      allp_3d(nump_3d)%small_header%end_dims(3) = kxp
      
      CALL total_pressure ( alls_3d(index_pp)%array , alls_2d(index_psc)%array , alls_1d(index_sigma)%array , &
                            ptop , imx , jmx , kxs , alls_3d(nums_3d)%array ) 

      !  Dot point 3d pressure.

      nums_3d = nums_3d + 1 
      nump_3d = nump_3d + 1
      index_pd = nums_3d
      ALLOCATE ( alls_3d(nums_3d)%array(imx,jmx,kxs) )
      ALLOCATE ( allp_3d(nump_3d)%array(imx,jmx,kxp) )
      alls_3d(nums_3d)%array = alls_3d(nums_3d-1)%array
      alls_3d(nums_3d)%small_header = alls_3d(nums_3d-1)%small_header
      alls_3d(nums_3d)%small_header%staggering(1:1) = 'D'
      alls_3d(nums_3d)%small_header%name(1:8) = 'P3DDOT  '
      alls_3d(nums_3d)%small_header%units(1:25) = 'Pa                       '
      alls_3d(nums_3d)%small_header%description(1:46) = 'Pressure on dot points                        '
      allp_3d(nump_3d)%small_header = alls_3d(nums_3d)%small_header
      allp_3d(nump_3d)%small_header%ordering(1:4) = 'YXP '
      allp_3d(nump_3d)%small_header%end_dims(3) = kxp

      DO k = 1 , kxs
         CALL crs2dot ( alls_3d(nums_3d)%array(:,:,k) , imx , jmx ) 
      END DO

      !  3d height diagnostic, and fill in the small header.

      nums_3d = nums_3d + 1
      nump_3d = nump_3d + 1
      index_z = nums_3d
      ALLOCATE ( alls_3d(nums_3d)%array(imx,jmx,kxs) )
      ALLOCATE ( allp_3d(nump_3d)%array(imx,jmx,kxp) )
      alls_3d(nums_3d)%small_header = alls_3d(nums_3d-1)%small_header
      alls_3d(nums_3d)%small_header%staggering(1:1) = 'C'
      alls_3d(nums_3d)%small_header%name(1:8) = 'H       '
      alls_3d(nums_3d)%small_header%units(1:25) = 'm                        '
      alls_3d(nums_3d)%small_header%description(1:46) = 'Geopotential height                           '
      allp_3d(nump_3d)%small_header = alls_3d(nums_3d)%small_header
      allp_3d(nump_3d)%small_header%ordering(1:4) = 'YXP '
      allp_3d(nump_3d)%small_header%end_dims(3) = kxp

      CALL height ( alls_2d(index_psc)%array , alls_1d(index_sigma)%array , &
                    p00 , ts0 , tlp , tiso , ptop , imx , jmx , kxs , &
                    alls_3d(nums_3d)%array ) 

      !  3d relative humidity on sigma surfaces.

      nums_3d = nums_3d + 1
      nump_3d = nump_3d + 1
      index_rh = nums_3d
      ALLOCATE ( alls_3d(nums_3d)%array(imx,jmx,kxs) )
      ALLOCATE ( allp_3d(nump_3d)%array(imx,jmx,kxp) )
      alls_3d(nums_3d)%small_header = alls_3d(nums_3d-1)%small_header
      alls_3d(nums_3d)%small_header%staggering(1:1) = 'C'
      alls_3d(nums_3d)%small_header%name(1:8) = 'RH      '
      alls_3d(nums_3d)%small_header%units(1:25) = '%                        '
      alls_3d(nums_3d)%small_header%description(1:46) = 'Relative humidity                             '
      allp_3d(nump_3d)%small_header = alls_3d(nums_3d)%small_header
      allp_3d(nump_3d)%small_header%ordering(1:4) = 'YXP '
      allp_3d(nump_3d)%small_header%end_dims(3) = kxp

      CALL qvtorh ( alls_3d(index_qv)%array , alls_3d(index_t)%array , alls_3d(index_pc)%array , & 
                    imx , jmx , kxs , alls_3d(nums_3d)%array ) 

      !  2d surface pressure

!     nums_2d = nums_2d + 1
!     index_psfc = nums_2d
!     ALLOCATE ( alls_2d(nums_2d)%array(imx,jmx) )
!     alls_2d(nums_2d)%small_header = alls_2d(nums_2d-1)%small_header
!     alls_2d(nums_2d)%small_header%staggering(1:1) = 'C'
!     alls_2d(nums_2d)%small_header%name(1:8) = 'PSFC    '
!     alls_2d(nums_2d)%small_header%units(1:25) = 'Pa                       '
!     alls_2d(nums_2d)%small_header%description(1:46) = 'Surface pressure                              '

!     CALL surface_pressure ( alls_3d(index_pp)%array , alls_2d(index_psc)%array , ptop , &
!                             imx , jmx , kxs , alls_2d(nums_2d)%array ) 

      nums_2d = nums_2d + 1
      index_psfc2 = nums_2d
      ALLOCATE ( alls_2d(nums_2d)%array(imx,jmx) )
      alls_2d(nums_2d)%small_header = alls_2d(nums_2d-1)%small_header
      alls_2d(nums_2d)%small_header%staggering(1:1) = 'C'
      alls_2d(nums_2d)%small_header%name(1:8) = 'PSFC    '
      alls_2d(nums_2d)%small_header%units(1:25) = 'Pa                       '
      alls_2d(nums_2d)%small_header%description(1:46) = 'Surface pressure assuming const rho pert      '

      CALL surface_pressure2 (  alls_3d(index_pp)%array , alls_3d(index_t)%array , alls_3d(index_qv)%array , &
                                alls_2d(index_ter)%array , alls_2d(index_psc)%array , &
                                alls_1d(index_sigma)%array , &
                                ptop , p00 , tlp , ts0 , &
                                imx , jmx , kxs , alls_2d(nums_2d)%array ) 

      !  2d sea level pressure

      nums_2d = nums_2d + 1
      index_slp = nums_2d
      ALLOCATE ( alls_2d(nums_2d)%array(imx,jmx) )
      alls_2d(nums_2d)%small_header = alls_2d(nums_2d-1)%small_header
      alls_2d(nums_2d)%small_header%staggering(1:1) = 'C'
      alls_2d(nums_2d)%small_header%name(1:8) = 'PSEALVLC'
      alls_2d(nums_2d)%small_header%units(1:25) = 'Pa                       '
      alls_2d(nums_2d)%small_header%description(1:46) = 'Sea level pressure                            '

!     CALL seapres ( alls_3d(index_z)%array , alls_3d(index_t)%array , alls_3d(index_qv)%array , alls_3d(index_pc)%array , &
!                    alls_2d(index_ter)%array , alls_2d(index_psfc)%array , imx , jmx , kxs , alls_2d(nums_2d)%array ) 
      CALL seapres ( alls_3d(index_z)%array , alls_3d(index_t)%array , alls_3d(index_qv)%array , alls_3d(index_pc)%array , &
                     alls_2d(index_ter)%array , alls_2d(index_psfc2)%array , imx , jmx , kxs , alls_2d(nums_2d)%array ) 

      !  2d sea level pressure, dot point

      nums_2d = nums_2d + 1
      index_slpd = nums_2d
      ALLOCATE ( alls_2d(nums_2d)%array(imx,jmx) )
      alls_2d(nums_2d)%small_header = alls_2d(nums_2d-1)%small_header
      alls_2d(nums_2d)%small_header%staggering(1:1) = 'D'
      alls_2d(nums_2d)%small_header%name(1:8) = 'PSEALVLD'
      alls_2d(nums_2d)%small_header%units(1:25) = 'Pa                       '
      alls_2d(nums_2d)%small_header%description(1:46) = 'Sea level pressure                            '

      alls_2d(index_slpd)%array = alls_2d(index_slp)%array
      CALL crs2dot ( alls_2d(index_slpd)%array , imx , jmx )

      !  2d computational lower pressure: MAX(slp+100, 100100) Pa

      nums_2d = nums_2d + 2
      index_slpc1 = nums_2d-1
      ALLOCATE ( alls_2d(nums_2d-1)%array(imx,jmx) )
      alls_2d(nums_2d-1)%small_header = alls_2d(nums_2d-2)%small_header
      alls_2d(nums_2d-1)%small_header%staggering(1:1) = 'C'
      alls_2d(nums_2d-1)%small_header%name(1:8) = 'SLPC+1  '
      alls_2d(nums_2d-1)%small_header%units(1:25) = 'Pa                       '
      alls_2d(nums_2d-1)%small_header%description(1:46) = 'MAX(( 100100. , slp(i,j) + 100. )      PAT1001'

      index_slpd1 = nums_2d
      ALLOCATE ( alls_2d(nums_2d)%array(imx,jmx) )
      alls_2d(nums_2d)%small_header = alls_2d(nums_2d-2)%small_header
      alls_2d(nums_2d)%small_header%staggering(1:1) = 'D'
      alls_2d(nums_2d)%small_header%name(1:8) = 'SLPD+1  '
      alls_2d(nums_2d)%small_header%units(1:25) = 'Pa                       '
      alls_2d(nums_2d)%small_header%description(1:46) = 'MAX(( 100100. , slp(i,j) + 100. )      PAT1001'

      CALL pat1001 ( alls_2d(index_slp)%array , imx , jmx , alls_2d(nums_2d-1)%array , alls_2d(nums_2d)%array ) 

      !  2d dot point latitude

      nums_2d = nums_2d + 1
      index_latd = nums_2d
      ALLOCATE ( alls_2d(nums_2d)%array(imx,jmx) )
      alls_2d(nums_2d)%small_header = alls_2d(nums_2d-1)%small_header
      alls_2d(nums_2d)%small_header%staggering(1:1) = 'D'
      alls_2d(nums_2d)%small_header%name(1:8) = 'LATITDOT'
      alls_2d(nums_2d)%small_header%units(1:25) = 'degrees                  '
      alls_2d(nums_2d)%small_header%description(1:46) = 'Latitude, + ==> north; - ==> south            '

      !  2d dot point longitude

      nums_2d = nums_2d + 1
      index_lond = nums_2d
      ALLOCATE ( alls_2d(nums_2d)%array(imx,jmx) )
      alls_2d(nums_2d)%small_header = alls_2d(nums_2d-1)%small_header
      alls_2d(nums_2d)%small_header%staggering(1:1) = 'D'
      alls_2d(nums_2d)%small_header%name(1:8) = 'LONGIDOT'
      alls_2d(nums_2d)%small_header%units(1:25) = 'degrees                  '
      alls_2d(nums_2d)%small_header%description(1:46) = 'Longitude, + ==> east; - ==> west             '

      !  We need the lower left corner of the latitude and longitude on cross points
      !  to compute the dot point arrays of latitude and longitude.

      DO loop = 1 , nums_2d
         IF      ( alls_2d(loop)%small_header%name(1:8) .EQ. 'LATITCRS' ) THEN
            lat1c = alls_2d(loop)%array(1,1)
         ELSE IF ( alls_2d(loop)%small_header%name(1:8) .EQ. 'LONGICRS' ) THEN
            lon1c = alls_2d(loop)%array(1,1)
         END IF
      END DO

      CALL lat_lon_dot ( alls_2d(index_latd)%array , alls_2d(index_lond)%array , imx , jmx , &
                         bhr(9,1)/1000. , bhr(3,1) , bhr(5,1) , bhr(6,1) , &
                         bhi(7,1) , lat1c , lon1c )

      nums_2d = nums_2d + 1
      index_fat1001 = nums_2d
      ALLOCATE ( alls_2d(nums_2d)%array(imx,jmx) )

      !  Vertically interpolate each 3d variable.  Later we decide if we really need them.

      DO loop = 1 , nums_3d

         IF      ( alls_3d(loop)%small_header%staggering(1:1) .EQ. 'C' ) THEN
            crsdot = 1
            index_pat1001 = index_slpc1
            index_p = index_pc
         ELSE IF ( alls_3d(loop)%small_header%staggering(1:1) .EQ. 'D' ) THEN
            crsdot = 0
            index_pat1001 = index_slpd1
            index_p = index_pd
         END IF

         IF      ( ( alls_3d(loop)%small_header%name(1:8) .EQ. 'U       ' ) .OR. &
                   ( alls_3d(loop)%small_header%name(1:8) .EQ. 'V       ' ) ) THEN
            IF      ( uv_interp_method .EQ. 'linear in p         ' ) THEN
               int_method = linear
            ELSE IF ( uv_interp_method .EQ. 'linear in log p     ' ) THEN
               int_method = log
            ELSE IF ( uv_interp_method .EQ. 'linear in p**kappa  ' ) THEN
               int_method = kappa
            END IF

            IF      ( uv_extrap_up .EQ. 'constant            ' ) THEN
               extrap_up = constant
            ELSE IF ( uv_extrap_up .EQ. 'extrapolate         ' ) THEN
               extrap_up = extrapolate
            END IF

            IF      ( uv_extrap_low .EQ. 'constant            ' ) THEN
               extrap_low = constant
            ELSE IF ( uv_extrap_low .EQ. 'extrapolate         ' ) THEN
               extrap_low = extrapolate
            END IF
         ELSE IF   ( alls_3d(loop)%small_header%name(1:8) .EQ. 'T       ' ) THEN
            IF      ( t_interp_method .EQ. 'linear in p         ' ) THEN
               int_method = linear
            ELSE IF ( t_interp_method .EQ. 'linear in log p     ' ) THEN
               int_method = log
            ELSE IF ( t_interp_method .EQ. 'linear in p**kappa  ' ) THEN
               int_method = kappa
            END IF

            IF      ( t_extrap_up .EQ. 'constant            ' ) THEN
               extrap_up = constant
            ELSE IF ( t_extrap_up .EQ. 'extrapolate         ' ) THEN
               extrap_up = extrapolate
            END IF

            IF      ( t_extrap_low .EQ. 'constant            ' ) THEN
               extrap_low = constant
            ELSE IF ( t_extrap_low .EQ. 'extrapolate         ' ) THEN
               extrap_low = extrapolate
            END IF
         ELSE IF   ( alls_3d(loop)%small_header%name(1:8) .EQ. 'H       ' ) THEN
            IF      ( height_interp_method .EQ. 'linear in p         ' ) THEN
               int_method = linear
            ELSE IF ( height_interp_method .EQ. 'linear in log p     ' ) THEN
               int_method = log
            ELSE IF ( height_interp_method .EQ. 'linear in p**kappa  ' ) THEN
               int_method = kappa
            END IF

            IF      ( height_extrap_up .EQ. 'constant            ' ) THEN
               extrap_up = constant
            ELSE IF ( height_extrap_up .EQ. 'extrapolate         ' ) THEN
               extrap_up = extrapolate
            END IF

            IF      ( height_extrap_low .EQ. 'constant            ' ) THEN
               extrap_low = constant
            ELSE IF ( height_extrap_low .EQ. 'extrapolate         ' ) THEN
               extrap_low = extrapolate
            END IF
         ELSE IF ( ( alls_3d(loop)%small_header%name(1:8) .EQ. 'Q       ' ) .OR. &
                   ( alls_3d(loop)%small_header%name(1:8) .EQ. 'CLW     ' ) .OR. &
                   ( alls_3d(loop)%small_header%name(1:8) .EQ. 'RNW     ' ) .OR. &
                   ( alls_3d(loop)%small_header%name(1:8) .EQ. 'SNOW    ' ) .OR. &
                   ( alls_3d(loop)%small_header%name(1:8) .EQ. 'GRAUPEL ' ) .OR. &
                   ( alls_3d(loop)%small_header%name(1:8) .EQ. 'ICE     ' ) ) THEN
            IF      ( moist_interp_method .EQ. 'linear in p         ' ) THEN
               int_method = linear
            ELSE IF ( moist_interp_method .EQ. 'linear in log p     ' ) THEN
               int_method = log
            ELSE IF ( moist_interp_method .EQ. 'linear in p**kappa  ' ) THEN
               int_method = kappa
            END IF

            IF      ( moist_extrap_up .EQ. 'constant            ' ) THEN
               extrap_up = constant
            ELSE IF ( moist_extrap_up .EQ. 'extrapolate         ' ) THEN
               extrap_up = extrapolate
            END IF

            IF      ( moist_extrap_low .EQ. 'constant            ' ) THEN
               extrap_low = constant
            ELSE IF ( moist_extrap_low .EQ. 'extrapolate         ' ) THEN
               extrap_low = extrapolate
            END IF
         ELSE IF   ( alls_3d(loop)%small_header%name(1:1) .EQ. 'P'        ) THEN
            IF      ( p_interp_method .EQ. 'linear in p         ' ) THEN
               int_method = linear
            ELSE IF ( p_interp_method .EQ. 'linear in log p     ' ) THEN
               int_method = log
            ELSE IF ( p_interp_method .EQ. 'linear in p**kappa  ' ) THEN
               int_method = kappa
            END IF

            IF      ( p_extrap_up .EQ. 'constant            ' ) THEN
               extrap_up = constant
            ELSE IF ( p_extrap_up .EQ. 'extrapolate         ' ) THEN
               extrap_up = extrapolate
            END IF

            IF      ( p_extrap_low .EQ. 'constant            ' ) THEN
               extrap_low = constant
            ELSE IF ( p_extrap_low .EQ. 'extrapolate         ' ) THEN
               extrap_low = extrapolate
            END IF
         ELSE
            IF      ( else_interp_method .EQ. 'linear in p         ' ) THEN
               int_method = linear
            ELSE IF ( else_interp_method .EQ. 'linear in log p     ' ) THEN
               int_method = log
            ELSE IF ( else_interp_method .EQ. 'linear in p**kappa  ' ) THEN
               int_method = kappa
            END IF

            IF      ( else_extrap_up .EQ. 'constant            ' ) THEN
               extrap_up = constant
            ELSE IF ( else_extrap_up .EQ. 'extrapolate         ' ) THEN
               extrap_up = extrapolate
            END IF

            IF      ( else_extrap_low .EQ. 'constant            ' ) THEN
               extrap_low = constant
            ELSE IF ( else_extrap_low .EQ. 'extrapolate         ' ) THEN
               extrap_low = extrapolate
            END IF
         END IF

         alls_2d(index_fat1001)%small_header = alls_3d(loop)%small_header
         alls_2d(index_fat1001)%small_header%end_dims(3) = 1
         alls_2d(index_fat1001)%small_header%description(40:46) = 'FAT1001'

         CALL fat1001 ( alls_3d(loop)%array(:,:,kxs) , alls_2d(index_pat1001)%array , alls_2d(index_ter)%array , &
                        alls_3d(loop)%small_header%name(1:8) , imx , jmx , crsdot , alls_2d(index_fat1001)%array )

         IF ( print_info ) THEN
            PRINT '(A,A,A,A,3I3,A)','Interpolating 3d field ',alls_3d(loop)%small_header%name(1:9), & 
                                     alls_3d(loop)%small_header%units(1:16),' ',int_method , extrap_up , extrap_low,'.'
         END IF
         CALL vinterp ( alls_3d(loop)%array , imx , jmx , kxs , crsdot , alls_3d(index_p)%array , &
                        alls_2d(index_pat1001)%array , alls_2d(index_fat1001)%array , &
                        allp_3d(loop)%array , kxp , int_method , extrap_up , extrap_low , print_info , newcoord )
       
         !  Reset lowest level (surface) height to terrain.

         IF ( alls_3d(loop)%small_header%name(1:8) .EQ. 'H       ' ) THEN
            allp_3d(loop)%array(:,:,1) = alls_2d(index_ter)%array
         END IF

      END DO

      !  Output the data in the V3 format, interpolated model output.

      CALL out_mmoutp       ( current_date(1:19) , imx , jmx , kxp , output_mmoutp       , &
                              itimes + 1 , interval , &
                              start_year , start_month , start_day , start_hour , &
                              start_minute , start_second , start_frac , &
                              newcoord , print_info )

      !  Output the data in the regridder V3 format.

      CALL out_regridder    ( current_date(1:19) , imx , jmx , kxp , output_regridder    , &
                              itimes + 1 , interval , ptop , &
                              start_year , start_month , start_day , start_hour , &
                              start_minute , start_second , start_frac , &
                              newcoord , print_info )

      !  Output the data in the pregrid intermediate format.

      CALL out_intermediate ( current_date(1:19) , imx , jmx , kxp , output_intermediate , &
                              itimes + 1 , interval , &
                              start_year , start_month , start_day , start_hour , &
                              start_minute , start_second , start_frac , &
                              newcoord , print_info )

      !  Clean up to do every time loop.

      DO loop = 1 , nums_3d
         DEALLOCATE(alls_3d(loop)%array)
      END DO
      DO loop = 1 , nums_2d
         DEALLOCATE(alls_2d(loop)%array)
      END DO
      DO loop = 1 , nums_1d
         DEALLOCATE(alls_1d(loop)%array)
      END DO

      DO loop = 1 , nump_3d
         DEALLOCATE(allp_3d(loop)%array)
      END DO
      DO loop = 1 , nump_2d
         DEALLOCATE(allp_2d(loop)%array)
      END DO
      DO loop = 1 , nump_1d
         DEALLOCATE(allp_1d(loop)%array)
      END DO

      !  Increment the times.

      itimes = itimes + 1

      IF ( current_date(1:19) .EQ. end_date(1:19) ) THEN
         CLOSE ( output_mmoutp    )
         CLOSE ( output_regridder )
         PRINT '(A)','Completed all requested dates.  Exiting program.'
         STOP 'Successful completion 99999'
      END IF

      !  Get the next time to process.

      CALL geth_newdate ( new_date , current_date , interval )
      IF ( print_info ) THEN
         PRINT '(A,A,A,A,A)','Finished with processing ',current_date,', starting with ',new_date,'.'
      END IF
      current_date = new_date

#ifdef F95
      CALL cpu_time ( end_time_loop )

      IF ( print_info ) THEN
         PRINT '(A,F6.2,A)','Time to process this loop is ',end_time_loop - start_time_loop,' s.'
      END IF
#endif

   END DO time_loop

   !  This is the "do it once section".

   CLOSE ( output_mmoutp    )
   CLOSE ( output_regridder )

END PROGRAM interpb
