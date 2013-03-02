MODULE all_io

   USE date_pack
   USE header_data
   USE util

   TYPE input_fields_3d
      REAL , POINTER , DIMENSION(:,:,:) :: array
      TYPE(sh)                          :: small_header
   END TYPE input_fields_3d

   TYPE input_fields_2d
      REAL , POINTER , DIMENSION(:,:)   :: array
      TYPE(sh)                          :: small_header
   END TYPE input_fields_2d

   TYPE input_fields_1d
      REAL , POINTER , DIMENSION(:)     :: array
      TYPE(sh)                          :: small_header
   END TYPE input_fields_1d

   TYPE(input_fields_3d) , ALLOCATABLE , DIMENSION(:) :: alls_3d
   TYPE(input_fields_2d) , ALLOCATABLE , DIMENSION(:) :: alls_2d
   TYPE(input_fields_1d) , ALLOCATABLE , DIMENSION(:) :: alls_1d

   TYPE(input_fields_3d) , ALLOCATABLE , DIMENSION(:) :: allp_3d
   TYPE(input_fields_2d) , ALLOCATABLE , DIMENSION(:) :: allp_2d
   TYPE(input_fields_1d) , ALLOCATABLE , DIMENSION(:) :: allp_1d

   INTEGER , PARAMETER :: max3d =  30 , &
                          max2d = 100 , &
                          max1d =  10

   INTEGER :: nums_3d , nums_2d , nums_1d
   INTEGER :: nump_3d , nump_2d , nump_1d

   INTEGER :: iprog_no_old

   CHARACTER(LEN=19) :: sh_date

   REAL , DIMENSION(:,:) , ALLOCATABLE :: tseasfc_sum , sfct_sum

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE assign_index ( index_pp , index_t , index_qv , &
                             index_psc , index_ter , &
                             index_sigma )

      IMPLICIT NONE

      INTEGER , INTENT(OUT) :: index_pp , index_t , index_qv , &
                               index_psc , index_ter ,  &
                               index_sigma

      INTEGER :: loop

      index_pp = -1
      index_t  = -1
      index_qv = -1

      DO loop = 1 , nums_3d
         IF ( alls_3d(loop)%small_header%name(1:8) .EQ. 'PP      ' ) THEN
            index_pp = loop
            EXIT
         END IF
      END DO 

      DO loop = 1 , nums_3d
         IF ( alls_3d(loop)%small_header%name(1:8) .EQ. 'T       ' ) THEN
            index_t  = loop
            EXIT
         END IF
      END DO 

      DO loop = 1 , nums_3d
         IF ( alls_3d(loop)%small_header%name(1:8) .EQ. 'Q       ' ) THEN
            index_qv = loop
            EXIT
         END IF
      END DO 

      index_psc = -1
      index_ter = -1

      DO loop = 1 , nums_2d
         IF ( alls_2d(loop)%small_header%name(1:8) .EQ. 'PSTARCRS' ) THEN
            index_psc = loop
            EXIT
         END IF
      END DO 

      DO loop = 1 , nums_2d
         IF ( alls_2d(loop)%small_header%name(1:8) .EQ. 'TERRAIN ' ) THEN
            index_ter = loop
            EXIT
         END IF
      END DO 

      index_sigma = -1

      DO loop = 1 , nums_1d
         IF ( alls_1d(loop)%small_header%name(1:8) .EQ. 'SIGMAH  ' ) THEN
            index_sigma = loop
            EXIT
         END IF
      END DO 

   END SUBROUTINE assign_index

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE read_data ( unit_fg , file_fg , count_fg , kx , newcoord , interval )

      IMPLICIT NONE

      !  Input variables.

      INTEGER , INTENT(IN) :: unit_fg , kx , interval
      INTEGER , INTENT(INOUT) :: count_fg 
      CHARACTER(LEN=132) , DIMENSION(100) , INTENT(IN) :: file_fg
      REAL , DIMENSION(kx) :: newcoord

      !  Local variables.

      INTEGER :: loop
      INTEGER :: loop3 , loop2 , loop1
      INTEGER :: flag , ok , bhi_i1
      REAL    :: dummy 
      INTEGER :: dts

      !  Loop over all of the variables.

      nums_3d = 0
      nums_2d = 0
      nums_1d = 0
      nump_3d = 0
      nump_2d = 0
      nump_1d = 0

      small_header%current_date = '1970-01-01_00:00:00.0000'
      var_loop : DO

         !  Start off with a flag, as ususal.  We are here after a big header
         !  flag, so we can just look for the small header without trying to 
         !  do too many error tests for which flag is expected.
   
         READ ( unit_fg , IOSTAT = ok ) flag
   
         IF      ( ( ok .LT. 0 ) .AND. ( count_fg .EQ. 100 ) ) THEN
            PRINT '(A,I8,A)','Unexpected EOF in forecast data: ',ok,'.'
            PRINT '(A)','Exhausted list of input file names in the namelist.input file.'
            PRINT '(A)','You may have asked for a time period that is not available.'
            STOP 'EOF_FG_read_all_files'
         ELSE IF ( ( ok .LT. 0 ) .AND. ( LEN_TRIM ( file_fg(count_fg+1) ) .EQ. 0 ) ) THEN
            PRINT '(A,I8,A)','Unexpected EOF in forecast data: ',ok,'.'
            PRINT '(A)','No subsequent file specified for input in the namelist.input file.'
            PRINT '(A)','You may have asked for a time period that is not available.'
            STOP 'EOF_in_FG_data'
         ELSE IF ( ( ok .LT. 0 ) .AND. ( LEN_TRIM ( file_fg(count_fg+1) ) .NE. 0 ) ) THEN
            PRINT '(A,I8,A)','Unexpected EOF in forecast data: ',ok,'.'
            PRINT '(A,A,A)','Moving to next file in the namelist.input file: ',TRIM ( file_fg(count_fg+1) ) ,'.'
            CLOSE ( unit_fg )
            count_fg = count_fg + 1
            CALL get_fg_file ( file_fg(count_fg) ,  unit_fg )
            DO loop = 1 , nums_3d
               DEALLOCATE ( alls_3d(loop)%array )
            END DO
            DO loop = 1 , nums_2d
               DEALLOCATE ( alls_2d(loop)%array )
            END DO
            DO loop = 1 , nums_1d
               DEALLOCATE ( alls_1d(loop)%array )
            END DO
            DO loop = 1 , nump_3d
               DEALLOCATE ( allp_3d(loop)%array )
            END DO
            DO loop = 1 , nump_2d
               DEALLOCATE ( allp_2d(loop)%array )
            END DO
            DO loop = 1 , nump_1d
               DEALLOCATE ( allp_1d(loop)%array )
            END DO
            nums_3d = 0
            nums_2d = 0
            nums_1d = 0
            nump_3d = 0
            nump_2d = 0
            nump_1d = 0
            CYCLE var_loop
         END IF

         IF ( flag .EQ. bh_flag ) THEN
            PRINT '(A,I1,A)','Wrong flag, wanted small header flag: ',flag,'.'
            PRINT '(A)','Assuming that this is an instance of multiple input files or cat''ed first-guess input files.'
            READ ( unit_fg ) bhi_i1
            IF ( ( bhi_i1 .EQ. 5 ) .OR. ( bhi_i1 .EQ. 11 ) ) THEN
               DO loop = 1 , nums_3d
                  DEALLOCATE ( alls_3d(loop)%array )
               END DO
               DO loop = 1 , nums_2d
                  DEALLOCATE ( alls_2d(loop)%array )
               END DO
               DO loop = 1 , nums_1d
                  DEALLOCATE ( alls_1d(loop)%array )
               END DO
               DO loop = 1 , nump_3d
                  DEALLOCATE ( allp_3d(loop)%array )
               END DO
               DO loop = 1 , nump_2d
                  DEALLOCATE ( allp_2d(loop)%array )
               END DO
               DO loop = 1 , nump_1d
                  DEALLOCATE ( allp_1d(loop)%array )
               END DO
               nums_3d = 0
               nums_2d = 0
               nums_1d = 0
               nump_3d = 0
               nump_2d = 0
               nump_1d = 0
               CYCLE var_loop
            ELSE
               PRINT '(A,I8,A)','This did not turn out to be sigma-level data, program ID = ',bhi_i1,'.'
               STOP 'Messed_up_multiple_input_files'
            END IF
         ELSE IF ( ( flag .EQ. eot_flag ) .AND. ( current_date(1:19) .EQ. small_header%current_date(1:19) ) ) THEN
            PRINT '(A,A,A)','Found the end of the time period that we want: ',current_date(1:19),'.'
            EXIT var_loop
         ELSE IF  ( flag .EQ. eot_flag ) THEN
            CALL geth_idts (  current_date(1:19) ,  small_header%current_date(1:19) , dts ) 
            IF ( ( bhi(1,1) .EQ. 11 ) .AND. ( ABS(dts) .LE. 2.*bhr(2,12) ) ) THEN
               PRINT '(A)','Found the approximate end of the time period that we want.'
               PRINT '(A,A,A)','Requested date = ',current_date(1:19),'.'
               PRINT '(A,A,A)','Data date      = ',small_header%current_date(1:19),'.'
               EXIT var_loop
            ELSE IF ( ( flag .EQ. eot_flag ) .AND. ( current_date(1:19) .GT. small_header%current_date(1:19) ) ) THEN
               PRINT '(A)','Found the end of the time period, but we want one after this.  Looping on.'
               PRINT '(A,A,A,A,A)','Target date: ',current_date(1:19),', Data date: ',small_header%current_date(1:19),'.'
               DO loop = 1 , nums_3d
                  DEALLOCATE ( alls_3d(loop)%array )
               END DO
               DO loop = 1 , nums_2d
                  DEALLOCATE ( alls_2d(loop)%array )
               END DO
               DO loop = 1 , nums_1d
                  DEALLOCATE ( alls_1d(loop)%array )
               END DO
               DO loop = 1 , nump_3d
                  DEALLOCATE ( allp_3d(loop)%array )
               END DO
               DO loop = 1 , nump_2d
                  DEALLOCATE ( allp_2d(loop)%array )
               END DO
               DO loop = 1 , nump_1d
                  DEALLOCATE ( allp_1d(loop)%array )
               END DO
               nums_3d = 0
               nums_2d = 0
               nums_1d = 0
               nump_3d = 0
               nump_2d = 0
               nump_1d = 0
               CYCLE var_loop
            ELSE IF ( ( flag .EQ. eot_flag ) .AND. ( current_date(1:19) .LT. small_header%current_date(1:19) ) ) THEN
               PRINT '(A)','Found the end of the time period, but we want one before this.'
               PRINT '(A)','Maybe you should look at the dates specified in the namelist and the dates I have found so far.'
               PRINT '(A,A,A)','Selected date from namelist = ',current_date(1:19),'.'
               PRINT '(A,A,A)','Date from input file        = ',small_header%current_date(1:19),'.'
               STOP 'No_can_do_find_the_date'
            END IF
         END IF
   
         !  We have the right flag, so get the header.
   
         READ ( unit_fg ) small_header%num_dims , small_header%start_dims , small_header%end_dims , &
                        small_header%xtime , small_header%staggering , small_header%ordering , &
                        small_header%current_date , small_header%name , small_header%units , &
                        small_header%description
   
         !  Let's allocate space for this one array.  Then read it in.

         IF      ( small_header%num_dims .eq. 3 ) THEN
            
            !  Increment the count of the fields.
          
            nums_3d = nums_3d + 1
            nump_3d = nump_3d + 1
    
            !  Allocate space for the data in the input array area.

            ALLOCATE ( alls_3d(nums_3d)%array(small_header%end_dims(1),small_header%end_dims(2),small_header%end_dims(3)) )
            ALLOCATE ( allp_3d(nump_3d)%array(small_header%end_dims(1),small_header%end_dims(2),kx                      ) )

            !  Use this space to input the data.

            READ ( unit_fg ) alls_3d(nums_3d)%array

            !  Assign the small_header data to this variable's small_header storage.

            alls_3d(nums_3d)%small_header = small_header
            allp_3d(nump_3d)%small_header = small_header

            !  Small fixes to header since this is the p-level data.
            
            allp_3d(nump_3d)%small_header%end_dims(3) = kx
            allp_3d(nump_3d)%small_header%ordering(3:3) = 'P'
   
         ELSE IF ( ( small_header%num_dims .eq. 2 ) .AND. (  small_header%ordering .EQ. 'YX  ' ) ) THEN
            
            !  Increment the count of the fields.
          
            nums_2d = nums_2d + 1
            nump_2d = nump_2d + 1
    
            !  Allocate space for the data in the input array area.

            ALLOCATE ( alls_2d(nums_2d)%array(small_header%end_dims(1),small_header%end_dims(2)) )
            ALLOCATE ( allp_2d(nump_2d)%array(small_header%end_dims(1),small_header%end_dims(2)) )

            !  Use this space to input the data.

            READ ( unit_fg ) alls_2d(nums_2d)%array
            allp_2d(nump_2d)%array = alls_2d(nums_2d)%array

            !  Assign the small_header data to this variable's small_header storage.

            alls_2d(nums_2d)%small_header = small_header
            allp_2d(nump_2d)%small_header = small_header
           
         ELSE IF ( ( small_header%num_dims .eq. 2 ) .AND. (  small_header%ordering .NE. 'YX  ' ) ) THEN

            READ ( unit_fg ) dummy

         ELSE IF ( small_header%num_dims .eq. 1 ) THEN
            
            !  Increment the count of the fields.
          
            nums_1d = nums_1d + 1
    
            !  Allocate space for the data in the input array area.

            ALLOCATE ( alls_1d(nums_1d)%array(small_header%end_dims(1)) )

            !  Use this space to input the data.

            READ ( unit_fg ) alls_1d(nums_1d)%array

            !  Assign the small_header data to this variable's small_header storage.

            alls_1d(nums_1d)%small_header = small_header
   
         END IF
   
      END DO var_loop

      !  Stick in the pressure levels in the allp_1d.

      nump_1d = nump_1d + 1

      ALLOCATE ( allp_1d(nump_1d)%array(kx) )
      allp_1d(nump_1d)%array = newcoord
      allp_1d(nump_1d)%small_header%num_dims    = 1
      allp_1d(nump_1d)%small_header%start_dims  = (/ 1 , 1 , 1 , 1 /)
      allp_1d(nump_1d)%small_header%end_dims    = (/kx , 1 , 1 , 1 /)
      allp_1d(nump_1d)%small_header%xtime       = small_header%xtime
      allp_1d(nump_1d)%small_header%staggering  = 'P   '
      allp_1d(nump_1d)%small_header%ordering    = 'P   '
      allp_1d(nump_1d)%small_header%current_date= small_header%current_date
      allp_1d(nump_1d)%small_header%name        = 'PRESSURE '
      allp_1d(nump_1d)%small_header%units       = 'Pa                       '
      allp_1d(nump_1d)%small_header%description = 'Pressure levels for vertical coordinate       '

      !  We want the date, we need the date, garoyle, bring many more dates.

      sh_date = small_header%current_date(1:19)

   END SUBROUTINE read_data

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE do_namelist ( input_file ,  &
                            start_year , start_month , start_day , start_hour , &
                            start_minute , start_second , start_frac , &
                            end_year ,   end_month ,   end_day ,   end_hour , &
                            end_minute ,   end_second ,   end_frac , &
                            interval ,  &
                            pressure_bu_no_sfc_Pa , &
                            print_info , &
                            uv_interp_method , t_interp_method , moist_interp_method , &
                            height_interp_method , p_interp_method , else_interp_method ,  &
                            uv_extrap_up , t_extrap_up , moist_extrap_up , height_extrap_up , p_extrap_up , else_extrap_up ,  &
                            uv_extrap_low , t_extrap_low , moist_extrap_low , height_extrap_low , p_extrap_low , else_extrap_low )

      IMPLICIT NONE

      !  Input variables.
   
      !  RECORD0
      
      CHARACTER(LEN=132) , DIMENSION(100) :: input_file
      NAMELIST /RECORD0/ input_file
      
      !  RECORD1
      
      INTEGER :: start_year , start_month , start_day , start_hour , &
                 start_minute , start_second , start_frac
      INTEGER ::   end_year ,   end_month ,   end_day ,   end_hour , &
                 end_minute ,   end_second ,   end_frac
      INTEGER :: interval

      NAMELIST /RECORD1/ start_year , start_month , start_day ,   &
                         start_hour ,  start_minute , start_second , start_frac ,  &
                         end_year ,   end_month ,   end_day ,     &
                         end_hour ,  end_minute ,   end_second ,   end_frac ,  &
                         interval
      
      !  RECORD2
      
      REAL , DIMENSION(1000) :: pressure_bu_no_sfc_Pa

      NAMELIST /RECORD2/ pressure_bu_no_sfc_Pa

      !  RECORD3

      LOGICAL :: print_info

      NAMELIST /RECORD3/ print_info

      !  RECORD4

      CHARACTER (LEN=20) :: uv_interp_method , t_interp_method , moist_interp_method , &
                            height_interp_method , p_interp_method , else_interp_method ,  &
                            uv_extrap_up , t_extrap_up , moist_extrap_up , height_extrap_up , p_extrap_up , else_extrap_up ,  &
                            uv_extrap_low , t_extrap_low , moist_extrap_low , height_extrap_low , p_extrap_low , else_extrap_low

      NAMELIST /RECORD4/    uv_interp_method , t_interp_method , moist_interp_method , &
                            height_interp_method , p_interp_method , else_interp_method ,  &
                            uv_extrap_up , t_extrap_up , moist_extrap_up , height_extrap_up , p_extrap_up , else_extrap_up ,  &
                            uv_extrap_low , t_extrap_low , moist_extrap_low , height_extrap_low , p_extrap_low , else_extrap_low

      !  Local variables.

      LOGICAL :: is_it_there = .FALSE.
      INTEGER , PARAMETER :: unit_nml=10

      !  Initializizations.

      start_minute  = 0 
      start_second  = 0
      start_frac    = 0
      end_minute    = 0
      end_second    = 0
      end_frac      = 0
      input_file = '                                                                                ' // &
                   '                                                    '
      print_info = .FALSE.
      uv_interp_method      = 'linear in p         '
      t_interp_method       = 'linear in log p     '
      moist_interp_method   = 'linear in p         '
      height_interp_method  = 'linear in log p     '
      p_interp_method       = 'linear in p         '
      else_interp_method    = 'linear in p         '

      uv_extrap_up          = 'constant            '
      t_extrap_up           = 'constant            '
      moist_extrap_up       = 'constant            '
      height_extrap_up      = 'extrapolate         '
      p_extrap_up           = 'extrapolate         '
      else_extrap_up        = 'constant            '

      uv_extrap_low         = 'constant            '
      t_extrap_low          = 'extrapolate         '
      moist_extrap_low      = 'constant            '
      height_extrap_low     = 'extrapolate         '
      p_extrap_low          = 'extrapolate         '
      else_extrap_low       = 'constant            ' 

!    &record4
!     uv_interp_method      = 'linear in p         '
!     t_interp_method       = 'linear in log p     '
!     moist_interp_method   = 'linear in p         '
!     height_interp_method  = 'linear in log p     '
!     p_interp_method       = 'linear in p         '
!     else_interp_method    = 'linear in p         '
!     uv_extrap_up          = 'constant            '
!     t_extrap_up           = 'constant            '
!     moist_extrap_up       = 'constant            '
!     height_extrap_up      = 'extrapolate         '
!     p_extrap_up           = 'extrapolate         '
!     else_extrap_up        = 'constant            '
!     uv_extrap_low         = 'constant            '
!     t_extrap_low          = 'extrapolate         '
!     moist_extrap_low      = 'constant            '
!     height_extrap_low     = 'extrapolate         '
!     p_extrap_low          = 'extrapolate         '
!     else_extrap_low       = 'constant            ' /
   
      !     Does the file exist?
      
      INQUIRE ( FILE = 'namelist.input' , EXIST = is_it_there )
      
      IF ( is_it_there ) THEN
      
         !  The file exists, get a unit number.
      
         OPEN ( FILE   = 'namelist.input' , &
                UNIT   =  unit_nml        , &
                STATUS = 'OLD'            , &
                FORM   = 'FORMATTED'      , &
                ACTION = 'READ'           , &
                ACCESS = 'SEQUENTIAL'     )
   
         !  File is opened, so read it.
   
         READ (unit_nml , NML = RECORD0 )
         READ (unit_nml , NML = RECORD1 )
         READ (unit_nml , NML = RECORD2 )
         READ (unit_nml , NML = RECORD3 )
!        READ (unit_nml , NML = RECORD4 )

         !  Tell the people what we wore.

         IF ( print_info ) THEN
!           WRITE (6    , NML = RECORD0 )
            WRITE (6    , NML = RECORD1 )
#ifndef IBM
!           WRITE (6    , NML = RECORD2 )
#endif
            WRITE (6    , NML = RECORD3 )
!           WRITE (6    , NML = RECORD4 )
         END IF
   
      ELSE
         PRINT '(A)','Could not find the namelist: "namelist.input".'
         STOP 'No_namelist_found'
      END IF
   
      CLOSE ( unit_nml )
   
   END SUBROUTINE do_namelist

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
   SUBROUTINE get_fg_file ( input_file , unit_fg )
   
      IMPLICIT NONE
   
      !  Input variables.

      CHARACTER(LEN=132) , INTENT(IN) :: input_file

      INTEGER , INTENT(IN) :: unit_fg

      !  Local variables.

      LOGICAL :: is_it_there = .FALSE.
   
      !  Does the file exist?
   
      INQUIRE ( FILE = TRIM(input_file) , EXIST = is_it_there )
   
      IF ( is_it_there ) THEN
         OPEN ( FILE   = TRIM(input_file) , &
                UNIT   =  unit_fg         ,  &
                STATUS = 'OLD'            , &
                FORM   = 'UNFORMATTED'    , &
                ACCESS = 'SEQUENTIAL'     )
         PRINT '(A)','Opened sigma level data for input.'
      
      ELSE
         PRINT '(A,A,A)','Could not find file ',TRIM(input_file),'.'
         STOP 'No_fg_file_found'
      
      END IF
   
   END SUBROUTINE get_fg_file
   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
   SUBROUTINE open_out_file ( mmoutp_unit , regridder_unit , domain_id )
   
      IMPLICIT NONE
      
      !  Input variables.

      INTEGER , INTENT(IN) :: mmoutp_unit , regridder_unit , domain_id
      
      !  Local variables

      CHARACTER(LEN=132) :: mmoutp_name , regridder_name

      WRITE ( mmoutp_name , '("MMOUTP_DOMAIN",i1)' ) domain_id
      
      OPEN ( FILE   =  TRIM(mmoutp_name)      , &
             UNIT   =  mmoutp_unit            , &
             STATUS = 'UNKNOWN'               , &
             FORM   = 'UNFORMATTED'           , &
             ACCESS = 'SEQUENTIAL'              )
             PRINT '(A)','Opened MMOUTP file for output.'

      WRITE ( regridder_name , '("REGRID_DOMAIN",i1)' ) domain_id
      
      OPEN ( FILE   =  TRIM(regridder_name)   , &
             UNIT   =  regridder_unit         , &
             STATUS = 'UNKNOWN'               , &
             FORM   = 'UNFORMATTED'           , &
             ACCESS = 'SEQUENTIAL'              )
             PRINT '(A)','Opened REGRID file for output.'

   END SUBROUTINE open_out_file

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE out_mmoutp ( current_date19 , &
                           imx , jmx , kxp , output_mmoutp , itimes , interval , &
                           start_year , start_month , start_day , start_hour , &
                           start_minute , start_second , start_frac , &
                           pressure , print_info )

      IMPLICIT NONE

      !  Input data.

      CHARACTER (LEN=19) , INTENT(IN) :: current_date19

      INTEGER , INTENT(IN)  :: imx , jmx , kxp , output_mmoutp , itimes , interval , &
                               start_year , start_month , start_day , start_hour , &
                               start_minute , start_second , start_frac

      REAL , DIMENSION(:) , INTENT(IN) :: pressure

      LOGICAL , INTENT(IN) :: print_info

      !  Local variables.

      INTEGER :: loop

      !  Variables for expanded domain cut down for the pass-through arrays.

      REAL , DIMENSION(imx,jmx) :: dum2d
      REAL , DIMENSION(imx,jmx,kxp) :: dum3d

      !  Data for small header.

      INTEGER                :: num_dims
      INTEGER , DIMENSION(4) :: start_dims
      INTEGER , DIMENSION(4) :: end_dims
      REAL                   :: xtime
      CHARACTER (LEN= 4)     :: staggering
      CHARACTER (LEN= 4)     :: ordering
      CHARACTER (LEN=24)     :: date
      CHARACTER (LEN= 9)     :: name
      CHARACTER (LEN=25)     :: units
      CHARACTER (LEN=46)     :: description
   
      INTEGER , PARAMETER :: index = 8

      !  If this is the first time in here, then we need to build and output 
      !  the big header.

      IF ( itimes .EQ. 1 ) THEN

         !  First, the flag, then the data.

         WRITE ( output_mmoutp ) bh_flag

         !  Just in case no one else cleared the big header, we should do it.  All
         !  of the header values for the model input index can be set to -999.

         IF ( index .NE. 3 ) THEN
            bhi (:,index) = -999
            bhr (:,index) = -999
            bhic(:,index) = '                                                                                '
            bhrc(:,index) = '                                                                                '
         END IF

         !  Identify the data as model input.

         bhi( 1,1) = index ; bhic( 1,1) = 'First guess fields interpolated from model output                               '
 
         !  A few lines about data format and version numbers.

         bhi( 2,index) = 1 ; bhic( 2,index) = 'INTERPB Version 3 MM5 System Format Edition Number                              '
         bhi( 3,index) = 3 ; bhic( 3,index) = 'INTERPB Program Version Number                                                  '
         bhi( 4,index) = 0 ; bhic( 4,index) = 'INTERPB Program Minor Revision Number                                           '
         
         !  Starting date of the INTERPB data, this was specified in the namelist.

         bhi ( 5,index) = start_year   
         bhic( 5,index) = 'Four-digit year of start time                                                   '
         bhi ( 6,index) = start_month
         bhic( 6,index) = 'Month of the year of the start time (1-12)                                      '
         bhi ( 7,index) = start_day
         bhic( 7,index) = 'Day of the month of the start time (1-31)                                       '
         bhi ( 8,index) = start_hour
         bhic( 8,index) = 'Hour of the day of the start time (0-23)                                        '
         bhi ( 9,index) = start_minute
         bhic( 9,index) = 'Minute of the start time (0-59)                                                 '
         bhi (10,index) = start_second
         bhic(10,index) = 'Second of the start time (0-59)                                                 '
         bhi (11,index) = start_frac
         bhic(11,index) = 'Ten thousandths of a second of the start time (0-9999)                          '

         !  How many pressure levels to expect.

         bhi (12,index) = kxp
         bhic(12,index) = 'Number of pressure levels in the output, bottom up, including surface           '

         !  Time interval between output files.

         bhr ( 1,index) = interval
         bhrc( 1,index) = 'Time difference (seconds) between isobaric interpolated model output files.     '

         !  Alrighty, we can output it now.
  
         WRITE ( output_mmoutp ) bhi , bhr , bhic , bhrc
         
      END IF

      !  Triplets of data follow: 1) small head flag, 2) small header, 3) data.

      DO loop = 1 , nump_3d
         IF ( allp_3d(loop)%small_header%name(1:3) .NE. 'P3D' ) THEN
            WRITE ( output_mmoutp ) sh_flag
            num_dims        = 3
            start_dims      = (/ 1, 1, 1, 1 /)
            end_dims        = (/ imx, jmx, kxp, 1 /)
            xtime           = allp_3d(loop)%small_header%xtime
            staggering      = allp_3d(loop)%small_header%staggering
            ordering        = allp_3d(loop)%small_header%ordering
            ordering(3:3)   = 'P'
            date            = current_date19 // '.0000'
            name            = allp_3d(loop)%small_header%name
            units           = allp_3d(loop)%small_header%units
            description     = allp_3d(loop)%small_header%description
            WRITE ( output_mmoutp ) num_dims , start_dims , end_dims , xtime , &
                             staggering , ordering , date , name , units , description
            dum3d = allp_3d(loop)%array
            CALL fill3d ( dum3d , imx , jmx , kxp , staggering )
            WRITE ( output_mmoutp ) dum3d 
            IF ( print_info ) THEN
               PRINT '(A,A,A)','Outputting 3D data ',name(1:8),'.'
            END IF
         END IF
      END DO

      !  Now the 2D data.

      DO loop = 1 , nums_2d
         IF ( ( alls_2d(loop)%small_header%ordering(1:2) .EQ. 'YX'  ) .AND. &
              ( alls_2d(loop)%small_header%name(1:3)     .NE. 'SLP' ) ) THEN
            WRITE ( output_mmoutp ) sh_flag
            num_dims        = 2
            start_dims      = (/ 1, 1, 1, 1 /)
            end_dims        = (/ imx, jmx, 1, 1 /)
            xtime           = alls_2d(loop)%small_header%xtime
            staggering      = alls_2d(loop)%small_header%staggering
            ordering        = 'YX  '
            date            = current_date19 // '.0000'
            name            = alls_2d(loop)%small_header%name
            IF ( name(1:3) .EQ. 'RH ' ) THEN
               name(1:8) = 'RH SFC  '
            END IF
            units           = alls_2d(loop)%small_header%units
            description     = alls_2d(loop)%small_header%description
            WRITE ( output_mmoutp ) num_dims , start_dims , end_dims , xtime , &
                             staggering , ordering , date , name , units , description
            dum2d = alls_2d(loop)%array
            CALL fill2d ( dum2d , imx , jmx , staggering )
            WRITE ( output_mmoutp ) dum2d 
            IF ( print_info ) THEN
               PRINT '(A,A,A)','Outputting 2D data ',name(1:8),'.'
            END IF
         END IF
      END DO

      !  Pressure

      WRITE ( output_mmoutp ) sh_flag
      num_dims        = 1
      start_dims      = (/ 1, 1, 1, 1 /)
      end_dims        = (/ kxp , 1 , 1 , 1 /)
      xtime           = 0
      staggering      = 'P   '
      ordering        = 'P   '
      date            = current_date19 // '.0000'
      name            = 'PRESSURE '
      units           = 'Pa                       '
      description     = 'Bottom-up pressure levels, includes surface   '
      WRITE ( output_mmoutp ) num_dims , start_dims , end_dims , xtime , &
                       staggering , ordering , date , name , units , description
      WRITE ( output_mmoutp ) pressure
      IF ( print_info ) THEN
         PRINT '(A,A,A)','Outputting 1D data ',name(1:8),'.'
      END IF
      
      !  The last thing you do before going to sleep at night is to put out the end of time flag.

      WRITE ( output_mmoutp ) eot_flag

   END SUBROUTINE out_mmoutp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE out_regridder ( current_date19 , &
                              imx , jmx , kxp , output_regridder , itimes , interval , ptop , &
                              start_year , start_month , start_day , start_hour , &
                              start_minute , start_second , start_frac , &
                              pressure , print_info )

      IMPLICIT NONE

      !  Input data.

      CHARACTER (LEN=19) , INTENT(IN) :: current_date19

      INTEGER , INTENT(IN)  :: imx , jmx , kxp , output_regridder , itimes , interval , &
                               start_year , start_month , start_day , start_hour , &
                               start_minute , start_second , start_frac

      REAL , INTENT(IN)  :: ptop

      REAL , DIMENSION(:) , INTENT(IN) :: pressure

      LOGICAL , INTENT(IN) :: print_info

      !  Local variables.

      INTEGER :: loop

      !  Variables for expanded domain cut down for the pass-through arrays.

      REAL , DIMENSION(imx,jmx) :: dum2d
      REAL , DIMENSION(imx,jmx,kxp) :: dum3d

      !  Data for small header.

      INTEGER                :: num_dims
      INTEGER , DIMENSION(4) :: start_dims
      INTEGER , DIMENSION(4) :: end_dims
      REAL                   :: xtime
      CHARACTER (LEN= 4)     :: staggering
      CHARACTER (LEN= 4)     :: ordering
      CHARACTER (LEN=24)     :: date
      CHARACTER (LEN= 9)     :: name
      CHARACTER (LEN=25)     :: units
      CHARACTER (LEN=46)     :: description
   
      INTEGER , PARAMETER :: index = 2

      LOGICAL :: found_landmask
      INTEGER :: index_landuse , i , j

      !  If this is the first time in here, then we need to build and output 
      !  the big header.

      IF ( itimes .EQ. 1 ) THEN

         !  First, the flag, then the data.

         WRITE ( output_regridder ) bh_flag

         !  Just in case no one else cleared the big header, we should do it.  All
         !  of the header values for the model input index can be set to -999.

         IF ( index .NE. 3 ) THEN
            bhi (:,index) = -999
            bhr (:,index) = -999
            bhic(:,index) = '                                                                                '
            bhrc(:,index) = '                                                                                '
         END IF

         !  Identify the data as model input.

         bhi( 1,1) = index ; bhic( 1,1) = 'First guess fields interpolated from model output                               '
 
         !  Change expanded domain option to non-expanded domain option in terrain header

         bhi( 8,1) = 0
         if ( bhi(13,1) == 1 ) then
              bhi( 9,1) = bhi( 5,1) 
              bhi(10,1) = bhi( 6,1) 
         endif

         !  A few lines about data format and version numbers.

         bhi( 2,index) = 1 ; bhic( 2,index) = 'INTERPB Version 3 MM5 System Format Edition Number                              '
         bhi( 3,index) = 2 ; bhic( 3,index) = 'INTERPB Program Version Number                                                  '
         bhi( 4,index) = 0 ; bhic( 4,index) = 'INTERPB Program Minor Revision Number                                           '
         
         !  Starting date of the INTERPB data, this was specified in the namelist.

         bhi ( 5,index) = start_year   
         bhic( 5,index) = 'Four-digit year of start time                                                   '
         bhi ( 6,index) = start_month
         bhic( 6,index) = 'Month of the year of the start time (1-12)                                      '
         bhi ( 7,index) = start_day
         bhic( 7,index) = 'Day of the month of the start time (1-31)                                       '
         bhi ( 8,index) = start_hour
         bhic( 8,index) = 'Hour of the day of the start time (0-23)                                        '
         bhi ( 9,index) = start_minute
         bhic( 9,index) = 'Minute of the start time (0-59)                                                 '
         bhi (10,index) = start_second
         bhic(10,index) = 'Second of the start time (0-59)                                                 '
         bhi (11,index) = start_frac
         bhic(11,index) = 'Ten thousandths of a second of the start time (0-9999)                          '

         !  How many pressure levels to expect.

         bhi (12,index) = kxp
         bhic(12,index) = 'Number of pressure levels in the output, bottom up, including surface           '

         !  Time interval between output files.

         bhr ( 1,index) = interval
         bhrc( 1,index) = 'Time difference (seconds) between isobaric interpolated model output files.     '
         bhr ( 2,index) = ptop
         bhrc( 2,index) = 'Top pressure used in analysis, pressure defining model lid (Pa)                 '

         !  Alrighty, we can output it now.
  
         WRITE ( output_regridder ) bhi , bhr , bhic , bhrc
         
      END IF

      !  Triplets of data follow: 1) small head flag, 2) small header, 3) data.

      DO loop = 1 , nump_3d
         IF ( ( allp_3d(loop)%small_header%name(1:8) .EQ. 'U       ' ) .OR. &
              ( allp_3d(loop)%small_header%name(1:8) .EQ. 'V       ' ) .OR. &
              ( allp_3d(loop)%small_header%name(1:8) .EQ. 'T       ' ) .OR. &
              ( allp_3d(loop)%small_header%name(1:8) .EQ. 'RH      ' ) .OR. &
              ( allp_3d(loop)%small_header%name(1:8) .EQ. 'H       ' ) ) THEN
            WRITE ( output_regridder ) sh_flag
            num_dims        = 3
            start_dims      = (/ 1, 1, 1, 1 /)
            end_dims        = (/ imx, jmx, kxp, 1 /)
            xtime           = allp_3d(loop)%small_header%xtime
            staggering      = allp_3d(loop)%small_header%staggering
            ordering        = allp_3d(loop)%small_header%ordering
            ordering(3:3)   = 'P'
            date            = current_date19 // '.0000'
            name            = allp_3d(loop)%small_header%name
            units           = allp_3d(loop)%small_header%units
            description     = allp_3d(loop)%small_header%description
            WRITE ( output_regridder ) num_dims , start_dims , end_dims , xtime , &
                             staggering , ordering , date , name , units , description
            dum3d = allp_3d(loop)%array
            CALL fill3d ( dum3d , imx , jmx , kxp , staggering )
            WRITE ( output_regridder ) dum3d 
            IF ( print_info ) THEN
               PRINT '(A,A,A)','Outputting 3D data ',name(1:8),'.'
            END IF
         ELSE
            IF ( print_info ) THEN
               PRINT '(A,A,A,I1,A)','Skipping output of 3d array ',allp_3d(loop)%small_header%name(1:8), &
                                  ' for REGRID_DOMAIN',bhi(13,1),' file.'
            END IF
         END IF
      END DO

      !  If there is no LANDMASK data

      found_landmask = .FALSE.
      DO loop = 1 , nums_2d
         IF      ( alls_2d(loop)%small_header%name(1:8) .EQ. 'LANDMASK' ) THEN
            found_landmask = .TRUE.
         ELSE IF ( alls_2d(loop)%small_header%name(1:8) .EQ. 'LAND USE' ) THEN
            index_landuse = loop
         ENDIF
      END DO
      IF ( .NOT. found_landmask ) THEN
         IF ( print_info ) THEN
            PRINT '(A,I2,A)','Creating LANDMASK field from LAND USE, using category ',bhi(23,1),'.'
         END IF
      END IF

      !  Now the 2D data.

      DO loop = 1 , nums_2d
         IF ( alls_2d(loop)%small_header%ordering(1:2) .EQ. 'YX' ) THEN
            IF ( ( alls_2d(loop)%small_header%name(1:8) .EQ. 'TERRAIN ' ) .OR. &
                 ( alls_2d(loop)%small_header%name(1:8) .EQ. 'MAPFACCR' ) .OR. &
                 ( alls_2d(loop)%small_header%name(1:8) .EQ. 'MAPFACDT' ) .OR. &
                 ( alls_2d(loop)%small_header%name(1:8) .EQ. 'CORIOLIS' ) .OR. &
                 ( alls_2d(loop)%small_header%name(1:8) .EQ. 'LATITCRS' ) .OR. &
                 ( alls_2d(loop)%small_header%name(1:8) .EQ. 'LONGICRS' ) .OR. &
                 ( alls_2d(loop)%small_header%name(1:8) .EQ. 'LAND USE' ) .OR. &
                 ( alls_2d(loop)%small_header%name(1:8) .EQ. 'LANDMASK' ) .OR. &
                 ( alls_2d(loop)%small_header%name(1:8) .EQ. 'SNOWCOVR' ) .OR. &
                 ( alls_2d(loop)%small_header%name(1:8) .EQ. 'TSEASFC ' ) .OR. &
                 ( alls_2d(loop)%small_header%name(1:8) .EQ. 'PSEALVLC' ) .OR. &
                 ( alls_2d(loop)%small_header%name(1:8) .EQ. 'PSEALVLD' ) .OR. &
                 ( alls_2d(loop)%small_header%name(1:8) .EQ. 'LATITDOT' ) .OR. &
                 ( alls_2d(loop)%small_header%name(1:8) .EQ. 'LONGIDOT' ) ) THEN
               WRITE ( output_regridder ) sh_flag
               num_dims        = 2
               start_dims      = (/ 1, 1, 1, 1 /)
               end_dims        = (/ imx, jmx, 1, 1 /)
               xtime           = alls_2d(loop)%small_header%xtime
               staggering      = alls_2d(loop)%small_header%staggering
               ordering        = alls_2d(loop)%small_header%ordering
               date            = current_date19 // '.0000'
               name            = alls_2d(loop)%small_header%name
               IF ( name(1:3) .EQ. 'RH ' ) THEN
                  name(1:8) = 'RH SFC  '
               END IF
               units           = alls_2d(loop)%small_header%units
               description     = alls_2d(loop)%small_header%description
               WRITE ( output_regridder ) num_dims , start_dims , end_dims , xtime , &
                                staggering , ordering , date , name , units , description
               dum2d = alls_2d(loop)%array
               CALL fill2d ( dum2d , imx , jmx , staggering )
               WRITE ( output_regridder ) dum2d 
               IF ( print_info ) THEN
                  PRINT '(A,A,A)','Outputting 2D data ',name(1:8),'.'
               END IF

               IF ( ( alls_2d(loop)%small_header%name(1:8) .EQ. 'LAND USE' ) .AND. ( .NOT. found_landmask ) ) THEN
                  WRITE ( output_regridder ) sh_flag
                  num_dims        = 2
                  start_dims      = (/ 1, 1, 1, 1 /)
                  end_dims        = (/ imx, jmx, 1, 1 /)
                  xtime           = alls_2d(loop)%small_header%xtime
                  staggering      = alls_2d(loop)%small_header%staggering
                  ordering        = alls_2d(loop)%small_header%ordering
                  date            = current_date19 // '.0000'
                  name(1:9)       = 'LANDMASK '
                  units           = 'FLAG: 0=water, 1=land    '
                  description     = 'Land mask: 0=water point, 1=land point        '
                  WRITE ( output_regridder ) num_dims , start_dims , end_dims , xtime , &
                                   staggering , ordering , date , name , units , description
                  DO j = 1 , jmx -1
                     DO i = 1 , imx -1
                        IF ( ABS( REAL(bhi(23,1)) - alls_2d(loop)%array(i,j) ) .LT. 0.1 ) THEN
                           dum2d(i,j) = 0.
                        ELSE
                           dum2d(i,j) = 1.
                        END IF
                     END DO
                  END DO
                  CALL fill2d ( dum2d , imx , jmx , staggering )
                  WRITE ( output_regridder ) dum2d 
                  IF ( print_info ) THEN
                     PRINT '(A,A,A)','Outputting 2D data ',name(1:8),'.'
                  END IF
               END IF
            ELSE
               IF ( print_info ) THEN
                  PRINT '(A,A,A,I1,A)','Skipping output of 2d array ',alls_2d(loop)%small_header%name(1:8), &
                                       ' for REGRID_DOMAIN',bhi(13,1),' file.'
               END IF
            END IF
         END IF
      END DO

      !  Pressure

      WRITE ( output_regridder ) sh_flag
      num_dims        = 1
      start_dims      = (/ 1, 1, 1, 1 /)
      end_dims        = (/ kxp , 1 , 1 , 1 /)
      xtime           = 0
      staggering      = 'P   '
      ordering        = 'P   '
      date            = current_date19 // '.0000'
      name            = 'PRESSURE '
      units           = 'Pa                       '
      description     = 'Bottom-up pressure levels, includes surface   '
      WRITE ( output_regridder ) num_dims , start_dims , end_dims , xtime , &
                       staggering , ordering , date , name , units , description
      WRITE ( output_regridder ) pressure
      IF ( print_info ) THEN
         PRINT '(A,A,A)','Outputting 1D data ',name(1:8),'.'
      END IF
      
      !  The last thing you do before going to sleep at night is to put out the end of time flag.

      WRITE ( output_regridder ) eot_flag

   END SUBROUTINE out_regridder

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE out_intermediate ( current_date19 , &
                                 imx , jmx , kxp , output_intermediate , itimes , interval , &
                                 start_year , start_month , start_day , start_hour , &
                                 start_minute , start_second , start_frac , &
                                 pressure , print_info )

      IMPLICIT NONE

      !  Input data.

      CHARACTER (LEN=19) , INTENT(IN) :: current_date19

      INTEGER , INTENT(IN)  :: imx , jmx , kxp , output_intermediate , itimes , interval , &
                               start_year , start_month , start_day , start_hour , &
                               start_minute , start_second , start_frac

      REAL , DIMENSION(:) , INTENT(IN) :: pressure

      LOGICAL , INTENT(IN) :: print_info

      !  Local variables.

      INTEGER :: i , j , k , loop

      !  The slabs of data to be output must be the correct size with respect
      !  to the dot/cross staggering.  No extraneous rows/columns are allowed
      !  to be floating around.

      REAL , DIMENSION(jmx-1,imx-1) :: dum2dc
      REAL , DIMENSION(jmx  ,imx  ) :: dum2dd

      !  Data from small header.

      REAL                   :: xtime
      CHARACTER (LEN= 4)     :: staggering
      CHARACTER (LEN=24)     :: date
      CHARACTER (LEN= 9)     :: name
      CHARACTER (LEN=25)     :: units
      CHARACTER (LEN=46)     :: description

      !  The I/O flag on the open statement.

      INTEGER :: open_status

      !  Some little ditties that are required when we output the data in the intermediate 
      !  format.  There are several different headers through which we must wade.

      INTEGER :: version = 3
      INTEGER :: llflag

      REAL :: level

      REAL :: lat1 , lat1d , lat1c , lon1 , lon1d , lon1c

      !  Open the output file.  Each time inside this routine (each time period) is a
      !  separate file.  The file names are created with a date/time stamp, consistent
      !  with the expectations of the regridder program.

      OPEN ( FILE   = 'FILE_MMOUTP:' // current_date19(1:13) , &
             UNIT   =  output_intermediate                   , &
             STATUS = 'UNKNOWN'                              , &
             FORM   = 'UNFORMATTED'                          , &
             ACCESS = 'SEQUENTIAL'                           , &
             IOSTAT =  open_status                             )

      IF ( open_status .NE. 0 ) THEN
         PRINT '(A,A,A)' , 'Error OPENing intermediate output file FILE_MMOUTP:',current_date19(1:13),'.'
         STOP 'Error_open_intermediate_for_output'
      ELSE
         PRINT '(A,A,A)','Opened file FILE_MMOUTP:',current_date19(1:13),' for intermediate format output.'
      END IF
 
      !  We need the lower left lat/lon for the cross and dot point.

      DO loop = 1 , nums_2d
         IF      ( alls_2d(loop)%small_header%name(1:8) .EQ. 'LATITDOT' ) THEN
            lat1d = alls_2d(loop)%array(1,1)
         ELSE IF ( alls_2d(loop)%small_header%name(1:8) .EQ. 'LATITCRS' ) THEN
            lat1c = alls_2d(loop)%array(1,1)
         ELSE IF ( alls_2d(loop)%small_header%name(1:8) .EQ. 'LONGIDOT' ) THEN
            lon1d = alls_2d(loop)%array(1,1)
         ELSE IF ( alls_2d(loop)%small_header%name(1:8) .EQ. 'LONGICRS' ) THEN
            lon1c = alls_2d(loop)%array(1,1)
         END IF
      END DO

      !  Projection constants for this data set.  We have to take the MM5 system projection
      !  info (1=Lambert, 2=polar, 3=Mercator) and map to to the intermediate format's
      !  request for said info (3=Lambert, 5=polar, 1=Mercator).  Could we have made this
      !  a bit less intuitive?

      IF      ( bhi(7,1) .EQ. 1 ) THEN
         llflag = 3
      ELSE IF ( bhi(7,1) .EQ. 2 ) THEN
         llflag = 5
      ELSE IF ( bhi(7,1) .EQ. 3 ) THEN
         llflag = 1
      END IF

      !  Loop over the 3d data, output the 2d slabs.

      all_3d : DO loop = 1 , nump_3d

         IF ( allp_3d(loop)%small_header%name(1:3) .NE. 'P3D' ) THEN

            !  The data that we need, which is contained in the small header, is as follows and I quote:
   
            xtime           = allp_3d(loop)%small_header%xtime
            staggering      = allp_3d(loop)%small_header%staggering
            date            = current_date19 // '.0000'
            name            = allp_3d(loop)%small_header%name
            units           = allp_3d(loop)%small_header%units
            description     = allp_3d(loop)%small_header%description
            IF ( name(1:2) .EQ. 'H ' ) THEN
               name(1:8) = 'HGT     '
            END IF
   
            !  Loop over each pressure level - this is the 3d section, ya know.  The intermediate
            !  format handles 2d slabs ONLY.
   
            all_press : DO k = 1 , kxp
   
               !  First output is the version number.  This is set to "3" right now.
      
               WRITE ( output_intermediate ) version
   
               !  State the pressure level on which we choose to be currently working.
   
               level = pressure(k)
             
               !  The regridder program expects a peculiar pressure level for the surface flag.
   
               IF ( ABS(level-100100.) .LT. 1. ) THEN
                  level = 200100
               END IF
   
               !  Many, many details that are a function of the staggering and a concern for us.
   
               IF      ( staggering(1:1) .EQ. 'C' ) THEN
   
                  !  The size of the array has to match the cross/dot staggering.  There can be
                  !  no extra rows/columns.  Write the next part of the header.  Note that this
                  !  array is only the size of the cross point data.
   
                  WRITE ( output_intermediate ) &
                  date, xtime, name, units, description, level, jmx-1, imx-1, llflag
           
                  !  Assign the staggering corrected lat/lon at grid point (1,1).
   
                  lat1 = lat1c
                  lon1 = lon1c
      
                  !  Load up the data into the 2d array for output.  From MM5 format, the data 
                  !  must be tranposed.  Note the cross point indices.
   
                  DO j = 1 , jmx - 1
                     DO i = 1 , imx - 1
                        dum2dc(j,i) = allp_3d(loop)%array(i,j,k)
                     END DO
                  END DO
   
               ELSE IF ( staggering(1:1) .EQ. 'D' ) THEN
   
                  !  The size of the array has to match the cross/dot staggering.  There can be
                  !  no extra rows/columns.  Write the next part of the header.  Note this array is
                  !  the full dot point size.
   
                  WRITE ( output_intermediate ) &
                  date, xtime, name, units, description, level, jmx  , imx  , llflag
           
                  !  Assign the staggering corrected lat/lon at grid point (1,1).
   
                  lat1 = lat1d
                  lon1 = lon1d
      
                  !  Load up the data into the 2d array for output.  From MM5 format, the data 
                  !  must be tranposed.  Note the dot point indices.
   
                  DO j = 1 , jmx
                     DO i = 1 , imx
                        dum2dd(j,i) = allp_3d(loop)%array(i,j,k)
                     END DO
                  END DO
   
               END IF
   
               !  Projection specific info.
   
               IF      ( bhi(7,1) .EQ. 1 ) THEN
                  WRITE ( output_intermediate ) lat1, lon1, bhr(9,1)/1000., bhr(9,1)/1000., &
                                                bhr(3,1), bhr(5,1), bhr(6,1)
               ELSE IF ( bhi(7,1) .EQ. 2 ) THEN
                  WRITE ( output_intermediate ) lat1, lon1, bhr(9,1)/1000., bhr(9,1)/1000., &
                                                bhr(3,1), bhr(5,1)
               ELSE IF ( bhi(7,1) .EQ. 3 ) THEN
                  WRITE ( output_intermediate ) lat1, lon1, bhr(9,1)/1000., bhr(9,1)/1000., &
                                                bhr(5,1)
               END IF
   
               !  Now we output the slice of data.  Either the dum2dc (cross point) or the dum2dd (dot
               !  point) array was filled.  This array is the correct size for the staggering.
   
               IF      ( staggering(1:1) .EQ. 'C' ) THEN
                  WRITE ( output_intermediate ) dum2dc
               ELSE IF ( staggering(1:1) .EQ. 'D' ) THEN
                  WRITE ( output_intermediate ) dum2dd
               END IF
   
               IF ( print_info ) THEN
                  PRINT '(A,A,A,I4,A)','Outputting 3D data ',name(1:8),', level = ',NINT(level/100.),' mb.'
               END IF
   
            END DO all_press

         END IF

      END DO all_3d

      !  Loop over the 2d data, output the 2d slabs.  There are precious few 2d fields that we
      !  want.  The trouble would arise if we had a field name identically to that which the
      !  terrain program would provide to regridder.

      all_2d : DO loop = 1 , nums_2d

         IF ( ( alls_2d(loop)%small_header%name(1:8) .EQ. 'SNOWCOVR' ) .OR. &
!             ( alls_2d(loop)%small_header%name(1:8) .EQ. 'TSEASFC ' ) .OR. &
              ( alls_2d(loop)%small_header%name(1:8) .EQ. 'GROUND T' ) .OR. &
              ( alls_2d(loop)%small_header%name(1:8) .EQ. 'PSEALVLC' ) ) THEN

            !  First output is the version number.  This is set to "3" right now.
   
            WRITE ( output_intermediate ) version
   
            !  The data that we need, which is contained in the small header, is as follows and I quote:
   
            xtime           = alls_2d(loop)%small_header%xtime
            staggering      = alls_2d(loop)%small_header%staggering
            date            = current_date19 // '.0000'
            name            = alls_2d(loop)%small_header%name
            IF      ( name(1:3) .EQ. 'RH '      ) THEN
               name(1:8) = 'RH SFC  '
            ELSE IF ( name(1:8) .EQ. 'TSEASFC ' ) THEN
               name(1:8) = 'SST     '
            ELSE IF ( name(1:8) .EQ. 'GROUND T' ) THEN
               name(1:8) = 'SST     '
            ELSE IF ( name(1:8) .EQ. 'PSEALVLC' ) THEN
               name(1:8) = 'PMSL    '
            END IF
            units           = alls_2d(loop)%small_header%units
            description     = alls_2d(loop)%small_header%description
   
            level = 200100
   
            !  Many, many details that are a function of the staggering and a concern for us.
   
            IF      ( staggering(1:1) .EQ. 'C' ) THEN
   
               !  The size of the array has to match the cross/dot staggering.  There can be
               !  no extra rows/columns.  Write the next part of the header.  Note that this
               !  array is only the size of the cross point data.
   
               WRITE ( output_intermediate ) &
               date, xtime, name, units, description, level, jmx-1, imx-1, llflag
        
               !  Assign the staggering corrected lat/lon at grid point (1,1).
   
               lat1 = lat1c
               lon1 = lon1c
   
               !  Load up the data into the 2d array for output.  From MM5 format, the data 
               !  must be tranposed.  Note the cross point indices.
   
               DO j = 1 , jmx - 1
                  DO i = 1 , imx - 1
                     dum2dc(j,i) = alls_2d(loop)%array(i,j)
                  END DO
               END DO
   
            ELSE IF ( staggering(1:1) .EQ. 'D' ) THEN
   
               !  The size of the array has to match the cross/dot staggering.  There can be
               !  no extra rows/columns.  Write the next part of the header.  Note this array is
               !  the full dot point size.
   
               WRITE ( output_intermediate ) &
               date, xtime, name, units, description, level, jmx  , imx  , llflag
        
               !  Assign the staggering corrected lat/lon at grid point (1,1).
   
               lat1 = lat1d
               lon1 = lon1d
   
               !  Load up the data into the 2d array for output.  From MM5 format, the data 
               !  must be tranposed.  Note the dot point indices.
   
               DO j = 1 , jmx
                  DO i = 1 , imx
                     dum2dd(j,i) = alls_2d(loop)%array(i,j)
                  END DO
               END DO
   
            END IF
   
            !  Projection specific info.
   
            IF      ( bhi(7,1) .EQ. 1 ) THEN
               WRITE ( output_intermediate ) lat1, lon1, bhr(9,1)/1000., bhr(9,1)/1000., &
                                             bhr(3,1), bhr(5,1), bhr(6,1)
            ELSE IF ( bhi(7,1) .EQ. 2 ) THEN
               WRITE ( output_intermediate ) lat1, lon1, bhr(9,1)/1000., bhr(9,1)/1000., &
                                             bhr(3,1), bhr(5,1)
            ELSE IF ( bhi(7,1) .EQ. 3 ) THEN
               WRITE ( output_intermediate ) lat1, lon1, bhr(9,1)/1000., bhr(9,1)/1000., &
                                             bhr(5,1)
            END IF
   
            !  Now we output the slice of data.  Either the dum2dc (cross point) or the dum2dd (dot
            !  point) array was filled.  This array is the correct size for the staggering.
   
            IF      ( staggering(1:1) .EQ. 'C' ) THEN
               WRITE ( output_intermediate ) dum2dc
            ELSE IF ( staggering(1:1) .EQ. 'D' ) THEN
               WRITE ( output_intermediate ) dum2dd
            END IF
   
            IF ( print_info ) THEN
               PRINT '(A,A,A,I4,A)','Outputting 2D data ',name(1:8),', level = ',NINT(level/100.),' mb.'
            END IF

         END IF

      END DO all_2d

      !  Oh, what a beautiful morning, oh what a beautiful day.  Rogers and Hammerstein would 
      !  certainly wax ecstatic and lyrical to a successful file output.

      IF ( print_info ) THEN
         PRINT '(A,A,A)','Completed outputting data into intermediate format for file FILE_MMOUTP:', &
         current_date19(1:13),'.'
      END IF

      !  Since this file that we output is a "once per time period" deal, we need
      !  to shut it down for the next time that we enter here.

      CLOSE ( output_intermediate ) 

   END SUBROUTINE out_intermediate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE add_qv

      IMPLICIT NONE

      INTEGER :: i , j , ks , kp

      !  The purpose for this routine is to artificially fill in an array of zero
      !  for the missing Qv array.  Why is it missing?  This may be model output that
      !  was run as a dry case.  Actually, I don't think there is any other way.

      !  How large to allocate the 3d arrays for the sigma and pressure Qv.  They are the
      !  same size as the other 3d sigma and pressure arrays.

      i  =  alls_3d(nums_3d)%small_header%end_dims(1)
      j  =  alls_3d(nums_3d)%small_header%end_dims(2)
      ks =  alls_3d(nums_3d)%small_header%end_dims(3)
      kp =  allp_3d(nump_3d)%small_header%end_dims(3)

      !  ALLOCATE the two 3d arrays.

      ALLOCATE ( alls_3d(nums_3d+1)%array(i,j,ks) ) 
      ALLOCATE ( allp_3d(nump_3d+1)%array(i,j,kp) ) 

      !  Fill the sigma level Qv array with zeros.  The pressure level guy will be
      !  interpolated a bit later.  All we need to do is have the sigma level data 
      !  available.

       alls_3d(nums_3d+1)%array = 0.

      !  Fill in all of the header things.

      alls_3d(nums_3d+1)%small_header = alls_3d(nums_3d)%small_header
      allp_3d(nump_3d+1)%small_header = allp_3d(nump_3d)%small_header

      !  We have to be cross point staggered.

      alls_3d(nums_3d+1)%small_header%staggering(1:4) = 'C   '
      allp_3d(nump_3d+1)%small_header%staggering(1:4) = 'C   '

      !  Give this lovely data a new name.

      alls_3d(nums_3d+1)%small_header%name(1:8) = 'Q       '
      allp_3d(nump_3d+1)%small_header%name(1:8) = 'Q       '

      !  Since we have added an array for each of the sigma and pressure 3d array holders,
      !  we need to bump their respective counters.

      nums_3d = nums_3d + 1
      nump_3d = nump_3d + 1
       
   END SUBROUTINE add_qv

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE add_snowcovr

      IMPLICIT NONE

      INTEGER :: i , j , ks , kp

      !  The purpose for this routine is to artificially fill in an array of zero
      !  for the missing SNOWCOVR array.  Why is it missing?  This may be model output that
      !  was run with SLAB and IFSNOW=0.

      !  How large to allocate the 2d arrays for the sigma and pressure snowcovr.  They are the
      !  same size as the other 2d sigma and pressure arrays.

      i  =  alls_2d(nums_2d)%small_header%end_dims(1)
      j  =  alls_2d(nums_2d)%small_header%end_dims(2)

      !  ALLOCATE the two 2d arrays.

      ALLOCATE ( alls_2d(nums_2d+1)%array(i,j) ) 
      ALLOCATE ( allp_2d(nump_2d+1)%array(i,j) ) 

      !  Fill the sigma level Qv array with zeros.  The pressure level guy will be
      !  interpolated a bit later.  All we need to do is have the sigma level data 
      !  available.

       alls_2d(nums_2d+1)%array = 0.

      !  Fill in all of the header things.

      alls_2d(nums_2d+1)%small_header = alls_2d(nums_2d)%small_header
      allp_2d(nump_2d+1)%small_header = allp_2d(nump_2d)%small_header

      !  We have to be cross point staggered.

      alls_2d(nums_2d+1)%small_header%staggering(1:4) = 'C   '
      allp_2d(nump_2d+1)%small_header%staggering(1:4) = 'C   '

      !  Give this lovely data a new name.

      alls_2d(nums_2d+1)%small_header%name(1:8) = 'SNOWCOVR'
      allp_2d(nump_2d+1)%small_header%name(1:8) = 'SNOWCOVR'

      !  Since we have added an array for each of the sigma and pressure 2d array holders,
      !  we need to bump their respective counters.

      nums_2d = nums_2d + 1
      nump_2d = nump_2d + 1
       
   END SUBROUTINE add_snowcovr

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END MODULE all_io
