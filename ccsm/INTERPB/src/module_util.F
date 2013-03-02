MODULE util

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
   SUBROUTINE check (kxp, newcoord)
     
   !  Check that the number of pressure layers is correct
   !  and in numerically descending order, given the namelist values.
     
      IMPLICIT NONE

      !  Input variables.
    
      INTEGER                      :: kxp
    
      REAL                         :: newcoord   ( : )
    
      !  Local variables.
    
      INTEGER                      :: k
      INTEGER                      :: kill_it = 0
    
      !  Check first pressure value newcoord(1) should be 100000
     
      IF ( newcoord(2) .LT. 100000 ) THEN
        PRINT '(A,F6.4)','Pressure value at bottom is ', newcoord(2),' Pa.'
        kill_it = kill_it + 1
      ENDIF
     
      !  Descending input pressure value order check.
    
      DO k = 1, kxp-1
         IF (newcoord(k+1) .GE. newcoord(k)) THEN
            PRINT '(A)','Pressure values are not in descending order.'
            kill_it = kill_it + 1
         ENDIF
      END DO
     
      !  Negative check.
    
      DO k = 1, kxp
         IF (newcoord(k) .LT. 0 ) THEN
            PRINT '(A)','Negative pressure requests are traditionally not obliged.'
            kill_it = kill_it + 1
         ENDIF
      END DO
     
      !  If any test was bad, stop.
     
      IF ( kill_it .GT. 0 ) THEN 
         PRINT '(A)','Troubles with input pressure data (pressure_bu_no_sfc_Pa).  Fix the namelist.input file.'
         STOP 'pressure_woes'
      END IF
    
   END SUBROUTINE CHECK

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE crs2dot(field,dim1,dim2)
   
      IMPLICIT NONE

      INTEGER :: dim1 , dim2
      REAL , DIMENSION(dim1,dim2) :: field,dummy
      
      dummy(2:dim1-1,2:dim2-1)           = ( field(1:dim1-2,1:dim2-2) + &
                                             field(1:dim1-2,2:dim2-1) + &
                                             field(2:dim1-1,1:dim2-2) + &
                                             field(2:dim1-1,2:dim2-1) ) * 0.25
   
      dummy(2:dim1-1,1:dim2:dim2-1)      = ( field(1:dim1-2,1:dim2-1:dim2-2) + &
                                             field(2:dim1-1,1:dim2-1:dim2-2) ) * 0.5
   
      dummy(1:dim1:dim1-1,2:dim2-1)      = ( field(1:dim1-1:dim1-2,1:dim2-2) + &
                                             field(1:dim1-1:dim1-2,2:dim2-1) ) * 0.5
   
      dummy(1:dim1:dim1-1,1:dim2:dim2-1) =   field(1:dim1-1:dim1-2,1:dim2-1:dim2-2)
   
      field                              =   dummy
   
   END SUBROUTINE crs2dot

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE dot2crs(field,dim1,dim2)

      IMPLICIT NONE

      INTEGER :: dim1 , dim2
      REAL , DIMENSION(dim1,dim2) :: field
      INTEGER :: i , j 
      DO j = 1 , dim2 - 1
         DO i = 1 , dim1 - 1
            field(i,j) = ( field(i  ,j  ) + & 
                           field(i+1,j  ) + & 
                           field(i  ,j+1) + & 
                           field(i+1,j+1) ) * 0.25
         END DO
      END DO

   END SUBROUTINE dot2crs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE clean_rh ( rh , iew , jns , rh_min , rh_max )

      IMPLICIT NONE
   
      INTEGER                        :: iew , jns
      REAL , DIMENSION ( iew , jns ) :: rh
      REAL                           :: rh_min , rh_max
   
      WHERE ( rh .GT. rh_max ) rh = rh_max
      WHERE ( rh .LT. rh_min ) rh = rh_min
   
   END SUBROUTINE clean_rh

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE fill (f, ix, jx, imx, jmx, ifirst, ilast, jfirst, jlast)

      IMPLICIT NONE

      INTEGER                     :: I
      INTEGER                     :: IFIRST
      INTEGER                     :: ILAST
      INTEGER                     :: IMX
      INTEGER                     :: IX
      INTEGER                     :: J
      INTEGER                     :: JFIRST
      INTEGER                     :: JLAST
      INTEGER                     :: JMX
      INTEGER                     :: JX

      REAL                        :: F           ( : , : )

      DO j = jfirst, jlast
         DO i = 1, ifirst - 1
            f(i,j) = f(ifirst,j)
         END DO
         DO i = ilast + 1, imx
            f(i,j) = f(ilast,j)
         END DO
      END DO

      DO j = 1, jfirst - 1
         f(:,j) = f(:,jfirst)
      END DO
      DO j = jlast + 1, jmx
         f(:,j) = f(:,jlast)
      END DO

   END SUBROUTINE fill

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE fill3d ( dum3d , imx , jmx , kxp , staggering )

      IMPLICIT NONE

      INTEGER , INTENT(IN) :: imx , jmx , kxp
      REAL , DIMENSION(imx,jmx,kxp) , INTENT(INOUT) :: dum3d
      CHARACTER (LEN=4) , INTENT(IN) :: staggering

      INTEGER :: i , j , k

      IF ( staggering(1:1) .EQ. 'C' ) THEN

         DO k = 1 , kxp
            DO j = 1 , jmx-1
               dum3d(imx,j,k) = dum3d(imx-1,j,k) 
            END DO 
            DO i = 1 , imx
               dum3d(i,jmx,k) = dum3d(i,jmx-1,k) 
            END DO 
         END DO

      END IF

   END SUBROUTINE fill3d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE fill2d ( dum2d , imx , jmx , staggering )

      IMPLICIT NONE

      INTEGER , INTENT(IN) :: imx , jmx 
      REAL , DIMENSION(imx,jmx) , INTENT(INOUT) :: dum2d
      CHARACTER (LEN=4) , INTENT(IN) :: staggering

      INTEGER :: i , j

      IF ( staggering(1:1) .EQ. 'C' ) THEN

         DO j = 1 , jmx-1
            dum2d(imx,j) = dum2d(imx-1,j) 
         END DO 
         DO i = 1 , imx
            dum2d(i,jmx) = dum2d(i,jmx-1) 
         END DO 

      END IF

   END SUBROUTINE fill2d


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE invert2 (dum3d, imx, jmx, kx)

      IMPLICIT NONE

      INTEGER                     :: I
      INTEGER                     :: IMX
      INTEGER                     :: J
      INTEGER                     :: JMX
      INTEGER                     :: K
      INTEGER                     :: KK
      INTEGER                     :: KX

      REAL                        :: DUM3D       ( : , : , : )
      REAL                        :: OUT

      DO K = 1, KX/2
         kk = kx - k + 1
         DO j = 1 , jmx
            DO i = 1 , imx
               out = dum3d(i,j,k)
               dum3d(i,j,k) = dum3d(i,j,kk)
               dum3d(i,j,kk) = out
            END DO
         END DO
      END DO

   END SUBROUTINE invert2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE smoother_desmoother ( slab , imx , jmx , passes , crsdot )
   
      IMPLICIT NONE
   
      INTEGER                        :: imx , jmx , passes , crsdot
      REAL , DIMENSION ( imx , jmx ) :: slab , & 
                                        slabnew
   
      REAL , DIMENSION ( 2 )         :: xnu
      INTEGER                        :: i , j , loop , n 
   
      xnu  =  (/ 0.50 , -0.52 /)
   
      !  The odd number passes of this are the "smoother", the even
      !  number passes are the "de-smoother" (note the differnt signs on xnu).
   
      smoothing_passes : DO loop = 1 , passes * 2
   
         n  =  2 - MOD ( loop , 2 )
    
         DO i = 2 , imx - 1 - crsdot
            DO j = 2 , jmx - 1 - crsdot
               slabnew(i,j) = slab(i,j) + xnu(n) *  & 
               ((slab(i,j+1) + slab(i,j-1)) * 0.5-slab(i,j))
            END DO
         END DO
    
         DO i = 2 , imx - 1 - crsdot
            DO j = 2 , jmx - 1 - crsdot
               slab(i,j) = slabnew(i,j)
            END DO
         END DO
    
         DO j = 2 , jmx - 1 - crsdot
            DO i = 2 , imx - 1 - crsdot
               slabnew(i,j) = slab(i,j) + xnu(n) *  &
               ((slab(i+1,j) + slab(i-1,j)) * 0.5-slab(i,j))
            END DO
         END DO
    
         DO i = 2 , imx - 1 - crsdot
            DO j = 2 , jmx - 1 - crsdot
               slab(i,j) = slabnew(i,j)
            END DO
         END DO
    
      END DO smoothing_passes
   
   END SUBROUTINE smoother_desmoother

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE smooth_5 ( field , iew , jns , passes , crsdot )

      IMPLICIT NONE
   
      INTEGER                        :: iew , jns , &
                                        passes    , &
                                        crsdot
      REAL , DIMENSION ( iew , jns ) :: field
   
      REAL , DIMENSION ( iew , jns ) :: temp
      INTEGER                        :: i , j , num_passes
   
      !  How may passes of this smoother are we using.
   
      smoothing_passes : DO num_passes = 1 , passes
   
         !  Apply 5-point stencil smoother on interior of the domain.
      
         DO j = 2 , jns - 1 - crsdot
            DO i = 2 , iew - 1 - crsdot
               temp(i,j) = ( field(i  ,j  ) * 4. +  & 
                             field(i+1,j  )      +  & 
                             field(i-1,j  )      +  & 
                             field(i  ,j+1)      +  & 
                             field(i  ,j-1)      )  * 0.125
            END DO
         END DO
   
         !  Apply 3-point stencil smoother on the boundaries.
      
         i = 1
         DO j = 2 , jns - 1 - crsdot
            temp(i,j) = ( field(i  ,j  ) * 2. +  & 
                          field(i  ,j+1)      +  & 
                          field(i  ,j-1)      )  * 0.25
         END DO
   
         i = iew - crsdot
         DO j = 2 , jns - 1 - crsdot
            temp(i,j) = ( field(i  ,j  ) * 2. +  & 
                          field(i  ,j+1)      +  & 
                          field(i  ,j-1)      )  * 0.25
         END DO
      
         j = 1
         DO i = 2 , iew - 1 - crsdot
            temp(i,j) = ( field(i  ,j  ) * 2. +  & 
                          field(i+1,j  )      +  & 
                          field(i-1,j  )      ) * 0.25
         END DO
      
         j = jns - crsdot
         DO i = 2 , iew - 1 - crsdot
            temp(i,j) = ( field(i  ,j  ) * 2. +  & 
                          field(i+1,j  )      +  & 
                          field(i-1,j  )      ) * 0.25
         END DO
      
         !  Store smoothed field back into original array.
      
         DO j = 2 , jns - 1 - crsdot
            DO i = 2 , iew - 1 - crsdot
               field(i,j) = temp(i,j)
            END DO
         END DO
      
         !  Store smoothed boundary field back into original array.
      
         DO j = 2 , jns - 1 - crsdot
            field(1         ,j) = temp(1         ,j)
            field(iew-crsdot,j) = temp(iew-crsdot,j)
         END DO
      
         DO i = 2 , iew - 1 - crsdot
            field(i,1         ) = temp(i,1         )
            field(i,jns-crsdot) = temp(i,jns-crsdot)
         END DO
   
      END DO smoothing_passes
   
   END SUBROUTINE smooth_5

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END MODULE util
