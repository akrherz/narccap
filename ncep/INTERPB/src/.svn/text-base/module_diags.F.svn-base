MODULE diags

   USE phys_consts
   USE util

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
   SUBROUTINE qvtorh ( q , t , p , imx , jmx , kxs , rh ) 

      IMPLICIT NONE

      INTEGER , INTENT(IN) :: imx , jmx , kxs
      REAL , INTENT(IN) , DIMENSION(:,:,:) :: q , t , p

      REAL , INTENT(OUT) , DIMENSION(:,:,:) :: rh

      !  Local variables.

      INTEGER :: i , j , k 
      REAL :: es , qs

      DO k = 1 , kxs 
         DO j = 1 , jmx - 1 
            DO i = 1 , imx - 1 
               es = svp1*10.*EXP(svp2*(t(i,j,k)-svpt0)/(t(i,j,k)-svp3))
               qs = eps*es/((0.01 *p(i,j,k))-es)
               rh(i,j,k) = MIN ( 100. , MAX ( 100.*q(i,j,k)/qs , 1. ) )
            END DO
         END DO
      END DO

   END SUBROUTINE qvtorh

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE surface_pressure ( pp , pstar , ptop , imx , jmx , kxs , psfc )

      IMPLICIT NONE

      !  Input variables.

      INTEGER , INTENT(IN)                    :: imx , jmx , kxs
      REAL    , INTENT(IN)                    :: ptop
      REAL    , INTENT(IN) , DIMENSION(:,:  ) :: pstar
      REAL    , INTENT(IN) , DIMENSION(:,:,:) :: pp

      !  Output variables.

      REAL    , INTENT(OUT) , DIMENSION(:,:)  :: psfc

      !  Local variables.

      INTEGER :: i , j

      DO j = 1 , jmx - 1 
         DO i = 1 , imx - 1 
            psfc(i,j) = pp(i,j,kxs) + pstar(i,j) + ptop
         END DO
      END DO

   END SUBROUTINE surface_pressure

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE surface_pressure2 ( pp , t , q , ter , pstar , sigmah , ptop , p0 , tlp , ts0 , &
                                  imx , jmx , kxs , psfc )

      IMPLICIT NONE

      !  Input variables.

      INTEGER , INTENT(IN)                    :: imx , jmx , kxs
      REAL    , INTENT(IN)                    :: ptop , p0 , tlp , ts0
      REAL    , INTENT(IN) , DIMENSION(:    ) :: sigmah
      REAL    , INTENT(IN) , DIMENSION(:,:  ) :: ter , pstar
      REAL    , INTENT(IN) , DIMENSION(:,:,:) :: pp , t , q

      !  Output variables.

      REAL    , INTENT(OUT) , DIMENSION(:,:)  :: psfc

      !  Local variables.

      INTEGER :: i , j
      REAL , DIMENSION(imx,jmx) :: tv , t0 , pr0

      !  Base state / reference values at the lowest half sigma level, oh, and let's pick
      !  up a couple boxes of virtual temperature at the store in case we need some.

      DO j = 1, jmx - 1
         DO i = 1, imx - 1
            pr0(i,j) = pstar(i,j) * sigmah(kxs) + ptop
            t0(i,j) = ts0 + tlp * LOG ( pr0(i,j) / p0 )
            tv(i,j) = t(i,j,kxs) * ( 1 + 0.608 * q(i,j,kxs) ) 
         END DO
      END DO

      !  Assuming density perturbation remains constant in lowest half sigma layer, compute
      !  surface pressure as inverse of pp diagnostic as originally performed in INTERPF.

      DO j = 1 , jmx - 1 
         DO i = 1 , imx - 1 
            psfc(i,j) = pp(i,j,kxs) * ( pstar(i,j) + ptop ) / pr0(i,j) + &
                       ( pstar(i,j) + ptop )                           - &
                       ( tv(i,j) - t0(i,j) ) / t(i,j,kxs) * ( pstar(i,j) + ptop - pr0(i,j) ) 
         END DO
      END DO

   END SUBROUTINE surface_pressure2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE fat1001 ( fsfc , psfc1001 , ter , name , imx , jmx , crsdot , fsfc1001 ) 

      IMPLICIT NONE

      INTEGER , INTENT(IN) :: imx , jmx , crsdot
      CHARACTER (LEN=8) , INTENT(IN) :: name
      REAL , INTENT(IN) , DIMENSION(:,:) :: ter
      REAL , INTENT(IN) , DIMENSION(:,:) :: psfc1001
      REAL , INTENT(IN) , DIMENSION(imx,jmx) :: fsfc

      REAL , INTENT(OUT) , DIMENSION(:,:) :: fsfc1001

      !  Local variables.

      INTEGER :: i , j

      IF      ( ( name(1:8) .EQ. 'U       ' ) .OR. ( name(1:8) .EQ. 'V       ' ) ) THEN
         DO j = 1 , jmx 
            DO i = 1 , imx
               fsfc1001(i,j) = fsfc(i,j)
            END DO
         END DO
      ELSE IF   ( name(1:8) .EQ. 'T       ' ) THEN
         DO j = 1 , jmx - 1
            DO i = 1 , imx - 1
               fsfc1001(i,j) = fsfc(i,j) + 0.0065 * ter(i,j)
            END DO
         END DO
      ELSE IF   ( name(1:8) .EQ. 'H       ' ) THEN
         DO j = 1 , jmx - 1
            DO i = 1 , imx - 1
               fsfc1001(i,j) = 0.
            END DO
         END DO
      ELSE IF   ( name(1:3) .EQ. 'P3D'      ) THEN
         DO j = 1 , jmx
            DO i = 1 , imx
               fsfc1001(i,j) = psfc1001(i,j)
            END DO
         END DO
      ELSE
         DO j = 1 , jmx - 1 
            DO i = 1 , imx - 1 
               fsfc1001(i,j) = fsfc(i,j)
            END DO
         END DO
      END IF

   END SUBROUTINE fat1001

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE pat1001 ( slp , imx , jmx , slp1001c , slp1001d ) 

      IMPLICIT NONE

      INTEGER , INTENT(IN) :: imx , jmx
      REAL , INTENT(IN) , DIMENSION(:,:) :: slp

      REAL , INTENT(OUT) , DIMENSION(:,:) :: slp1001c , slp1001d

      !  Local variables.

      INTEGER :: i , j

      DO j = 1 , jmx - 1
         DO i = 1 , imx - 1 
            slp1001c(i,j) = MIN ( 100100. , slp(i,j) + 100. ) 
            slp1001d(i,j) = slp1001c(i,j)
         END DO
      END DO

     
      CALL crs2dot ( slp1001d , imx , jmx ) 

   END SUBROUTINE pat1001

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE seapres  ( h , t , q , p , ter , psfc , imx , jmx , kxs , slp )

      IMPLICIT NONE

      INTEGER , INTENT(IN) :: imx , jmx , kxs
      REAL, INTENT(IN), DIMENSION(:,:) :: psfc , ter
      REAL, INTENT(IN), DIMENSION(:,:,:) :: h , t , q , p

      REAL, INTENT(OUT), DIMENSION(:,:) :: slp

      !  Local variables.

      INTEGER :: i , j , k
      INTEGER , DIMENSION(imx,jmx) ::  klevel

      REAL , DIMENSION(imx,jmx) :: tsfc , t0
      REAL :: p100mbup , t100mbup , h100mbup

      REAL , PARAMETER :: xterm  = gamma * r / g
      REAL , PARAMETER :: pconst = 1.0e4

      LOGICAL :: l1 , l2 , l3

      !  Initialize location for 100 mb above surface.

      DO j = 1 , jmx - 1
         DO i = 1 , imx - 1
            klevel(i,j) = 0
         END DO
      END DO

     !  Find k level in sigma space that traps 100 mn above the surface.

      DO j = 1 , jmx - 1
         DO i = 1 , imx - 1
            DO k = kxs , 2 , -1
              IF ( ( p(i,j,k  ) .GE. psfc(i,j)-pconst ) .AND. &
                   ( p(i,j,k-1) .LT. psfc(i,j)-pconst ) ) THEN
                 klevel(i,j) = k
              END IF
            END DO
         END DO
      END DO

      !  Did we find what we needed?

      DO j = 1 , jmx - 1
         DO i = 1 , imx - 1
            IF ( klevel(i,j) .LT. 1 ) THEN
               PRINT '(A)','Error in finding level 100 mb above surface.'
               PRINT '(A,3I4,A)','(i,j,k)=',i,j,k,'.'
               PRINT *,'pressure string = ',p(i,j,kxs:1:-1)
               PRINT '(A,F9.1,A)','Surface pressure = ',psfc(i,j),'hPa.'
               STOP 'Error finding 100 mb above for slp'
            END IF 
         END DO
      END DO

      !  Compute surface temp (tsfc)and sea level temp (t0).

      DO j = 1 , jmx - 1
         DO i = 1 , imx - 1
            k = klevel(i,j)
            p100mbup = psfc(i,j)-pconst
            t100mbup = t(i,j,k-1)*(1+0.608*q(i,j,k-1)) - &
                       (t(i,j,k-1)*(1+0.608*q(i,j,k-1))-t(i,j,k)*(1+0.608*q(i,j,k  )))* &
                       LOG(p100mbup/p(i,j,k-1))/LOG(p(i,j,k)/p(i,j,k-1))
            tsfc(i,j) = t100mbup * ( psfc(i,j) /  p100mbup ) ** xterm
!           h100mbup = ( h(i,j,k  )*(LOG(p100mbup)-LOG(p(i,j,k-1))) + &
!                        h(i,j,k-1)*(LOG(p(i,j,k))-LOG(p100mbup)) ) / &
!                        ( LOG(p(i,j,k))-LOG(p(i,j,k-1)))
            h100mbup = ter(i,j) - r/g*LOG(p100mbup/psfc(i,j))*(0.5*(t100mbup + tsfc(i,j)))
            t0(i,j) = t100mbup + gamma * h100mbup
         END DO
      END DO

      !  A correction to the sea-level temperature, in case it is too warm.

      DO j = 1 , jmx - 1
         DO i = 1 , imx - 1
            l1 = t0(i,j)   .LT. tc
            l2 = tsfc(i,j) .LE. tc
            l3 = .NOT. l1
            IF ( .NOT. ( l1 .AND. l2 ) ) THEN
               IF      ( l2 .AND. l3 ) THEN
                  t0(i,j) = tc
               ELSE IF ( ( .NOT. l2 ) .AND. l3 ) THEN
                  t0(i,j) = tc - 0.005 * ( tsfc(i,j) - tc ) ** 2
               END IF
            END IF
         END DO
      END DO

      !  Compute slp.

      DO j = 1 , jmx - 1
         DO i = 1 , imx - 1
            slp(i,j) = psfc(i,j) * EXP ( 2.*g*ter(i,j)/(r*(tsfc(i,j)+t0(i,j))))
         END DO
      END DO

   END SUBROUTINE seapres


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE height ( psc , sigma , p00 , ts0 , tlp , tiso , ptop , imx , jmx , kxs , z ) 

      INTEGER , INTENT(IN) :: imx , jmx , kxs
      REAL , INTENT(IN) , DIMENSION(:) :: sigma
      REAL , INTENT(IN) , DIMENSION(:,:) :: psc
      REAL , INTENT(IN) :: p00 , ts0 , tlp , tiso , ptop 

      REAL , INTENT(OUT) , DIMENSION(:,:,:) :: z

      !  Local variables.

      INTEGER :: i , j , k 
      REAL :: piso , pr0 , alnp , ziso

!     DO k = 1 , kxs
!        DO j = 1 , jmx - 1
!           DO i = 1 , imx - 1 
!              z(i,j,k) = - ( r * tlp / ( 2. * g ) * ( LOG ( ( psc(i,j) * sigma(k) + ptop ) / p00 ) ) **2 + &
!                             r * ts0 /        g   *   LOG ( ( psc(i,j) * sigma(k) + ptop ) / p00 )         )
!           END DO
!        END DO
!     END DO

      piso = p00 * EXP ( ( tiso - ts0 ) / tlp )

      DO k = 1, kxs
         DO j = 1, jmx - 1
            DO i = 1, imx - 1
               pr0 = psc(i,j) * sigma(k) + ptop
               IF ( pr0 .LT. piso ) THEN
                  alnp     = LOG ( piso / p00 )
                  ziso     = - ( r * tlp / ( 2. * g ) * alnp * alnp + r * ts0 / g * alnp )
                  z(i,j,k) = ziso - r * tiso / g * LOG ( pr0 / piso )
               ELSE
                  alnp     = LOG ( pr0 / p00 )
                  z(i,j,k) = - ( r * tlp / ( 2. * g ) * alnp * alnp + r * ts0 / g * alnp )
               END IF
            END DO
         END DO
      END DO

   END SUBROUTINE height

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE total_pressure ( pp , pstar , sigma , ptop , ix , jx , kx , p )

      IMPLICIT NONE

      !  Input variables.

      INTEGER , INTENT(IN)                    :: ix , jx , kx
      REAL    , INTENT(IN)                    :: ptop
      REAL    , INTENT(IN) , DIMENSION(:    ) :: sigma
      REAL    , INTENT(IN) , DIMENSION(:,:  ) :: pstar
      REAL    , INTENT(IN) , DIMENSION(:,:,:) :: pp

      !  Output variables.

      REAL    , INTENT(OUT) , DIMENSION(:,:,:)  :: p

      !  Local variables.

      INTEGER :: i , j , k

      DO k = 1 , kx
         DO j = 1 , jx - 1 
            DO i = 1 , ix - 1
               p(i,j,k) = pp(i,j,k) + pstar(i,j)*sigma(k) + ptop
            END DO
         END DO
      END DO

   END SUBROUTINE total_pressure

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END MODULE diags
