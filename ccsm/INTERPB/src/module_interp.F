MODULE interp

   USE phys_consts

CONTAINS

   SUBROUTINE vinterp ( sdata , imx , jmx , kxs , crsdot , p , &
                        pat1001 , fat1001 , pdata , kxp , int_meth , extrap_up , extrap_low , print_info , pres )

      IMPLICIT NONE

      INTEGER , INTENT(IN) :: imx , jmx , kxs , kxp , crsdot , int_meth , extrap_up , extrap_low
      REAL , INTENT(IN) , DIMENSION(imx,jmx,kxs) :: sdata , p
      REAL , INTENT(IN) , DIMENSION(imx,jmx) :: pat1001 , fat1001
      REAL , INTENT(IN) , DIMENSION(kxp) :: pres 
      LOGICAL , INTENT(IN) :: print_info

      REAL , INTENT(OUT) , DIMENSION(imx,jmx,kxp) :: pdata

      !  Local variables.

      INTEGER :: i , j , ks , kp 
      INTEGER , DIMENSION(imx,jmx) :: kstart , kend

      !  Find valid starting location.  The kstart in each column is the index
      !  of the pressure level data which has the largest pressure (lowest level)
      !  that is trapped by sigma levels.

      DO j = 1 , jmx - crsdot
         DO i = 1 , imx - crsdot
            kstart(i,j) = 1
         END DO
      END DO

      DO j = 1 , jmx - crsdot
         DO i = 1 , imx - crsdot
            pressure_search_start : DO kp = 2 , kxp
               each_sigma_start : DO ks = kxs , 2 , -1
                  IF ( ( p(i,j,ks  ) .GE. pres(kp) ) .AND. &
                       ( p(i,j,ks-1) .LT. pres(kp) ) ) THEN
                     kstart(i,j) = kp
                     EXIT pressure_search_start
                  ELSE IF ( p(i,j,ks  ) .LT. pres(kp) ) THEN
                     EXIT each_sigma_start
                  END IF
               END DO each_sigma_start
            END DO pressure_search_start
         END DO
      END DO

      !  Find valid ending location.  The kend in each column is the index
      !  of the pressure level data which has the smallest pressure (highest level)
      !  that is trapped by sigma levels.

      DO j = 1 , jmx - crsdot
         DO i = 1 , imx - crsdot
            kend(i,j) = kxp
         END DO
      END DO

      DO j = 1 , jmx - crsdot
         DO i = 1 , imx - crsdot
            pressure_search_end : DO kp = kxp , 2 , -1
               each_sigma_end : DO ks = 2 , kxs
                  IF ( ( p(i,j,ks  ) .GE. pres(kp) ) .AND.  &
                       ( p(i,j,ks-1) .LT. pres(kp) ) ) THEN
                     kend(i,j) = kp
                     EXIT pressure_search_end
                  ELSE IF ( p(i,j,ks-1) .GT. pres(kp) ) THEN
                     EXIT each_sigma_end
                  END IF
               END DO each_sigma_end
            END DO pressure_search_end
         END DO
      END DO
         
      !  Compute trapped interpolations.
      
      IF      ( int_meth .EQ. 1 ) THEN
         DO j = 1 , jmx - crsdot
            DO i = 1 , imx - crsdot
               DO kp = kstart(i,j) , kend(i,j)
                  sigma_search_1 : DO ks = kxs , 2 , -1
                     IF ( ( p(i,j,ks  ) .GE. pres(kp) ) .AND. &
                          ( p(i,j,ks-1) .LT. pres(kp) ) ) THEN
                        pdata(i,j,kp) = ( sdata(i,j,ks-1)*( p(i,j,ks  ) - pres(kp) ) + &
                                          sdata(i,j,ks  )*( pres(kp)  - p(i,j,ks-1) ) )  / &
                                          ( p(i,j,ks  ) -  p(i,j,ks-1) )
                        EXIT sigma_search_1
                     END IF
                  END DO sigma_search_1
               END DO
            END DO
         END DO
      ELSE IF ( int_meth .EQ. 2 ) THEN
         DO j = 1 , jmx - crsdot
            DO i = 1 , imx - crsdot
               DO kp = kstart(i,j) , kend(i,j)
                  sigma_search_2 : DO ks = kxs , 2 , -1
                     IF ( ( p(i,j,ks  ) .GE. pres(kp) ) .AND. &
                          ( p(i,j,ks-1) .LT. pres(kp) ) ) THEN
                        pdata(i,j,kp) = ( sdata(i,j,ks-1)*( LOG(p(i,j,ks  )/pres(kp)) ) + &
                                          sdata(i,j,ks  )*( LOG(pres(kp)/p(i,j,ks-1)) ) )  / &
                                          ( LOG(p(i,j,ks  )/p(i,j,ks-1)) )
                        EXIT sigma_search_2
                     END IF
                  END DO sigma_search_2
               END DO
            END DO
         END DO
      ELSE IF ( int_meth .EQ. 3 ) THEN
         DO j = 1 , jmx - crsdot
            DO i = 1 , imx - crsdot
               DO kp = kstart(i,j) , kend(i,j)
                  sigma_search_3 : DO ks = kxs , 2 , -1
                     IF ( ( p(i,j,ks  ) .GE. pres(kp) ) .AND. &
                          ( p(i,j,ks-1) .LT. pres(kp) ) ) THEN
                        pdata(i,j,kp) = ( sdata(i,j,ks-1)*( p(i,j,ks  )**rovcp - pres(kp)**rovcp ) + &
                                          sdata(i,j,ks  )*( pres(kp)**rovcp  - p(i,j,ks-1)**rovcp ) )  / &
                                          ( p(i,j,ks  )**rovcp -  p(i,j,ks-1)**rovcp )
                        EXIT sigma_search_3
                     END IF
                  END DO sigma_search_3
               END DO
            END DO
         END DO
      END IF
         
      !  Do extrapolations from bottom up to the level below the kstart.
      
      IF      ( ( int_meth .EQ. 1 ) .AND. ( extrap_low .EQ. 1 ) ) THEN
         DO j = 1 , jmx - crsdot
            DO i = 1 , imx - crsdot
               DO kp = 2 , kstart(i,j) - 1
                  pdata(i,j,kp) = sdata(i,j,kxs)
               END DO
            END DO
         END DO
      ELSE IF ( ( int_meth .EQ. 1 ) .AND. ( extrap_low .EQ. 2 ) ) THEN
         DO j = 1 , jmx - crsdot
            DO i = 1 , imx - crsdot
               DO kp = 2 , kstart(i,j) - 1
!                 pdata(i,j,kp) = ( sdata(i,j,kxs-1)*( p(i,j,kxs  ) - pres(kp) ) + &
!                                   sdata(i,j,kxs  )*( pres(kp)  - p(i,j,kxs-1) ) )  / &
!                                   ( p(i,j,kxs  ) -  p(i,j,kxs-1) )
                  pdata(i,j,kp) = ( fat1001(i,j)    *( p(i,j,kxs  ) - pres(kp) ) + &
                                    sdata(i,j,kxs  )*( pres(kp)  - pat1001(i,j)) ) / &
                                    ( p(i,j,kxs  ) -  pat1001(i,j) )
               END DO
            END DO
         END DO
      ELSE IF ( ( int_meth .EQ. 2 ) .AND. ( extrap_low .EQ. 1 ) ) THEN
         DO j = 1 , jmx - crsdot
            DO i = 1 , imx - crsdot
               DO kp = 2 , kstart(i,j) - 1
                  pdata(i,j,kp) = sdata(i,j,kxs)
               END DO
            END DO
         END DO
      ELSE IF ( ( int_meth .EQ. 2 ) .AND. ( extrap_low .EQ. 2 ) ) THEN
         DO j = 1 , jmx - crsdot
            DO i = 1 , imx - crsdot
               DO kp = 2 , kstart(i,j) - 1
                  pdata(i,j,kp) = ( fat1001(i,j)    *( LOG(p(i,j,kxs  )/pres(kp)) ) + &
                                    sdata(i,j,kxs  )*( LOG(pres(kp)/pat1001(i,j)) ) )  / &
                                    ( LOG(p(i,j,kxs  )/pat1001(i,j)) )
               END DO
            END DO
         END DO
      ELSE IF ( ( int_meth .EQ. 3 ) .AND. ( extrap_low .EQ. 1 ) ) THEN
         DO j = 1 , jmx - crsdot
            DO i = 1 , imx - crsdot
               DO kp = 2 , kstart(i,j) - 1
                  pdata(i,j,kp) = sdata(i,j,kxs)
               END DO
            END DO
         END DO
      ELSE IF ( ( int_meth .EQ. 3 ) .AND. ( extrap_low .EQ. 2 ) ) THEN
         DO j = 1 , jmx - crsdot
            DO i = 1 , imx - crsdot
               DO kp = 2 , kstart(i,j) - 1
                  pdata(i,j,kp) = ( fat1001(i,j)    *( p(i,j,kxs  )**rovcp - pres(kp)**rovcp ) + &
                                    sdata(i,j,kxs  )*( pres(kp)**rovcp  - pat1001(i,j)**rovcp) )  / &
                                    ( p(i,j,kxs  )**rovcp -  pat1001(i,j)**rovcp )
               END DO
            END DO
         END DO
      END IF
         
      !  Do extrapolations from the level above kend to the top requested pressure level.
      
      IF      ( ( int_meth .EQ. 1 ) .AND. ( extrap_up .EQ. 1 ) ) THEN
         DO j = 1 , jmx - crsdot
            DO i = 1 , imx - crsdot
               DO kp = kend(i,j) + 1 , kxp
                  pdata(i,j,kp) = sdata(i,j,1)
               END DO
            END DO
         END DO
      ELSE IF ( ( int_meth .EQ. 1 ) .AND. ( extrap_up .EQ. 2 ) ) THEN
         DO j = 1 , jmx - crsdot
            DO i = 1 , imx - crsdot
               DO kp = kend(i,j) + 1 , kxp
                  pdata(i,j,kp) = ( sdata(i,j,1)*( p(i,j,2) - pres(kp) ) + &
                                    sdata(i,j,2)*( pres(kp)  - p(i,j,1) ) )  / &
                                    ( p(i,j,2) -  p(i,j,1) )
               END DO
            END DO
         END DO
      ELSE IF ( ( int_meth .EQ. 2 ) .AND. ( extrap_up .EQ. 1 ) ) THEN
         DO j = 1 , jmx - crsdot
            DO i = 1 , imx - crsdot
               DO kp = kend(i,j) + 1 , kxp
                  pdata(i,j,kp) = sdata(i,j,1)
               END DO
            END DO
         END DO
      ELSE IF ( ( int_meth .EQ. 2 ) .AND. ( extrap_up .EQ. 2 ) ) THEN
         DO j = 1 , jmx - crsdot
            DO i = 1 , imx - crsdot
               DO kp = kend(i,j) + 1 , kxp
                  pdata(i,j,kp) = ( sdata(i,j,1)*( LOG(p(i,j,2)/pres(kp)) ) + &
                                    sdata(i,j,2)*( LOG(pres(kp)/p(i,j,1)) ) )  / &
                                    ( LOG(p(i,j,2)/p(i,j,1)) )
               END DO
            END DO
         END DO
      ELSE IF ( ( int_meth .EQ. 3 ) .AND. ( extrap_up .EQ. 1 ) ) THEN
         DO j = 1 , jmx - crsdot
            DO i = 1 , imx - crsdot
               DO kp = kend(i,j) + 1 , kxp
                  pdata(i,j,kp) = sdata(i,j,1)
               END DO
            END DO
         END DO
      ELSE IF ( ( int_meth .EQ. 3 ) .AND. ( extrap_up .EQ. 2 ) ) THEN
         DO j = 1 , jmx - crsdot
            DO i = 1 , imx - crsdot
               DO kp = kend(i,j) + 1 , kxp
                  pdata(i,j,kp) = ( sdata(i,j,1)*( p(i,j,2)**rovcp - pres(kp)**rovcp ) + &
                                    sdata(i,j,2)*( pres(kp)**rovcp  - p(i,j,1)**rovcp ) )  / &
                                    ( p(i,j,2)**rovcp -  p(i,j,1)**rovcp )
               END DO
            END DO
         END DO
      END IF

      IF ( print_info ) THEN
         PRINT '(A)','  Pressure           Value '
         PRINT '(A)','    (Pa)'
         PRINT '(A)','------------     ------------'
         DO kp = 2 , kxp
            PRINT '(F12.5,5X,G12.5)',pres(kp),pdata(imx/2,jmx/2,kp)
         END DO
      END IF

      !  Now assign the lowest level sigma data to the "surface" location in the pressure level array.
      !  This is a bad assumption for the height field, but we fix it in the main program.

      DO j = 1 , jmx - crsdot
         DO i = 1 , imx - crsdot
            pdata(i,j,1) = sdata(i,j,kxs)
         END DO
      END DO

   END SUBROUTINE vinterp

END MODULE interp
