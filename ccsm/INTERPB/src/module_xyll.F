MODULE xyll

   USE map_utils

   TYPE(proj_info) :: projx

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE lat_lon_dot ( lat , lon , imx , jmx , &
                            ds , xlonc , truelat1 , truelat2 , &
                            projection , lat1c , lon1c ) 

      IMPLICIT NONE

      !  Input variables.

      REAL , INTENT(IN) :: ds , xlonc , truelat1 , truelat2
      REAL , INTENT(IN) :: lat1c , lon1c
      INTEGER , INTENT(IN) :: imx , jmx , projection

      !  Output variables.

      REAL , INTENT(OUT) , DIMENSION(:,:) :: lat , lon

      !  Local variables.

      INTEGER :: i , j

      IF      ( projection .EQ. 1 ) THEN
         CALL map_set(proj_lc  ,lat1c,lon1c,ds*1000.,xlonc,truelat1,truelat2,jmx-1,imx-1,projx)
         DO j = 1 , jmx
            DO i = 1 , imx
               CALL ij_to_latlon ( projx , REAL(j)-0.5 , REAL(i)-0.5 , lat(i,j) , lon(i,j) )
            END DO
         END DO
      ELSE IF ( projection .EQ. 2 ) THEN
         CALL map_set(proj_ps  ,lat1c,lon1c,ds*1000.,xlonc,truelat1,truelat2,jmx-1,imx-1,projx)
         DO j = 1 , jmx
            DO i = 1 , imx
               CALL ij_to_latlon ( projx , REAL(j)-0.5 , REAL(i)-0.5 , lat(i,j) , lon(i,j) )
            END DO
         END DO
      ELSE IF ( projection .EQ. 3 ) THEN
         CALL map_set(proj_merc,lat1c,lon1c,ds*1000.,xlonc,truelat1,truelat2,jmx-1,imx-1,projx)
         DO j = 1 , jmx
            DO i = 1 , imx
               CALL ij_to_latlon ( projx , REAL(j)-0.5 , REAL(i)-0.5 , lat(i,j) , lon(i,j) )
            END DO
         END DO
      END IF

   END SUBROUTINE lat_lon_dot

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END MODULE xyll
