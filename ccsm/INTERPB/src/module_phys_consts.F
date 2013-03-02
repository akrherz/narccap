MODULE phys_consts

      REAL , PARAMETER :: cp          =  1004.0
      REAL , PARAMETER :: g           =     9.81
      REAL , PARAMETER :: r           =   287.04
      REAL , PARAMETER :: rovcp       =  r / cp

      REAL , PARAMETER :: gamma       =     6.5e-3

      REAL , PARAMETER :: tc          =   290.65  ! 273.15 + 17.5

      REAL , PARAMETER :: eps         =     0.622
      REAL , PARAMETER :: xlvc1       =     3.1484e6
      REAL , PARAMETER :: xlvc2       =     2.37e3

      REAL , PARAMETER :: svp1        =     0.6112
      REAL , PARAMETER :: svp2        =    17.67
      REAL , PARAMETER :: svp3        =    29.65
      REAL , PARAMETER :: svpt0       =   273.15

END MODULE phys_consts
