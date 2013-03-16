narccap
=======

These codes help with postprocessing our NARCCAP MM5 Runs, they may or may not
be useful to you and are provided without warranty!

	i_dot = 130 ;
	j_dot = 155 ;
	i_cross = 129 ;
	j_cross = 154 ;
	pressure = 26 ;

Here's where the local processing is being done

CCSM
----
    howard  -- /mnt/howard/narccap/exp2a.1/output/
    osgood  -- /mnt/osgood/narccap/exp2a.2/output

HADCM3
------
    thumper -- /tera10/akrherz/narccap/hadcm3/Run.contemporary.deepsoil_off
    comet -- /mnt/comet/narccap/hadcm3/Run.scenario.deepsoil_off
        
NCEP REANALYSIS
---------------
    stanley -- /stanley/narccap/ncep/Run.NCEP

landtyp Land-Cover Type	-- MMOUTP has land_use(i_cross, j_cross)
mrsofc	Capacity of Soil to Store Water	kg m-2	 
orog    Surface Altitude -- MMOUTP has terrain(i_cross, j_cross)
rootd	Root Depth	m	 
sftlf	Land Area Fraction	1 -- REGRID has LANDMASK
 