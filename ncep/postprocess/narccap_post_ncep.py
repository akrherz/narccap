"""
Uber script to generate NARCCAP archive specification files
This consumes the resulf of running the preprocess.py script
"""
import netCDF4
import mx.DateTime
import numpy
import sys
import os

if len(sys.argv) < 3:
    print 'Usage: python narccap_post.py RUNID VNAME'
    print '  RUNID is S: scenario  C: contempory'
    print '  VNAME is the standard NARCCAP Variable Name'
    print 'Example: python narccap_post.py C pr'
    sys.exit()

RUNID = sys.argv[1]
VNAME = sys.argv[2]
if len(sys.argv) == 4:
    PLEVEL = int(sys.argv[3])
else:
    PLEVEL = None
META  = {}
if RUNID == 'C':
    DATADIR = "../Run.NCEP"
    META['title'] = 'ISU MM5 model output prepared for NARCCAP present-day climate using NCEP/DOE Reanalysis'
    META['prefix'] = 'MM5I'
    META['experiment_id'] = 'present-day climate using NCEP/DOE Reanalysis'
    TIMES = [
    mx.DateTime.DateTime(1979,1,1),
    mx.DateTime.DateTime(1981,1,1),
    mx.DateTime.DateTime(1986,1,1),
    mx.DateTime.DateTime(1991,1,1),
    mx.DateTime.DateTime(1996,1,1),
    mx.DateTime.DateTime(2001,1,1),
    mx.DateTime.DateTime(2004,11,30)
    ]

HOURLY3, DAILY = (1,2)

LOOP1 = None
LOOP2 = None

VARS = {
 'prw'  : {'units'  : 'kg m-2',
            'source' : 'MMOUTP',
          'table' : 3,
          'interval' : HOURLY3,
            'coordinates'  : "lon lat",
          'cell_methods'  : 'time: Instantaneous',
            'long_name'  : 'Precipitable Water',
            'standard_name'  : 'precipitable_water'},
 'zmla'  : {'units'  : 'm',
            'source' : 'MMOUTP',
            'ncsource'  : 'pbl_hgt',
          'table' : 3,
          'interval' : HOURLY3,
            'coordinates'  : "lon lat",
          'cell_methods'  : 'time: Instantaneous',
            'long_name'  : 'Atmospheric Boundary Layer Thickness',
            'standard_name'  : 'atmosphere_boundary_layer_thickness'},
 'tauv'  : {'units'  : 'Pa',
            'source' : 'NCOUT',
            'ncsource'  : 'vmomflux',
          'table' : 3,
          'interval' : HOURLY3,
            'coordinates'  : "lon lat",
            'positive': 'down',
            'cell_methods'  : 'time: average (interval: 3 hours)',
            'long_name'  : 'Surface Downward Flux of Northward Momentum',
            'standard_name'  : 'surface_downward_northward_flux'},
 'tauu'  : {'units'  : 'Pa',
            'source' : 'NCOUT',
            'ncsource'  : 'umomflux',
          'table' : 3,
          'interval' : HOURLY3,
            'coordinates'  : "lon lat",
            'positive': 'down',
            'cell_methods'  : 'time: average (interval: 3 hours)',
            'long_name'  : 'Surface Downward Flux of Eastern Momentum',
            'standard_name'  : 'surface_downward_eastward_flux'},
 'snm'  : {'units'  : 'kg m-2 s-1',
            'source' : 'NCOUT',
            'ncsource'  : 'snowmelt',
          'table' : 3,
          'interval' : HOURLY3,
            'coordinates'  : "lon lat",
            'cell_methods'  : 'time: average (interval: 3 hours)',
            'long_name'  : 'Snow Melt',
            'standard_name'  : 'surface_snow_melt_flux'},
 'snd'  : {'units'  : 'm',
            'source' : 'MMOUTP',
            'ncsource'  : 'snowh',
          'table' : 3,
          'interval' : HOURLY3,
            'coordinates'  : "lon lat",
          'cell_methods'  : 'time: Instantaneous',
            'long_name'  : 'Snow Depth',
            'standard_name'  : 'surface_snow_thickness'},
 'swe'  : {'units'  : 'mm',
            'source' : 'MMOUTP',
            'ncsource'  : 'weasd',
          'table' : 3,
          'interval' : HOURLY3,
            'coordinates'  : "lon lat",
          'cell_methods'  : 'time: Instantaneous',
            'long_name'  : 'Snow Water Equivalent',
            'standard_name'  : 'lwe_thickness_of_surface_snow_amount'},
 'rsut'  : {'units'  : 'W m-2',
            'source' : 'NCOUT',
            'ncsource'  : 'osw1',
             'positive': 'up',
          'table' : 3,
          'interval' : HOURLY3,
            'coordinates'  : "lon lat",
            'cell_methods'  : 'time: average (interval: 3 hours)',
            'long_name'  : 'TOA Reflected Shortwave Radiation',
            'standard_name'  : 'toa_outgoing_shortwave_flux'},
 'rsus'  : {'units'  : 'W m-2',
            'source' : 'NCOUT',
            'ncsource'  : 'osw1', 
             'positive': 'up',
          'table' : 3,
          'interval' : HOURLY3,
            'coordinates'  : "lon lat",
            'cell_methods'  : 'time: average (interval: 3 hours)',
            'long_name'  : 'Surface Upwelling Shortwave Radiation',
            'standard_name'  : 'surface_upwelling_shortwave_flux_in_air'},
 'rsdt'  : {'units'  : 'W m-2',
            'source' : 'NCOUT',
            'ncsource'  : 'toasw', 
            'quo' : 1.0/18.0,
             'positive': 'down',
          'table' : 3,
          'interval' : HOURLY3,
            'coordinates'  : "lon lat",
            'cell_methods'  : 'time: average (interval: 3 hours)',
            'long_name'  : 'TOA Incident Shortwave Radiation',
            'standard_name'  : 'toa_incoming_shortwave_flux'},
 'rlut'  : {'units'  : 'W m-2',
            'source' : 'NCOUT',
            'ncsource'  : 'olw1',
             'positive': 'up',
          'table' : 3,
          'interval' : HOURLY3,
            'coordinates'  : "lon lat",
            'cell_methods'  : 'time: average (interval: 3 hours)',
            'long_name'  : 'Outgoing Longwave Radiation',
            'standard_name'  : 'toa_outgoing_longwave_flux'},
 'rlus'  : {'units'  : 'W m-2',
            'source' : 'NCOUT',
            'ncsource'  : 'ulw1',
             'positive': 'up',
          'table' : 3,
          'interval' : HOURLY3,
            'coordinates'  : "lon lat",
            'cell_methods'  : 'time: average (interval: 3 hours)',
            'long_name'  : 'Surface Upwelling Longwave Radiation',
            'standard_name'  : 'surface_upwelling_longwave_flux_in_air'},
 'rlds'  : {'units'  : 'W m-2',
            'source' : 'NCOUT',
            'ncsource'  : 'glw1',
             'positive': 'down',
          'table' : 3,
          'interval' : HOURLY3,
            'coordinates'  : "lon lat",
            'cell_methods'  : 'time: average (interval: 3 hours)',
            'long_name'  : 'Surface Downwelling Longwave Radiation',
            'standard_name'  : 'surface_downwelling_longwave_flux_in_air'},
 'psl' : {'source' : 'MMOUTP',
          'ncsource'  : 'psealvlc',
          'units'  : 'Pa',
          'table' : 3,
          'interval' : HOURLY3,
          'coordinates' : "lon lat",
          'cell_methods'  : 'time: Instantaneous',
          'long_name'  : 'Sea Level Pressure',
          'standard_name'  : 'air_pressure_at_sea_level'},
 'mrso'  : {'units'  : 'kg m-2',
                'table' : 3,
            'source' : 'MMOUTP',
            'ncsource'  : None,
            'coordinates'  : "lon lat",
          'interval' : HOURLY3,
          'cell_methods'  : 'time: Instantaneous',
            'long_name'  : 'Total Soil Moisture Content',
            'standard_name'  : 'soil_moisture_content'},
 'mrros'  : {'units'  : 'kg m-2 s-1',
		'table' : 3,
            'source' : 'MMOUTP',
            'ncsource'  : None,
            'quo': 10800,
            'coordinates'  : "lon lat",
          'interval' : HOURLY3,
            'cell_methods'  : 'time: average (interval: 3 hours)',
            'long_name'  : 'Surface Runoff',
            'standard_name'  : 'surface_runoff_flux'},
 'mrro'  : {'units'  : 'kg m-2 s-1',
            'source' : 'MMOUTP',
            'ncsource'  : None,
            'quo': 10800,
            'coordinates'  : "lon lat",
          'interval' : HOURLY3,
            'cell_methods'  : 'time: average (interval: 3 hours)',
            'long_name'  : 'Surface and Subsurface Runoff',
            'standard_name'  : 'runoff_flux'},
 'mrfso'  : {'units'  : 'kg m-2',
            'source' : 'MMOUTP',
            'ncsource'  : None,
            'coordinates'  : "lon lat",
          'interval' : HOURLY3,
          'cell_methods'  : 'time: Instantaneous',
            'long_name'  : 'Soil Frozen Water Content',
            'standard_name'  : 'soil_frozen_water_content'},
 'evps'  : {'units'  : 'kg m-2 s-1',
            'source' : 'NCOUT',
            'ncsource'  : 'surfevap',
            'coordinates'  : "lon lat",
          'interval' : HOURLY3,
            'cell_methods'  : 'time: average (interval: 3 hours)',
            'long_name'  : 'Surface Evaporation of Condensed Water',
            'standard_name'  : 'water_evaporation_flux'},
 'clt'  : {'units'  : '1',
            'source' : 'NCOUT',
            'ncsource'  : 'totcldavG',
            'coordinates'  : "lon lat",
          'interval' : HOURLY3,
            'cell_methods'  : 'time: average (interval: 3 hours)',
            'long_name'  : 'Total Cloud Fraction',
            'standard_name'  : 'cloud_area_fraction'},
 'hfss'  : {'units'  : 'W m-2',
            'source' : 'NCOUT',
            'ncsource'  : 'hfx1',
            'coordinates'  : "lon lat",
             'positive': 'up',
          'interval' : HOURLY3,
            'cell_methods'  : 'time: average (interval: 3 hours)',
            'long_name'  : 'Surface Sensible Heat Flux',
            'standard_name'  : 'surface_upward_sensible_heat_flux'},
 'hfls'  : {'units'  : 'W m-2',
            'source' : 'NCOUT',
            'ncsource'  : 'qfx1',
            'coordinates'  : "lon lat",
             'positive': 'up',
          'interval' : HOURLY3,
            'cell_methods'  : 'time: average (interval: 3 hours)',
            'long_name'  : 'Surface Latent Heat Flux',
            'standard_name'  : 'surface_upward_latent_heat_flux'},
 'vas' : {'source' : 'MMOUTP',
          'ncsource'  : 'v10',
          'units'  : 'm s-1',
          'interval' : HOURLY3,
          'coordinates'  : 'lon lat height',
          'cell_methods'  : 'time: Instantaneous',
          'long_name'  : 'Meridional Surface Wind Speed',
          'standard_name'  : 'northward_wind'},
 'uas' : {'source' : 'MMOUTP',
          'ncsource'  : 'u10',
          'units'  : 'm s-1',
          'interval' : HOURLY3,
          'coordinates'  : 'lon lat height',
          'cell_methods'  : 'time: Instantaneous',
          'long_name'  : 'Zonal Surface Wind Speed',
          'standard_name'  : 'eastward_wind'},
 'ts' : {
   'ncsource': 'ground_t', # was tseasfc
   'interval' : HOURLY3,
   'source' : 'MMOUTP',
   'table': 3,
   'long_name' : 'Surface (skin) Temperature',
   'units' : 'K',
   'standard_name' : 'surface_temperature',
   'cell_methods' : 'time: instantaneous',
   'coordinates': 'lon lat',
 },
 'rsds' : {
   'ncsource': 'gsw1',
   'interval' : HOURLY3,
   'source' : 'NCOUT',
   'long_name' : 'Surface Downwelling Shortwave Radiation',
   'units' : 'W m-2',
   'standard_name' : 'surface_downwelling_shortwave_flux_in_air',
   'cell_methods' : 'time: mean (interval: 3 hours)',
   'coordinates': 'lon lat',
 },
 'huss' : {
           'ncsource' : 'q2',
   'interval' : HOURLY3,
   'source' : 'MMOUTP',
   'long_name' : 'Surface Specific Humidity',
   'units' : 'kg kg-1',
   'standard_name' : 'specific_humidity',
   'cell_methods' : 'time: instantaneous',
   'coordinates': 'lon lat height',
 },
 'ta' : {
   '3d': True,
    'ncsource' : 't',
   'interval' : HOURLY3,
   'source' : 'MMOUTP',
   'long_name' : 'Temperature',
   'units' : 'K',
   'standard_name' : 'air_temperature',
   'cell_methods' : 'time: instantaneous',
   'coordinates': 'lon lat plev',
 },
 'ua' : {
   '3d': True,
    'ncsource' : 'u',
   'interval' : HOURLY3,
   'source' : 'MMOUTP',
   'long_name' : 'Zonal Wind Component',
   'units' : 'm s-1',
   'standard_name' : 'eastward_wind',
   'cell_methods' : 'time: instantaneous',
   'coordinates': 'lon lat plev',
 },
 'wa' : {
   '3d': True,
    'ncsource' : 'w',
   'interval' : HOURLY3,
   'source' : 'MMOUTP',
   'long_name' : 'Vertical Wind Component',
   'units' : 'm s-1',
   'standard_name' : 'upward_air_velocity',
   'cell_methods' : 'time: instantaneous',
   'coordinates': 'lon lat plev',
 },
 'va' : {                                                                       
   '3d': True,                                                                  
    'ncsource' : 'v',                                                           
   'interval' : HOURLY3,                                                        
   'source' : 'MMOUTP',                                                         
   'long_name' : 'Meridional Wind Component',                                   
   'units' : 'm s-1',                                                           
   'standard_name' : 'northward_wind',                                          
   'cell_methods' : 'time: instantaneous',                                         'coordinates': 'lon lat plev',                                               
 },  
 'zg' : {
   '3d': True,
    'ncsource' : 'h',
   'table': 3,
   'interval' : HOURLY3,
   'source' : 'MMOUTP',
   'long_name' : '500 hPa Geopotential Height',
   'units' : 'm',
   'standard_name' : 'geopotential_height',
   'cell_methods' : 'time: instantaneous',
   'coordinates': 'lon lat plev',
 },
 'clw' : {
   '3d': True,
    'ncsource' : 'clw',
   'interval' : HOURLY3,
   'source' : 'MMOUTP',
   'long_name' : 'Cloud Liquid Water Fraction of Layer',
   'units' : '1',
   'standard_name' : 'mass_fraction_of_cloud_liquid_water_in_air',
   'cell_methods' : 'time: instantaneous',
   'coordinates': 'lon lat plev',
 },
 'cli' : {
   '3d': True,
    'ncsource' : 'ice',
   'interval' : HOURLY3,
   'source' : 'MMOUTP',
   'long_name' : 'Cloud Ice Fraction of Layer',
   'units' : '1',
   'standard_name' : 'mass_fraction_of_cloud_ice_in_air',
   'cell_methods' : 'time: instantaneous',
   'coordinates': 'lon lat plev',
 },
 'hus' : {
   '3d': True,
    'ncsource' : 'q',
   'interval' : HOURLY3,
   'source' : 'MMOUTP',
   'long_name' : 'Specific Humidity',
   'units' : 'kg kg-1',
   'standard_name' : 'specific_humidity',
   'cell_methods' : 'time: instantaneous',
   'coordinates': 'lon lat plev',
 },
 'prc' : {
   'interval' : HOURLY3,
   'source' : 'MMOUTP',
   'ncsource' : 'rain_con',
   'table': 3,
   'quo' : 1080.0, # cm to mm to /s
   'long_name' : 'Convective Precipitation',
   'units' : 'kg m-2 s-1',
   'standard_name' : 'convective_precipitation_flux',
   'cell_methods' : 'time: mean (interval: 3 hours)',
   'coordinates': 'lon lat',
 },
 'pr' : {
   'interval' : HOURLY3,
   'source' : 'MMOUTP',
   'long_name' : 'Precipitation',
   'units' : 'kg m-2 s-1',
   'quo' : 1080.0, # cm to mm to /s
   'standard_name' : 'precipitation_flux',
   'cell_methods' : 'time: mean (interval: 3 hours)',
   'coordinates': 'lon lat',
 },
 'soilt4' : {
    'ncsource': 'soil_t_4',
   'interval' : HOURLY3,
   'source' : 'MMOUTP',
   'long_name' : 'Soil Temperature at Layer 4',
   'units' : 'K',
   'standard_name' : 'unknown',
   'cell_methods' : 'time: instantaneous',
   'coordinates': 'lon lat height',
 },
 'tas' : {
    'ncsource': 't2',
   'interval' : HOURLY3,
   'source' : 'MMOUTP',
   'long_name' : 'Surface Air Temperature',
   'units' : 'K',
   'standard_name' : 'air_temperature',
   'cell_methods' : 'time: instantaneous',
   'coordinates': 'lon lat height',
 },
 'ps' : {
    'ncsource': 'psfc',
   'interval' : HOURLY3,
   'source' : 'MMOUTP',
   'long_name' : 'Surface Pressure',
   'units' : 'Pa',
   'standard_name' : 'surface_air_pressure',
   'cell_methods' : 'time: instantaneous',
   'coordinates': 'lon lat',
 },
 'soilt1' : {
    'ncsource': 'soil_t_1',
   'interval' : HOURLY3,
   'source' : 'MMOUTP',
   'long_name' : 'Soil Temperature at Layer 1',
   'units' : 'K',
   'standard_name' : 'unknown',
   'cell_methods' : 'time: instantaneous',
   'coordinates': 'lon lat height',
 },
 'soilm3' : {
    'ncsource': 'soil_m_3',
   'interval' : HOURLY3,
   'source' : 'MMOUTP',
   'long_name' : 'Soil Moisture at Layer 3',
   'units' : 'm3 m-3',
   'standard_name' : 'unknown',
   'cell_methods' : 'time: instantaneous',
   'coordinates': 'lon lat height',
 },
 'tasmax' : {
    'ncsource': 'tmax',
   'interval' : DAILY,
   'source' : 'NCOUT',
   'long_name' : 'Maximum Daily Surface Air Temperature',
   'units' : 'K',
   'standard_name' : 'air_temperature',
   'cell_methods' : 'time: maximum (interval: 24 hours)',
   'npfunc': numpy.max,
   'npfunc2': numpy.maximum,
   'coordinates': 'lon lat height',
 },
 'tasmin' : {
    'ncsource': 'tmin',
   'interval' : DAILY,
   'source' : 'NCOUT',
   'long_name' : 'Minimum Daily Surface Air Temperature',
   'units' : 'K',
   'standard_name' : 'air_temperature',
   'cell_methods' : 'time: minimum (interval: 24 hours)',
   'npfunc': numpy.min,
   'npfunc2': numpy.minimum,
   'coordinates': 'lon lat height',
 },
 'spdmax' : {
    'ncsource': 'maxwnd',
   'interval' : DAILY,
   'source' : 'NCOUT',
   'long_name' : 'Maximum Daily 10-Meter Wind Speed',
   'units' : 'm s-1',
   'standard_name' : 'wind_speed_of_gust',
   'cell_methods' : 'time: maximum (interval: 24 hours)',
   'npfunc': numpy.max,
   'npfunc2': numpy.maximum,
   'coordinates': 'lon lat height',
 },
 'sic' : {
    'ncsource': 'seaice',
   'interval' : DAILY,
   'source' : 'MMOUTP',
   'long_name' : 'Daily Average Sea-ice Fraction',
   'units' : 'fraction in [0,1]',
   'standard_name' : 'sea_ice_fraction',
   'cell_methods' : 'time: instantaneous',
   'npfunc': numpy.average,
   'coordinates': 'lon lat',
 }
}



def create_file(VNAME, ts0, ts1):
    """
    Create a CF compliant file for NARCCAP
    """
    fp = '../final/%s_%s_ncep_%s03.nc' % (VNAME, META['prefix'], ts0.strftime("%Y%m%d"))
    if PLEVEL is not None:
        fp = '../final/%s_%s_ncep_p%03i_%s.nc' % (VNAME, META['prefix'],
                                        PLEVEL, ts0.strftime("%Y%m%d%H"))
    nc = netCDF4.Dataset(fp, 'w', format='NETCDF3_CLASSIC')
    nc.Conventions = 'CF-1.0'
    nc.title = META['title']
    nc.history = 'Subsetted by archiver and python netCDF routines, %s' % (mx.DateTime.now().strftime("%Y-%m-%d %H:%M"),)
    nc.contact1 = 'Daryl Herzmann (akrherz@iastate.edu)'
    nc.contact2 = '3015 Agronomy Hall, Iowa State Univ.,Ames, Iowa, USA'
    nc.realization = '1'
    nc.experiment_id = META['experiment_id']
    # Should be table 1 for daily ?
    nc.table_id = 'Table %s' % (VARS[VNAME].get('table', 1),)
    nc.project_id = 'NARCCAP'
    nc.source = 'MM5(2002): atmosphere: MM5v3.6.3 non-hydrostatic; sst/sea ice: AMIPII; land: Noah;  Convection: Kain-Fritsch 2; Radiation: RRTM; PBL: MRF; Explicit Moisture: Reisner Mixed-Phase; Buffer: 15 point exponential; Horizontal Resolution: 50km; Vertical Levels: 24'
    nc.institution = 'ISU (Iowa State University, Ames, Iowa, USA)'

    tsteps = int((ts1 - ts0).days)
    if VARS[VNAME]['interval'] == HOURLY3:
        tsteps *= 8
    print ' + Created NetCDF File %s has %s time steps' % (fp, tsteps)
    nc.createDimension('time', 0)
    nc.createDimension('bnds', 2)
    if VNAME in ['ua','va']:
        nc.createDimension('xc', 124)
        nc.createDimension('yc', 100)
        latgrid = 'latitdot'
        longrid = 'longidot'
        xgend = 139
        ygend = 115
    else:
        nc.createDimension('xc', 123)
        nc.createDimension('yc', 99)
        latgrid = 'latitcrs'
        longrid = 'longicrs'
        xgend = 138
        ygend = 114

    # Create Time Dimension
    tm = nc.createVariable('time', 'd', ('time',))
    tm.long_name = 'time'
    tm.standard_name = 'time'
    tm.axis = 'T'
    tm.calendar = 'gregorian'
    tm.units = 'days since %s 00:00:0.0' % (TIMES[0].strftime("%Y-%m-%d"),)
    tm.bounds = 'time_bnds'

    tmb = nc.createVariable('time_bnds', 'd', ('time', 'bnds'))

    xc = nc.createVariable('xc', 'd', ('xc',))
    xc.long_name = 'x-coordinate Cartesian system'
    xc.standard_name = 'projection_x_coordinate'
    xc.axis = 'X'
    xc.units = 'm'

    yc = nc.createVariable('yc', 'd', ('yc',))
    yc.long_name = 'y-coordinate Cartesian system'
    yc.standard_name = 'projection_y_coordinate'
    yc.axis = 'Y'
    yc.units = 'm'

    lat = nc.createVariable('lat', 'd', ('yc', 'xc'))
    lat.long_name = 'latitude'
    lat.standard_name = 'latitude'
    lat.axis = 'Y'
    lat.units = 'degrees_north'
    
    lon = nc.createVariable('lon', 'd', ('yc', 'xc'))
    lon.long_name = 'longitude'
    lon.standard_name = 'longitude'
    lon.axis = 'X'
    lon.units = 'degrees_east'

    p = nc.createVariable('Lambert_Conformal', 'c', ())
    p.grid_mapping_name = "lambert_conformal_conic"
    p.false_easting = 3825000.
    p.false_northing = 3187500.

    v = nc.createVariable(VNAME, 'f', ('time','yc','xc'), fill_value=1e20)
    v.units = VARS[VNAME]['units'] 
    v.standard_name = VARS[VNAME]['standard_name'] 
    v.long_name = VARS[VNAME]['long_name'] 
    v.cell_methods = VARS[VNAME]['cell_methods']
    v.missing_value = numpy.array(1e20, v.dtype)
    v.coordinates = VARS[VNAME]['coordinates']
    v.grid_mapping = 'Lambert_Conformal'
    if VARS[VNAME].has_key('positive'):
        v.positive = VARS[VNAME]['positive']


    # write tm
    offset = int((ts0 - TIMES[0]).days)
    if VARS[VNAME]['interval'] == HOURLY3:
        tm[:] = offset + numpy.arange(0.125, (tsteps/8) + 0.125, 0.125)
        # write tmb
        tmb[:,0] = offset + numpy.arange(0., (tsteps/8), 0.125)
        tmb[:,1] = offset + numpy.arange(0.125, (tsteps/8)+0.125, 0.125)
    else:
        tm[:] = offset + numpy.arange(0.25, tsteps , 1.)
        # write tmb
        tmb[:,0] = offset + numpy.arange(0.25, tsteps, 1.)
        tmb[:,1] = offset + numpy.arange(1.25, tsteps+1., 1.)
   
    nc2 = netCDF4.Dataset('%s/%s_DOMAIN1_0001.nc' % (
                                DATADIR, VARS[VNAME]['source']), 'r')
    nc3 = netCDF4.Dataset('%s/MMOUTP_DOMAIN1_0001.nc' % (
                                DATADIR), 'r')
    # write lat
    # write lat
    lat[:] = nc3.variables[latgrid][15:-15,15:-16]
    lon[:] = nc3.variables[longrid][15:-15,15:-16]
    nc3.close()
    xc[:] = numpy.arange(15,xgend) * nc2.variables['grid_ds'][:] * 1000.0
    yc[:] = numpy.arange(15,ygend) * nc2.variables['grid_ds'][:] * 1000.0
    p.standard_parallel = [nc2.variables['stdlat_2'][:], nc2.variables['stdlat_1'][:]]
    p.longitude_of_central_meridian = nc2.variables['coarse_cenlon'][:]
    p.latitude_of_projection_origin = nc2.variables['coarse_cenlat'][:]
    nc2.close()

    # Generate singletons
    if VNAME in ['huss','tas','tasmax', 'tasmin','uas','vas']:
        height = nc.createVariable('height', 'd')
        height.long_name = "height"
        height.standard_name = "height"
        height.units = "m"
        height.positive = "up"
        height.axis = "Z"
        if VNAME in ['huss','tas','tasmax','tasmin']:
            height[:] = 2.
        elif VNAME in ['uas','vas']:
            height[:] = 10.

    if PLEVEL is not None:
        plev = nc.createVariable('plev', 'd')
        plev.long_name = "pressure"
        plev.standard_name = "air_pressure"
        plev.units = "Pa"
        plev.positive = "down"
        plev.axis = "Z"
        plev[:] = PLEVEL
        

    nc.close()
    return fp

def compute1d(VNAME, fp, ts0, ts1):
    """
    6z to 6z is our local day, so we need 9z to 6z inclusive.  The file
    starts at 3z, so offset 2 to 10
    """
    nc = netCDF4.Dataset(fp, 'a')
    ncv = nc.variables[VNAME]

    lookfor = ts0.strftime("minutes since %Y-%m-%d")
    # Figure out when our data begins!
    for i in range(1,1225):
        fp2 = '%s/MMOUTP_DOMAIN1_%04i.nc' % (DATADIR, i)
        nc2 = netCDF4.Dataset(fp2, 'r')
        # minutes since 1983-11-01 03:00:16
        if nc2.variables['time'].units.find(lookfor) == 0:
           nc2.close()
           break
        nc2.close()

    print 'For timestamp %s We found file %s' % (ts0, fp2)
    now = ts0 + mx.DateTime.RelativeDateTime(hours=6)
    oneday = mx.DateTime.RelativeDateTime(days=1)
    cnter = 0
    while now < ts1:
        #if now.month == 2 and now.day == 29:
        #    now += oneday
        #    continue
        fp2 = '%s/%s_DOMAIN1_%04i.nc' % (DATADIR, VARS[VNAME]['source'], i,)
        nc2 = netCDF4.Dataset(fp2, 'r')
        fp3 = '%s/MMOUTP_DOMAIN1_%04i.nc' % (DATADIR,  i,)
        nc3 = netCDF4.Dataset(fp3, 'r')
        # Figure out timestamp base
        # Figure out timestamp base
        tsbase = mx.DateTime.strptime(nc3.variables['time'].units[14:27], 
		'%Y-%m-%d %H')
        nc3.close()
        tsteps = len(nc2.variables['time'][:])

        # Okay, we need to go from 6z to 6z of the next day
        offset1 = int((now - tsbase).hours / 3)
        offset2 = int(((now + oneday) - tsbase).hours / 3)
        if offset2 > tsteps:
            i += 1
            offset2 = tsteps
        data = VARS[VNAME]['npfunc'](nc2.variables[VARS[VNAME]['ncsource']][offset1:offset2,15:-15,15:-16], axis=0)
        nc2.close()


        if offset2 == tsteps: # Need to step ahead and get next file
            fp2 = '%s/%s_DOMAIN1_%04i.nc' % (DATADIR, VARS[VNAME]['source'], i,)
            nc2 = netCDF4.Dataset(fp2, 'r')
            data2 = nc2.variables[VARS[VNAME]['ncsource']][:1,15:-15,15:-16]
            nc2.close()
            #print 'data IN', numpy.average(data)
            #print 'data2 IN', numpy.average(data2)
            if numpy.max(data2) != 0 and VNAME not in ['sic',]:
                data = VARS[VNAME]['npfunc2'](data, data2)
            else:
                print 'Skipping TS2 Computation'
        print '%s %s %02i %02i Avg: %.3f' % (now.strftime("%Y%m%d"), fp2, offset1, offset2, numpy.average(data))

        ncv[cnter] = data.astype('f')
        cnter += 1
        now += mx.DateTime.RelativeDateTime(days=1)
    nc.close()


def compute3h(VNAME, fp, ts0, ts1):
    """
    This is just a straight dumping of data from the NC or MMOUTP files
    to the resulting netCDF file.
    """
    lookfor = ts0.strftime("minutes since %Y-%m-%d")                            
    # Figure out when our data begins!                                          
    for i in range(1,1500):                                                     
        fp2 = '%s/%s_DOMAIN1_%04i.nc' % (DATADIR, VARS[VNAME]['source'], i,)    
        nc2 = netCDF4.Dataset(fp2, 'r')                                         
        # minutes since 1983-11-01 03:00:16                                     
        if nc2.variables['time'].units.find(lookfor) == 0:                      
           nc2.close()                                                          
           break                                                                
        nc2.close()                                                             
                                                                                
    print 'For timestamp %s We found file %s' % (ts0, fp2)       

    nc = netCDF4.Dataset(fp, 'a')
    ncv = nc.variables[VNAME]
    total = len(nc.variables['time'][:])
    now = ts0
    v = 0
    while v < total:
        fp2 = '%s/%s_DOMAIN1_%04i.nc' % (DATADIR, VARS[VNAME]['source'], i,)
        nc2 = netCDF4.Dataset(fp2, 'r')
        tsbase = mx.DateTime.strptime(nc2.variables['time'].units[14:27], 
		'%Y-%m-%d %H')
        tsteps = len(nc2.variables['time'][:])

        if VNAME == 'mrfso':
            data = ((nc2.variables['soil_m_1'][:,15:-15,15:-16] - nc2.variables['soil_w_1'][:,15:-15,15:-16]) * 0.10  + (nc2.variables['soil_m_2'][:,15:-15,15:-16] - nc2.variables['soil_w_2'][:,15:-15,15:-16]) * 0.30  + (nc2.variables['soil_m_3'][:,15:-15,15:-16] - nc2.variables['soil_w_3'][:,15:-15,15:-16]) * 0.60  + (nc2.variables['soil_m_4'][:,15:-15,15:-16] - nc2.variables['soil_w_4'][:,15:-15,15:-16]) * 1.00  ) * 1000.0

        elif VNAME == 'mrso':
            data = ((nc2.variables['soil_m_1'][:,15:-15,15:-16]) * 0.10  + (nc2.variables['soil_m_2'][:,15:-15,15:-16] ) * 0.30  + (nc2.variables['soil_m_3'][:,15:-15,15:-16] ) * 0.60  + (nc2.variables['soil_m_4'][:,15:-15,15:-16] ) * 1.00  ) * 1000.0

        elif VNAME == 'prw': # Precipitable Water
            # http://docs.lib.noaa.gov/rescue/mwr/067/mwr-067-04-0100.pdf
            bogus = nc2.variables['soil_m_1'][:,15:-15,15:-16]
            data = numpy.zeros( bogus.shape )
            plevels = nc2.variables['pressure'][:]
            for i in range(1,len(plevels)-1):
                dp = plevels[i] - plevels[i+1]
                q1 = nc2.variables['q'][:,i,15:-15,15:-16]
                q2 = nc2.variables['q'][:,i+1,15:-15,15:-16]
                # inches, convert to mm , which is kg m-2
                data += .0002 * (dp) * ((q1 + q2) * 1000.0) * 25.4

        elif VNAME == 'mrro':
            surface = nc2.variables['sfcrnoff'][:,15:-15,15:-16]
            subsurface = nc2.variables['ugdrnoff'][:,15:-15,15:-16]
            data = numpy.zeros( surface.shape )
            tsteps = surface.shape[0]
            for i in range(tsteps):
                if LOOP1 is None:
                    s0 = surface[i]
                    ss0 = subsurface[i]
                else:
                    s0 = surface[i] - LOOP1
                    ss0 = subsurface[i] - LOOP2
                LOOP1 = surface[i]
                LOOP2 = subsurface[i]
                data[i] = s0 + ss0

        elif VNAME == 'mrros':
            surface = nc2.variables['sfcrnoff'][:,15:-15,15:-16]
            data = numpy.zeros( surface.shape )
            tsteps = surface.shape[0]
            for i in range(tsteps):
                if LOOP1 is None:
                    s0 = surface[i]
                else:
                    s0 = surface[i] - LOOP1
                LOOP1 = surface[i]
                data[i] = s0 

        elif VNAME == 'huss':
            q2 = nc2.variables['q2'][:,15:-15,15:-16]
            data = q2 / (1.0 + q2)
        elif VNAME == 'pr':
            data = nc2.variables['rain_con'][:,15:-15,15:-16] + nc2.variables['rain_non'][:,15:-15,15:-16]
        else:
            ncs = VARS[VNAME]['ncsource']
            if PLEVEL is not None:
                l = list(nc2.variables['pressure'][:]).index(PLEVEL)
                data = nc2.variables[ncs][:,l,15:-15,15:-16]
            else:
                data = nc2.variables[ncs][:,15:-15,15:-16]
        # Okay, we have the data var loaded up
        v2 = v + tsteps
        ed = tsteps
        if v2 > total:
            ed -= (v2 - total)
            v2 = total
        
        print "i=%4i tsteps=%2i %5i %5i/%5i %.5f" % (i, tsteps, v, v2, total, 
           numpy.max(data) / VARS[VNAME].get('quo', 1.0))
        ncv[v:v2] = data[0:ed] / VARS[VNAME].get('quo', 1.0)
        nc2.close()
        v = v2
        i += 1
    nc.close()

for i in range(len(TIMES)-1):
    ts0 = TIMES[i]
    ts1 = TIMES[i+1]
    fp = create_file(VNAME, ts0, ts1 )
    if VARS[VNAME]['interval'] == HOURLY3:
        compute3h(VNAME, fp, ts0, ts1)
    else:
        compute1d(VNAME, fp, ts0, ts1)
