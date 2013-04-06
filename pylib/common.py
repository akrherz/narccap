"""
 Common stuff for the NARCCAP POST scripts
"""
import numpy
HOURLY3, DAILY = (1,2)

VARS = {
 'swe'  : {'units'  : 'mm',
           'source' : 'MMOUTP',
          'ncsource'  : 'weasd',
          'table' : 3,
          'interval' : HOURLY3,
          'coordinates'  : "lon lat",
          'cell_methods'  : 'time: Instantaneous',
          'long_name'  : 'Snow Water Equivalent',
          'standard_name'  : 'lwe_thickness_of_surface_snow_amount'},
 'prw'  : {'units'  : 'kg m-2',
            'source' : 'MMOUTP',
          'table' : 3,
          'interval' : HOURLY3,
            'coordinates'  : "lon lat",
            'long_name'  : 'Precipitable Water',
            'standard_name'  : 'atmosphere_water_vapor_content'},
 'zmla'  : {'units'  : 'm',
            'source' : 'MMOUTP',
            'ncsource'  : 'pbl_hgt',
          'table' : 3,
          'interval' : HOURLY3,
            'coordinates'  : "lon lat",
            'long_name'  : 'Atmospheric Boundary Layer Thickness',
            'standard_name'  : 'atmosphere_boundary_layer_thickness'},
 'tauv'  : {'units'  : 'Pa',
            'source' : 'NCOUT',
            'ncsource'  : 'vmomflux',
          'table' : 3,
          'interval' : HOURLY3,
            'coordinates'  : "lon lat",
            'positive': 'down',
            'cell_methods'  : 'time: mean (interval: 3 hours)',
            'long_name'  : 'Surface Downward Flux of Northward Momentum',
            'standard_name'  : 'surface_downward_northward_stress'},
 'tauu'  : {'units'  : 'Pa',
            'source' : 'NCOUT',
            'ncsource'  : 'umomflux',
          'table' : 3,
          'interval' : HOURLY3,
            'coordinates'  : "lon lat",
            'positive': 'down',
            'cell_methods'  : 'time: mean (interval: 3 hours)',
            'long_name'  : 'Surface Downward Flux of Eastern Momentum',
            'standard_name'  : 'surface_downward_eastward_stress'},
 'snm'  : {'units'  : 'kg m-2 s-1',
            'source' : 'NCOUT',
            'ncsource'  : 'snowmelt',
          'table' : 3,
          'interval' : HOURLY3,
            'coordinates'  : "lon lat",
            'cell_methods'  : 'time: mean (interval: 3 hours)',
            'long_name'  : 'Snow Melt',
            'standard_name'  : 'surface_snow_melt_flux'},
 'snd'  : {'units'  : 'm',
            'source' : 'MMOUTP',
            'ncsource'  : 'snowh',
          'table' : 3,
          'interval' : HOURLY3,
            'coordinates'  : "lon lat",
            'long_name'  : 'Snow Depth',
            'standard_name'  : 'surface_snow_thickness'},
 'rsut'  : {'units'  : 'W m-2',
            'source' : 'NCOUT',
            'ncsource'  : 'osw1',
             'positive': 'up',
          'table' : 3,
          'interval' : HOURLY3,
            'coordinates'  : "lon lat",
            'cell_methods'  : 'time: mean (interval: 3 hours)',
            'long_name'  : 'TOA Reflected Shortwave Radiation',
            'standard_name'  : 'toa_outgoing_shortwave_flux'},
 'rsus'  : {'units'  : 'W m-2',
            'source' : 'NCOUT',
            'ncsource'  : 'osw1', 
             'positive': 'up',
          'table' : 3,
          'interval' : HOURLY3,
            'coordinates'  : "lon lat",
            'cell_methods'  : 'time: mean (interval: 3 hours)',
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
            'cell_methods'  : 'time: mean (interval: 3 hours)',
            'long_name'  : 'TOA Incident Shortwave Radiation',
            'standard_name'  : 'toa_incoming_shortwave_flux'},
 'rlut'  : {'units'  : 'W m-2',
            'source' : 'NCOUT',
            'ncsource'  : 'olw1',
             'positive': 'up',
          'table' : 3,
          'interval' : HOURLY3,
            'coordinates'  : "lon lat",
            'cell_methods'  : 'time: mean (interval: 3 hours)',
            'long_name'  : 'Outgoing Longwave Radiation',
            'standard_name'  : 'toa_outgoing_longwave_flux'},
 'rlus'  : {'units'  : 'W m-2',
            'source' : 'NCOUT',
            'ncsource'  : 'ulw1',
             'positive': 'up',
          'table' : 3,
          'interval' : HOURLY3,
            'coordinates'  : "lon lat",
            'cell_methods'  : 'time: mean (interval: 3 hours)',
            'long_name'  : 'Surface Upwelling Longwave Radiation',
            'standard_name'  : 'surface_upwelling_longwave_flux_in_air'},
 'rlds'  : {'units'  : 'W m-2',
            'source' : 'NCOUT',
            'ncsource'  : 'glw1',
             'positive': 'down',
          'table' : 3,
          'interval' : HOURLY3,
            'coordinates'  : "lon lat",
            'cell_methods'  : 'time: mean (interval: 3 hours)',
            'long_name'  : 'Surface Downwelling Longwave Radiation',
            'standard_name'  : 'surface_downwelling_longwave_flux_in_air'},
 'psl' : {'source' : 'MMOUTP',
          'ncsource'  : 'psealvlc',
          'units'  : 'Pa',
          'table' : 3,
          'interval' : HOURLY3,
          'coordinates' : "lon lat",
          'long_name'  : 'Sea Level Pressure',
          'standard_name'  : 'air_pressure_at_sea_level'},
 'mrso'  : {'units'  : 'kg m-2',
                'table' : 3,
            'source' : 'MMOUTP',
            'ncsource'  : None,
            'coordinates'  : "lon lat",
          'interval' : HOURLY3,
            'long_name'  : 'Total Soil Moisture Content',
            'standard_name'  : 'soil_moisture_content'},
 'mrros'  : {'units'  : 'kg m-2 s-1',
        'table' : 3,
            'source' : 'MMOUTP',
            'ncsource'  : None,
            'quo': 10800,
            'coordinates'  : "lon lat",
          'interval' : HOURLY3,
            'cell_methods'  : 'time: mean (interval: 3 hours)',
            'long_name'  : 'Surface Runoff',
            'standard_name'  : 'surface_runoff_flux'},
 'mrro'  : {'units'  : 'kg m-2 s-1',
            'source' : 'MMOUTP',
            'ncsource'  : None,
            'quo': 10800,
            'table': 3,
            'coordinates'  : "lon lat",
          'interval' : HOURLY3,
            'cell_methods'  : 'time: mean (interval: 3 hours)',
            'long_name'  : 'Surface and Subsurface Runoff',
            'standard_name'  : 'runoff_flux'},
 'mrfso'  : {'units'  : 'kg m-2',
            'source' : 'MMOUTP',
            'table': 3,
            'ncsource'  : None,
            'coordinates'  : "lon lat",
          'interval' : HOURLY3,
            'long_name'  : 'Soil Frozen Water Content',
            'standard_name'  : 'soil_frozen_water_content'},
 'evps'  : {'units'  : 'kg m-2 s-1',
            'source' : 'NCOUT',
            'table': 3,
            'ncsource'  : 'surfevap',
            'coordinates'  : "lon lat",
          'interval' : HOURLY3,
            'cell_methods'  : 'time: mean (interval: 3 hours)',
            'long_name'  : 'Surface Evaporation of Condensed Water',
            'standard_name'  : 'water_evaporation_flux'},
 'clt'  : {'units'  : '1',
            'source' : 'NCOUT',
            'ncsource'  : 'totcldavG',
            'coordinates'  : "lon lat",
            'table': 3,
          'interval' : HOURLY3,
            'cell_methods'  : 'time: mean (interval: 3 hours)',
            'long_name'  : 'Total Cloud Fraction',
            'standard_name'  : 'cloud_area_fraction'},
 'hfss'  : {'units'  : 'W m-2',
            'source' : 'NCOUT',
            'ncsource'  : 'hfx1',
            'coordinates'  : "lon lat",
            'table': 3,
             'positive': 'up',
          'interval' : HOURLY3,
            'cell_methods'  : 'time: mean (interval: 3 hours)',
            'long_name'  : 'Surface Sensible Heat Flux',
            'standard_name'  : 'surface_upward_sensible_heat_flux'},
 'hfls'  : {'units'  : 'W m-2',
            'source' : 'NCOUT',
            'ncsource'  : 'qfx1',
            'table': 3,
            'coordinates'  : "lon lat",
             'positive': 'up',
          'interval' : HOURLY3,
            'cell_methods'  : 'time: mean (interval: 3 hours)',
            'long_name'  : 'Surface Latent Heat Flux',
            'standard_name'  : 'surface_upward_latent_heat_flux'},
 'vas' : {'source' : 'MMOUTP',
          'ncsource'  : 'v10',
          'units'  : 'm s-1',
            'table': 2,
          'interval' : HOURLY3,
          'coordinates'  : 'lon lat level',
          'long_name'  : 'Meridional Surface Wind Speed',
          'standard_name'  : 'northward_wind'},
 'uas' : {'source' : 'MMOUTP',
          'ncsource'  : 'u10',
            'table': 2,
          'units'  : 'm s-1',
          'interval' : HOURLY3,
          'coordinates'  : 'lon lat level',
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
   'coordinates': 'lon lat',
 },
 'rsds' : {
   'ncsource': 'gsw1',
   'interval' : HOURLY3,
   'source' : 'NCOUT',  
   'table': 2,
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
    'table': 2,
   'units' : 'kg kg-1',
   'standard_name' : 'specific_humidity',
   'coordinates': 'lon lat level',
 },
 'ta' : {
   '3d': True,
    'ncsource' : 't',
   'interval' : HOURLY3,
   'source' : 'MMOUTP',
   'long_name' : 'Temperature',
    'table': 5,
   'units' : 'K',
   'standard_name' : 'air_temperature',
   'coordinates': 'lon lat level',
 },
 'ua' : {
   '3d': True,
    'ncsource' : 'u',
   'interval' : HOURLY3,
   'source' : 'MMOUTP',
   'long_name' : 'Zonal Wind Component',
            'table': 5,
   'units' : 'm s-1',
   'standard_name' : 'eastward_wind',
   'coordinates': 'lon lat level',
 },
 'wa' : {
   '3d': True,
    'ncsource' : 'w',
   'interval' : HOURLY3,
   'source' : 'MMOUTP',
   'long_name' : 'Vertical Wind Component',
            'table': 5,
   'units' : 'm s-1',
   'standard_name' : 'upward_air_velocity',
   'coordinates': 'lon lat level',
 },
 'va' : {                                                                       
   '3d': True,                                                                  
    'ncsource' : 'v',                                                           
   'interval' : HOURLY3,                                                        
   'source' : 'MMOUTP',                                                         
   'long_name' : 'Meridional Wind Component',                                   
   'units' : 'm s-1', 
   'table': 5,                                 
   'standard_name' : 'northward_wind',   
   'coordinates': 'lon lat level',                                       
 },  
 'zg500' : {
   '3d': True,
    'ncsource' : 'h',
   'table': 3,
   'interval' : HOURLY3,
   'source' : 'MMOUTP',
   'long_name' : '500 hPa Geopotential Height',
   'units' : 'm',
   'standard_name' : 'geopotential_height',
   'coordinates': 'lon lat level',
 },
 'clw' : {
   '3d': True,
    'ncsource' : 'clw',
   'interval' : HOURLY3,
            'table': 5,
   'source' : 'MMOUTP',
   'long_name' : 'Cloud Liquid Water Fraction of Layer',
   'units' : '1',
   'standard_name' : 'mass_fraction_of_cloud_liquid_water_in_air',
   'coordinates': 'lon lat level',
 },
 'cli' : {
   '3d': True,
    'ncsource' : 'ice',
   'interval' : HOURLY3,
   'source' : 'MMOUTP',
            'table': 5,
   'long_name' : 'Cloud Ice Fraction of Layer',
   'units' : '1',
   'standard_name' : 'mass_fraction_of_cloud_ice_in_air',
   'coordinates': 'lon lat level',
 },
 'hus' : {
   '3d': True,
    'ncsource' : 'q',
   'interval' : HOURLY3,
   'source' : 'MMOUTP',
   'long_name' : 'Specific Humidity',
   'units' : 'kg kg-1',
            'table': 5,
   'standard_name' : 'specific_humidity',
   'coordinates': 'lon lat level',
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
            'table': 2,
   'quo' : 1080.0, # cm to mm to /s
   'standard_name' : 'precipitation_flux',
   'cell_methods' : 'time: mean (interval: 3 hours)',
   'coordinates': 'lon lat',
 },
# NOT REQUIRED FOR NARCCAP
 'soilt4' : {
    'ncsource': 'soil_t_4',
   'interval' : HOURLY3,
   'source' : 'MMOUTP',
   'long_name' : 'Soil Temperature at Layer 4',
   'units' : 'K',
   'standard_name' : 'unknown',
   'coordinates': 'lon lat level',
 },
 'soilt1' : {
    'ncsource': 'soil_t_1',
   'interval' : HOURLY3,
   'source' : 'MMOUTP',
   'long_name' : 'Soil Temperature at Layer 1',
   'units' : 'K',
   'standard_name' : 'unknown',
   'coordinates': 'lon lat level',
 },
 'soilm3' : {
    'ncsource': 'soil_m_3',
   'interval' : HOURLY3,
   'source' : 'MMOUTP',
   'long_name' : 'Soil Moisture at Layer 3',
   'units' : 'm3 m-3',
   'standard_name' : 'unknown',
   'coordinates': 'lon lat level',
 },
 'tas' : {
    'ncsource': 't2',
   'interval' : HOURLY3,
   'source' : 'MMOUTP',
            'table': 2,
   'long_name' : 'Surface Air Temperature',
   'units' : 'K',
   'standard_name' : 'air_temperature',
   'coordinates': 'lon lat level',
 },
 'ps' : {
    'ncsource': 'psfc',
   'interval' : HOURLY3,
   'source' : 'MMOUTP',
            'table': 2,
   'long_name' : 'Surface Pressure',
   'units' : 'Pa',
   'standard_name' : 'surface_air_pressure',
   'coordinates': 'lon lat',
 },
 'tasmax' : {
    'ncsource': 'tmax',
   'interval' : DAILY,
   'source' : 'NCOUT',
   'long_name' : 'Maximum Daily Surface Air Temperature',
   'units' : 'K',
            'table': 1,
   'standard_name' : 'air_temperature',
   'cell_methods' : 'time: maximum (interval: 24 hours)',
   'npfunc': numpy.max,
   'npfunc2': numpy.maximum,
   'coordinates': 'lon lat level',
 },
 'tasmin' : {
    'ncsource': 'tmin',
   'interval' : DAILY,
   'source' : 'NCOUT',
   'long_name' : 'Minimum Daily Surface Air Temperature',
            'table': 1,
   'units' : 'K',
   'standard_name' : 'air_temperature',
   'cell_methods' : 'time: minimum (interval: 24 hours)',
   'npfunc': numpy.min,
   'npfunc2': numpy.minimum,
   'coordinates': 'lon lat level',
 },
 'spdmax' : {
    'ncsource': 'maxwnd',
   'interval' : DAILY,
   'source' : 'NCOUT',
   'long_name' : 'Maximum Daily 10-Meter Wind Speed',
            'table': 1,
   'units' : 'm s-1',
   'standard_name' : 'wind_speed_of_gust',
   'cell_methods' : 'time: maximum (interval: 24 hours)',
   'npfunc': numpy.max,
   'npfunc2': numpy.maximum,
   'coordinates': 'lon lat level',
 },
 'sic' : {
    'ncsource': 'seaice',
   'interval' : DAILY,
   'source' : 'MMOUTP',
   'long_name' : 'Daily Average Sea-ice Fraction',
            'table': 1,
   'units' : '1',
   'standard_name' : 'sea_ice_area_fraction',
   'cell_methods' : 'time: mean (interval: 24 hours)',
   'npfunc': numpy.average,
   'coordinates': 'lon lat',
 }

} # end of VARS

def do_singleton(nc, VNAME, PLEVEL):
    # Generate singletons
    if VNAME in ['huss', 'tas', 'tasmax', 'tasmin', 'uas', 'vas', 'spdmax', 
                 'zg500']:
        level = nc.createVariable('level', 'd')
        level.long_name = "height"
        level.standard_name = "height"
        if VNAME == 'zg500':
            level.units = 'hPa'
        else:
            level.units = 'm'
        level.positive = "up"
        level.axis = "Z"
        if VNAME in ['huss','tas','tasmax','tasmin']:
            level[:] = 2.
        elif VNAME in ['uas','vas', 'spdmax']:
            level[:] = 10.
        elif VNAME in ['zg500',]:
            level[:] = 500.

    if PLEVEL is not None:
        level = nc.createVariable('level', 'd')
        level.long_name = "pressure"
        level.standard_name = "air_pressure"
        level.units = "Pa"
        level.positive = "down"
        level.axis = "Z"
        level[:] = PLEVEL