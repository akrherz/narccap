"""
 Common stuff for the NARCCAP POST scripts
"""
import numpy
HOURLY3, DAILY = (1,2)

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

} # end of VARS