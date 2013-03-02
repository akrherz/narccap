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
    print '  RUNID in S: scenario  C: contempory'
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
if RUNID == 'S':
    DATADIR = "data.FUTURE"
    META['title'] = 'ISU MM5 model output prepared for NARCCAP scenario climate using HADCM3'
    META['experiment_id'] = 'scenario climate using HADCM3'
    TIMES = [
    mx.DateTime.DateTime(2038,1,1),
    mx.DateTime.DateTime(2041,1,1),
    mx.DateTime.DateTime(2046,1,1),
    mx.DateTime.DateTime(2051,1,1),
    mx.DateTime.DateTime(2056,1,1),
    mx.DateTime.DateTime(2061,1,1),
    mx.DateTime.DateTime(2066,1,1),
    mx.DateTime.DateTime(2071,1,1)
    ]
if RUNID == 'C':
    DATADIR = "data"
    META['title'] = 'ISU MM5 model output prepared for NARCCAP contemporary climate using HADCM3'
    META['experiment_id'] = 'contemporary climate using HADCM3'
    TIMES = [
    mx.DateTime.DateTime(1968,1,1),
    mx.DateTime.DateTime(1971,1,1),
    mx.DateTime.DateTime(1976,1,1),
    mx.DateTime.DateTime(1981,1,1),
    mx.DateTime.DateTime(1986,1,1),
    mx.DateTime.DateTime(1991,1,1),
    mx.DateTime.DateTime(1996,1,1),
    mx.DateTime.DateTime(2001,1,1)
    ]

HOURLY3, DAILY = (1,2)

PLEVELS = [1050, 1000, 850, 700, 500, 300]

VARS = {
 'ts' : {
   'ncsource': 'tseasfc',
   'interval' : HOURLY3,
   'source' : 'MMOUTP',
   'long_name' : 'Surface (skin) Temperature',
   'units' : 'K',
   'standard_name' : 'surface_temperature',
   'cell_methods' : 'time: instantaneous',
   'coordinates': 'yc xc',
 },
 'rsds' : {
   'ncsource': 'gsw1',
   'interval' : HOURLY3,
   'source' : 'NCOUT',
   'long_name' : 'Surface Downwelling Shortwave Radiation',
   'units' : 'W m-2',
   'standard_name' : 'surface_downwelling_shortwave_flux_in_air',
   'cell_methods' : 'time: mean (interval: 3 hours)',
   'coordinates': 'yc xc',
 },
 'huss' : {
           'ncsource' : 'q2',
   'interval' : HOURLY3,
   'source' : 'MMOUTP',
   'long_name' : 'Surface Specific Humidity',
   'units' : 'kg kg-1',
   'standard_name' : 'specific_humidity',
   'cell_methods' : 'time: instantaneous',
   'coordinates': 'yc xc height',
 },
        # cli and clw are a bit concerning, MM5 has it as kg/kg, but not
        # certain if that is against dry air or cloud ice.  
'cli' : {
   '3d': True,
    'ncsource' : 'ice',
   'interval' : HOURLY3,
   'source' : 'MMOUTP',
   'long_name' : 'Cloud Ice Fraction of Layer',
   'units' : '1',
   'standard_name' : 'mass_fraction_of_cloud_ice_in_air',
   'cell_methods' : 'time: instantaneous',
   'coordinates': 'yc xc plev',
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
   'coordinates': 'yc xc plev',
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
   'coordinates': 'yc xc plev',
 },
 'va' : {
   '3d': True,
    'ncsource' : 'v',
   'interval' : HOURLY3,
   'source' : 'MMOUTP',
   'long_name' : 'Meridional Wind Component',
   'units' : 'm s-1',
   'standard_name' : 'northward_wind',
   'cell_methods' : 'time: instantaneous',
   'coordinates': 'yc xc plev',
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
   'coordinates': 'yc xc plev',
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
   'coordinates': 'yc xc plev',
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
   'coordinates': 'yc xc plev',
 },
 'pr' : {
   'interval' : HOURLY3,
   'source' : 'MMOUTP',
   'long_name' : 'Precipitation',
   'units' : 'kg m-2 s-1',
   'standard_name' : 'precipitation_flux',
   'cell_methods' : 'time: mean (interval: 3 hours)',
   'coordinates': 'yc xc',
 },
 'soilt4' : {
    'ncsource': 'soil_t_4',
   'interval' : HOURLY3,
   'source' : 'MMOUTP',
   'long_name' : 'Soil Temperature at Layer 4',
   'units' : 'K',
   'standard_name' : 'unknown',
   'cell_methods' : 'time: instantaneous',
 },
 'tas' : {
    'ncsource': 't2',
   'interval' : HOURLY3,
   'source' : 'MMOUTP',
   'long_name' : 'Surface Air Temperature',
   'units' : 'K',
   'standard_name' : 'air_temperature',
   'cell_methods' : 'time: instantaneous',
 },
 'soilt1' : {
    'ncsource': 'soil_t_1',
   'interval' : HOURLY3,
   'source' : 'MMOUTP',
   'long_name' : 'Soil Temperature at Layer 1',
   'units' : 'K',
   'standard_name' : 'unknown',
   'cell_methods' : 'time: instantaneous',
 },
 'tasmax' : {
    'ncsource': 'tmax',
   'interval' : DAILY,
   'source' : 'NCOUT',
   'long_name' : 'Maximum Daily Surface Air Temperature',
   'units' : 'K',
   'standard_name' : 'air_temperature',
   'cell_methods' : 'time: maximum within days',
   'npfunc': numpy.max,
 },
 'tasmin' : {
    'ncsource': 'tmin',
   'interval' : DAILY,
   'source' : 'NCOUT',
   'long_name' : 'Minimum Daily Surface Air Temperature',
   'units' : 'K',
   'standard_name' : 'air_temperature',
   'cell_methods' : 'time: minimum within days',
   'npfunc': numpy.min,
 },
 'spdmax' : {
    'ncsource': 'maxwnd',
   'interval' : DAILY,
   'source' : 'NCOUT',
   'long_name' : 'Maximum Daily 10-Meter Wind Speed',
   'units' : 'm s-1',
   'standard_name' : 'wind_speed_of_gust',
   'cell_methods' : 'time: maximum within days',
   'npfunc': numpy.max,
 },
 'sic' : {
    'ncsource': 'seaice',
   'interval' : DAILY,
   'source' : 'MMOUTP',
   'long_name' : 'Daily Average Sea-ice Fraction',
   'units' : 'fraction in [0,1]',
   'standard_name' : 'sea_ice_fraction',
   'npfunc': numpy.average,
 }
}



def create_file(VNAME, ts0, ts1):
    """
    Create a CF compliant file for NARCCAP
    """
    fp = '../final/%s_MM5I_hadcm3_%s.nc' % (VNAME, ts0.strftime("%Y%m%d%H"))
    if PLEVEL is not None:
        fp = '../final/%s_MM5I_hadcm3_p%03i_%s.nc' % (VNAME, 
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
    nc.table_id = 'Table 2'
    nc.project_id = 'NARCCAP'
    nc.source = 'MM5(2002): atmosphere: MM5v3.6.3 non-hydrostatic; sst/sea ice: AMIPII; land: Noah;  Convection: Kain-Fritsch 2; Radiation: RRTM; PBL: MRF; Explicit Moisture: Reisner Mixed-Phase; Buffer: 15 point exponential; Horizontal Resolution: 50km; Vertical Levels: 24'
    nc.institution = 'ISU (Iowa State University, Ames, Iowa, USA)'

    tsteps = (ts1.year - ts0.year) * 360
    if VARS[VNAME]['interval'] == HOURLY3:
        tsteps *= 8
    print ' + Created NetCDF File %s has %s time steps' % (fp, tsteps)
    nc.createDimension('time', 0)
    nc.createDimension('bnds', 2)
    nc.createDimension('xc', 123)
    nc.createDimension('yc', 99)

    # Create Time Dimension
    tm = nc.createVariable('time', 'd', ('time',))
    tm.long_name = 'time'
    tm.standard_name = 'time'
    tm.axis = 'T'
    tm.calendar = '360_day'
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

    # write tm
    offset = (ts0.year - TIMES[0].year) * 360
    if VARS[VNAME]['interval'] == HOURLY3:
        tm[:] = offset + numpy.arange(0.125, (tsteps/8) + 0.125, 0.125)
        # write tmb
        tmb[:,0] = offset + numpy.arange(0., (tsteps/8), 0.125)
        tmb[:,1] = offset + numpy.arange(0.125, (tsteps/8)+0.125, 0.125)
    else:
        tm[:] = offset + numpy.arange(0., tsteps , 1.)
        # write tmb
        tmb[:,0] = offset + numpy.arange(0., tsteps, 1.)
        tmb[:,1] = offset + numpy.arange(1., tsteps+1., 1.)
   
    nc2 = netCDF4.Dataset('%s/%s_DOMAIN1_0001.nc' % (
                                DATADIR, VARS[VNAME]['source']), 'r')
    # write lat
    lat[:] = nc2.variables['latitcrs'][15:-15,15:-16]
    lon[:] = nc2.variables['longicrs'][15:-15,15:-16]
    xc[:] = numpy.arange(15,138) * nc2.variables['grid_ds'][:] * 1000.0
    yc[:] = numpy.arange(15,114) * nc2.variables['grid_ds'][:] * 1000.0
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

    # Now we dance
    cnter = 0
    for yr in range(ts0.year, ts1.year):
        for mo in range(1,13):
            for dy in range(1,31):
                i = (yr - TIMES[0].year) * 36 + 1
                i += (mo -1) * 3
                i += ((dy-1)/ 10)
                offset1 = ((dy - 1) * 8 + 2) % 80
                offset2 = ((dy) * 8 + 2) % 80
                if offset2 < offset1:
                    offset2 = -1
                fp2 = '%s/%s_DOMAIN1_%04i.nc' % (DATADIR, VARS[VNAME]['source'], i,)
                if not os.path.isfile(fp2):
                    print 'Missing File: %s, continuing' % (fp2,)
                    continue
                nc2 = netCDF4.Dataset(fp2, 'r')
                data = VARS[VNAME]['npfunc'](nc2.variables[VARS[VNAME]['ncsource']][offset1:offset2,15:-15,15:-16], axis=0)
                nc2.close()
                # Figure out which files we need!
                if dy in [10,20,30]: # Uh oh, need two files!
                    fp2 = '%s/%s_DOMAIN1_%04i.nc' % (DATADIR, VARS[VNAME]['source'], i+1)
                    nc2 = netCDF4.Dataset(fp2, 'r')
                    data2 = VARS[VNAME]['npfunc'](nc2.variables[VARS[VNAME]['ncsource']][:2,15:-15,15:-16], axis=0)
                    data = numpy.where( data2 > data, data2, data )
                    nc2.close()

                ncv[cnter] = data.astype('f')
                print '%s-%s-%s OFF1: %s OFF2: %s I: %s' % (yr, mo, dy, 
                    offset1, offset2, i)
                cnter += 1
    nc.close()


def compute3h(VNAME, fp, ts0, ts1):
    """
    This is just a straight dumping of data from the NC or MMOUTP files
    to the resulting netCDF file.
    """
    offset = (ts0.year - TIMES[0].year) * 36 + 1
    offset2 = (ts1.year - TIMES[0].year) * 36 + 1
    nc = netCDF4.Dataset(fp, 'a')
    ncv = nc.variables[VNAME]
    # Instead of writing data one 10 day chunk at a time, load it all up
    # first into a variable and then do a bulk write, I think this is 
    # much faster, maybe not...
    #result = numpy.zeros( numpy.shape(ncv), 'f')
    v = 0
    for i in range(offset,offset2):
        fp2 = '%s/%s_DOMAIN1_%04i.nc' % (DATADIR, VARS[VNAME]['source'], i,)
        if not os.path.isfile(fp2):
            print 'Missing File: %s, continuing' % (fp2,)
            v += 80
            continue
        nc2 = netCDF4.Dataset(fp2, 'r')
        if VNAME == 'pr':
            data = nc2.variables['rain_con'][:,15:-15,15:-16] + nc2.variables['rain_non'][:,15:-15,15:-16]
            # cm to mm to kg/m2/s
            data = (data * 10.0 / 10800.0).astype('f')
        else:
            ncs = VARS[VNAME]['ncsource']
            if PLEVEL is not None:
                l = PLEVELS.index(PLEVEL)
                data = nc2.variables[ncs][:,l,15:-15,15:-16]
            else:
                data = nc2.variables[ncs][:,15:-15,15:-16]
        ncv[v:v+80] = data
        nc2.close()
        v += 80
    #ncv[:] = result
    nc.close()

for i in range(len(TIMES)-1):
    ts0 = TIMES[i]
    ts1 = TIMES[i+1]
    fp = create_file(VNAME, ts0, ts1 )
    if VARS[VNAME]['interval'] == HOURLY3:
        compute3h(VNAME, fp, ts0, ts1)
    else:
        compute1d(VNAME, fp, ts0, ts1)
