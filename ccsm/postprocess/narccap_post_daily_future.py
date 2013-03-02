import netCDF3
import mx.DateTime
import numpy

BASETS = mx.DateTime.DateTime(2038,1,1,0)

daily = {
 'tasmax' : {
   'long_name' : 'Maximum Daily Surface Air Temperature',
   'units' : 'K',
   'standard_name' : 'air_temperature',
   'cell_methods' : 'time: maximum within days',
 },
 'tasmin' : {
   'long_name' : 'Minimum Daily Surface Air Temperature',
   'units' : 'K',
   'standard_name' : 'air_temperature',
   'cell_methods' : 'time: minimum within days',
 },
 'spdmax' : {
   'long_name' : 'Maximum Daily 10-Meter Wind Speed',
   'units' : 'm s-1',
   'standard_name' : 'wind_speed_of_gust',
   'cell_methods' : 'time: maximum within days',
 },
 'sic' : {
   'long_name' : 'Daily Average Sea-ice Fraction',
   'units' : 'fraction in [0,1]',
   'standard_name' : 'sea_ice_fraction',
 }
}

times = [
 mx.DateTime.DateTime(2038,1,1),
 mx.DateTime.DateTime(2041,1,1),
 mx.DateTime.DateTime(2046,1,1),
 mx.DateTime.DateTime(2051,1,1),
 mx.DateTime.DateTime(2056,1,1),
 mx.DateTime.DateTime(2061,1,1),
 mx.DateTime.DateTime(2066,1,1),
 mx.DateTime.DateTime(2071,1,1)
]

def create_file(vname, ts0, ts1):
    """
    Create a CF compliant file for NARCCAP
    """
    fp = '%s_MM5I_hadcm3_%s.nc' % (vname, ts0.strftime("%Y%m%d%H"))
    nc = netCDF3.Dataset(fp, 'w')
    nc.Conventions = 'CF-1.0'
    nc.title = 'ISU MM5 model output prepared for NARCCAP future climate using HADCM3'
    nc.history = 'Subsetted by archiver and python netCDF routines'
    nc.contact1 = 'Daryl Herzmann (akrherz@iastate.edu)'
    nc.contact2 = '3015 Agronomy Hall, Iowa State Univ.,Ames, Iowa, USA'
    nc.realization = '1'
    nc.experiment_id = 'future climate using HADCM3'
    # Should be table 1 for daily ?
    nc.table_id = 'Table 2'
    nc.project_id = 'NARCCAP'
    nc.source = 'MM5(2002): atmosphere: MM5v3.6.3 non-hydrostatic; sst/sea ice: AMIPII; land: Noah;  Convection: Kain-Fritsch 2; Radiation: RRTM; PBL: MRF; Explicit Moisture: Reisner Mixed-Phase; Buffer: 15 point exponential; Horizontal Resolution: 50km; Vertical Levels: 24'
    nc.institution = 'ISU (Iowa State University, Ames, Iowa, USA)'


    nc.createDimension('time', 0)
    nc.createDimension('bnds', 2)
    nc.createDimension('xc', 154)
    nc.createDimension('yc', 129)

    # Create Time Dimension
    tm = nc.createVariable('time', 'd', ('time',))
    tm.long_name = 'time'
    tm.standard_name = 'time'
    tm.axis = 'T'
    tm.calendar = '360_day'
    tm.units = 'days since %s 00:00:0.0' % (ts0.strftime("%Y-%m-%d"),)
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
    p.false_northing = 3200000.

    v = nc.createVariable(vname, 'f', ('time','yc','xc'))
    v.units = daily[vname]['units'] 
    v.standard_name = daily[vname]['standard_name'] 
    v.long_name = daily[vname]['long_name'] 
    v.cell_methods = daily[vname]['cell_methods']
    v.missing_value = numpy.array(1e20, v.dtype)
    v._FillValue = numpy.array(1e20, v.dtype)
    v.coordinates = 'lon lat'
    v.grid_mapping = 'Lambert_Conformal'

    # write tm
    steps  = (ts1.year - ts0.year) * 360
    offset = (ts0.year - BASETS.year) * 360
    tm[:] = offset + numpy.arange(0, steps, 1.0)

    # write tmb
    tmb[:,0] = numpy.arange(0., steps, 1.0)
    tmb[:,1] = numpy.arange(1., steps+1., 1.0)


    nc2 = netCDF3.Dataset('data.FUTURE/MMOUTP_DOMAIN1_0001.nc', 'r')
    # write lat
    lat[:] = nc2.variables['latitcrs'][:]
    lon[:] = nc2.variables['longicrs'][:]
    xc[:] = numpy.arange(0,154) * nc2.variables['grid_ds'][:] * 1000.0
    yc[:] = numpy.arange(0,129) * nc2.variables['grid_ds'][:] * 1000.0
    p.standard_parallel = [nc2.variables['stdlat_2'][:], nc2.variables['stdlat_1'][:]]
    p.longitude_of_central_meridian = nc2.variables['coarse_cenlon'][:]
    p.latitude_of_projection_origin = nc2.variables['coarse_cenlat'][:]
    nc2.close()

    nc.close()
    return fp

def compute(vname, fp, ts0, ts1):
    """
    6z to 6z is our local day, so we need 9z to 6z inclusive.  The file
    starts at 3z, so offset 2 to 10
    """
    nc = netCDF3.Dataset(fp, 'a')
    ncv = nc.variables[vname]

    # Now we dance
    cnter = 0
    for yr in range(ts0.year, ts1.year):
        for mo in range(1,13):
            for dy in range(1,31):
                i = (yr - BASETS.year) * 36 + 1
                i += (mo -1) * 3
                i += ((dy-1)/ 10)
                offset1 = ((dy - 1) * 8 + 2) % 80
                offset2 = ((dy) * 8 + 2) % 80
                if offset2 < offset1:
                    offset2 = -1
                nc2 = netCDF3.Dataset('data.FUTURE/NCOUT_DOMAIN1_%04i.nc' % (i,), 'r')
                data = numpy.min(nc2.variables['tmin'][offset1:offset2], axis=0)
                nc2.close()
                # Figure out which files we need!
                if dy in [10,20,30]: # Uh oh, need two files!
                    nc2 = netCDF3.Dataset('data.FUTURE/NCOUT_DOMAIN1_%04i.nc' % (i+1,), 'r')
                    data2 = numpy.min(nc2.variables['tmin'][:2], axis=0)
                    data = numpy.where( data2 > data, data2, data )
                    nc2.close()

                ncv[cnter] = data.astype('f')
                print '%s-%s-%s OFF1: %s OFF2: %s I: %s' % (yr, mo, dy, 
                    offset1, offset2, i)
                cnter += 1
    nc.close()

for i in range(len(times)-1):
    ts0 = times[i]
    ts1 = times[i+1]
    fp = create_file('tasmin', ts0, ts1 )
    compute('tasmin', fp, ts0, ts1)
