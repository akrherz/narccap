"""
Generate the table 4 variable lat
"""
import netCDF4
import datetime
import numpy

fn = 'lat_MM5I.nc'
nc = netCDF4.Dataset(fn, 'w', format='NETCDF3_CLASSIC')
nc.Conventions = 'CF-1.0'
nc.title = 'ISU MM5 model output prepared for NARCCAP'
nc.history = 'python netCDF routines, %s' % (datetime.datetime.now().strftime("%Y-%m-%d %H:%M"),)
nc.contact1 = 'Daryl Herzmann (akrherz@iastate.edu)'
nc.contact2 = '3015 Agronomy Hall, Iowa State Univ.,Ames, Iowa, USA'
nc.realization = '1'
nc.experiment_id = 'all'
nc.table_id = 'Table 4'
nc.project_id = 'NARCCAP'
nc.source = 'MM5(2002): atmosphere: MM5v3.6.3 non-hydrostatic; sst/sea ice: AMIPII; land: Noah;  Convection: Kain-Fritsch 2; Radiation: RRTM; PBL: MRF; Explicit Moisture: Reisner Mixed-Phase; Buffer: 15 point exponential; Horizontal Resolution: 50km; Vertical Levels: 24'
nc.institution = 'ISU (Iowa State University, Ames, Iowa, USA)'

nc.createDimension('xc', 124)
nc.createDimension('yc', 99)

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
lat.units = 'degrees_north'

p = nc.createVariable('Lambert_Conformal', 'c', ())
p.grid_mapping_name = "lambert_conformal_conic"
p.false_easting = 3825000.
p.false_northing = 3187500.

nc2 = netCDF4.Dataset('MMOUTP_DOMAIN1_1002.nc')

lat[:] = nc2.variables['latitcrs'][15:-15,15:-15]
xc[:] = numpy.arange(15,139) * nc2.variables['grid_ds'][:] * 1000.0
yc[:] = numpy.arange(15,114) * nc2.variables['grid_ds'][:] * 1000.0
p.standard_parallel = [nc2.variables['stdlat_2'][:], nc2.variables['stdlat_1'][:]]
p.longitude_of_central_meridian = nc2.variables['coarse_cenlon'][:]
p.latitude_of_projection_origin = nc2.variables['coarse_cenlat'][:]

nc2.close()
nc.close()
