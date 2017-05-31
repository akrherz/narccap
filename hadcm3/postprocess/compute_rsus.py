"""compute rsus by rsds * albedo

This is likely not strictly accurate, but close enough for horseshoes
"""
from __future__ import print_function
import sys
import glob
import datetime

import numpy as np
import netCDF4

RUNID = sys.argv[1]
META = {}
if RUNID == 'S':
    DATADIR = "../Run.scenario.deepsoil_off"
    META['title'] = ('ISU MM5 model output prepared for NARCCAP '
                     'scenario climate using HADCM3')
    META['prefix'] = 'MM5I'
    META['experiment_id'] = 'scenario climate using HADCM3'
    TIMES = [
        datetime.datetime(2038, 1, 1),
        datetime.datetime(2041, 1, 1),
        datetime.datetime(2046, 1, 1),
        datetime.datetime(2051, 1, 1),
        datetime.datetime(2056, 1, 1),
        datetime.datetime(2061, 1, 1),
        datetime.datetime(2066, 1, 1),
        datetime.datetime(2071, 1, 1)
    ]
if RUNID == 'C':
    DATADIR = "../Run.contemporary.deepsoil_off"
    META['title'] = ('ISU MM5 model output prepared for NARCCAP '
                     'contemporary climate using HADCM3')
    META['prefix'] = 'MM5I'
    META['experiment_id'] = 'contemporary climate using HADCM3'
    TIMES = [
        datetime.datetime(1968, 1, 1),
        datetime.datetime(1971, 1, 1),
        datetime.datetime(1976, 1, 1),
        datetime.datetime(1981, 1, 1),
        datetime.datetime(1986, 1, 1),
        datetime.datetime(1991, 1, 1),
        datetime.datetime(1996, 1, 1),
        datetime.datetime(2001, 1, 1)
    ]


def create_file(ts0, ts1):
    """
    Create a CF compliant file for NARCCAP
    """
    ncfn = '../final/rsus_MM5I_hadcm3_%s03.nc' % (ts0.strftime("%Y%m%d"),)
    nc = netCDF4.Dataset(ncfn, 'w', format='NETCDF3_CLASSIC')
    nc.Conventions = 'CF-1.0'
    nc.title = META['title']
    nc.history = ('rsds * albedo, computed on %s'
                  ) % (datetime.datetime.now().strftime("%Y-%m-%d %H:%M"),)
    nc.contact1 = 'Daryl Herzmann (akrherz@iastate.edu)'
    nc.contact2 = '3015 Agronomy Hall, Iowa State Univ.,Ames, Iowa, USA'
    nc.realization = '1'
    nc.experiment_id = META['experiment_id']
    nc.table_id = 'Table XX'
    nc.project_id = 'NARCCAP'
    nc.source = ('MM5(2002): atmosphere: MM5v3.6.3 non-hydrostatic; '
                 'sst/sea ice: AMIPII; land: Noah;  Convection: '
                 'Kain-Fritsch 2; Radiation: RRTM; PBL: MRF; Explicit '
                 'Moisture: Reisner Mixed-Phase; Buffer: '
                 '15 point exponential; Horizontal Resolution: 50km; '
                 'Vertical Levels: 24; Deepsoil: False')
    nc.institution = 'ISU (Iowa State University, Ames, Iowa, USA)'

    tsteps = int((ts1.year - ts0.year) * 360) * 8
    print(' + Created NetCDF File %s has %s time steps' % (ncfn, tsteps))
    nc.createDimension('time', 0)
    nc.createDimension('bnds', 2)
    nc.createDimension('xc', 124)
    nc.createDimension('yc', 99)
    latgrid = 'latitcrs'
    longrid = 'longicrs'

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
    lat.units = 'degrees_north'

    lon = nc.createVariable('lon', 'd', ('yc', 'xc'))
    lon.long_name = 'longitude'
    lon.standard_name = 'longitude'
    lon.units = 'degrees_east'

    p = nc.createVariable('Lambert_Conformal', 'c', ())
    p.grid_mapping_name = "lambert_conformal_conic"
    p.false_easting = 3825000.
    p.false_northing = 3187500.

    v = nc.createVariable('rsus', 'f', ('time', 'yc', 'xc'), fill_value=1e20)
    v.units = 'W m-2'
    v.standard_name = 'surface_upwelling_shortwave_flux_in_air'
    v.long_name = 'Surface Upwelling Shortwave Radiation'
    v.cell_methods = 'time: mean (interval: 3 hours)'
    v.missing_value = np.array(1e20, v.dtype)
    v.coordinates = "lon lat"
    v.grid_mapping = 'Lambert_Conformal'
    v.positive = "up"

    # write tm
    offset = (ts0.year - TIMES[0].year) * 360
    tm[:] = offset + np.arange(0.125, (tsteps/8) + 0.125, 0.125)
    tmb[:, 0] = offset + np.arange(0., (tsteps/8), 0.125)
    tmb[:, 1] = offset + np.arange(0.125, (tsteps/8)+0.125, 0.125)

    nc2 = netCDF4.Dataset(('%s/NCOUT_DOMAIN1_0001.nc'
                           ) % (DATADIR, ), 'r')
    # write lat
    lat[:] = nc2.variables[latgrid][15:-15, 15:-15]
    lon[:] = nc2.variables[longrid][15:-15, 15:-15] + 360.0
    xc[:] = np.arange(15, 139) * nc2.variables['grid_ds'][:] * 1000.0
    yc[:] = np.arange(15, 114) * nc2.variables['grid_ds'][:] * 1000.0
    p.standard_parallel = [nc2.variables['stdlat_2'][:],
                           nc2.variables['stdlat_1'][:]]
    p.longitude_of_central_meridian = nc2.variables['coarse_cenlon'][:]
    p.latitude_of_projection_origin = nc2.variables['coarse_cenlat'][:]
    nc2.close()

    nc.close()
    return ncfn


def process(i, ts):
    """Do this rsds file"""
    rsdsnc = netCDF4.Dataset(
        "../final/rsds_MM5I_hadcm3_%s03.nc" % (ts.strftime("%Y%m%d"),))
    albnc = netCDF4.Dataset(
        "../final/albsrfc_MM5I_hadcm3_%s03.nc" % (ts.strftime("%Y%m%d"),))
    rsusncfn = create_file(ts, TIMES[i+1])
    rsusnc = netCDF4.Dataset(rsusncfn, 'a')
    albvar = albnc.variables["albsrfc"]
    rsdsvar = rsdsnc.variables['rsds']
    rsusvar = rsusnc.variables['rsus']
    for tstep in range(rsusvar.shape[0]):
        val = albvar[tstep, :, :] * rsdsvar[tstep, :, :]
        print("step=%s/%s/%s max=%.4f" % (tstep, rsusvar.shape[0],
                                          albvar.shape[0], np.max(val)))
        rsusvar[tstep, :, :] = val
    rsusnc.close()


def main(argv):
    """Go Main Go"""
    for i, ts in enumerate(TIMES[:-1]):
        process(i, ts)


if __name__ == '__main__':
    main(sys.argv)
