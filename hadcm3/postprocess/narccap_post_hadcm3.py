"""
Uber script to generate NARCCAP archive specification files
This consumes the result of running the preprocess.py script
"""
from __future__ import print_function
import os
import datetime
import sys

import netCDF4
import numpy as np

sys.path.insert(0, "../../pylib")
import common

if len(sys.argv) < 4:
    print("Usage: python narccap_post.py RUNID VNAME")
    print("  RUNID is S: scenario  C: contempory")
    print("  DSON  is T: true/on  F: false/off")
    print("  VNAME is the standard NARCCAP Variable Name")
    print("Example: python narccap_post.py C pr")
    sys.exit()

RUNID = sys.argv[1]
DSON = sys.argv[2]
VNAME = sys.argv[3]
if len(sys.argv) == 5:
    PLEVEL = int(sys.argv[4])
else:
    PLEVEL = None
META = {}
if RUNID == "S":
    DATADIR = "../Run.scenario"
    META["title"] = (
        "ISU MM5 model output prepared for NARCCAP "
        "scenario climate using HADCM3"
    )
    META["prefix"] = "MM5I"
    META["experiment_id"] = "scenario climate using HADCM3"
    TIMES = [
        datetime.datetime(2038, 1, 1),
        datetime.datetime(2041, 1, 1),
        datetime.datetime(2046, 1, 1),
        datetime.datetime(2051, 1, 1),
        datetime.datetime(2056, 1, 1),
        datetime.datetime(2061, 1, 1),
        datetime.datetime(2066, 1, 1),
        datetime.datetime(2071, 1, 1),
    ]
if RUNID == "C":
    DATADIR = "../Run.contemporary"
    META["title"] = (
        "ISU MM5 model output prepared for NARCCAP "
        "contemporary climate using HADCM3"
    )
    META["prefix"] = "MM5I"
    META["experiment_id"] = "contemporary climate using HADCM3"
    TIMES = [
        datetime.datetime(1968, 1, 1),
        datetime.datetime(1971, 1, 1),
        datetime.datetime(1976, 1, 1),
        datetime.datetime(1981, 1, 1),
        datetime.datetime(1986, 1, 1),
        datetime.datetime(1991, 1, 1),
        datetime.datetime(1996, 1, 1),
        datetime.datetime(2001, 1, 1),
    ]
if RUNID == "T":
    DATADIR = "data.TEST"
    META["title"] = "ISU MM5 HADCM3 TEST"
    META["prefix"] = "TEST"
    META["experiment_id"] = "TEST 7 JUL 2011"
    TIMES = [
        datetime.datetime(1968, 1, 1),
        datetime.datetime(1971, 1, 1),
        datetime.datetime(1976, 1, 1),
        datetime.datetime(1981, 1, 1),
        datetime.datetime(1986, 1, 1),
        datetime.datetime(1991, 1, 1),
        datetime.datetime(1996, 1, 1),
        datetime.datetime(2001, 1, 1),
    ]
if DSON == "T":
    DATADIR += ".deepsoil_on"
else:
    DATADIR += ".deepsoil_off"

HOURLY3, DAILY = (1, 2)

LOOP1 = None
LOOP2 = None


def create_file(VNAME, ts0, ts1):
    """
    Create a CF compliant file for NARCCAP
    """
    fp = "../final/%s_%s_hadcm3_%s03.nc" % (
        VNAME,
        META["prefix"],
        ts0.strftime("%Y%m%d"),
    )
    if PLEVEL is not None:
        fp = ("../final/%s_%s_hadcm3_p%03i_%s.nc") % (
            VNAME,
            META["prefix"],
            PLEVEL,
            ts0.strftime("%Y%m%d%H"),
        )
    nc = netCDF4.Dataset(fp, "w", format="NETCDF3_CLASSIC")
    nc.Conventions = "CF-1.0"
    nc.title = META["title"]
    nc.history = ("Subsetted by archiver and python netCDF routines, %s") % (
        datetime.datetime.now().strftime("%Y-%m-%d %H:%M"),
    )
    nc.contact1 = "Daryl Herzmann (akrherz@iastate.edu)"
    nc.contact2 = "3015 Agronomy Hall, Iowa State Univ.,Ames, Iowa, USA"
    nc.realization = "1"
    nc.experiment_id = META["experiment_id"]
    nc.table_id = "Table %s" % (common.VARS[VNAME]["table"],)
    nc.project_id = "NARCCAP"
    nc.source = (
        "MM5(2002): atmosphere: MM5v3.6.3 non-hydrostatic; "
        "sst/sea ice: AMIPII; land: Noah;  Convection: "
        "Kain-Fritsch 2; Radiation: RRTM; PBL: MRF; Explicit "
        "Moisture: Reisner Mixed-Phase; Buffer: "
        "15 point exponential; Horizontal Resolution: 50km; "
        "Vertical Levels: 24; Deepsoil: %s"
    ) % (DSON,)
    nc.institution = "ISU (Iowa State University, Ames, Iowa, USA)"

    tsteps = (ts1.year - ts0.year) * 360
    if common.VARS[VNAME]["interval"] == HOURLY3:
        tsteps *= 8
    print(" + Created NetCDF File %s has %s time steps" % (fp, tsteps))
    nc.createDimension("time", 0)
    if "cell_methods" in common.VARS[VNAME]:
        nc.createDimension("bnds", 2)
    nc.createDimension("xc", 124)
    nc.createDimension("yc", 99)
    latgrid = "latitcrs"
    longrid = "longicrs"

    # Create Time Dimension
    tm = nc.createVariable("time", "d", ("time",))
    tm.long_name = "time"
    tm.standard_name = "time"
    tm.axis = "T"
    tm.calendar = "360_day"
    tm.units = "days since %s 00:00:0.0" % (TIMES[0].strftime("%Y-%m-%d"),)
    if "cell_methods" in common.VARS[VNAME]:
        tm.bounds = "time_bnds"

        tmb = nc.createVariable("time_bnds", "d", ("time", "bnds"))

    xc = nc.createVariable("xc", "d", ("xc",))
    xc.long_name = "x-coordinate Cartesian system"
    xc.standard_name = "projection_x_coordinate"
    xc.axis = "X"
    xc.units = "m"

    yc = nc.createVariable("yc", "d", ("yc",))
    yc.long_name = "y-coordinate Cartesian system"
    yc.standard_name = "projection_y_coordinate"
    yc.axis = "Y"
    yc.units = "m"

    lat = nc.createVariable("lat", "d", ("yc", "xc"))
    lat.long_name = "latitude"
    lat.standard_name = "latitude"
    lat.units = "degrees_north"

    lon = nc.createVariable("lon", "d", ("yc", "xc"))
    lon.long_name = "longitude"
    lon.standard_name = "longitude"
    lon.units = "degrees_east"

    p = nc.createVariable("Lambert_Conformal", "c", ())
    p.grid_mapping_name = "lambert_conformal_conic"
    p.false_easting = 3825000.0
    p.false_northing = 3187500.0

    v = nc.createVariable(VNAME, "f", ("time", "yc", "xc"), fill_value=1e20)
    v.units = common.VARS[VNAME]["units"]
    v.standard_name = common.VARS[VNAME]["standard_name"]
    v.long_name = common.VARS[VNAME]["long_name"]
    if "cell_methods" in common.VARS[VNAME]:
        v.cell_methods = common.VARS[VNAME]["cell_methods"]
    v.missing_value = np.array(1e20, v.dtype)
    v.coordinates = common.VARS[VNAME]["coordinates"]
    v.grid_mapping = "Lambert_Conformal"
    if "positive" in common.VARS[VNAME]:
        v.positive = common.VARS[VNAME]["positive"]

    # write tm
    offset = (ts0.year - TIMES[0].year) * 360
    if common.VARS[VNAME]["interval"] == HOURLY3:
        tm[:] = offset + np.arange(0.125, (tsteps / 8) + 0.125, 0.125)
        # write tmb
        if "cell_methods" in common.VARS[VNAME]:
            tmb[:, 0] = offset + np.arange(0.0, (tsteps / 8), 0.125)
            tmb[:, 1] = offset + np.arange(0.125, (tsteps / 8) + 0.125, 0.125)
    else:
        tm[:] = offset + np.arange(0.25, tsteps, 1.0)
        # write tmb
        if "cell_methods" in common.VARS[VNAME]:
            tmb[:, 0] = offset + np.arange(0.25, tsteps, 1.0)
            tmb[:, 1] = offset + np.arange(1.25, tsteps + 1.0, 1.0)

    nc2 = netCDF4.Dataset(
        "%s/%s_DOMAIN1_0001.nc" % (DATADIR, common.VARS[VNAME]["source"]), "r"
    )
    # write lat
    lat[:] = nc2.variables[latgrid][15:-15, 15:-15]
    lon[:] = nc2.variables[longrid][15:-15, 15:-15] + 360.0
    xc[:] = np.arange(15, 139) * nc2.variables["grid_ds"][:] * 1000.0
    yc[:] = np.arange(15, 114) * nc2.variables["grid_ds"][:] * 1000.0
    p.standard_parallel = [
        nc2.variables["stdlat_2"][:],
        nc2.variables["stdlat_1"][:],
    ]
    p.longitude_of_central_meridian = nc2.variables["coarse_cenlon"][:]
    p.latitude_of_projection_origin = nc2.variables["coarse_cenlat"][:]
    nc2.close()

    # Generate singletons
    common.do_singleton(nc, VNAME, PLEVEL)

    nc.close()
    return fp


def compute1d(VNAME, fp, ts0, ts1):
    """
    6z to 6z is our local day, so we need 9z to 6z inclusive.  The file
    starts at 3z, so offset 2 to 10
    """
    nc = netCDF4.Dataset(fp, "a")
    ncv = nc.variables[VNAME]

    # Now we dance
    cnter = 0
    for yr in range(ts0.year, ts1.year):
        for mo in range(1, 13):
            for dy in range(1, 31):
                i = (yr - TIMES[0].year) * 36 + 1
                i += (mo - 1) * 3
                i += (dy - 1) / 10
                offset1 = ((dy - 1) * 8 + 2) % 80
                offset2 = ((dy) * 8 + 2) % 80
                if offset2 < offset1:
                    offset2 = -1
                fp2 = "%s/%s_DOMAIN1_%04i.nc" % (
                    DATADIR,
                    common.VARS[VNAME]["source"],
                    i,
                )
                if not os.path.isfile(fp2):
                    print("Missing File: %s, continuing" % (fp2,))
                    continue
                nc2 = netCDF4.Dataset(fp2, "r")
                data = common.VARS[VNAME]["npfunc"](
                    nc2.variables[common.VARS[VNAME]["ncsource"]][
                        offset1:offset2, 15:-15, 15:-15
                    ],
                    axis=0,
                )
                nc2.close()
                # Figure out which files we need!
                if dy in [10, 20, 30]:  # Uh oh, need two files!
                    fp2 = ("%s/%s_DOMAIN1_%04i.nc") % (
                        DATADIR,
                        common.VARS[VNAME]["source"],
                        i + 1,
                    )
                    if not os.path.isfile(fp2):
                        print("Missing File: %s, continuing" % (fp2,))
                        continue
                    nc2 = netCDF4.Dataset(fp2, "r")
                    data2 = common.VARS[VNAME]["npfunc"](
                        nc2.variables[common.VARS[VNAME]["ncsource"]][
                            :2, 15:-15, 15:-15
                        ],
                        axis=0,
                    )
                    data = np.where(data2 > data, data2, data)
                    nc2.close()

                ncv[cnter] = data.astype("f")
                print(
                    ("%s-%s-%s OFF1: %s OFF2: %s I: %s")
                    % (yr, mo, dy, offset1, offset2, i)
                )
                cnter += 1
    nc.close()


def compute3h(VNAME, fp, ts0, ts1):
    """
    This is just a straight dumping of data from the NC or MMOUTP files
    to the resulting netCDF file.
    """
    global LOOP1, LOOP2
    offset = (ts0.year - TIMES[0].year) * 36 + 1
    offset2 = (ts1.year - TIMES[0].year) * 36 + 1
    nc = netCDF4.Dataset(fp, "a")
    ncv = nc.variables[VNAME]
    # Instead of writing data one 10 day chunk at a time, load it all up
    # first into a variable and then do a bulk write, I think this is
    # much faster, maybe not...
    # result = np.zeros( np.shape(ncv), 'f')
    v = 0
    for i in range(offset, offset2):
        fp2 = "%s/%s_DOMAIN1_%04i.nc" % (
            DATADIR,
            common.VARS[VNAME]["source"],
            i,
        )
        if not os.path.isfile(fp2):
            print("Missing File: %s, continuing" % (fp2,))
            v += 80
            continue
        nc2 = netCDF4.Dataset(fp2, "r")
        if VNAME == "mrfso":
            data = (
                (
                    nc2.variables["soil_m_1"][:, 15:-15, 15:-15]
                    - nc2.variables["soil_w_1"][:, 15:-15, 15:-15]
                )
                * 0.10
                + (
                    nc2.variables["soil_m_2"][:, 15:-15, 15:-15]
                    - nc2.variables["soil_w_2"][:, 15:-15, 15:-15]
                )
                * 0.30
                + (
                    nc2.variables["soil_m_3"][:, 15:-15, 15:-15]
                    - nc2.variables["soil_w_3"][:, 15:-15, 15:-15]
                )
                * 0.60
                + (
                    nc2.variables["soil_m_4"][:, 15:-15, 15:-15]
                    - nc2.variables["soil_w_4"][:, 15:-15, 15:-15]
                )
                * 1.00
            ) * 1000.0

        elif VNAME == "mrso":
            data = (
                (nc2.variables["soil_m_1"][:, 15:-15, 15:-15]) * 0.10
                + (nc2.variables["soil_m_2"][:, 15:-15, 15:-15]) * 0.30
                + (nc2.variables["soil_m_3"][:, 15:-15, 15:-15]) * 0.60
                + (nc2.variables["soil_m_4"][:, 15:-15, 15:-15]) * 1.00
            ) * 1000.0

        elif VNAME == "prw":  # Precipitable Water
            # http://docs.lib.noaa.gov/rescue/mwr/067/mwr-067-04-0100.pdf
            bogus = nc2.variables["soil_m_1"][:, 15:-15, 15:-15]
            data = np.zeros(bogus.shape)
            plevels = nc2.variables["pressure"][:]
            for j in range(1, len(plevels) - 1):
                dp = plevels[j] - plevels[j + 1]
                q1 = nc2.variables["q"][:, j, 15:-15, 15:-15]
                q2 = nc2.variables["q"][:, j + 1, 15:-15, 15:-15]
                # inches, convert to mm , which is kg m-2
                data += 0.0002 * (dp) * ((q1 + q2) * 1000.0) * 25.4

        elif VNAME == "mrro":
            surface = nc2.variables["sfcrnoff"][:, 15:-15, 15:-15]
            subsurface = nc2.variables["ugdrnoff"][:, 15:-15, 15:-15]
            data = np.zeros(surface.shape)
            tsteps = surface.shape[0]
            for j in range(tsteps):
                if LOOP1 is None:
                    s0 = surface[j]
                    ss0 = subsurface[j]
                else:
                    s0 = surface[j] - LOOP1
                    ss0 = subsurface[j] - LOOP2
                LOOP1 = surface[j]
                LOOP2 = subsurface[j]
                data[j] = s0 + ss0

        elif VNAME == "mrros":
            surface = nc2.variables["sfcrnoff"][:, 15:-15, 15:-15]
            data = np.zeros(surface.shape)
            tsteps = surface.shape[0]
            for j in range(tsteps):
                if LOOP1 is None:
                    s0 = surface[j]
                else:
                    s0 = surface[j] - LOOP1
                LOOP1 = surface[j]
                data[j] = s0

        elif VNAME == "huss":
            q2 = nc2.variables["q2"][:, 15:-15, 15:-15]
            data = q2 / (1.0 + q2)

        elif VNAME in ["ua", "va"]:
            ncs = common.VARS[VNAME]["ncsource"]
            plev = list(nc2.variables["pressure"][:]).index(PLEVEL)
            # dot grid surrounds cross making this striding easy
            proj = (
                nc2.variables[ncs][:, plev, :-1, :-1]
                + nc2.variables[ncs][:, plev, 1:, 1:]
                + nc2.variables[ncs][:, plev, :-1, 1:]
                + nc2.variables[ncs][:, plev, 1:, :-1]
            ) / 4.0

            data = proj[:, 15:-15, 15:-15]

        elif VNAME == "pr":
            data = (
                nc2.variables["rain_con"][:, 15:-15, 15:-15]
                + nc2.variables["rain_non"][:, 15:-15, 15:-15]
            )
            # Floor data at zero, prevent small negative numbers
            data = np.where(data > 0, data, 0)
        elif VNAME == "prc":
            data = nc2.variables["rain_con"][:, 15:-15, 15:-15]
            # Floor data at zero, prevent small negative numbers
            data = np.where(data > 0, data, 0)
        elif VNAME == "zg500":
            plev = list(nc2.variables["pressure"][:]).index(500)
            data = nc2.variables["h"][:, plev, 15:-15, 15:-15]

        else:
            ncs = common.VARS[VNAME]["ncsource"]
            if PLEVEL is not None:
                plev = list(nc2.variables["pressure"][:]).index(PLEVEL)
                data = nc2.variables[ncs][:, plev, 15:-15, 15:-15]
            else:
                data = nc2.variables[ncs][:, 15:-15, 15:-15]
        print(
            ("%s %5i %5i %.5f")
            % (
                fp2,
                v,
                v + 80,
                np.max(data) / common.VARS[VNAME].get("quo", 1.0),
            )
        )
        ncv[v : v + 80] = data / common.VARS[VNAME].get("quo", 1.0)
        nc2.close()
        v += 80
    # ncv[:] = result
    nc.close()


def main():
    """Go Main Go"""
    for tidx in range(len(TIMES) - 1):
        ts0 = TIMES[tidx]
        ts1 = TIMES[tidx + 1]
        fp = create_file(VNAME, ts0, ts1)
        if common.VARS[VNAME]["interval"] == HOURLY3:
            compute3h(VNAME, fp, ts0, ts1)
        else:
            compute1d(VNAME, fp, ts0, ts1)


if __name__ == "__main__":
    main()
