"""
Uber script to generate NARCCAP archive specification files
This consumes the result of running the preprocess.py script
"""

from __future__ import print_statement
import netCDF4
import mx.DateTime
import numpy
import sys

sys.path.insert(0, "../../pylib")
import common

if len(sys.argv) < 3:
    print("Usage: python narccap_post.py RUNID VNAME")
    print("  RUNID is S: scenario  C: contempory")
    print("  VNAME is the standard NARCCAP Variable Name")
    print("Example: python narccap_post.py C pr")
    sys.exit()

RUNID = sys.argv[1]
VNAME = sys.argv[2]
if len(sys.argv) == 4:
    PLEVEL = int(sys.argv[3])
else:
    PLEVEL = None
META = {}
if RUNID == "C":
    DATADIR = "../Run.NCEP"
    META["title"] = (
        "ISU MM5 model output prepared for NARCCAP present-day climate using NCEP/DOE Reanalysis"
    )
    META["prefix"] = "MM5I"
    META["experiment_id"] = "present-day climate using NCEP/DOE Reanalysis"
    TIMES = [
        mx.DateTime.DateTime(1979, 1, 1),
        mx.DateTime.DateTime(1981, 1, 1),
        mx.DateTime.DateTime(1986, 1, 1),
        mx.DateTime.DateTime(1991, 1, 1),
        mx.DateTime.DateTime(1996, 1, 1),
        mx.DateTime.DateTime(2001, 1, 1),
        mx.DateTime.DateTime(2004, 11, 30),
    ]

HOURLY3, DAILY = (1, 2)

LOOP1 = None
LOOP2 = None


def create_file(VNAME, ts0, ts1):
    """
    Create a CF compliant file for NARCCAP
    """
    fp = "../final/%s_%s_ncep_%s03.nc" % (
        VNAME,
        META["prefix"],
        ts0.strftime("%Y%m%d"),
    )
    if PLEVEL is not None:
        fp = "../final/%s_%s_ncep_p%03i_%s.nc" % (
            VNAME,
            META["prefix"],
            PLEVEL,
            ts0.strftime("%Y%m%d%H"),
        )
    nc = netCDF4.Dataset(fp, "w", format="NETCDF3_CLASSIC")
    nc.Conventions = "CF-1.0"
    nc.title = META["title"]
    nc.history = "Subsetted by archiver and python netCDF routines, %s" % (
        mx.DateTime.now().strftime("%Y-%m-%d %H:%M"),
    )
    nc.contact1 = "Daryl Herzmann (akrherz@iastate.edu)"
    nc.contact2 = "3015 Agronomy Hall, Iowa State Univ.,Ames, Iowa, USA"
    nc.realization = "1"
    nc.experiment_id = META["experiment_id"]
    nc.table_id = "Table %s" % (common.VARS[VNAME]["table"],)
    nc.project_id = "NARCCAP"
    nc.source = "MM5(2002): atmosphere: MM5v3.6.3 non-hydrostatic; sst/sea ice: AMIPII; land: Noah;  Convection: Kain-Fritsch 2; Radiation: RRTM; PBL: MRF; Explicit Moisture: Reisner Mixed-Phase; Buffer: 15 point exponential; Horizontal Resolution: 50km; Vertical Levels: 24"
    nc.institution = "ISU (Iowa State University, Ames, Iowa, USA)"

    tsteps = int((ts1 - ts0).days)
    if common.VARS[VNAME]["interval"] == HOURLY3:
        tsteps *= 8
    print(" + Created NetCDF File %s has %s time steps" % (fp, tsteps))
    nc.createDimension("time", 0)
    if common.VARS[VNAME].has_key("cell_methods"):
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
    tm.calendar = "gregorian"
    tm.units = "days since %s 00:00:0.0" % (TIMES[0].strftime("%Y-%m-%d"),)
    if common.VARS[VNAME].has_key("cell_methods"):
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
    if common.VARS[VNAME].has_key("cell_methods"):
        v.cell_methods = common.VARS[VNAME]["cell_methods"]
    v.missing_value = numpy.array(1e20, v.dtype)
    v.coordinates = common.VARS[VNAME]["coordinates"]
    v.grid_mapping = "Lambert_Conformal"
    if common.VARS[VNAME].has_key("positive"):
        v.positive = common.VARS[VNAME]["positive"]

    # write tm
    offset = int((ts0 - TIMES[0]).days)
    if common.VARS[VNAME]["interval"] == HOURLY3:
        tm[:] = offset + numpy.arange(0.125, (tsteps / 8) + 0.125, 0.125)
        # write tmb
        if common.VARS[VNAME].has_key("cell_methods"):
            tmb[:, 0] = offset + numpy.arange(0.0, (tsteps / 8), 0.125)
            tmb[:, 1] = offset + numpy.arange(
                0.125, (tsteps / 8) + 0.125, 0.125
            )
    else:
        tm[:] = offset + numpy.arange(0.25, tsteps, 1.0)
        # write tmb
        if common.VARS[VNAME].has_key("cell_methods"):
            tmb[:, 0] = offset + numpy.arange(0.25, tsteps, 1.0)
            tmb[:, 1] = offset + numpy.arange(1.25, tsteps + 1.0, 1.0)

    nc2 = netCDF4.Dataset(
        "%s/%s_DOMAIN1_0001.nc" % (DATADIR, common.VARS[VNAME]["source"]), "r"
    )
    nc3 = netCDF4.Dataset("%s/MMOUTP_DOMAIN1_0001.nc" % (DATADIR), "r")
    # write lat
    # write lat
    lat[:] = nc3.variables[latgrid][15:-15, 15:-15]
    lon[:] = nc3.variables[longrid][15:-15, 15:-15] + 360.0
    nc3.close()
    xc[:] = numpy.arange(15, 139) * nc2.variables["grid_ds"][:] * 1000.0
    yc[:] = numpy.arange(15, 114) * nc2.variables["grid_ds"][:] * 1000.0
    p.standard_parallel = [
        nc2.variables["stdlat_2"][:],
        nc2.variables["stdlat_1"][:],
    ]
    p.longitude_of_central_meridian = nc2.variables["coarse_cenlon"][:]
    p.latitude_of_projection_origin = nc2.variables["coarse_cenlat"][:]
    nc2.close()

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

    lookfor = ts0.strftime("minutes since %Y-%m-%d")
    # Figure out when our data begins!
    for i in range(1, 1225):
        fp2 = "%s/MMOUTP_DOMAIN1_%04i.nc" % (DATADIR, i)
        nc2 = netCDF4.Dataset(fp2, "r")
        # minutes since 1983-11-01 03:00:16
        if nc2.variables["time"].units.find(lookfor) == 0:
            nc2.close()
            break
        nc2.close()

    print("For timestamp %s We found file %s" % (ts0, fp2))
    now = ts0 + mx.DateTime.RelativeDateTime(hours=6)
    oneday = mx.DateTime.RelativeDateTime(days=1)
    cnter = 0
    while now < ts1:
        # if now.month == 2 and now.day == 29:
        #    now += oneday
        #    continue
        fp2 = "%s/%s_DOMAIN1_%04i.nc" % (
            DATADIR,
            common.VARS[VNAME]["source"],
            i,
        )
        nc2 = netCDF4.Dataset(fp2, "r")
        fp3 = "%s/MMOUTP_DOMAIN1_%04i.nc" % (DATADIR, i)
        nc3 = netCDF4.Dataset(fp3, "r")
        # Figure out timestamp base
        # Figure out timestamp base
        tsbase = mx.DateTime.strptime(
            nc3.variables["time"].units[14:27], "%Y-%m-%d %H"
        )
        nc3.close()
        tsteps = len(nc2.variables["time"][:])

        # Okay, we need to go from 6z to 6z of the next day
        offset1 = int((now - tsbase).hours / 3)
        offset2 = int(((now + oneday) - tsbase).hours / 3)
        if offset2 > tsteps:
            i += 1
            offset2 = tsteps
        data = common.VARS[VNAME]["npfunc"](
            nc2.variables[common.VARS[VNAME]["ncsource"]][
                offset1:offset2, 15:-15, 15:-15
            ],
            axis=0,
        )
        nc2.close()

        if offset2 == tsteps:  # Need to step ahead and get next file
            fp2 = "%s/%s_DOMAIN1_%04i.nc" % (
                DATADIR,
                common.VARS[VNAME]["source"],
                i,
            )
            nc2 = netCDF4.Dataset(fp2, "r")
            data2 = nc2.variables[common.VARS[VNAME]["ncsource"]][
                :1, 15:-15, 15:-15
            ]
            nc2.close()
            # print 'data IN', numpy.average(data)
            # print 'data2 IN', numpy.average(data2)
            if numpy.max(data2) != 0 and VNAME not in ["sic"]:
                data = common.VARS[VNAME]["npfunc2"](data, data2)
            else:
                print("Skipping TS2 Computation")
        print(
            "%s %s %02i %02i Avg: %.3f"
            % (
                now.strftime("%Y%m%d"),
                fp2,
                offset1,
                offset2,
                numpy.average(data),
            )
        )

        ncv[cnter] = data.astype("f")
        cnter += 1
        now += mx.DateTime.RelativeDateTime(days=1)
    nc.close()


def compute3h(VNAME, fp, ts0, ts1, running):
    """
    This is just a straight dumping of data from the NC or MMOUTP files
    to the resulting netCDF file.
    """
    global LOOP1, LOOP2
    lookfor = ts0.strftime("minutes since %Y-%m-%d")
    # Figure out when our data begins!
    for i in range(1, 1226):
        fp2 = "%s/MMOUTP_DOMAIN1_%04i.nc" % (DATADIR, i)
        nc2 = netCDF4.Dataset(fp2, "r")
        # minutes since 1983-11-01 03:00:16
        if nc2.variables["time"].units.find(lookfor) == 0:
            nc2.close()
            break
        nc2.close()

    print("For timestamp %s We found file %s" % (ts0, fp2))

    nc = netCDF4.Dataset(fp, "a")
    ncv = nc.variables[VNAME]
    total = len(nc.variables["time"][:])
    now = ts0
    v = 0
    while v < total:
        fp2 = "%s/%s_DOMAIN1_%04i.nc" % (
            DATADIR,
            common.VARS[VNAME]["source"],
            i,
        )
        nc2 = netCDF4.Dataset(fp2, "r")
        tsbase = mx.DateTime.strptime(
            nc2.variables["time"].units[14:27], "%Y-%m-%d %H"
        )
        tsteps = len(nc2.variables["time"][:])

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
            data = numpy.zeros(bogus.shape)
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
            data = numpy.zeros(surface.shape)
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
            data = numpy.zeros(surface.shape)
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
            l = list(nc2.variables["pressure"][:]).index(PLEVEL)
            # dot grid surrounds cross making this striding easy
            proj = (
                nc2.variables[ncs][:, l, :-1, :-1]
                + nc2.variables[ncs][:, l, 1:, 1:]
                + nc2.variables[ncs][:, l, :-1, 1:]
                + nc2.variables[ncs][:, l, 1:, :-1]
            ) / 4.0

            data = proj[:, 15:-15, 15:-15]

        elif VNAME in ["pr", "prc"]:  # Its acumulated
            if VNAME == "pr":
                pr = (
                    nc2.variables["rain_con"][:, 15:-15, 15:-15]
                    + nc2.variables["rain_non"][:, 15:-15, 15:-15]
                )
            else:
                pr = nc2.variables["rain_con"][:, 15:-15, 15:-15]
            data = numpy.zeros(pr.shape)
            for j in range(tsteps):
                if running is None:
                    running = numpy.zeros(pr[0].shape)
                data[j] = pr[j] - running
                running = pr[j]
            # Floor data at zero, prevent small negative numbers
            data = numpy.where(data > 0, data, 0)
        elif VNAME == "zg500":
            l = list(nc2.variables["pressure"][:]).index(500)
            data = nc2.variables["h"][:, l, 15:-15, 15:-15]
        else:
            ncs = common.VARS[VNAME]["ncsource"]
            if PLEVEL is not None:
                l = list(nc2.variables["pressure"][:]).index(PLEVEL)
                data = nc2.variables[ncs][:, l, 15:-15, 15:-15]
            else:
                data = nc2.variables[ncs][:, 15:-15, 15:-15]
        # Okay, we have the data var loaded up
        v2 = v + tsteps
        ed = tsteps
        if v2 > total:
            ed -= v2 - total
            v2 = total

        print(
            "i=%4i tsteps=%2i %5i %5i/%5i %.5f %s"
            % (
                i,
                tsteps,
                v,
                v2,
                total,
                numpy.max(data) / common.VARS[VNAME].get("quo", 1.0),
                numpy.shape(data),
            )
        )
        ncv[v:v2] = data[0:ed] / common.VARS[VNAME].get("quo", 1.0)
        nc2.close()
        v = v2
        i += 1
    nc.close()
    return running


if __name__ == "__main__":
    running = None
    for i in range(len(TIMES) - 1):
        ts0 = TIMES[i]
        ts1 = TIMES[i + 1]
        fp = create_file(VNAME, ts0, ts1)
        if common.VARS[VNAME]["interval"] == HOURLY3:
            running = compute3h(VNAME, fp, ts0, ts1, running)
        else:
            compute1d(VNAME, fp, ts0, ts1)
