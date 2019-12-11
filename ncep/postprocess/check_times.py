import netCDF3
import mx.DateTime
import numpy
import os

lts = 0
for i in range(1, 1192):
    fp = "../Run.NCEP/MMOUTP_DOMAIN1_%04i.nc" % (i,)
    nc = netCDF3.Dataset(fp)
    tm = nc.variables["time"]
    tstr = tm.units.replace("minutes since ", " ")
    ts = mx.DateTime.strptime(tstr, "%Y-%m-%d %H:%M:%S")
    ts2 = ts + mx.DateTime.RelativeDateTime(minutes=tm[-1])
    if (ts.day >= 19 and ts.month == 2) or len(tm[:]) % 2 != 0:
        print "%s %s %s %s" % (
            fp,
            ts.strftime("%Y%m%d%H"),
            ts2.strftime("%Y%m%d%H"),
            numpy.shape(nc.variables["rain_con"][:]),
        )
    s = numpy.shape(nc.variables["rain_con"][:])
    lts = ts2
    nc.close()
