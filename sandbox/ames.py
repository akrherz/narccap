import netCDF4
import numpy as np
import mx.DateTime
from pyIEM import mesonet
import iemdb

coop = iemdb.connect("coop", bypass=True)
ccursor = coop.cursor()

times = [
    mx.DateTime.DateTime(1968, 1, 1),
    mx.DateTime.DateTime(1970, 1, 1),
    mx.DateTime.DateTime(1975, 1, 1),
    mx.DateTime.DateTime(1980, 1, 1),
    mx.DateTime.DateTime(1985, 1, 1),
    mx.DateTime.DateTime(1990, 1, 1),
    mx.DateTime.DateTime(1995, 1, 1),
    mx.DateTime.DateTime(2000, 1, 1),
    mx.DateTime.DateTime(2001, 1, 1),
]
times = [
    mx.DateTime.DateTime(2038, 1, 1),
    mx.DateTime.DateTime(2041, 1, 1),
    mx.DateTime.DateTime(2046, 1, 1),
    mx.DateTime.DateTime(2051, 1, 1),
    mx.DateTime.DateTime(2056, 1, 1),
    mx.DateTime.DateTime(2061, 1, 1),
    mx.DateTime.DateTime(2066, 1, 1),
    mx.DateTime.DateTime(2071, 1, 1),
]


def findIJ(lon, lat):
    nc = netCDF4.Dataset("tasmax_MM5I_hadcm3_1970010100.nc", "r")
    lats = nc.variables["lat"][:]
    lons = nc.variables["lon"][:]
    (y, x) = np.shape(lats)
    mindist = 1000.0
    for j in range(y):
        for i in range(x):
            dist = ((lats[j, i] - lat) ** 2 + (lons[j, i] - lon) ** 2) ** 0.5
            if dist < mindist:
                myi = i
                myj = j
                mindist = dist
    nc.close()
    return myi, myj


results = np.zeros((33, 12), np.float)
models = np.zeros((33, 12), np.float)

i, j = findIJ(-93.62, 41.99)
for k in range(len(times) - 1):
    ts0 = times[k]
    ts1 = times[k + 1]
    nc = netCDF4.Dataset(
        "tasmax_MM5I_hadcm3_%s00.nc" % (ts0.strftime("%Y%m%d"),), "r"
    )
    data = nc.variables["tasmax"][:, j, i]
    # print np.shape(nc.variables['tasmax'][:])
    # print nc.variables['lat'][j,i]
    # print nc.variables['lon'][j,i]
    cnter = 0
    for yr in range(ts0.year, ts1.year):
        for mo in range(1, 13):
            offset = (yr - ts0.year) * 360
            offset += (mo - 1) * 30
            offset2 = offset + 30
            model = mesonet.k2f(np.average(data[offset:offset2]))
            # sql = "SELECT avg(high) from alldata where stationid = 'ia0200' and year = %s and month = %s" % (yr, mo)
            # ccursor.execute(sql)
            # row = ccursor.fetchone()
            # print yr, mo, model, row[0],  model - float(row[0])
            # results[yr-1968,mo-1] = model - float(row[0])
            models[yr - 2038, mo - 1] = model
    nc.close()

avg = np.average(results, axis=0)
x = np.max(results, axis=0)
n = np.min(results, axis=0)
std = np.std(results, axis=0)

"""
print "%4s %6s %6s %6s %6s" % ("MON", "MIN", "AVG", "MAX", "STDDEV") 
for mo in range(1,13):
  print "%4s %6.2f %6.2f %6.2f %6.2f" % (mo, n[mo-1], avg[mo-1], x[mo-1], std[mo-1])

print "%4s %6.2f %6.2f %6.2f %6.2f" % ("ALL", np.min(results), np.average(results), np.max(results), np.std(results))

for yr in range(1968,2001):
  print "%4s %6.2f %6.2f %6.2f %6.2f" % (yr, np.min(results[yr-1968,:]), np.average(results[yr-1968,:]), np.max(results[yr-1968,:]), np.std(results[yr-1968,:]))
"""
for yr in range(2038, 2071, 1):
    print "%4s" % (yr,),
    for mo in range(1, 13):
        print "%5.1f" % (models[yr - 2038, mo - 1],),
    print
