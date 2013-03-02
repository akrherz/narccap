import netCDF4
import numpy as np
import mx.DateTime
from pyIEM import mesonet
import iemdb
import os

times = [
 mx.DateTime.DateTime(1968,1,1),
 mx.DateTime.DateTime(1971,1,1),
 mx.DateTime.DateTime(1976,1,1),
 mx.DateTime.DateTime(1981,1,1),
 mx.DateTime.DateTime(1986,1,1),
 mx.DateTime.DateTime(1991,1,1),
 mx.DateTime.DateTime(1996,1,1),
 mx.DateTime.DateTime(2001,1,1)
]

def findIJ(lon,lat):
  nc = netCDF4.Dataset('../final/pr_MM5I_hadcm3_1968010100.nc', 'r')
  lats = nc.variables['lat'][:]
  lons = nc.variables['lon'][:]
  (y,x) = np.shape(lats)
  mindist = 1000.
  for j in range(y):
     for i in range(x):
        dist = ((lats[j,i] - lat) ** 2 + (lons[j,i] - lon) ** 2 ) ** .5
        if dist < mindist:
          myi = i
          myj = j
          mindist = dist
  nc.close()
  return myi,myj


i,j = findIJ(-93.62, 41.99)
print i,j
sys.exit()
for k in range(len(times)-1):
    ts0 = times[k]
    ts1 = times[k+1]
    fp = '../final/tasmax_MM5I_hadcm3_%s00.nc' % (ts0.strftime('%Y%m%d'),)
    if not os.path.isfile(fp):
        continue
    nc = netCDF4.Dataset(fp, 'r')
    data = nc.variables['tasmax'][:,j,i] 
    #print np.shape(nc.variables['tasmax'][:])
    #print nc.variables['lat'][j,i]
    #print nc.variables['lon'][j,i]
    cnter = 0
    for yr in range(ts0.year,ts1.year):
      for mo in range(1,13):
        for dy in range(1,31):
            print yr, mo, dy, mesonet.k2f( data[cnter])
            cnter += 1
    nc.close()
