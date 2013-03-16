import netCDF4
import numpy as np
import mx.DateTime
import os

times2 = [
 mx.DateTime.DateTime(1968,1,1),
 mx.DateTime.DateTime(1971,1,1),
 mx.DateTime.DateTime(1976,1,1),
 mx.DateTime.DateTime(1981,1,1),
 mx.DateTime.DateTime(1986,1,1),
 mx.DateTime.DateTime(1991,1,1),
 mx.DateTime.DateTime(1996,1,1),
 mx.DateTime.DateTime(2001,1,1)
]
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

def findIJ(lon,lat):
  nc = netCDF4.Dataset('../final/swe_MM5I_hadcm3_2038010103.nc', 'r')
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

results = np.zeros( (33,12), np.float)
models = np.zeros( (33,12), np.float)

i,j = findIJ(-93.62, 41.99)
for k in range(len(times)-1):
    ts0 = times[k]
    ts1 = times[k+1]
    fp = '../final/swe_MM5I_hadcm3_%s03.nc' % (ts0.strftime('%Y%m%d'),)
    if not os.path.isfile(fp):
        continue
    nc = netCDF4.Dataset(fp, 'r')
    data = nc.variables['swe'][:,j,i] 
    #print np.shape(nc.variables['tasmax'][:])
    #print nc.variables['lat'][j,i]
    #print nc.variables['lon'][j,i]
    cnter = 0
    for yr in range(ts0.year,ts1.year):
      for mo in range(1,13):
        offset = (yr-ts0.year)*360*8
        offset += (mo-1)*30*8
        offset2 = offset + (30*8)
        model = np.average(data[offset:offset2] )
        #sql = "SELECT avg(high) from alldata where stationid = 'ia0200' and year = %s and month = %s" % (yr, mo)
        #ccursor.execute(sql)
        #row = ccursor.fetchone()
        #print yr, mo, model, row[0],  model - float(row[0])
        #results[yr-1968,mo-1] = model - float(row[0])
        models[yr-times[0].year,mo-1] = model 
    nc.close()

avg = np.average(results,axis=0)
x = np.max(results,axis=0)
n = np.min(results,axis=0)
std = np.std(results,axis=0)

"""
print "%4s %6s %6s %6s %6s" % ("MON", "MIN", "AVG", "MAX", "STDDEV") 
for mo in range(1,13):
  print "%4s %6.2f %6.2f %6.2f %6.2f" % (mo, n[mo-1], avg[mo-1], x[mo-1], std[mo-1])

print "%4s %6.2f %6.2f %6.2f %6.2f" % ("ALL", np.min(results), np.average(results), np.max(results), np.std(results))

for yr in range(1968,2001):
  print "%4s %6.2f %6.2f %6.2f %6.2f" % (yr, np.min(results[yr-1968,:]), np.average(results[yr-1968,:]), np.max(results[yr-1968,:]), np.std(results[yr-1968,:]))
"""
for yr in range(times[0].year,times[-1].year,1):
  print "%4s" % (yr,),
  for mo in range(1,13):
    print "%7.1f" % ( models[yr-times[0].year,mo-1],),
  print "%7.1f" % ( np.sum(models[yr-times[0].year,:]),),
  print
