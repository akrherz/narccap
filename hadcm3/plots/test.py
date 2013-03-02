import iemplot
import netCDF4
import numpy

nc = netCDF4.Dataset("/tmp/snd_MM5I_hadcm3_2066010103.nc", 'r')

snd = nc.variables['snd']
lat = nc.variables['lat'][:]
lon = nc.variables['lon'][:]

data = numpy.max( snd[-2400:,:,:], 0)
print data.shape

import iemplot
cfg = {
 '_conus': True,
 'cnLevels' : [0,0.1,0.2,0.5,1,2,3,5,10,20,50,100,150],
 'cnLevelSelectionMode' : 'ExplicitLevels',
 'cnFillMode'	: 'CellFill',
}

tmpfp = iemplot.simple_grid_fill(lon, lat, data, cfg)
iemplot.postprocess(tmpfp, None)


nc.close()
