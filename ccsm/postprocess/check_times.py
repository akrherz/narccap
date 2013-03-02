import netCDF3
import mx.DateTime
import numpy
import os

lts = 0
for i in range(1,1192):
  fp = "Run.contemporary.deepsoil_off/MMOUTP_DOMAIN1_%04i.nc" % (i,)
  nc = netCDF3.Dataset(fp)
  tm = nc.variables['time']
  tstr = tm.units.replace("minutes since ", ' ')
  ts = mx.DateTime.strptime(tstr, '%Y-%m-%d %H:%M:%S')
  print fp, ts,  numpy.shape( nc.variables['rain_con'][:] )
  s = numpy.shape( nc.variables['rain_con'][:] )
  if lts == ts:
    print 'DUP', fp
  lts = ts
  if s[0] != 80:
    print 'BAD', fp
    #cmd = "bbcp root@metl26:/mnt/sdd1/narccap/output.FUTURE/MMOUT_DOMAIN1_%03i.gz data.FUTURE/" % (i,)
    #os.system(cmd)
    #os.rename("data.FUTURE/MMOUT_DOMAIN1_%03i.gz" % (i,), "data.FUTURE/MMOUT_DOMAIN1_%04i.gz" % (i,))
    #os.unlink(fp)
  nc.close()
