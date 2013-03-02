# Do some diagnostics on the time dimension

import sys
import glob
import os
import lib
import mx.DateTime

ts0 = mx.DateTime.DateTime(1968,1,1,0)

datadir = sys.argv[1]
os.chdir( datadir )
# Look for any MMOUT and NCOUT files
files = glob.glob('MMOUT_DOMAIN1_[0-9][0-9][0-9]')
files.sort()
sz = 0
for file in files:
    if file ==  "MMOUTP_DOMAIN1_000":
        continue
    ts = lib.extract_times( file )
    print '%s, %s - %s' % (file, ts[0].strftime("%Y %m %d %H"),
                    ts[-1].strftime("%Y %m %d %H"))
    sz += len(ts)
    ts = ts[-1] - mx.DateTime.RelativeDateTime(hours=((sz-1)*3))
    if ts != ts0:
        print 'OFF!', ts
print 'Total dimension size is %s' % (sz,)
