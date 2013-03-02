# Make all file names have 3 char IDs

import os
print 'Hi Logan'
os.chdir("data")
for i in range(73,1196):
	for p in ['MMOUTP','NCOUT']:
		ofp = "%s_DOMAIN1_%04i.nc" % (p, i+1)
		nfp = "%s_DOMAIN1_%04i.nc" % (p, i)
                if os.path.isfile(nfp):
                  print nfp
                  sys.exit()
		os.rename(ofp,nfp)
