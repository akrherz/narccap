# Make all file names have 3 char IDs

import os
print 'Hi Logan'
os.chdir("data")
for i in range(1,1000):
	for p in ['MMOUTP','NCOUT']:
		ofp = "%s_DOMAIN1_%03i.nc" % (p, i)
		nfp = "%s_DOMAIN1_%04i.nc" % (p, i)
		os.rename(ofp,nfp)
