# Make all file names have 3 char IDs

import os
os.chdir("Run.scenario.deepsoil_off")
for i in range(1,1000):
	for p in ['MMOUT','NCOUT']:
		ofp = "%s_DOMAIN1_%02i" % (p, i)
		nfp = "%s_DOMAIN1_%04i" % (p, i)
                if os.path.isfile(ofp):
                  print ofp
		  os.rename(ofp,nfp)
