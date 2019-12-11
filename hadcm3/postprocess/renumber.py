# Make all file names have 3 char IDs

import os

print "Hi Logan"
os.chdir("Run.scenario.deepsoil_off")
for i in range(1139, 1199):
    for p in ["MMOUTP", "NCOUT"]:
        for suffix in ["nc", "gz"]:
            if suffix == "gz" and p == "MMOUTP":
                p = "MMOUT"
            ofp = "%s_DOMAIN1_%04i.%s" % (p, i + 1, suffix)
            nfp = "%s_DOMAIN1_%04i.%s" % (p, i, suffix)
            if os.path.isfile(nfp):
                print nfp
                sys.exit()
            print ofp, nfp
            os.rename(ofp, nfp)
