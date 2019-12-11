# Make all file names have 3 char IDs

import os

print "Hi Logan"
os.chdir("Run.NCEP")
for i in range(1224, 1188, -1):
    for p in ["MMOUTP", "NCOUT"]:
        for suffix in [".nc", ""]:
            if suffix == "" and p == "MMOUTP":
                p = "MMOUT"
            ofp = "%s_DOMAIN1_%04i%s" % (p, i, suffix)
            nfp = "%s_DOMAIN1_%04i%s" % (p, i + 1, suffix)
            if os.path.isfile(nfp):
                print nfp
                sys.exit()
            print ofp, nfp
            os.rename(ofp, nfp)
