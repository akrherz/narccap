import os
import glob

os.chdir("Run.NCEP")

for year in range(1979,2005):
  for mo in range(1,13):
    dir = '%s%02i' % (year, mo)
    files = glob.glob(dir+'/*OUT*')
    for file in files:
       tick = file.split("_")[-1]
       newfn = '%s_%04i' % (file.split("/")[-1].replace("_"+tick, ''), float(tick))
       #print("%s %s" % (file, newfn))
       os.system("mv %s %s" % (file, newfn))
