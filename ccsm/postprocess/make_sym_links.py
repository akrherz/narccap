import os
import glob

os.chdir("Run.scenario")

for year in range(2038,2071):
  for mo in range(1,13):
    dir = '/mnt/osgood/narccap/exp2a.2/output/%s%02i' % (year, mo)
    files = glob.glob(dir+'/*OUT*')
    for file in files:
       tick = file.split("_")[-1]
       newfn = '%s_%04i' % (file.split("/")[-1].replace("_"+tick, ''), float(tick))
       os.system("ln -s %s %s" % (file, newfn))
