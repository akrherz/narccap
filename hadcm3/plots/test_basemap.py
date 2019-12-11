import iemplot
import netCDF4
import numpy
import numpy.ma
from mpl_toolkits.basemap import Basemap, cm
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.colors import ListedColormap, BoundaryNorm
import Image
import mx.DateTime

ncfp = "soilt4_MM5I_hadcm3_2066010103.nc"
nc = netCDF4.Dataset("/tmp/" + ncfp, "r")
ncfp = "soilt4_MM5I_hadcm3_1996010103.nc"
nc2 = netCDF4.Dataset("/tmp/" + ncfp, "r")

# snd = nc.variables['snd']
data = nc.variables["soilt4"][:, :, :] - 273.15
data2 = nc2.variables["soilt4"][:, :, :] - 273.15
data = numpy.average(data, 0) - numpy.average(data2, 0)
lats = nc.variables["lat"][:]
lons = nc.variables["lon"][:]
xc = nc.variables["xc"][:]
yc = nc.variables["yc"][:]


# data = numpy.ma.max( snd[:,:,:], 0)
# data.mask = numpy.where( data < 0.001, True, False)

fig = plt.figure(num=None, figsize=(10.24, 7.68))
ax = plt.axes([0.01, 0.1, 0.9, 0.8], axisbg=(0.4471, 0.6235, 0.8117))

lcc = nc.variables["Lambert_Conformal"]
# 		Lambert_Conformal:false_easting = 3825000. ;
# 		Lambert_Conformal:false_northing = 3187500. ;

# proj4string = "+lon_0=-97.0 +y_0=3187500.0 +R=6371200.0 +proj=lcc +x_0=3825000.0 +units=m +lat_2=60.0 +lat_1=47.5 +lat_0=30.0"
map = Basemap(
    projection="lcc",
    # urcrnrlat=numpy.max(lats),
    # llcrnrlat=numpy.min(lats),
    # urcrnrlon=numpy.max(lons),
    # llcrnrlon=numpy.min(lons),
    height=numpy.max(yc) - numpy.min(yc),
    width=numpy.max(xc) - numpy.min(xc),
    lat_1=lcc.standard_parallel[0],
    lat_2=lcc.standard_parallel[1],
    lat_0=lcc.latitude_of_projection_origin,
    lon_0=lcc.longitude_of_central_meridian,
    rsphere=6371200.0,
    resolution="l",
    area_thresh=10000,
    ax=ax,
)

map.fillcontinents(color="0.7", zorder=0)

# values = [0., 0.001, 0.005, 0.01, 0.02, 0.05, 0.10, 0.20, 0.25,
#          0.5, 0.75, 0.90, 1., 1.5, 2., 5., 10., 20., 50., 100.]

maxV = numpy.max(data) + 2
minV = numpy.min(data) - 2
# values = numpy.arange(minV, maxV, (maxV-minV)/20.0)
values = numpy.arange(-5.0, 5.0, 0.5)

colors = [
    (0, 0, 0),
    (24, 24, 112),
    (16, 78, 139),
    (23, 116, 205),
    (72, 118, 255),
    (91, 172, 237),
    (173, 215, 230),
    (209, 237, 237),
    (229, 239, 249),
    (242, 255, 255),
    (253, 245, 230),
    (255, 228, 180),
    (243, 164, 96),
    (237, 118, 0),
    (205, 102, 29),
    (224, 49, 15),
    (237, 0, 0),
    (205, 0, 0),
    (139, 0, 0),
]
colors = numpy.array(colors) / 255.0

# colors = ["#FFFFFF", "#7FFF00", "#00CD00", "#008B00", "#104E8B", "#1E90FF",
#             "#00B2EE", "#00EEEE", "#8968CD", "#912CEE", "#8B008B", "#8B0000",
#             "#CD0000", "#EE4000", "#FF7F00", "#CD8500", "#FFD700", "#EEEE00",
#             "#FFFF00", "#7FFF00","#000000"]
cmap = ListedColormap(colors[:])
cmap.set_over(color=colors[-1])
cmap.set_under(color=colors[0])
normer = BoundaryNorm(values, cmap.N, clip=False)


# data_proj,x,y = map.transform_scalar(data,lons,lats,nx,ny,returnxy=True)
cax = plt.axes([0.902, 0.2, 0.03, 0.6])
res = map.imshow(data, interpolation="nearest", cmap=cmap, norm=normer, ax=ax)
clr = plt.colorbar(res, cax=cax, format="%g")
clr.set_ticks(values)
map.drawcoastlines(ax=ax)
map.drawstates(ax=ax)
map.drawcountries(ax=ax)

logo = Image.open("NARCCAP.png")
ax3 = plt.axes(
    [0.05, 0.89, 0.2, 0.1],
    frameon=False,
    axisbg=(0.4471, 0.6235, 0.8117),
    yticks=[],
    xticks=[],
)
ax3.imshow(logo, origin="lower")

cax.text(
    -0.5,
    0.5,
    "Temperature Difference [C]",
    transform=cax.transAxes,
    size="medium",
    color="k",
    horizontalalignment="center",
    verticalalignment="center",
    rotation="vertical",
)

# ax.text(0.27, 1.05, 'File: %s' % (ncfp,), transform=ax.transAxes,
#     size=13, horizontalalignment='left', verticalalignment='center')
# ax.text(0.27, 1.02, 'Experiment: %s' % (nc.experiment_id,), transform=ax.transAxes,
#     size=13, horizontalalignment='left', verticalalignment='center')

ax.text(
    0.27,
    1.08,
    "Soil Layer 4 Temperature [soil_t_4]\nScenario Avg (2066-2070) minus Contempory Avg (1996-2000)",
    transform=ax.transAxes,
    size=15,
    horizontalalignment="left",
    verticalalignment="center",
)

ax.text(
    0.01,
    -0.03,
    "Plot Generated: %s" % (mx.DateTime.gmt().strftime("%d %b %Y %H%MZ"),),
    transform=ax.transAxes,
    size="small",
    horizontalalignment="left",
    verticalalignment="bottom",
)

# x,y = map(lons,lats)
# clevs = [0,0.01,0.05,0.1,1,2,10]
# cs = map.contourf(x,y,data, clevs)
fig.savefig("test.png")
