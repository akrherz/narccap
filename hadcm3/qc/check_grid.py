import pyproj
import netCDF4

lcc = pyproj.Proj(
    proj="lcc", lat_0=47.5, lat_1=60, lat_2=30, lon_0=-97, ellps="WGS84"
)

nc = netCDF4.Dataset("../postprocess/data/NCOUT_DOMAIN1_0001.nc", "r")
lats = nc.variables["latitcrs"][:]
lons = nc.variables["longicrs"][:]

(lat0, lon0) = lats[0, 0], lons[0, 0]

print "Cross", lcc(lon0, lat0)

lats = nc.variables["latitdot"][:]
lons = nc.variables["longidot"][:]

(lat0, lon0) = lats[0, 0], lons[0, 0]

print "Dot", lcc(lon0, lat0)

nc.close()
