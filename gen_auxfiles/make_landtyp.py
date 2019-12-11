"""
Generate the table 4 variable landtype
"""
import netCDF4
import datetime
import numpy

fn = "landtyp_MM5I.nc"
nc = netCDF4.Dataset(fn, "w", format="NETCDF3_CLASSIC")
nc.Conventions = "CF-1.0"
nc.title = "ISU MM5 model output prepared for NARCCAP"
nc.history = "python netCDF routines, %s" % (
    datetime.datetime.now().strftime("%Y-%m-%d %H:%M"),
)
nc.contact1 = "Daryl Herzmann (akrherz@iastate.edu)"
nc.contact2 = "3015 Agronomy Hall, Iowa State Univ.,Ames, Iowa, USA"
nc.realization = "1"
nc.experiment_id = "all"
nc.table_id = "Table 4"
nc.project_id = "NARCCAP"
nc.source = "MM5(2002): atmosphere: MM5v3.6.3 non-hydrostatic; sst/sea ice: AMIPII; land: Noah;  Convection: Kain-Fritsch 2; Radiation: RRTM; PBL: MRF; Explicit Moisture: Reisner Mixed-Phase; Buffer: 15 point exponential; Horizontal Resolution: 50km; Vertical Levels: 24"
nc.institution = "ISU (Iowa State University, Ames, Iowa, USA)"

nc.createDimension("xc", 124)
nc.createDimension("yc", 99)

xc = nc.createVariable("xc", "d", ("xc",))
xc.long_name = "x-coordinate Cartesian system"
xc.standard_name = "projection_x_coordinate"
xc.axis = "X"
xc.units = "m"

yc = nc.createVariable("yc", "d", ("yc",))
yc.long_name = "y-coordinate Cartesian system"
yc.standard_name = "projection_y_coordinate"
yc.axis = "Y"
yc.units = "m"

lat = nc.createVariable("lat", "d", ("yc", "xc"))
lat.long_name = "latitude"
lat.standard_name = "latitude"
lat.units = "degrees_north"

lon = nc.createVariable("lon", "d", ("yc", "xc"))
lon.long_name = "longitude"
lon.standard_name = "longitude"
lon.units = "degrees_east"

p = nc.createVariable("Lambert_Conformal", "c", ())
p.grid_mapping_name = "lambert_conformal_conic"
p.false_easting = 3825000.0
p.false_northing = 3187500.0

v = nc.createVariable("landtyp", "i", ("yc", "xc"))
v.units = "1"
v.standard_name = "land_cover"
# v.missing_value = numpy.array(1e20, v.dtype)
v.flag_values = (
    "1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24"
)
v.flag_meanings = "Urban_and_Built-Up_Land Dryland_Cropland_and_Pasture Irrigated_Cropland_and_Pasture Mixed_Dryland/Irrigated_Cropland_and_Pasture Cropland/Grassland_Mosaic Cropland/Woodland_Mosaic Grassland Shrubland Mixed_Shrubland/Grassland Savanna Deciduous_Broadleaf_Forest Deciduous_Needleleaf_Forest Evergreen_Broadleaf_Forest Evergreen_Needleleaf_Forest Mixed_Forest Water_Bodies Herbaceous_Wetland Wooded_Wetland Barren_or_Sparsely_Vegetated Herbaceous_Tundra Wooded_Tundra Mixed_Tundra Bare_Ground_Tundra Snow_or_Ice"
v.coordinates = "lon lat"
v.grid_mapping = "Lambert_Conformal"

nc2 = netCDF4.Dataset("MMOUTP_DOMAIN1_1002.nc")

lat[:] = nc2.variables["latitcrs"][15:-15, 15:-15]
lon[:] = nc2.variables["longicrs"][15:-15, 15:-15] + 360.0
xc[:] = numpy.arange(15, 139) * nc2.variables["grid_ds"][:] * 1000.0
yc[:] = numpy.arange(15, 114) * nc2.variables["grid_ds"][:] * 1000.0
p.standard_parallel = [
    nc2.variables["stdlat_2"][:],
    nc2.variables["stdlat_1"][:],
]
p.longitude_of_central_meridian = nc2.variables["coarse_cenlon"][:]
p.latitude_of_projection_origin = nc2.variables["coarse_cenlat"][:]


data = nc2.variables["land_use"]
v[:] = data[15:-15, 15:-15]

nc2.close()
nc.close()
