
from pprint import pprint

import geopandas

from qygsf.profiles.geoprofiles import GeoProfileSet
from qygsf.georeferenced.rasters import *
from qygsf.io.plot import plot_grid
from qygsf.io.gdal.vector import read_linestring_geometries
from qygsf.io.plot import plot_line
from qygsf.profiles.geoprofiles import GeoProfile
from qygsf.profiles.profilers import SegmentProfiler
from qygsf.profiles.plot import plot
from qygsf.io.geology import try_extract_flat_georeferenced_attitudes
from qygsf.profiles.profilers import georef_attitudes_3d_from_grid
from qygsf.io.gdal.raster import try_read_raster_band
from qygsf.profiles.profilers import ParallelProfilers
from qygsf.io.gdal.vector import reading_line_shapefile


# Geologic profiles with qygsf

## 1 - Mount Alpi (Lucania, Southern Italy)

### DEM input

source_data = "/home/mauro/Documents/projects/gsf/example_data/mt_alpi/malpi_aster_w4u3.tif"

success, result = try_read_raster_band(raster_source=source_data)

print(success)

geotransform, projection, band_params, data = result

for info in (geotransform, projection, band_params, data):
    print(f"info type: {type(info)}:\n\n{info}\n\n")

geo_array = GeoArray(
    inGeotransform=geotransform,
    epsg_code=32633,
    inLevels=[data]
)

fig = plot_grid(geo_array)

### Source profile

src_profile_shapefile_pth = "/home/mauro/Documents/projects/gsf/example_data/mt_alpi/profile.shp"
profiles = read_linestring_geometries(src_profile_shapefile_pth)

print(profiles)

line = profiles.line()
plot_line(fig, line)

### Initializing a geoprofile

geoprofile = GeoProfile()

### Creating a linear profiler

profiler = SegmentProfiler(
    segment2d=Segment2D(
        line.start_pt(),
        line.end_pt()
    ),
    densify_distance=5,  # meters sampling distance along profile
    epsg_cd=32633
)

profiler.num_pts()

### Adding a topographic profile

topo_profile = profiler.profile_grid(geo_array)

print(topo_profile)
print(f"Type of 'topo_profile': {type(topo_profile)}")

geoprofile.topo_profiles = topo_profile

fig = plot(geoprofile)
fig.show()

### Plotting geological attitudes

attitudes_shape = "/home/mauro/Documents/projects/gsf/example_data/mt_alpi/attitudes.shp"

attitudes = geopandas.read_file(attitudes_shape)

print(attitudes)

print(attitudes.crs)

ax = attitudes.plot()



success, result = try_extract_flat_georeferenced_attitudes(
    geodataframe=attitudes,
    azim_fldnm="dip_dir",
    dip_ang_fldnm="dip_ang"
)

print(success)

if not success:

    msg = result
    print(msg)

else:

    attitudes = result

print("Number of attitudes: {}".format(len(attitudes)))

pprint(attitudes)


attitudes_3d = georef_attitudes_3d_from_grid(
    structural_data=attitudes,
    height_source=geo_array,
)

for att3d in attitudes_3d:
    print(att3d)

mapping_method = {}
mapping_method['method'] = 'nearest'

projected_attitudes = profiler.map_georef_attitudes_to_section(
    attitudes_3d=attitudes_3d,
    mapping_method=mapping_method)

for att in projected_attitudes:
    print("\n" + str(att))

geoprofile.profile_attitudes = projected_attitudes

fig = plot(geoprofile)

## 2 - Timpa San Lorenzo area (Calabria, Southern Italy) (in preparation)

src_dem = "/home/mauro/Documents/projects/gsf/example_data/timpa_san_lorenzo/tsl_tinitaly_w84u32.tif"
src_profile = "/home/mauro/Documents/projects/gsf/example_data/timpa_san_lorenzo/profile.shp"


success, result = try_read_raster_band(raster_source=src_dem)

print(success, result)

geotransform, projection, band_params, data = result

geoarray = GeoArray(
    inGeotransform=geotransform,
    epsg_code=32633,
    inLevels=[data]
)

fig = plot_grid(geoarray)
fig.show()

profiles = read_linestring_geometries(src_profile)
line = profiles.line()
plot_line(fig, line)


geoprofiles = GeoProfileSet()

base_profiler = SegmentProfiler(
    start_pt=line.start_pt(),
    end_pt=line.end_pt(),
    densify_distance=5,
    epsg_cd=32633
)

base_profiler

base_profiler.num_pts()


multiple_profilers = ParallelProfilers.fromBaseProfiler(
    base_profiler=base_profiler,
    profs_num=15,
    profs_offset=500,
    profs_arr="central"
)

multiple_profilers

topo_profiles = multiple_profilers.profile_grid(geoarray)

print(type(topo_profiles))

geoprofiles.topo_profiles_set = topo_profiles

fig = plot(geoprofiles)

### Line intersections


faults_shape = "/home/mauro/Documents/projects/gsf/example_data/timpa_san_lorenzo/faults.shp"
faults = geopandas.read_file(faults_shape)
faults

type(faults)

type(faults['geometry'][0])


success, answer = reading_line_shapefile(
    shp_path=faults_shape,
    flds=['fid', 'name'],
    read_z=False)

print(success)

results = answer

pprint(results)

line_types = set(map(lambda val: type(val[0]), results))
#print(line_types)

tsl_faults = [rec[0] for rec in results if rec[1][1] == 'Timpa San Lorenzo thrust flat']

pprint(tsl_faults)

## Test of intersection with single line profile

base_profiler

topo_profile = base_profiler.profile_grid(geoarray)

print(topo_profile)

#sys.exit()

geoprofile = GeoProfile()

geoprofile.topo_profiles = topo_profile
fig = plot(geoprofile)
fig.show()

# pprint(tsl_faults)

line_intersections = base_profiler.intersect_lines(tsl_faults)

pprint(line_intersections)
print(f"type {type(line_intersections)}")
print("Printed line intersections")

profile_intersections = base_profiler.parse_profile_intersections(line_intersections)

pprint(profile_intersections)
print(f"type {type(profile_intersections)}")
print("Printed profile_intersections")

geoprofile.lines_intersections = profile_intersections

fig = plot(geoprofile)
fig.show()

print("Completed successfully")
