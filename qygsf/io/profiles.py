import xml.dom
from enum import Enum, auto
from math import degrees, asin, cos, radians, ceil
from typing import Tuple, Union, List

import numpy as np

from ..profiles.geoprofiles import GeoProfile_
from ..profiles.sets import ProfileElevations
from ..geometries.shapes.space2d import Line2D
from ..geometries.shapes.space3d import Point3D, Line3D
from ..geometries.shapes.space4d import Line4D
from ..utils.qgis_utils.rasters import get_min_dem_resolution
from ..utils.time import standard_gpstime_to_datetime
from ..utils.qgis_utils.project import projectCrs
from ..utils.qgis_utils.rasters import interpolate_z
from ..utils.qgis_utils.points import TrackPointGPX


def topoline_from_dem(
    resampled_trace2d,
    dem,
    dem_params
):

    project_crs = projectCrs()
    if dem.crs() != projectCrs():
        trace2d_in_dem_crs = resampled_trace2d.crs_project(
            project_crs,
            dem.crs()
        )
    else:
        trace2d_in_dem_crs = resampled_trace2d

    ln3dtProfile = Line3D()
    for trace_pt2d_dem_crs, trace_pt2d_project_crs in zip(trace2d_in_dem_crs.pts, resampled_trace2d.pts):
        fInterpolatedZVal = interpolate_z(dem, dem_params, trace_pt2d_dem_crs)
        pt3dtPoint = Point3D(trace_pt2d_project_crs.x,
                             trace_pt2d_project_crs.y,
                             fInterpolatedZVal)
        ln3dtProfile.add_pt(pt3dtPoint)

    return ln3dtProfile


def topoprofiles_from_gpxfile(
        source_gpx_path: str,
        invert_profile: bool
) -> ProfileElevations:

    doc = xml.dom.minidom.parse(source_gpx_path)

    # define track name
    try:
        trkname = doc.getElementsByTagName('trk')[0].getElementsByTagName('name')[0].firstChild.data
    except:
        trkname = ''

    # get raw track point values (lat, lon, elev, time)
    track_raw_data = []
    for trk_node in doc.getElementsByTagName('trk'):
        for trksegment in trk_node.getElementsByTagName('trkseg'):
            for tkr_pt in trksegment.getElementsByTagName('trkpt'):
                track_raw_data.append((tkr_pt.getAttribute("lat"),
                                       tkr_pt.getAttribute("lon"),
                                       tkr_pt.getElementsByTagName("ele")[0].childNodes[0].data,
                                       tkr_pt.getElementsByTagName("time")[0].childNodes[0].data))

    # reverse profile orientation if requested
    if invert_profile:
        track_data = track_raw_data[::-1]
    else:
        track_data = track_raw_data

    # create list of TrackPointGPX elements
    track_points = []
    for val in track_data:
        gpx_trackpoint = TrackPointGPX(*val)
        track_points.append(gpx_trackpoint)

    # check for the presence of track points
    if len(track_points) == 0:
        raise Exception("No track point found in this file")

    # calculate delta elevations between consecutive points
    delta_elev_values = [np.nan]
    for ndx in range(1, len(track_points)):
        delta_elev_values.append(track_points[ndx].elev - track_points[ndx - 1].elev)

    # convert original values into ECEF values (x, y, z in ECEF global coordinate system)
    trk_ECEFpoints = [trck.as_pt3dt() for trck in track_points]

    # calculate 3D distances between consecutive points
    dist_3D_values = [np.nan]
    for ndx in range(1, len(trk_ECEFpoints)):
        dist_3D_values.append(trk_ECEFpoints[ndx].distance(trk_ECEFpoints[ndx - 1]))

    # calculate slope along track
    dir_slopes = []
    for delta_elev, dist_3D in zip(delta_elev_values, dist_3D_values):
        try:
            slope = degrees(asin(delta_elev / dist_3D))
        except:
            slope = 0.0
        dir_slopes.append(slope)

    # calculate horizontal distance along track
    horiz_dist_values = []
    for slope, dist_3D in zip(dir_slopes, dist_3D_values):
        try:
            horiz_dist_values.append(dist_3D * cos(radians(slope)))
        except:
            horiz_dist_values.append(np.nan)

    # defines the cumulative 2D distance values
    cum_distances_2D = [0.0]
    for ndx in range(1, len(horiz_dist_values)):
        cum_distances_2D.append(cum_distances_2D[-1] + horiz_dist_values[ndx])

    # defines the cumulative 3D distance values
    cum_distances_3D = [0.0]
    for ndx in range(1, len(dist_3D_values)):
        cum_distances_3D.append(cum_distances_3D[-1] + dist_3D_values[ndx])

    lat_values = [track.lat for track in track_points]
    lon_values = [track.lon for track in track_points]
    time_values = [track.time for track in track_points]
    elevations = [track.elev for track in track_points]

    topo_profiles = ProfileElevations()

    #topo_profiles.line_source = gpx_source
    topo_profiles.inverted = invert_profile

    topo_profiles.lons = np.asarray(lon_values)
    topo_profiles.lats = np.asarray(lat_values)
    topo_profiles.times = time_values
    topo_profiles.surface_names = [trkname]  # [] required for compatibility with DEM case
    topo_profiles.incr_len_2d = np.asarray(cum_distances_2D)
    topo_profiles.incr_len_3d = [np.asarray(cum_distances_3D)]  # [] required for compatibility with DEM case
    topo_profiles.z_array = [np.asarray(elevations)]  # [] required for compatibility with DEM case
    topo_profiles.dir_slopes = [np.asarray(dir_slopes)]  # [] required for compatibility with DEM case

    return topo_profiles


def try_extract_track_from_gpxfile(
        source_gpx_path: str,
        invert_profile: bool
) -> Tuple[bool, Union[str, Tuple[str, Line4D]]]:

    try:

        doc = xml.dom.minidom.parse(source_gpx_path)

        # define track name
        try:
            profile_name = doc.getElementsByTagName('trk')[0].getElementsByTagName('name')[0].firstChild.data
        except:
            profile_name = ''

        # get raw track point values (lat, lon, elev, time)
        track_raw_data = []
        for trk_node in doc.getElementsByTagName('trk'):
            for trksegment in trk_node.getElementsByTagName('trkseg'):
                for tkr_pt in trksegment.getElementsByTagName('trkpt'):
                    track_raw_data.append((tkr_pt.getAttribute("lat"),
                                           tkr_pt.getAttribute("lon"),
                                           tkr_pt.getElementsByTagName("ele")[0].childNodes[0].data,
                                           tkr_pt.getElementsByTagName("time")[0].childNodes[0].data))

        # create list of TrackPointGPX elements
        track_points = []
        for val in track_raw_data:
            lat, lon, ele, time = val
            time = standard_gpstime_to_datetime(time)
            gpx_trackpoint = TrackPointGPX(
                lat,
                lon,
                ele,
                time
            )
            track_points.append(gpx_trackpoint)

        # check for the presence of track points
        if len(track_points) == 0:
            return False, "No track point found in this file"

        # project track points to QGIS project CRS

        projected_pts = []
        for track_pt in track_points:

            projected_pt = track_pt.project(
                    dest_crs=projectCrs()
            )
            projected_pts.append(projected_pt)

        profile_line = Line4D(pts=projected_pts)

        if invert_profile:

            profile_line = profile_line.invert_direction()

        return True, (profile_name, profile_line)

    except Exception as e:

        return False, str(e)


class TrackSource(Enum):
    """
    The profile source type.
    """

    UNDEFINED  = auto()
    LINE_LAYER = auto()
    DIGITATION = auto()
    POINT_LIST = auto()
    GPX_FILE   = auto()


class GPXElevationUsage(Enum):
    """
    The profile source type.
    """

    NOT_USED = auto()
    USE_WITH_DEMS = auto()
    ONLY_TO_USE = auto()


def try_prepare_single_topo_profiles(
    profile_line: Union[Line2D, Line3D],
    track_source: TrackSource,
    gpx_elevation_usage: GPXElevationUsage,
    selected_dems: List,
    selected_dem_parameters: List,
    gpx_track_name: str
    ) -> Tuple[bool, Union[str, GeoProfile_]]:

    pt_num_threshold = 1e4

    try:

        # get DEMs resolutions in project CRS and choose the min value

        if track_source != TrackSource.GPX_FILE or \
           gpx_elevation_usage != GPXElevationUsage.ONLY_TO_USE:

            sample_distance = get_min_dem_resolution(
                selected_dems,
                selected_dem_parameters
            )

            # check total number of points in line(s) to create

            estimated_total_num_pts = 0

            profile_length = profile_line.length_2d
            profile_num_pts = profile_length / sample_distance
            estimated_total_num_pts += profile_num_pts

            estimated_total_num_pts = int(ceil(estimated_total_num_pts))

            if estimated_total_num_pts > pt_num_threshold:
                return False, f"There are {estimated_total_num_pts} estimated points (limit is {pt_num_threshold}) in profile(s) to create.\nTry increasing sample distance value"

            dem_named_3dlines = topo_lines_from_dems(
                source_profile_line=profile_line,
                sample_distance=sample_distance,
                selected_dems=selected_dems,
                selected_dem_parameters=selected_dem_parameters
            )

            if dem_named_3dlines is None:
                return False, "Debug: profile not created"

        else:

            dem_named_3dlines = None

        # Manage GPX data

        if track_source == TrackSource.GPX_FILE and \
           gpx_elevation_usage != GPXElevationUsage.NOT_USED:

            gpx_named_3dlines = [(gpx_track_name, profile_line)]

        else:

            gpx_named_3dlines = None

        print(f"DEBUG: profile line: {profile_line}")
        print(f"DEBUG: profile line num points: {profile_line.num_pts}")

        if dem_named_3dlines is None and gpx_named_3dlines is None:
            named_3dlines = None
        elif dem_named_3dlines is None:
            named_3dlines = gpx_named_3dlines
        elif gpx_named_3dlines is None:
            named_3dlines = dem_named_3dlines
        else:
            named_3dlines = dem_named_3dlines + gpx_named_3dlines

        if named_3dlines is None:
            return False, "Unable to create profiles"

        geoprofile = GeoProfile_()
        geoprofile.original_line = profile_line
        geoprofile.set_topo_profiles(named_3dlines)

        return True, geoprofile

    except Exception as e:

        return False, str(e)


def topo_lines_from_dems(
        source_profile_line,
        sample_distance,
        selected_dems,
        selected_dem_parameters
) -> List[Tuple[str, Line3D]]:

    resampled_line = source_profile_line.densify_2d_line(sample_distance)  # line resampled by sample distance

    # calculate 3D profiles from DEMs

    lines3d = []

    for dem, dem_params in zip(selected_dems, selected_dem_parameters):

        line3d = topoline_from_dem(
            resampled_line,
            dem,
            dem_params
        )

        lines3d.append((dem.name(), line3d))

    """
    # setup topoprofiles properties

    topo_profiles = ProfileElevations()

    topo_profiles.x_array = np.asarray(resampled_line.x_array)
    topo_profiles.y_array = np.asarray(resampled_line.y_array)
    topo_profiles.surface_names = [dem.name() for dem in selected_dems]
    topo_profiles.incremental_length_2d = np.asarray(resampled_line.incr_len_2d)
    topo_profiles.incremental_length_3d = [np.asarray(cl3dt.incr_len_3d) for cl3dt in lines3d]
    topo_profiles.z_array = [cl3dt.z_array for cl3dt in lines3d]
    topo_profiles.dir_slopes = [np.asarray(cl3dt.dir_slopes) for cl3dt in lines3d]
    topo_profiles.dem_params = [DEMParams(dem, params) for (dem, params) in
                                zip(selected_dems, selected_dem_parameters)]
    """

    return lines3d