import xml.dom.minidom

from qygsf.utils.qgis_utils.lines import *
from qygsf.geometries.shapes.features import *
from qygsf.georeferenced.geodetic import *
from qygsf.utils.qgis_utils.project import *
from qygsf.utils.qgis_utils.rasters import *


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


class GeoProfilesSet(object):
    """
    Represents a set of ProfileElements instances,
    stored as a list
    """

    def __init__(self, name=""):

        self._name = name
        self._geoprofiles = []  # list of GeoProfile instances
        self.profiles_created = False
        self.plot_params = None

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, new_name):

        self._name = new_name

    @property
    def geoprofiles(self):

        return self._geoprofiles

    @property
    def geoprofiles_num(self):

        return len(self._geoprofiles)

    def geoprofile(self, ndx):

        return self._geoprofiles[ndx]

    def append(self, geoprofile):

        self._geoprofiles.append(geoprofile)

    def insert(self, ndx, geoprofile):

        self._geoprofiles.insert(ndx, geoprofile)

    def move(self, ndx_init, ndx_final):

        geoprofile = self._geoprofiles.pop(ndx_init)
        self.insert(ndx_final, geoprofile)

    def move_up(self, ndx):

        self.move(ndx, ndx -1)

    def move_down(self, ndx):

        self.move(ndx, ndx + 1)

    def remove(self, ndx):

        _ = self._geoprofiles.pop(ndx)


class GeoProfile:
    """
    Class representing the topographic and geological elements
    embodying a single geological profile.
    """

    def __init__(self):

        """
        self.source_data_type = None
        self.original_line = None
        self.sample_distance = None  # max spacing along profile; float
        self.resampled_line = None
        """

        self.named_lines = []  # list of name and Line
        self.geoplane_attitudes = []
        self.geosurfaces = []
        self.geosurfaces_ids = []
        self.lineaments = []
        self.outcrops = []

    def set_topo_profiles(self,
        topo_profiles
                          ):

        self.named_lines = topo_profiles

    def add_topo_profiles(self, topo_profiles):

        self.named_lines += topo_profiles

    def add_intersections_pts(self, intersection_list):

        self.lineaments += intersection_list

    def add_intersections_lines(self,
        formation_list,
        intersection_line3d_list,
        intersection_polygon_s_list2
    ):

        self.outcrops = list(
            zip(
                formation_list,
                intersection_line3d_list,
                intersection_polygon_s_list2
            )
        )

    """
    def get_current_dem_names(self):

        return self.topo_profiles.name
    """

    def s_max(self):

        return max([line.length_2d for _, line in self.named_lines])

    def min_z_topo(self):

        return min([line.z_min for _, line in self.named_lines])

    def max_z_topo(self):

        return max([line.z_max for _, line in self.named_lines])

    def min_z_plane_attitudes(self):

        # TODO:  manage case for possible nan p_z values
        return min([plane_attitude.pt_3d.p_z for plane_attitude_set in self.geoplane_attitudes for plane_attitude in
                    plane_attitude_set if 0.0 <= plane_attitude.sign_hor_dist <= self.s_max()])

    def max_z_plane_attitudes(self):

        # TODO:  manage case for possible nan p_z values
        return max([plane_attitude.pt_3d.p_z for plane_attitude_set in self.geoplane_attitudes for plane_attitude in
                    plane_attitude_set if 0.0 <= plane_attitude.sign_hor_dist <= self.s_max()])

    def min_z_curves(self):

        return min([pt_2d.p_y for multiline_2d_list in self.geosurfaces for multiline_2d in multiline_2d_list for line_2d in
                    multiline_2d.lines for pt_2d in line_2d.pts if 0.0 <= pt_2d.p_x <= self.s_max()])

    def max_z_curves(self):

        return max([pt_2d.p_y for multiline_2d_list in self.geosurfaces for multiline_2d in multiline_2d_list for line_2d in
                    multiline_2d.lines for pt_2d in line_2d.pts if 0.0 <= pt_2d.p_x <= self.s_max()])

    def z_min(self):

        min_z = self.min_z_topo()

        if len(self.geoplane_attitudes) > 0:
            min_z = min([min_z, self.min_z_plane_attitudes()])

        if len(self.geosurfaces) > 0:
            min_z = min([min_z, self.min_z_curves()])

        return min_z

    def z_max(self):

        max_z = self.max_z_topo()

        if len(self.geoplane_attitudes) > 0:
            max_z = max([max_z, self.max_z_plane_attitudes()])

        if len(self.geosurfaces) > 0:
            max_z = max([max_z, self.max_z_curves()])

        return max_z

    def add_plane_attitudes(self, plane_attitudes):

        self.geoplane_attitudes.append(plane_attitudes)

    def add_curves(self, lMultilines, lIds):

        self.geosurfaces.append(lMultilines)
        self.geosurfaces_ids.append(lIds)


class ProfileElevations:

    def __init__(self):

        #self.line_source = None
        self.dem_params = []
        self.gpx_params = None

        self.x_array = None
        self.y_array = None
        self.times = None
        self.incr_len_2d = None

        self.surface_names = []

        self.incr_len_3d = []
        self.z_array = []
        self.dir_slopes = []

        self.inverted = None

        self.statistics_calculated = False
        self.profile_created = False

    def max_s(self):

        return self.incr_len_2d[-1]

    def min_z(self):

        return min(list(map(np.nanmin, self.z_array)))

    def max_z(self):

        return max(list(map(np.nanmax, self.z_array)))

    @property
    def absolute_slopes(self):

        return list(map(np.fabs, self.dir_slopes))

"""
class DEMParams(object):

    def __init__(self, layer, params):

        self.layer = layer
        self.params = params
"""

class PlaneAttitude(object):

    def __init__(self, rec_id, source_point_3d, source_geol_plane, point_3d, slope_rad, dwnwrd_sense, sign_hor_dist):

        self.id = rec_id
        self.src_pt_3d = source_point_3d
        self.src_geol_plane = source_geol_plane
        self.pt_3d = point_3d
        self.slope_rad = slope_rad
        self.dwnwrd_sense = dwnwrd_sense
        self.sign_hor_dist = sign_hor_dist


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

    ln3dtProfile = Line()
    for trace_pt2d_dem_crs, trace_pt2d_project_crs in zip(trace2d_in_dem_crs.pts, resampled_trace2d.pts):
        fInterpolatedZVal = interpolate_z(dem, dem_params, trace_pt2d_dem_crs)
        pt3dtPoint = Point(trace_pt2d_project_crs.x,
                           trace_pt2d_project_crs.y,
                           fInterpolatedZVal)
        ln3dtProfile.add_pt(pt3dtPoint)

    return ln3dtProfile


def topo_lines_from_dems(
        source_profile_line,
        sample_distance,
        selected_dems,
        selected_dem_parameters
) -> List[Tuple[str, Line]]:

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
        dist_3D_values.append(trk_ECEFpoints[ndx].dist_3d(trk_ECEFpoints[ndx - 1]))

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
) -> Tuple[bool, Union[str, Tuple[str, Line]]]:

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

        profile_line = Line(pts=projected_pts)

        if invert_profile:

            profile_line = profile_line.invert_direction()

        return True, (profile_name, profile_line)

    except Exception as e:

        return False, str(e)


def get_min_dem_resolution(
    selected_dems: List,
    selected_dem_parameters: List
) -> numbers.Real:

    dem_resolutions_prj_crs_list = []

    for dem, dem_params in zip(selected_dems, selected_dem_parameters):
        dem_resolutions_prj_crs_list.append(
            get_dem_resolution_in_prj_crs(
                dem,
                dem_params,
                projectCrs())
        )

    min_dem_resolution = min(dem_resolutions_prj_crs_list)

    if min_dem_resolution > 1:
        sample_distance = round(min_dem_resolution)
    else:
        sample_distance = min_dem_resolution

    return sample_distance


def try_prepare_single_topo_profiles(
    profile_line: Line,
    track_source: TrackSource,
    gpx_elevation_usage: GPXElevationUsage,
    selected_dems: List,
    selected_dem_parameters: List,
    gpx_track_name: str
    ) -> Tuple[bool, Union[str, GeoProfile]]:

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

        geoprofile = GeoProfile()
        geoprofile.original_line = profile_line
        geoprofile.set_topo_profiles(named_3dlines)

        return True, geoprofile

    except Exception as e:

        return False, str(e)


def intersect_with_dem(
        demLayer,
        demParams,
        project_crs,
        lIntersPts
):
    """
    
    :param demLayer: 
    :param demParams:
    :param project_crs: 
    :param lIntersPts: 
    :return: a list of Point instances
    """

    # project to Dem CRS
    if demParams.crs != project_crs:
        lQgsPoints = [QgsPointXY(pt.x, pt.y) for pt in lIntersPts]
        lDemCrsIntersQgsPoints = [project_qgs_point(qgsPt, project_crs, demParams.crs) for qgsPt in
                                               lQgsPoints]
        lDemCrsIntersPts = [Point(qgispt.x(), qgispt.y()) for qgispt in lDemCrsIntersQgsPoints]
    else:
        lDemCrsIntersPts = lIntersPts

    # interpolate z values from Dem
    lZVals = [interpolate_z(demLayer, demParams, pt) for pt in lDemCrsIntersPts]

    lXYZVals = [(pt2d.x, pt2d.y, z) for pt2d, z in zip(lIntersPts, lZVals)]

    return [Point(x, y, z) for x, y, z in lXYZVals]


def calculate_profile_lines_intersection(
        multilines2d_list,
        id_list,
        profile_line2d
):

    profile_segment2d_list = profile_line2d.as_segments()

    profile_segment2d = profile_segment2d_list[0]

    intersection_list = []
    for ndx, multiline2d in enumerate(multilines2d_list):
        if id_list is None:
            multiline_id = ''
        else:
            multiline_id = id_list[ndx]
        for line2d in multiline2d.lines:
            for line_segment2d in line2d.as_segments():
                try:
                    intersection_point2d = profile_segment2d.intersection_2d_pt(line_segment2d)
                except ZeroDivisionError:
                    continue
                if intersection_point2d is None:
                    continue
                if line_segment2d.contains_2d_pt(intersection_point2d) and \
                   profile_segment2d.contains_2d_pt(intersection_point2d):
                    intersection_list.append([intersection_point2d, multiline_id])

    return intersection_list


def intersection_distances_by_profile_start_list(profile_line, intersections):

    # convert the profile line
    # from a CartesianLine2DT to a CartesianSegment2DT
    profile_segment2d_list = profile_line.as_segments()
    # debug
    assert len(profile_segment2d_list) == 1
    profile_segment2d = profile_segment2d_list[0]

    # determine distances for each point in intersection list
    # creating a list of float values
    distance_from_profile_start_list = []
    for intersection in intersections:
        distance_from_profile_start_list.append(profile_segment2d.start_pt.dist_2d(intersection[0]))

    return distance_from_profile_start_list


def calculate_pts_in_projection(pts_in_orig_crs, srcCrs, destCrs):

    pts_in_prj_crs = []
    for pt in pts_in_orig_crs:
        qgs_pt = QgsPointXY(pt.x, pt.y)
        qgs_pt_prj_crs = project_qgs_point(qgs_pt, srcCrs, destCrs)
        pts_in_prj_crs.append(Point(qgs_pt_prj_crs.x(), qgs_pt_prj_crs.y()))
    return pts_in_prj_crs


def profile_polygon_intersection(profile_qgsgeometry, polygon_layer, inters_polygon_classifaction_field_ndx):

    intersection_polyline_polygon_crs_list = []

    if polygon_layer.selectedFeatureCount() > 0:
        features = polygon_layer.selectedFeatures()
    else:
        features = polygon_layer.getFeatures()

    for polygon_feature in features:
        # retrieve every (selected) feature with its geometry and attributes

        # fetch geometry
        poly_geom = polygon_feature.geometry()

        intersection_qgsgeometry = poly_geom.intersection(profile_qgsgeometry)

        try:
            if intersection_qgsgeometry.isEmpty():
                continue
        except:
            try:
                if intersection_qgsgeometry.isGeosEmpty():
                    continue
            except:
                return False, "Missing function for checking empty geometries.\nPlease upgrade QGIS"

        if inters_polygon_classifaction_field_ndx >= 0:
            attrs = polygon_feature.attributes()
            polygon_classification = attrs[inters_polygon_classifaction_field_ndx]
        else:
            polygon_classification = None

        if intersection_qgsgeometry.isMultipart():
            lines = intersection_qgsgeometry.asMultiPolyline()
        else:
            lines = [intersection_qgsgeometry.asPolyline()]

        for line in lines:
            intersection_polyline_polygon_crs_list.append([polygon_classification, line])

    return True, intersection_polyline_polygon_crs_list


def extract_multiline2d_list(
        structural_line_layer,
        project_crs
):

    line_orig_crs_geoms_attrs = line_geoms_attrs(structural_line_layer)

    line_orig_geom_list3 = [geom_data[0] for geom_data in line_orig_crs_geoms_attrs]
    line_orig_crs_MultiLine2D_list = [xytuple_l2_to_MultiLine(xy_list2) for xy_list2 in line_orig_geom_list3]
    line_orig_crs_clean_MultiLine2D_list = [multiline_2d.remove_coincident_points() for multiline_2d in
                                            line_orig_crs_MultiLine2D_list]

    # get CRS information
    structural_line_layer_crs = structural_line_layer.crs()

    # project input line layer to project CRS
    line_proj_crs_MultiLine2D_list = [multiline2d.crs_project(structural_line_layer_crs, project_crs) for
                                          multiline2d in line_orig_crs_clean_MultiLine2D_list]

    return line_proj_crs_MultiLine2D_list


def define_plot_structural_segment(structural_attitude, profile_length, vertical_exaggeration, segment_scale_factor=70.0):

    ve = float(vertical_exaggeration)
    intersection_point = structural_attitude.pt_3d
    z0 = intersection_point.z

    h_dist = structural_attitude.sign_hor_dist
    slope_rad = structural_attitude.slope_rad
    intersection_downward_sense = structural_attitude.dwnwrd_sense
    length = profile_length / segment_scale_factor

    s_slope = sin(float(slope_rad))
    c_slope = cos(float(slope_rad))

    if c_slope == 0.0:
        height_corr = length / ve
        structural_segment_s = [h_dist, h_dist]
        structural_segment_z = [z0 + height_corr, z0 - height_corr]
    else:
        t_slope = s_slope / c_slope
        width = length * c_slope

        length_exag = width * sqrt(1 + ve*ve * t_slope*t_slope)

        corr_width = width * length / length_exag
        corr_height = corr_width * t_slope

        structural_segment_s = [h_dist - corr_width, h_dist + corr_width]
        structural_segment_z = [z0 + corr_height, z0 - corr_height]

        if intersection_downward_sense == "left":
            structural_segment_z = [z0 - corr_height, z0 + corr_height]

    return structural_segment_s, structural_segment_z


def calculate_projected_3d_pts(
    struct_pts,
    structural_pts_crs,
    demObj
):

    demCrs = demObj.params.crs

    # check if on-the-fly-projection is set on
    project_crs = projectCrs()

    # set points in the project crs
    if structural_pts_crs != project_crs:
        struct_pts_in_prj_crs = calculate_pts_in_projection(struct_pts, structural_pts_crs, project_crs)
    else:
        struct_pts_in_prj_crs = copy.deepcopy(struct_pts)

        # project the source points from point layer crs to DEM crs
    # if the two crs are different
    if structural_pts_crs != demCrs:
        struct_pts_in_dem_crs = calculate_pts_in_projection(struct_pts, structural_pts_crs, demCrs)
    else:
        struct_pts_in_dem_crs = copy.deepcopy(struct_pts)

        # - 3D structural points, with x, y, and z extracted from the current DEM
    struct_pts_z = get_zs_from_dem(struct_pts_in_dem_crs, demObj)

    assert len(struct_pts_in_prj_crs) == len(struct_pts_z)

    return [Point(pt.x, pt.y, z) for (pt, z) in zip(struct_pts_in_prj_crs, struct_pts_z)]