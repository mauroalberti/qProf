from __future__ import division

import numpy as np


class Profile_Elements(object):
    def __init__(self):

        self.profile_source_type = None
        self.source_profile_line2dt = None
        self.sample_distance = None  # max spacing along profile; float
        self.topoline_colors = None
        self.plot_params = None

        self.resamp_src_line = None
        self.topo_profiles = None
        self.plane_attitudes = []
        self.curves = []
        self.curves_ids = []
        self.intersection_pts = []
        self.intersection_lines = []

    def set_topo_profiles(self, topo_profiles):

        self.topo_profiles = topo_profiles

    def add_intersections_pts(self, intersection_list):

        self.intersection_pts += intersection_list

    def add_intersections_lines(self, formation_list, intersection_line3d_list, intersection_polygon_s_list2):

        self.intersection_lines = zip(formation_list, intersection_line3d_list, intersection_polygon_s_list2)

    def get_current_dem_names(self):

        return self.topo_profiles.names

    def max_s(self):
        return self.topo_profiles.max_s()

    def min_z_topo(self):
        return self.topo_profiles.min_z()

    def max_z_topo(self):
        return self.topo_profiles.max_z()

    def min_z_plane_attitudes(self):

        # TODO:  manage case for possible nan p_z values
        return min([plane_attitude.pt_3d.p_z for plane_attitude_set in self.plane_attitudes for plane_attitude in
                    plane_attitude_set if 0.0 <= plane_attitude.sign_hor_dist <= self.get_max_s()])

    def max_z_plane_attitudes(self):

        # TODO:  manage case for possible nan p_z values
        return max([plane_attitude.pt_3d.p_z for plane_attitude_set in self.plane_attitudes for plane_attitude in
                    plane_attitude_set if 0.0 <= plane_attitude.sign_hor_dist <= self.get_max_s()])

    def min_z_curves(self):

        return min([pt_2d.p_y for multiline_2d_list in self.curves for multiline_2d in multiline_2d_list for line_2d in
                    multiline_2d.lines for pt_2d in line_2d.pts if 0.0 <= pt_2d.p_x <= self.get_max_s()])

    def max_z_curves(self):

        return max([pt_2d.p_y for multiline_2d_list in self.curves for multiline_2d in multiline_2d_list for line_2d in
                    multiline_2d.lines for pt_2d in line_2d.pts if 0.0 <= pt_2d.p_x <= self.get_max_s()])

    def min_z(self):

        min_z = self.min_z_topo()

        if len(self.plane_attitudes) > 0:
            min_z = min([min_z, self.min_z_plane_attitudes()])

        if len(self.curves) > 0:
            min_z = min([min_z, self.min_z_curves()])

        return min_z

    def max_z(self):

        max_z = self.max_z_topo()

        if len(self.plane_attitudes) > 0:
            max_z = max([max_z, self.max_z_plane_attitudes()])

        if len(self.curves) > 0:
            max_z = max([max_z, self.max_z_curves()])

        return max_z

    def add_plane_attitudes(self, plane_attitudes):

        self.plane_attitudes.append(plane_attitudes)

    def add_curves(self, multiline_2d_list, ids_list):

        self.curves.append(multiline_2d_list)
        self.curves_ids.append(ids_list)

class DEMParams(object):

    def __init__(self, layer, params):
        self.layer = layer
        self.params = params

class TopoLine3D(object):

    def __init__(self, name, line3d):
        self.name = name
        self.profile_3d = line3d  # class CartesianLine3DT, a list of CartesianPoint3DT

    def min_z(self):
        return self.profile_3d.z_min()

    def max_z(self):
        return self.profile_3d.z_max()

    def mean_z(self):
        return self.profile_3d.z_mean()

    def var_z(self):
        return self.profile_3d.z_var()

    def std_z(self):
        return self.profile_3d.z_std()

    def x_list(self):
        return [pt_3d.p_x for pt_3d in self.profile_3d.pts]

    def y_list(self):
        return [pt_3d.p_y for pt_3d in self.profile_3d.pts]

    def z_list(self):
        return [pt_3d.p_z for pt_3d in self.profile_3d.pts]

    def directional_slopes(self):
        return self.profile_3d.slopes_list()

    def length_2d(self):
        return self.profile_3d.length_2d()

    def get_increm_dist_3d(self):
        return self.profile_3d.incremental_length_3d()

    def get_increm_dist_2d(self):
        return self.profile_3d.incremental_length_2d()

class TopoProfiles(object):

    def __init__(self):

        self.xs = None
        self.ys = None
        self.lons = None
        self.lats = None
        self.times = None
        self.names = []
        self.s = None
        self.s3d = []
        self.elevs = []
        self.dir_slopes = []
        self.dem_params = []
        self.gpx_params = None
        self.colors = []
        self.statistics_defined = False
        self.profile_defined = False

    def max_s(self):
        return self.s[-1]

    def min_z(self):
        return min(map(np.nanmin, self.elevs))

    def max_z(self):
        return max(map(np.nanmax, self.elevs))

    @property
    def absolute_slopes(self):
        return map(np.fabs, self.dir_slopes)

class PlaneAttitude(object):
    def __init__(self, rec_id, source_point_3d, source_geol_plane, point_3d, slope_rad, dwnwrd_sense, sign_hor_dist):
        self.id = rec_id
        self.src_pt_3d = source_point_3d
        self.src_geol_plane = source_geol_plane
        self.pt_3d = point_3d
        self.slope_rad = slope_rad
        self.dwnwrd_sense = dwnwrd_sense
        self.sign_hor_dist = sign_hor_dist
