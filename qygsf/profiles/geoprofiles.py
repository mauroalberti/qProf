
from .sets import *


z_padding = 0.2


class GeoProfile_:
    """
    Class representing the topographic and geological elements
    representing a single geological profile.
    """

    def __init__(self):

        """
        self.source_data_type = None
        self.original_line = None
        self.sample_distance = None  # max spacing along profile; float
        self.resampled_line = None
        """

        self.named_topoprofiles = []  # list of name and Line

        self.profile_attitudes = []
        self.projected_lines = []
        self.projected_lines_ids = []
        self.line_intersections = []
        self.polygons_intersections = []
        self.lines_intersections = None

    def clear_topo_profile(self):
        """

        :return:
        """

        self.named_topoprofiles = None

    def set_topo_profiles(self,
        topo_profiles
    ):

        self.named_topoprofiles = topo_profiles

    @property
    def topoprofiles(self):
        """

        :return:
        """

        return self.named_topoprofiles

    def add_topo_profiles(self, topo_profiles):

        self.named_topoprofiles += topo_profiles

    def add_intersections_pts(self, intersection_list):

        self.line_intersections += intersection_list

    def clear_lines_intersections(self):
        """
        Clear line intersections content.

        :return:
        """

        self.lines_intersections = None

    def add_intersections_lines(self,
        formation_list,
        intersection_line3d_list,
        intersection_polygon_s_list2
    ):

        self.polygons_intersections = list(
            zip(
                formation_list,
                intersection_line3d_list,
                intersection_polygon_s_list2
            )
        )

    def s_min(self):
        """

        :return:
        """

        return self.topoprofiles.s_min()

    def s_max(self):

        return np.nanmax([line.length_2d() for _, line in self.named_topoprofiles])

    def min_z_topo(self):

        return np.nanmin([line.z_min() for _, line in self.named_topoprofiles])

    def max_z_topo(self):

        return np.nanmax([line.z_max() for _, line in self.named_topoprofiles])

    def min_z_plane_attitudes(self):

        # TODO:  manage case for possible nan p_z values
        return np.nanmin([plane_attitude.pt_3d.p_z for plane_attitude_set in self.profile_attitudes for plane_attitude in
                    plane_attitude_set if 0.0 <= plane_attitude.sign_hor_dist <= self.s_max()])

    def max_z_plane_attitudes(self):

        # TODO:  manage case for possible nan p_z values
        return np.nanmax([plane_attitude.pt_3d.p_z for plane_attitude_set in self.profile_attitudes for plane_attitude in
                    plane_attitude_set if 0.0 <= plane_attitude.sign_hor_dist <= self.s_max()])

    def min_z_curves(self):

        return np.nanmin([pt_2d.p_y for multiline_2d_list in self.projected_lines for multiline_2d in multiline_2d_list for line_2d in
                    multiline_2d.lines for pt_2d in line_2d.pts if 0.0 <= pt_2d.p_x <= self.s_max()])

    def max_z_curves(self):

        return np.nanmax([pt_2d.p_y for multiline_2d_list in self.projected_lines for multiline_2d in multiline_2d_list for line_2d in
                    multiline_2d.lines for pt_2d in line_2d.pts if 0.0 <= pt_2d.p_x <= self.s_max()])

    def z_min(self):

        min_z = self.min_z_topo()

        if len(self.profile_attitudes) > 0:
            min_z = np.nanmin([min_z, self.min_z_plane_attitudes()])

        if len(self.projected_lines) > 0:
            min_z = np.nanmin([min_z, self.min_z_curves()])

        return min_z

    def z_max(self):

        max_z = self.max_z_topo()

        if len(self.profile_attitudes) > 0:
            max_z = np.nanmax([max_z, self.max_z_plane_attitudes()])

        if len(self.projected_lines) > 0:
            max_z = np.nanmax([max_z, self.max_z_curves()])

        return max_z

    def clear_attitudes(self):
        """
        Clear projected _attitudes content.

        :return:
        """

        self.profile_attitudes = None

    def add_plane_attitudes(self, plane_attitudes):

        self.profile_attitudes.append(plane_attitudes)

    def add_curves(self, lMultilines, lIds):

        self.projected_lines.append(lMultilines)
        self.projected_lines_ids.append(lIds)


class GeoProfilesSet_:
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
    representing a single geological profile.
    """

    def __init__(self,
                 topo_profile: Optional[TopographicProfile] = None,
                 profile_attitudes: Optional[AttitudesProfile] = None,
                 lines_intersections: Optional[IntersectionsProfile] = None,
                 polygons_intersections: Optional[IntersectionsProfile] = None
                 ):

        if topo_profile:
            check_type(topo_profile, "Topographic profile", TopographicProfile)

        if profile_attitudes:
            check_type(profile_attitudes, "Attitudes", AttitudesProfile)

        if lines_intersections:
            check_type(lines_intersections, "Line intersections", IntersectionsProfile)

        if polygons_intersections:
            check_type(polygons_intersections, "Polygon intersections", IntersectionsProfile)

        self._topo_profile = topo_profile
        self._profile_attitudes = profile_attitudes
        self._lines_intersections = lines_intersections
        self._polygons_intersections = polygons_intersections

    @property
    def topo_profile(self):
        """

        :return:
        """

        return self._topo_profile

    @topo_profile.setter
    def topo_profile(self,
                     profile: TopographicProfile):
        """

        :param profile: the profile.
        :return:

        """

        check_type(profile, "Topographic profile", TopographicProfile)
        self._topo_profile = profile

    def clear_topo_profile(self):
        """

        :return:
        """

        self._topo_profile = None

    @property
    def profile_attitudes(self):
        """

        :return:
        """

        return self._profile_attitudes

    @profile_attitudes.setter
    def profile_attitudes(self,
                          prj_attitudes: AttitudesProfile):
        """
        Set the projected _attitudes content.

        :param prj_attitudes: projected _attitudes.
        """

        check_type(prj_attitudes, "Projected _attitudes", List)
        for el in prj_attitudes:
            check_type(el, "Projected attitude", ProfileAttitude)

        self._profile_attitudes = prj_attitudes

    def clear_attitudes(self):
        """
        Clear projected _attitudes content.

        :return:
        """

        self._profile_attitudes = None

    @property
    def lines_intersections(self):
        """

        :return:
        """

        return self._lines_intersections

    @lines_intersections.setter
    def lines_intersections(self,
                            lines_intersections: IntersectionsProfile):
        """
        Set the line intersections content.

        :param lines_intersections: line intersections.
        :type lines_intersections: LinesIntersections.
        :return:
        """

        check_type(lines_intersections, "Lines intersections", IntersectionsProfile)

        self._lines_intersections = lines_intersections

    def clear_lines_intersections(self):
        """
        Clear line intersections content.

        :return:
        """

        self._lines_intersections = None

    @property
    def polygons_intersections(self):
        """

        :return:
        """

        return self._polygons_intersections

    @polygons_intersections.setter
    def polygons_intersections(self,
                               polygons_intersections: IntersectionsProfile):
        """
        Set the polygons intersections content.

        :param polygons_intersections: polygons intersections.
        :return:
        """

        check_type(polygons_intersections, "Polygons intersections", IntersectionsProfile)

        self._polygons_intersections = polygons_intersections

    '''
    def plot_attitudes(self, color="red"):
        """

        :return:
        """

        self.fig = self._attitudes.plot(
            self.fig,
            self.length_2d(),
            color=color
        )
    '''

    '''
    def plot(self,
             topo_profile=True,
             attitudes=True,
             lines_intersections=True,
             polygon_intersections=True,
             topo_profile_color="blue",
             attitudes_color="red",
             line_intersections_color="orange",
             aspect: numbers.Real = 1,
             width: numbers.Real = 18.5,
             height: numbers.Real = 10.5,
             **params
             ):
        """

        :param topo_profile:
        :param attitudes:
        :param lines_intersections:
        :param polygon_intersections:
        :param line_projections:
        :return:
        """

        fig, ax = plt.subplots()
        fig.set_size_inches(width, height)

        ax.set_aspect(aspect)

        if 'plot_z_min' in params and 'plot_z_max' in params:
            plot_z_min = params['plot_z_min']
            plot_z_max = params['plot_z_max']
        else:

            z_range = self.z_max() - self.z_min()
            plot_z_min = self.z_min() - z_padding * z_range
            plot_z_max = self.z_max() + z_padding * z_range

        if topo_profile and self._topo_profile:

            ax.set_ylim([plot_z_min, plot_z_max])
            ax.plot(
                self._topo_profile.s(),
                self._topo_profile.z(),
                color=topo_profile_color
            )

        if attitudes and self._attitudes:

            self._attitudes.plot(
                fig,
                self.length_2d(),
                color=attitudes_color
            )

        if lines_intersections and self._lines_intersections:

            self._lines_intersections.plot(
                fig,
                self.length_2d(),
                color=line_intersections_color
            )

        if polygon_intersections and self._polygons_intersections:

            self._polygons_intersections.plot(
                fig,
                self.length_2d()
            )
    '''

    def s_min(self):
        """

        :return:
        """

        return self.topo_profile.s_min()

    def s_max(self):
        """

        :return:
        """

        return self.topo_profile.s_max()

    def z_min(self):
        """

        :return:
        """

        return self.topo_profile.z_min()

    def z_max(self):
        """

        :return:
        """

        return self.topo_profile.z_max()

    '''
    def add_intersections_pts(self, intersection_list):
        """

        :param intersection_list:
        :return:
        """

        self._lines_intersections += intersection_list

    def add_intersections_lines(self, formation_list, intersection_line3d_list, intersection_polygon_s_list2):
        """

        :param formation_list:
        :param intersection_line3d_list:
        :param intersection_polygon_s_list2:
        :return:
        """

        self._polygons_intersections = zip(formation_list, intersection_line3d_list, intersection_polygon_s_list2)
    '''

    def length_2d(self) -> numbers.Real:
        """
        Returns the 2D length of the profile.

        :return: the 2D profile length.
        :rtype: numbers.Real.
        """

        return self._topo_profile.profile_length()


class GeoProfileSet:
    """
    Represents a set of Geoprofile elements,
    stored as a list
    """

    def __init__(self,
                 topo_profiles: Optional[TopographicProfileSet] = None,
                 profile_attitudes: Optional[AttitudesSet] = None,
                 lines_intersections_set: Optional[PointSegmentCollectionsSet] = None,
                 polygons_intersections_set: Optional[PointSegmentCollectionsSet] = None
                 ):

        if topo_profiles:
            #print("Num. topographic profiles: {}".format(len(topo_profiles)))
            check_type(topo_profiles, "Topographic profiles set", TopographicProfileSet)

        if profile_attitudes:
            #print("Num. attitude profiles: {}".format(len(profile_attitudes)))
            check_type(profile_attitudes, "Attitudes set", AttitudesSet)

        if lines_intersections_set:
            check_type(lines_intersections_set, "Lines intersections set", PointSegmentCollectionsSet)

        if polygons_intersections_set:
            check_type(polygons_intersections_set, "Polygons_intersections set", PointSegmentCollectionsSet)

        self._topo_profiles_set = topo_profiles
        self._attitudes_set = profile_attitudes
        self._lines_intersections_set = lines_intersections_set
        self._polygons_intersections_set = polygons_intersections_set

    def parameters(self) -> List:
        """
        Returns all the attributes of the class.

        :return:
        """

        return [
            self._topo_profiles_set,
            self._attitudes_set,
            self._lines_intersections_set,
            self._polygons_intersections_set
        ]

    @property
    def topo_profiles_set(self):
        """

        :return:
        """

        return self._topo_profiles_set

    @topo_profiles_set.setter
    def topo_profiles_set(self,
        topo_profiles_set: TopographicProfileSet):
        """

        :param topo_profiles_set: the scalar values profiles.
        :type topo_profiles_set: TopographicProfile.
        :return:

        """

        check_type(topo_profiles_set, "Topographic profiles set", TopographicProfileSet)
        self._topo_profiles_set = topo_profiles_set

    @property
    def profile_attitudes(self):
        """

        :return:
        """

        return self._attitudes_set

    @profile_attitudes.setter
    def profile_attitudes(self,
                          attitudes_set: AttitudesSet):
        """

        :param attitudes_set: the attitudes set.
        :type attitudes_set: AttitudesSet.
        :return:

        """

        check_type(attitudes_set, "Attitudes set", AttitudesSet)
        self._attitudes_set = attitudes_set

    @property
    def lines_intersections_set(self):
        """

        :return:
        """

        return self._lines_intersections_set

    @lines_intersections_set.setter
    def lines_intersections_set(self,
                                lines_intersections_set: PointSegmentCollectionsSet):
        """

        :param lines_intersections_set: the lines intersections set.
        :type lines_intersections_set: LinesIntersectionsSet.
        :return:

        """

        check_type(lines_intersections_set, "Line intersections set", PointSegmentCollectionsSet)
        self._lines_intersections_set = lines_intersections_set

    @property
    def polygons_intersections_set(self):
        """

        :return:
        """

        return self._polygons_intersections_set

    @polygons_intersections_set.setter
    def polygons_intersections_set(self,
                                   polygons_intersections_set: PointSegmentCollectionsSet):
        """

        :param polygons_intersections_set: the polygons intersections set.
        :type polygons_intersections_set: PolygonsIntersectionsSet.
        :return:

        """

        check_type(polygons_intersections_set, "Polygons intersections set", PointSegmentCollectionsSet)
        self._polygons_intersections_set = polygons_intersections_set

    def num_profiles(self) -> numbers.Integral:
        """
        Returns the number of profiles in the geoprofile set.

        :return: number of profiles in the geoprofile set.
        :rtype: numbers.Integral.
        """

        return np.nanmax(list(map(lambda lst: len(lst) if lst else 0, self.parameters())))

    def extract_geoprofile(
            self,
            ndx: numbers.Integral
    ) -> GeoProfile_:
        """
        Returns a geoprofile referencing slices of stored data.

        :param ndx: the index of the geoprofile.
        :type ndx: numbers.Integral.
        :return: the extracted Geoprofile or None.
        :rtype: GeoProfile.
        :raise: Exception.
        """

        if ndx not in range(self.num_profiles()):
            raise Exception("Geoprofile set range is in 0-{} but {} got".format(self.num_profiles() - 1, ndx))

        return GeoProfile(
            topo_profile=self.topo_profiles_set[ndx] if self.topo_profiles_set and ndx < len(self.topo_profiles_set) else None,
            profile_attitudes=self.profile_attitudes[ndx] if self.profile_attitudes and ndx < len(self.profile_attitudes) else None,
            lines_intersections=self.lines_intersections_set[ndx] if self.lines_intersections_set and ndx < len(self.lines_intersections_set) else None,
            polygons_intersections=self.polygons_intersections_set[ndx] if self.polygons_intersections_set and ndx < len(self.polygons_intersections_set) else None
        )

    def s_min(self):
        """

        :return:
        """

        return self.topo_profiles_set.s_min()

    def s_max(self):
        """

        :return:
        """

        return self.topo_profiles_set.s_max()

    def z_min(self):
        """

        :return:
        """

        return self.topo_profiles_set.z_min()

    def z_max(self):
        """

        :return:
        """

        return self.topo_profiles_set.z_max()

    ## inherited - TO CHECK

    def profiles_svals(self) -> List[List[numbers.Real]]:
        """
        Returns the list of the s values for the profiles.

        :return: list of the s values.
        :rtype
        """

        return [topoprofile.profile_s() for topoprofile in self._topo_profiles_set]

    def profiles_zs(self) -> List[numbers.Real]:
        """
        Returns the elevations of the profiles.

        :return: the elevations.
        :rtype: list of numbers.Real values.
        """

        return [topoprofile.elevations() for topoprofile in self._topo_profiles_set]

    def profiles_lengths_2d(self) -> List[numbers.Real]:
        """
        Returns the 2D lengths of the profiles.

        :return: the 2D profiles lengths.
        :rtype: list of numbers.Real values.
        """

        return [topoprofile.profile_length_2d() for topoprofile in self._topo_profiles_set]

    def profiles_lengths_3d(self) -> List[numbers.Real]:
        """
        Returns the 3D lengths of the profiles.

        :return: the 3D profiles lengths.
        :rtype: list of numbers.Real values.
        """

        return [topoprofile.profile_length_3d() for topoprofile in self._topo_profiles_set]

    def max_length_2d(self) -> Optional[numbers.Real]:
        """
        Returns the maximum 2D length of profiles.

        :return: the maximum profiles lengths.
        :rtype: an optional numbers.Real value.
        """

        lengths = self.profiles_lengths_2d()

        if lengths:
            return np.nanmax(lengths)
        else:
            return None

    def add_plane_attitudes(self, plane_attitudes):
        """

        :param plane_attitudes:
        :return:
        """

        self._attitudes_set.append(plane_attitudes)


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


def intersection_distances_by_profile_start_list(
        profile_line,
        intersections
):

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


def define_plot_structural_segment(
        structural_attitude,
        profile_length,
        vertical_exaggeration,
        segment_scale_factor=70.0
):

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


