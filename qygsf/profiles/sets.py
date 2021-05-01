from .chains import *
from ..geometries.shapes.space2d import XYArrayPair
from ..georeferenced.geoshapes3d import GeoPointSegmentCollections3D


class TopographicProfileSet(list):
    """

    Class storing a set of topographic profiles.
    """

    def __init__(self, topo_profiles_set: List[XYArrayPair]):
        """
        Instantiates a topographic profile set.

        :param topo_profiles_set: the topographic profile set.
        :type topo_profiles_set: List[TopographicProfile].
        """

        check_type(topo_profiles_set, "Topographic profiles set", List)
        for el in topo_profiles_set:
            check_type(el, "Topographic profile", XYArrayPair)

        super(TopographicProfileSet, self).__init__(topo_profiles_set)

    def s_min(self) -> Optional[numbers.Real]:
        """
        Returns the minimum s value in the topographic profiles.

        :return: the minimum s value in the profiles.
        :rtype: optional numbers.Real.
        """

        return np.nanmin([prof.x_min() for prof in self])

    def s_max(self) -> Optional[numbers.Real]:
        """
        Returns the maximum s value in the topographic profiles.

        :return: the maximum s value in the profiles.
        :rtype: optional numbers.Real.
        """

        return np.nanmax([prof.x_max() for prof in self])

    def z_min(self) -> Optional[numbers.Real]:
        """
        Returns the minimum elevation value in the topographic profiles.

        :return: the minimum elevation value in the profiles.
        :rtype: optional numbers.Real.
        """

        return np.nanmin([prof.y_min() for prof in self])

    def z_max(self) -> Optional[numbers.Real]:
        """
        Returns the maximum elevation value in the topographic profiles.

        :return: the maximum elevation value in the profiles.
        :rtype: optional numbers.Real.
        """

        return np.nanmax([prof.y_max() for prof in self])

    def natural_elev_range(self) -> Tuple[numbers.Real, numbers.Real]:
        """
        Returns the elevation range of the profiles.

        :return: minimum and maximum values of the considered topographic profiles.
        :rtype: tuple of two floats.
        """

        return self.z_min(), self.z_max()

    def topoprofiles_params(self):
        """

        :return:
        """

        return self.s_min(), self.s_max(), self.z_min(), self.z_max()

    def elev_stats(self) -> List[Dict]:
        """
        Returns the elevation statistics for the topographic profiles.

        :return: elevation statistics for the topographic profiles.
        :rtype: list of dictionaries.
        """

        return [topoprofile.elev_stats() for topoprofile in self]

    def slopes(self) -> List[List[Optional[numbers.Real]]]:
        """
        Returns a list of the slopes of the topographic profiles in the geoprofile.

        :return: list of the slopes of the topographic profiles.
        :rtype: list of list of topographic slopes.
        """

        return [topoprofile.slopes() for topoprofile in self]

    def abs_slopes(self) -> List[List[Optional[numbers.Real]]]:
        """
        Returns a list of the absolute slopes of the topographic profiles in the geoprofile.

        :return: list of the absolute slopes of the topographic profiles.
        :rtype: list of list of topographic absolute slopes.
        """

        return [topoprofile.abs_slopes() for topoprofile in self]

    def slopes_stats(self) -> List[Dict]:
        """
        Returns the directional slopes statistics
        for each profile.

        :return: the list of the profiles statistics.
        :rtype: List of dictionaries.
        """

        return [topoprofile.slopes_stats() for topoprofile in self]

    def absslopes_stats(self) -> List[Dict]:
        """
        Returns the absolute slopes statistics
        for each profile.

        :return: the list of the profiles statistics.
        :rtype: List of dictionaries.
        """

        return [topoprofile.absslopes_stats() for topoprofile in self]


class AttitudesSet(list):
    """

    Class storing a set of topographic profiles.
    """

    def __init__(self,
        attitudes_set: Optional[List[ProfileAttitude]] = None
    ):
        """
        Instantiates an attitudes set.

        :param attitudes_set: the optional initial attitudes set.
        :type attitudes_set: List[Attitude].
        """

        if attitudes_set:
            check_type(attitudes_set, "Attitude set", List)
            for el in attitudes_set:
                check_type(el, "Attitude", ProfileAttitude)

        super(AttitudesSet, self).__init__()

        if attitudes_set:
            self.append(attitudes_set)


class PointSegmentCollectionsSet(list):
    """

    Class storing a set of point-segment collections.
    """

    def __init__(self, ptsegm_collects_set: List[GeoPointSegmentCollections3D]):
        """
        Instantiates an lines intersections set.

        :param ptsegm_collects_set: the PointSegmentCollections set.
        :type ptsegm_collects_set: List[PointSegmentCollections].
        """

        check_type(ptsegm_collects_set, "Point-segment collections set", List)
        for el in ptsegm_collects_set:
            check_type(el, "Point-segment collections", GeoPointSegmentCollections3D)

        super(PointSegmentCollectionsSet, self).__init__(ptsegm_collects_set)


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

        return np.nanmin(list(map(np.nanmin, self.z_array)))

    def max_z(self):

        return np.nanmax(list(map(np.nanmax, self.z_array)))

    @property
    def absolute_slopes(self):

        return list(map(np.fabs, self.dir_slopes))