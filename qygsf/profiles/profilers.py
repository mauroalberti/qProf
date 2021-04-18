
from typing import List, Iterable
from operator import attrgetter

from ..georeferenced.geoshapes2d import *
from ..georeferenced.rasters import *
from ..geology.base import *

from .sets import *
from ..orientations.orientations import Axis
from ..geometries.shapes.space2d import PointSegmentCollection2D


def georef_attitudes_3d_from_grid(
        structural_data: List[GeorefAttitude],
        height_source: GeoArray,
) -> List[GeorefAttitude]:
    """
    Create a set of 3D georeferenced attitudes, extracting heights from a grid.

    :param structural_data: the set of georeferenced attitudes
    :type structural_data: List[GeorefAttitude]
    :param height_source: the elevation source
    :type height_source: GeoArray
    :return: list of GeorefAttitude values
    :rtype: List[GeorefAttitude]
    :raise: Exception.
    """

    if not isinstance(height_source, GeoArray):

        raise Exception(
            "Height source should be GeoArray instead of {}".format(
                type(height_source)
            )
        )

    attitudes_3d = []

    for georef_attitude in structural_data:

        pt3d = height_source.interpolate_bilinear_point(
            pt=georef_attitude.posit
        )

        if pt3d:
            attitudes_3d.append(
                GeorefAttitude(
                    georef_attitude.id,
                    pt3d,
                    georef_attitude.attitude
                )
            )

    return attitudes_3d


class LinearProfiler:
    """
    Class storing a linear (straight) profile.
    It intends to represent a vertical profile.
    In a possible future implementations, it would be superseded by a
    plane profiler, not necessarily vertical.
    """

    def __init__(self,
                 start_pt: Union[Point2D, Point3D],
                 end_pt: Union[Point2D, Point3D],
                 densify_distance: numbers.Real,
                 epsg_cd: numbers.Integral
                 ):
        """
        Instantiates a 2D linear profile object.
        It is represented by two 2D points and by a densify distance.

        :param start_pt: the profile start point.
        :param end_pt: the profile end point.
        :param densify_distance: the distance with which to densify the segment profile.
        """

        check_type(start_pt, "Input start point", (Point2D, Point3D))

        check_type(end_pt, "Input end point", (Point2D, Point3D))

        '''
        if start_pt.georeferenced != end_pt.georeferenced:
            raise Exception("Both points must have same CRS")
        '''

        if start_pt.distance(end_pt) == 0.0:
            raise Exception("Input segment length cannot be zero")

        check_type(densify_distance, "Input densify distance", numbers.Real)

        if not isfinite(densify_distance):
            raise Exception("Input densify distance must be finite")

        if densify_distance <= 0.0:
            raise Exception("Input densify distance must be positive")

        self._start_pt = Point2D(
            x=start_pt.x,
            y=start_pt.y
        )

        self._end_pt = Point2D(
            x=end_pt.x,
            y=end_pt.y
        )

        self._crs = Crs(epsg_cd)

        self._densify_dist = float(densify_distance)

    def start_pt(self) -> Point2D:
        """
        Returns a copy of the segment start point.

        :return: start point copy.
        :rtype: Point.
        """

        return self._start_pt.clone()

    def end_pt(self) -> Point2D:
        """
        Returns a copy of the segment end point.

        :return: end point copy.
        :rtype: Point.
        """

        return self._end_pt.clone()

    def to_line(self) -> Line2D:
        """
        Convert to a line.

        """

        return Line2D(
            pts=[
                self.start_pt(),
                self.end_pt()
            ]
        )

    def densify_dist(self) -> numbers.Real:
        """
        Returns the densify distance of the profiler.

        :return: the densify distance of the profiler.
        :rtype: numbers.Real.
        """

        return self._densify_dist

    def __repr__(self):
        """
        Representation of a profile instance.

        :return: the textual representation of the instance.
        :rtype: str.
        """

        return "LinearProfiler(\n\tstart_pt = {},\n\tend_pt = {},\n\tdensify_distance = {})".format(
            self.start_pt(),
            self.end_pt(),
            self.densify_dist()
        )

    @property
    def crs(self) -> Crs:
        """
        Returns the CRS of the profile.

        :return: the CRS of the profile.
        :rtype: Crs.
        """

        return self._crs

    @property
    def epsg_code(self) -> numbers.Integral:
        """
        Returns the EPSG code of the profile.

        :return: the EPSG code of the profile.
        """

        return self.crs.epsg_code

    def clone(self) -> 'LinearProfiler':
        """
        Returns a deep copy of the current linear profiler.

        :return: a deep copy of the current linear profiler.
        :rtype: LinearProfiler.
        """

        return LinearProfiler(
            start_pt=self.start_pt().clone(),
            end_pt=self.end_pt().clone(),
            densify_distance=self.densify_dist(),
            epsg_cd=self.epsg_code
        )

    def segment(self) -> Segment2D:
        """
        Returns the horizontal segment representing the profile.

        :return: segment representing the profile.
        :rtype: Segment.
        """

        return Segment2D(
            start_pt=self._start_pt,
            end_pt=self._end_pt
        )

    def length(self) -> numbers.Real:
        """
        Returns the length of the profiler section.

        :return: length of the profiler section.
        :rtype: numbers.Real.
        """

        return self.segment().length()

    def vector(self) -> Vect2D:
        """
        Returns the horizontal vector representing the profile.

        :return: vector representing the profile.
        :rtype: Vect.
        """

        return Vect2D(
            x=self.segment().delta_x(),
            y=self.segment().delta_y()
        )

    def versor(self) -> Vect2D:
        """
        Returns the horizontal versor (unit vector) representing the profile.

        :return: vector representing the profile.
        :rtype: Vect.
        """

        return self.vector().versor()

    def densified_2d_steps(self) -> array:
        """
        Returns an array made up by the incremental steps (2D distances) along the profile.

        :return: array storing incremental steps values.
        :rtype: array.
        """

        return self.segment().densify2d_asSteps(self._densify_dist)

    def num_pts(self) -> numbers.Integral:
        """
        Returns the number of steps making up the profile.

        :return: number of steps making up the profile.
        :rtype: numbers.Integral.
        """

        return len(self.densified_2d_points())

    def densified_2d_points(self) -> List[Point2D]:
        """
        Returns the list of densified points.

        :return: list of densified points.
        """

        return self.segment().densify2d_asPts(densify_distance=self._densify_dist)

    def vertical_plane(self) -> CPlane3D:
        """
        Returns the vertical plane of the segment, as a Cartesian plane.

        :return: the vertical plane of the segment, as a Cartesian plane.
        """

        return Segment3D.from2D(self.segment()).vertical_plane()

    def normal_versor(self) -> Vect3D:
        """
        Returns the perpendicular (horizontal) versor to the profile (vertical) plane.

        :return: the perpendicular (horizontal) versor to the profile (vertical) plane.
        """

        return self.vertical_plane().normVersor()

    def left_norm_vers(self) -> Vect2D:
        """
        Returns the left horizontal normal versor.

        :return: the left horizontal normal versor.
        """

        return Vect2D(
            x= - self.versor().y,
            y= self.versor().x
        )

    def right_norm_vers(self) -> Vect2D:
        """
        Returns the right horizontal normal versor.

        :return: the right horizontal normal versor.
        """

        return Vect2D(
            x= self.versor().y,
            y= - self.versor().x
        )

    def vector_offset(self,
                      vect: Vect2D
                      ) -> 'LinearProfiler':
        """
        Returns a new LinearProfiler instance, horizontally offset by the
        provided vector horizontal components.
        """

        dx, dy = vect.x, vect.y

        shifted_start_pt = self.start_pt().shift(
            sx=dx,
            sy=dy
        )

        shifted_end_pt = self.end_pt().shift(
            sx=dx,
            sy=dy
        )

        return LinearProfiler(
            start_pt=shifted_start_pt,
            end_pt=shifted_end_pt,
            densify_distance=self.densify_dist(),
            epsg_cd=self.epsg_code
        )

    def right_parallel_offset(self,
                              offset: numbers.Real) -> 'LinearProfiler':
        """
        Returns a copy of the current linear profiler, offset to the right by the provided offset distance.

        :param offset: the lateral offset to apply to create the new LinearProfiler.
        :return: the offset linear profiler.
        """

        return self.vector_offset(vect=self.right_norm_vers().scale(offset))

    def left_parallel_offset(self,
                             offset: numbers.Real) -> 'LinearProfiler':
        """
        Returns a copy of the current linear profiler, offset to the left by the provided offset distance.

        :param offset: the lateral offset to apply to create the new LinearProfiler.
        :return: the offset linear profiler.
        """

        return self.vector_offset(vect=self.left_norm_vers().scale(offset))

    def point_in_profile(self, pt: Point2D) -> bool:
        """
        Checks whether a point lies in the profiler plane.

        :param pt: the point to check.
        :return: whether the point lie in the profiler plane.
        :raise; Exception.
        """

        check_type(pt, 'Point2D', Point2D)
        return self.vertical_plane().isPointInPlane(pt)

    def point_distance(self, pt: Point3D) -> numbers.Real:
        """
        Calculates the point distance from the profiler plane.

        :param pt: the point to check.
        :type pt: Point.
        :return: the point distance from the profiler plane.
        :rtype: numbers.Real.
        :raise; Exception.
        """

        return self.vertical_plane().pointDistance(pt)

    def sample_grid(
            self,
            grid: GeoArray) -> array:
        """
        Sample grid values along the profiler points.

        :param grid: the input grid
        :type grid: GeoArray.
        :return: array storing the z values sampled from the grid,
        :rtype: array.
        :raise: Exception
        """

        if not isinstance(grid, GeoArray):
            raise Exception("Input grid must be a GeoArray but is {}".format(type(grid)))

        if self.crs != grid.crs:
            raise Exception("Input grid EPSG code must be {} but is {}".format(self.epsg_code, grid.epsg_code))

        return array('d', [grid.interpolate_bilinear(pt_2d.x, pt_2d.y) for pt_2d in self.densified_2d_points()])

    def profile_grid(
            self,
            geoarray: GeoArray) -> TopographicProfile:
        """
        Create profile from one geoarray.

        :param geoarray: the source geoarray.
        :type geoarray: GeoArray.
        :return: the profile of the scalar variable stored in the geoarray.
        :rtype: TopographicProfile.
        :raise: Exception.
        """

        check_type(geoarray, "GeoArray", GeoArray)

        return TopographicProfile(
            s_array=self.densified_2d_steps(),
            z_array=self.sample_grid(geoarray))

    def profile_grids(self,
        *grids: Iterable[GeoArray]
    ) -> List[TopographicProfile]:
        """
        Create profiles of one or more grids.

        :param grids: a set of one or more grids.
        :return:
        """

        for ndx, grid in enumerate(grids):

            check_type(grid, "{} grid".format(ndx+1), GeoArray)

        for ndx, grid in enumerate(grids):

            check_crs(self, grid)

        topo_profiles = []

        for grid in grids:

            topo_profiles.append(
                TopographicProfile(
                    s_array=self.densified_2d_steps(),
                    z_array=self.sample_grid(grid)
                )
            )

        return topo_profiles

    def intersect_line(self,
                       mline: Union[Line2D, GeoMultiLine2D],
                       ) -> PointSegmentCollection2D:
        """
        Calculates the intersection with a line/multiline.
        Note: the intersections are intended flat (in a 2D plane, not 3D).

        :param mline: the line/multiline to intersect profile with
        :type mline: Union[Line, GeoMultiLine]
        :return: the possible intersections
        :rtype: PointSegmentCollection2D
        """

        return mline.intersectSegment(self.segment())

    def intersect_lines(self,
        mlines: Iterable[Union[Line2D, GeoMultiLine2D]],
) -> List[List[Optional[Union[Point2D, Segment2D]]]]:
        """
        Calculates the intersection with a set of lines/multilines.
        Note: the intersections are intended flat (in a 2D plane, not 3D).

        :param mlines: an iterable of Lines or MultiLines to intersect profile with
        :type mlines: Iterable[Union[Line, GeoMultiLine]]
        :return: the possible intersections
        :rtype: List[List[Optional[Point, Segment]]]
        """

        results = [self.intersect_line(mline) for mline in mlines]
        valid_results = [[ndx, GeoPointSegmentCollection2D(geoms=res, epsg_code=self.epsg_code)] for ndx, res in enumerate(results) if res]

        return GeoPointSegmentCollections2D(valid_results)

    def intersect_polygon(self,
        mpolygon: GeoMPolygon2D,
        ) -> GeoLines2D:
        """
        Calculates the intersection with a shapely polygon/multipolygon.
        Note: the intersections are considered flat, i.e., in a 2D plane, not 3D.

        :param mpolygon: the shapely polygon/multipolygon to intersect profile with
        :type mpolygon: qygsf.spatial.vectorial.polygons.MGeoPolygon
        :return: the possible intersections
        :rtype: GeoLines3D
        """

        check_type(
            mpolygon,
            "Polygon",
            GeoMPolygon2D
        )

        line_shapely, epsg_code = line2d_to_shapely(self.to_line())
        return mpolygon.intersect_line(line=line_shapely)

    def intersect_polygons(self,
       mpolygons: List[GeoMPolygon2D]
    ) -> List[GeoLines2D]:
        """
        Calculates the intersection with a set of shapely polygon/multipolygon.
        Note: the intersections are intended flat (in a 2D plane, not 3D).

        :param mpolygons: the shapely set of polygon/multipolygon to intersect profile with
        :type mpolygons: List[GeoMPolygon]
        :return: the possible intersections
        :rtype: List[GeoLines]
        """

        results = []

        for mpolygon in mpolygons:

            result = self.intersect_polygon(
                mpolygon=mpolygon
            )
            if result:
                results.append(result)

        return results

    def point_along_profile_signed_s(
            self,
            pt: Point2D) -> Optional[numbers.Real]:
        """
        Calculates the point along-profile signed distance (positive in the segment direction, negative otherwise)
        from the profile start.
        The projected point must already lay in the profile vertical plane, otherwise an exception is raised.

        The implementation assumes (and also verifies) that the point lies in the profile vertical plane.
        Given that case, it calculates the signed distance from the section start point,
        by using the triangle law of sines.

        :param pt: the point on the section.
        :return: the signed distance along the profile or None if outside the segment.
        """

        if not self.point_in_profile(pt):
            raise Exception(f"Projected point should lie in the profile plane but there is a distance of {self.point_distance(pt)} units")

        if pt.is_coincident(self.start_pt()):
            return 0.0

        # the vector starting at the profile start and ending at the given point

        projected_vector = Segment2D(self.start_pt(), pt).vector()

        # the angle between the profile vector and the previous vector
        cos_alpha = self.vector().cosine_of_angle(projected_vector)

        signed_distance = projected_vector.length * cos_alpha

        return signed_distance

    def segment_along_profile_signed_s_tuple(self,
                                             segment: Segment2D
                                             ) -> Tuple[Optional[numbers.Real], Optional[numbers.Real]]:
        """
        Calculates the segment distances from the profiles start.
        The segment must already lay in the profile vertical plane, otherwise None is returned.

        :param segment: the analysed segment
        :return: the segment vertices distances from the profile start
        """

        segment_start_distance = self.point_along_profile_signed_s(segment.start_pt)
        segment_end_distance = self.point_along_profile_signed_s(segment.end_pt)

        return segment_start_distance, segment_end_distance

    def pt_segm_along_profile_signed_s(self,
                                       geom: Union[Point2D, Segment2D]
                                       ) -> array:
        """
        Calculates the point or segment signed distances from the profiles start.

        :param geom: point or segment
        :return: the distance(s) from the profile start
        """

        if isinstance(geom, Point2D):
            return array('d', [self.point_along_profile_signed_s(geom)])
        elif isinstance(geom, Segment2D):
            return array('d', [*self.segment_along_profile_signed_s_tuple(geom)])
        else:
            return NotImplemented

    def get_intersection_slope(self,
        intersection_vector: Vect3D
    ) -> Tuple[numbers.Real, str]:
        """
        Calculates the slope (in radians) and the downward sense ('left', 'right' or 'vertical')
        for a profile-laying vector.

        :param intersection_vector: the profile-plane lying vector.
        :type intersection_vector: Vect,
        :return: the slope (in radians) and the downward sense.
        :rtype: Tuple[numbers.Real, str].
        :raise: Exception.
        """

        if not isinstance(intersection_vector, Vect3D):
            raise Exception("Input argument should be Vect but is {}".format(type(intersection_vector)))

        angle = degrees(acos(self.normal_versor().cosine_of_angle(intersection_vector)))
        if abs(90.0 - angle) > 1.0e-4:
            raise Exception("Input argument should lay in the profile plane")

        slope_radians = abs(radians(intersection_vector.slope_degr()))

        scalar_product_for_downward_sense = self.vector().dot_product(intersection_vector.downward())
        if scalar_product_for_downward_sense > 0.0:
            intersection_downward_sense = "right"
        elif scalar_product_for_downward_sense == 0.0:
            intersection_downward_sense = "vertical"
        else:
            intersection_downward_sense = "left"

        return slope_radians, intersection_downward_sense

    def calculate_axis_intersection(self,
                                    map_axis: Axis,
                                    structural_pt: Point2D) -> Optional[Point3D]:
        """
        Calculates the optional intersection point between an axis passing through a point
        and the profiler plane.

        :param map_axis: the projection axis.
        :type map_axis: Axis.
        :param structural_pt: the point through which the axis passes.
        :type structural_pt: Point.
        :return: the optional intersection point.
        :raise: Exception.
        """

        if not isinstance(map_axis, Axis):
            raise Exception("Map axis should be Axis but is {}".format(type(map_axis)))

        if not isinstance(structural_pt, Point2D):
            raise Exception("Structural point should be Point but is {}".format(type(structural_pt)))

        axis_versor = map_axis.as_direction().as_versor()

        l, m, n = axis_versor.x, axis_versor.y, axis_versor.z

        axis_param_line = ParamLine3D(structural_pt, l, m, n)

        return axis_param_line.intersect_cartes_plane(self.vertical_plane())

    def calculate_intersection_versor(
            self,
            attitude_plane: Plane,
            attitude_pt: Point3D) -> Optional[Vect3D]:
        """
        Calculate the intersection versor between the plane profiler and
        a geological plane with location defined by a Point.

        :param attitude_plane:
        :type attitude_plane: Plane,
        :param attitude_pt: the attitude point.
        :type attitude_pt: Point.
        :return:
        """

        if not isinstance(attitude_plane, Plane):
            raise Exception("Attitude plane should be Plane but is {}".format(type(attitude_plane)))

        if not isinstance(attitude_pt, Point3D):
            raise Exception("Attitude point should be Point3D but is {}".format(type(attitude_pt)))

        '''
        if self.crs != attitude_pt.crs:
            raise Exception("Attitude point should has EPSG {} but has {}".format(self.epsg_code, attitude_pt.epsg_code))
        '''

        putative_inters_versor = self.vertical_plane().intersVersor(
            CPlane3D.from_geological_plane(attitude_plane, attitude_pt))

        if not putative_inters_versor.is_valid:
            return None

        return putative_inters_versor

    def nearest_attitude_projection(
            self,
            georef_attitude: GeorefAttitude) -> Point3D:
        """
        Calculates the nearest projection of a given attitude on a vertical plane.

        :param georef_attitude: geological attitude.
        :type georef_attitude: GeorefAttitude
        :return: the nearest projected point on the vertical section.
        :rtype: qygsf.spatial.vectorial.geometries.Point.
        :raise: Exception.
        """

        if not isinstance(georef_attitude, GeorefAttitude):
            raise Exception("georef_attitude point should be GeorefAttitude but is {}".format(type(georef_attitude)))

        """
        if self.crs != georef_attitude.posit.crs:
            raise Exception("Attitude point should has EPSG {} but has {}".format(self.epsg_code, georef_attitude.posit.epsg_code))
        """

        attitude_cplane = CPlane3D.from_geological_plane(
            geol_plane=georef_attitude.attitude,
            pt=georef_attitude.posit)
        # attitude_cplane = georef_attitude.attitude.to_cartesian_plane(georef_attitude.posit)
        intersection_versor = self.vertical_plane().intersVersor(attitude_cplane)
        dummy_inters_pt = self.vertical_plane().intersPoint(attitude_cplane)
        dummy_structural_vect = Segment3D(dummy_inters_pt, georef_attitude.posit).vector()
        dummy_distance = dummy_structural_vect.dot_product(intersection_versor)
        offset_vector = intersection_versor.scale(dummy_distance)

        projected_pt = Point3D(
            x=dummy_inters_pt.x + offset_vector.x,
            y=dummy_inters_pt.y + offset_vector.y,
            z=dummy_inters_pt.z + offset_vector.z
        )

        return projected_pt

    def map_attitude_to_section(
            self,
            georef_attitude: GeorefAttitude,
            map_axis: Optional[Axis] = None,
            max_profile_distance: Optional[numbers.Real] = None
    ) -> Optional[ProfileAttitude]:
        """
        Project a georeferenced attitude to the section.

        :param georef_attitude: the georeferenced attitude.
        :type georef_attitude: GeorefAttitude.
        :param map_axis: the map axis.
        :type map_axis: Optional[Axis].
        :param max_profile_distance: the maximum projection distance between the attitude and the profile
        :type max_profile_distance: Optional[numbers.Real]
        :return: the optional planar attitude on the profiler vertical plane.
        :rtype: Optional[PlanarAttitude].
        """

        if not isinstance(georef_attitude, GeorefAttitude):
            raise Exception("Georef attitude should be GeorefAttitude but is {}".format(type(georef_attitude)))

        """
        if self.crs != georef_attitude.posit.crs:
            raise Exception("Attitude point should has EPSG {} but has {}".format(self.epsg_code, georef_attitude.posit.epsg_code))
        """

        if map_axis:
            if not isinstance(map_axis, Axis):
                raise Exception("Map axis should be Axis but is {}".format(type(map_axis)))

        # intersection versor

        intersection_versor = self.calculate_intersection_versor(
            attitude_plane=georef_attitude.attitude,
            attitude_pt=georef_attitude.posit
        )

        # calculate slope of geological plane onto section plane

        slope_radians, intersection_downward_sense = self.get_intersection_slope(intersection_versor)

        # intersection point

        if map_axis is None:
            intersection_point_3d = self.nearest_attitude_projection(
                georef_attitude=georef_attitude)
        else:
            intersection_point_3d = self.calculate_axis_intersection(
                map_axis=map_axis,
                structural_pt=georef_attitude.posit)

        if not intersection_point_3d:
            return None

        # distance along projection vector

        dist = georef_attitude.posit.distance(intersection_point_3d)

        if max_profile_distance is not None and dist > max_profile_distance:
            return None

        # horizontal spat_distance between projected structural point and profile start

        signed_distance_from_section_start = self.point_along_profile_signed_s(intersection_point_3d.to2d())

        # solution for current structural point

        profile_attitude = ProfileAttitude(
            rec_id=georef_attitude.id,
            s=signed_distance_from_section_start,
            z=intersection_point_3d.z,
            slope_degr=degrees(slope_radians),
            down_sense=intersection_downward_sense,
            dist=dist,
            src_dip_dir=georef_attitude.attitude.dd,
            src_dip_ang=georef_attitude.attitude.da
        )
        # print("Profile attitude: {}".format(profile_attitude))

        return profile_attitude

    def map_georef_attitudes_to_section(
        self,
        attitudes_3d: List[GeorefAttitude],
        mapping_method: dict,
        max_profile_distance: Optional[numbers.Real] = None
    ) -> Optional[List[ProfileAttitude]]:
        """
        Projects a set of georeferenced space3d attitudes onto the section profile.

        :param attitudes_3d: the set of georeferenced space3d attitudes to plot on the section.
        :type attitudes_3d: List[GeorefAttitude]
        :param mapping_method: the method to map the attitudes to the section.
        ;type mapping_method; Dict.
        :param max_profile_distance: the maximum projection distance between the attitude and the profile
        :type max_profile_distance: Optional[numbers.Real]
        :return: sorted list of ProfileAttitude values.
        :rtype: Optional[List[Attitude]]
        :raise: Exception.
        """

        if mapping_method['method'] == 'nearest':
            results = [self.map_attitude_to_section(georef_att, max_profile_distance=max_profile_distance) for georef_att in attitudes_3d]
        elif mapping_method['method'] == 'common axis':
            map_axis = Axis(mapping_method['trend'], mapping_method['plunge'])
            results = [self.map_attitude_to_section(georef_att, map_axis, max_profile_distance=max_profile_distance) for georef_att in attitudes_3d]
        elif mapping_method['method'] == 'individual axes':
            if len(mapping_method['individual_axes_values']) != len(attitudes_3d):
                raise Exception(
                    "Individual axes values are {} but attitudes are {}".format(
                        len(mapping_method['individual_axes_values']),
                        len(attitudes_3d)
                    )
                )

            results = []
            for georef_att, (trend, plunge) in zip(attitudes_3d, mapping_method['individual_axes_values']):
                try:
                    map_axis = Axis(trend, plunge)
                    results.append(self.map_attitude_to_section(georef_att, map_axis, max_profile_distance=max_profile_distance))
                except Exception as e:
                    print("Exception while processing individual axes values: {}".format(e))
                    continue
        else:
            return NotImplemented

        if results is None:
            return None

        results = list(filter(lambda res: res is not None, results))

        if len(results) == 0:
            return None

        return AttitudesProfile(sorted(results, key=attrgetter('s')))

    def parse_profile_intersections(
            self,
            intersections: GeoPointSegmentCollections2D
    ) -> IntersectionsProfile:
        """
        Parse the profile intersections for incorporation
        as elements in a geoprofile.

        :param intersections: the intersections
        :return:
        """

        parsed_intersections = []

        for line_id, inters_geoms in intersections:

            intersections_arrays = [self.pt_segm_along_profile_signed_s(geom) for geom in inters_geoms]

            parsed_intersections.append(ArrayList(line_id, intersections_arrays))

        return IntersectionsProfile(parsed_intersections)


class ParallelProfiler(list):
    """
    Parallel linear profiler.
    """

    def __init__(self,
                 parallel_profiler: List[LinearProfiler]):
        """
        :param parallel_profiler: a set of parallel profilers
        :return:
        """

        check_type(parallel_profiler, "Profilers", List)
        for el in parallel_profiler:
            check_type(el, "Profiler", LinearProfiler)
        template_profile = parallel_profiler[0]
        for el in parallel_profiler[1:]:
            check_crs(
                template_element=template_profile,
                checked_element=el
            )

        super(ParallelProfiler, self).__init__(parallel_profiler)

        self._crs = Crs(template_profile.epsg_code)

    @classmethod
    def fromBaseProfiler(cls,
                         base_profiler: LinearProfiler,
                         profs_num: numbers.Integral,
                         profs_offset: numbers.Real,
                         profs_arr: str = "central",  # one of: "left", "central", "right"
                         ):
        """
        Initialize the parallel linear profilers.

        :param base_profiler: the base profiler.
        :type base_profiler: LinearProfiler.
        :param profs_num: the number of profilers to create.
        :type profs_num: numbers.Integral.
        :param profs_offset: the lateral offset between profilers.
        :type profs_offset: numbers.Real.
        :param profs_arr: profiles arrangement: one of "left", "central", "right".
        :type: str.
        :return: the parallel linear profilers.
        :type: ParallelProfilers.
        :raise: Exception.

        """

        check_type(base_profiler, "Base profiler", LinearProfiler)

        check_type(profs_num, "Profilers number", numbers.Integral)
        if profs_num < 2:
            raise Exception("Profilers number must be >= 2")

        check_type(profs_arr, "Profilers arrangement", str)
        if profs_arr not in ["central", "left", "right"]:
            raise Exception("Profilers arrangement must be 'left', 'central' (default) or 'right'")

        if profs_arr == "central" and profs_num % 2 != 1:
            raise Exception("When profilers arrangement is 'central' profilers number must be odd")

        if profs_arr == "central":

            side_profs_num = profs_num // 2
            num_left_profs = num_right_profs = side_profs_num

        elif profs_arr == "left":

            num_left_profs = profs_num -1
            num_right_profs = 0

        else:

            num_right_profs = profs_num -1
            num_left_profs = 0

        profilers = []

        for i in range(num_left_profs, 0, -1):

            current_offset = profs_offset * i

            profilers.append(base_profiler.left_parallel_offset(offset=current_offset))

        profilers.append(base_profiler.clone())

        for i in range(1, num_right_profs + 1):

            current_offset = profs_offset * i

            profilers.append(base_profiler.right_parallel_offset(offset=current_offset))

        return cls(profilers)

    def __repr__(self) -> str:
        """
        Represents a parallel linear profilers set.

        :return: the textual representation of the parallel linear profiler set.
        :rtype: str.
        """

        inner_profilers = "\n".join([repr(profiler) for profiler in self])
        return "ParallelProfilers([\n{}]\n)".format(inner_profilers)

    @property
    def crs(self) -> Crs:
        """
        Returns the CRS of the profiles.

        :return: the CRS of the profiles.
        :rtype: Crs.
        """

        return self._crs

    @property
    def epsg_code(self) -> numbers.Integral:
        """
        Returns the EPSG code of the profiles.

        :return: the EPSG code of the profiles.
        :rtype: numbers.Real.
        """

        return self.crs.epsg_code

    def profile_grid(
            self,
            geoarray: GeoArray) -> TopographicProfileSet:
        """
        Create profile from one geoarray.

        :param geoarray: the source geoarray.
        :type geoarray: GeoArray.
        :return: list of profiles of the scalar variable stored in the geoarray.
        :rtype: TopographicProfileSet.
        :raise: Exception.
        """

        topo_profiles = []

        for profiler in self:

            topo_profiles.append(profiler.profile_grid(geoarray))

        return TopographicProfileSet(topo_profiles)

    def map_georef_attitudes_to_section(
            self,
            attitudes_3d: List[GeorefAttitude],
            mapping_method: dict,
            max_profile_distance: Optional[numbers.Real] = None
    ) -> AttitudesSet:
        """
        Projects a set of georeferenced space3d attitudes onto the section profile.

        Projects a set of georeferenced space3d attitudes onto the section profile,

        :param attitudes_3d: the set of georeferenced space3d attitudes to plot on the section.
        :type attitudes_3d: List[GeorefAttitude]
        :param mapping_method: the method to map the attitudes to the section.
        ;type mapping_method; Dict.
        :param max_profile_distance: the maximum projection distance between the attitude and the profile
        :type max_profile_distance: Optional[numbers.Real]
        :return: an attitudes set
        :rtype: AttitudesSet
        :raise: Exception.
        """

        attitudes_set = AttitudesSet()

        for profiler in self:

            profile_attitudes = profiler.map_georef_attitudes_to_section(
                attitudes_3d=attitudes_3d,
                mapping_method=mapping_method,
                max_profile_distance=max_profile_distance
            )

            attitudes_set.append(profile_attitudes)

        #print(attitudes_set)

        return attitudes_set
