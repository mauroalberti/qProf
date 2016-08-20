# -*- coding: utf-8 -*-

from __future__ import division

from math import sqrt, floor, ceil, sin, cos, tan, radians, asin, acos, atan, atan2, degrees, isnan
import numpy as np
import copy
from .qgs_tools import project_qgs_point, qgs_point_2d
from .utils import array_from_function, almost_zero
from array_utils import point_solution, formula_to_grid
from .errors import AnaliticSurfaceIOException, AnaliticSurfaceCalcException

MINIMUM_SEPARATION_THRESHOLD = 1e-10
MINIMUM_VECTOR_MAGNITUDE = 1e-10


class CartesianPoint2DT(object):
    def __init__(self, x=np.nan, y=np.nan, t=None):

        self._x = x
        self._y = y
        self._t = t

    @property
    def p_x(self):

        return self._x

    @property
    def p_y(self):

        return self._y

    @property
    def p_t(self):
        return self._t

    def clone(self):

        return CartesianPoint2DT(self.p_x, self.p_y, self.p_t)

    def spat_distance(self, another):

        return sqrt((self.p_x - another.p_x) ** 2 + (self.p_y - another.p_y) ** 2)

    def to_point3dt(self, z=np.nan, t=None):

        return CartesianPoint3DT(self.p_x, self.p_y, z, t)

    def translate_with_vector(self, displacement_vector):

        return CartesianPoint2DT(self.p_x + displacement_vector.x, self.p_y + displacement_vector.y, self.p_t)

    def spat_coincident_with(self, another, tolerance=MINIMUM_SEPARATION_THRESHOLD):

        if self.spat_distance(another) > tolerance:
            return False
        else:
            return True

    def crs_project_2d(self, srcCrs, destCrs):

        qgis_pt = qgs_point_2d(self.p_x, self.p_y)
        destCrs_qgis_pt = project_qgs_point(qgis_pt, srcCrs, destCrs)

        return CartesianPoint2DT(destCrs_qgis_pt.x(), destCrs_qgis_pt.y(), self.p_t)


class CartesianSegment2DT(object):
    def __init__(self, start_pt2dt, end_pt2dt):

        self._start_pt = start_pt2dt.clone()
        self._end_pt = end_pt2dt.clone()

    @property
    def start_pt(self):

        return self._start_pt

    @property
    def end_pt(self):

        return self._end_pt

    def clone(self):

        return CartesianSegment2DT(self.start_pt, self.end_pt)

    def vector_2d(self):

        return CartesianVector2D(self.end_pt.p_x - self.start_pt.p_x,
                                 self.end_pt.p_y - self.start_pt.p_y)

    def increasing_x(self):

        if self.end_pt.p_x < self.start_pt.p_x:
            return CartesianSegment2DT(self.end_pt, self.start_pt)
        else:
            return self.clone()

    def segment_m(self):

        return (self.end_pt.p_y - self.start_pt.p_y) / (self.end_pt.p_x - self.start_pt.p_x)

    def segment_p(self):

        return self.start_pt.p_y - self.segment_m() * self.start_pt.p_x

    def intersection_pt(self, another):

        assert self.length_2d > 0.0
        assert another.length_2d > 0.0

        # at least one segment vertical
        if self.start_pt.p_x == self.end_pt.p_x:
            x0 = self.start_pt.p_x
            try:
                m1, p1 = another.segment_m(), another.segment_p()
            except:
                return None
            y0 = m1 * x0 + p1
        elif another.start_pt.p_x == another.end_pt.p_x:
            x0 = another.start_pt.p_x
            try:
                m1, p1 = self.segment_m(), self.segment_p()
            except:
                return None
            y0 = m1 * x0 + p1
        else:
            m0, p0 = self.segment_m(), self.segment_p()
            m1, p1 = another.segment_m(), another.segment_p()
            x0 = (p1 - p0) / (m0 - m1)
            y0 = m0 * x0 + p0

        return CartesianPoint2DT(x0, y0)

    @property
    def segment_x_range(self):

        if self.start_pt.p_x < self.end_pt.p_x:
            return self.start_pt.p_x, self.end_pt.p_x
        else:
            return self.end_pt.p_x, self.start_pt.p_x

    @property
    def segment_y_range(self):

        if self.start_pt.p_y < self.end_pt.p_y:
            return self.start_pt.p_y, self.end_pt.p_y
        else:
            return self.end_pt.p_y, self.start_pt.p_y

    def fast_contains_pt(self, pt2d):
        """
        to work properly, requires that the pt lies on the line defined by the segment
        """

        range_x = self.segment_x_range
        range_y = self.segment_y_range

        if range_x[0] <= pt2d.p_x <= range_x[1] or \
                                range_y[0] <= pt2d.p_y <= range_y[1]:
            return True
        else:
            return False

    @property
    def length_2d(self):

        return self.start_pt.spat_distance(self.end_pt)

    def contains_pt(self, pt2d):

        segment_length = self.length_2d
        segmentstart_pt2d_distance = self.start_pt.spat_distance(pt2d)
        segmentend_pt2d_distance = self.end_pt.spat_distance(pt2d)

        if segmentstart_pt2d_distance > segment_length or \
                        segmentend_pt2d_distance > segment_length:
            return False
        else:
            return True

    @property
    def delta_x(self):

        return self.end_pt.p_x - self.start_pt.p_x

    @property
    def delta_y(self):

        return self.end_pt.p_y - self.start_pt.p_y

    def scale(self, scale_factor):

        delta_x = self.delta_x * scale_factor
        delta_y = self.delta_y * scale_factor

        return CartesianSegment2DT(self.start_pt, CartesianPoint2DT(self.start_pt.p_x + delta_x, self.start_pt.p_y + delta_y))

    def segment_3d(self):

        return CartesianSegment3DT(self.start_pt.to_point3dt(), self.end_pt.to_point3dt())

    def densify(self, densify_distance):

        assert densify_distance > 0.0

        segment_length = self.length_2d

        assert segment_length > 0.0

        generator_vector = self.vector_2d().versor_2d().scale(densify_distance)
        interpolated_line = CartesianLine2DT([self.start_pt])
        n = 0
        while True:
            n += 1
            new_pt = self._start_pt.translate_with_vector(generator_vector.scale(n))
            if self.start_pt.spat_distance(new_pt) >= segment_length:
                break
            interpolated_line = interpolated_line.add_pt(new_pt)
        interpolated_line = interpolated_line.add_pt(self.end_pt)

        return interpolated_line


class CartesianVector2D(object):
    def __init__(self, x=np.nan, y=np.nan):
        self._x = x
        self._y = y

    @property
    def x(self):
        return self._x

    @property
    def y(self):
        return self._y

    def clone(self):
        return CartesianVector2D(self.x, self.y)

    @property
    def length(self):
        return sqrt(self.x * self.x + self.y * self.y)

    def scale(self, scale_factor):
        return CartesianVector2D(self.x * scale_factor, self.y * scale_factor)

    def versor_2d(self):
        return self.scale(1.0 / self.length)

    def add(self, another):
        return CartesianVector2D(self.x + another.x, self.y + another.y)

    def minus(self, another):
        return self.add(another.scale(-1))

    def as_vector3d(self, z=0.0):
        return CartesianVector3D(self.x, self.y, z)


class CartesianLine2DT(object):
    def __init__(self, pts_2dt=None):

        if pts_2dt is None:
            pts_2dt = []

        self._pts = [pt_2d.clone() for pt_2d in pts_2dt]

    @property
    def pts(self):

        return self._pts

    @property
    def num_pts(self):

        return len(self.pts)

    def clone(self):

        return CartesianLine2DT(self.pts)

    def reverse_direction(self):

        return CartesianLine2DT(self.pts[::-1])

    def add_pt(self, pt_2dt):

        return CartesianLine2DT(self.pts + [pt_2dt])

    def add_pts(self, pts_2dt):

        return CartesianLine2DT(self.pts + pts_2dt)

    @property
    def num_points(self):

        return len(self.pts)

    @property
    def x_list(self):

        return [pt_2dt.p_x for pt_2dt in self.pts]

    @property
    def y_list(self):

        return [pt_2dt.p_y for pt_2dt in self.pts]

    def xy_lists(self):

        return self.x_list, self.y_list

    @property
    def x_min(self):

        return min([x for x in self.x_list if not isnan(x)])

    @property
    def x_max(self):

        return max([x for x in self.x_list if not isnan(x)])

    @property
    def y_min(self):

        return min([y for y in self.y_list if not isnan(y)])

    @property
    def y_max(self):

        return max([y for y in self.y_list if not isnan(y)])

    def remove_coincident_successive_points(self):

        assert self.num_points > 0

        new_line = CartesianLine2DT([self.pts[0]])
        for ndx in range(1, self.num_points):
            if not self.pts[ndx].spat_coincident_with(new_line.pts[-1]):
                new_line = new_line.add_pt(self.pts[ndx])

        return new_line

    def as_segments2dt(self):

        pts_pairs = zip(self.pts[:-1], self.pts[1:])

        return [CartesianSegment2DT(pt_a, pt_b) for (pt_a, pt_b) in pts_pairs]

    def densify(self, sample_distance):

        assert sample_distance > 0.0

        densified_line_list = [segment.densify(sample_distance) for segment in self.as_segments2dt()]

        assert len(densified_line_list) > 0

        return CartesianMultiLine2DT(densified_line_list).as_line2dt().remove_coincident_successive_points()

    @property
    def length(self):

        length = 0.0
        for ndx in range(self.num_points - 1):
            length += self.pts[ndx].spat_distance(self.pts[ndx + 1])

        return length

    @property
    def incremental_length(self):

        incremental_length_list = []
        length = 0.0
        incremental_length_list.append(length)
        for ndx in range(self.num_points - 1):
            length += self.pts[ndx].spat_distance(self.pts[ndx + 1])
            incremental_length_list.append(length)

        return incremental_length_list

    def crs_project(self, srcCrs, destCrs):

        points = []
        for point in self.pts:
            destCrs_point = point.crs_project_2d(srcCrs, destCrs)
            points.append(destCrs_point)

        return CartesianLine2DT(points)


class CartesianMultiLine2DT(object):
    # CartesianMultiLine2DT is a list of CartesianLine2DT objects

    def __init__(self, lines_list=None):

        if lines_list is None:
            lines_list = []
        self._lines = [line_2d.clone() for line_2d in lines_list]

    @property
    def lines(self):

        return self._lines

    def add(self, line):

        return CartesianMultiLine2DT(self.lines + [line])

    def clone(self):

        return CartesianMultiLine2DT(self.lines)

    @property
    def num_parts(self):

        return len(self.lines)

    @property
    def num_points(self):

        num_elements = map(lambda x: len(x.pts), self.lines)
        return reduce(lambda x, y: x + y, num_elements)

    def is_continuous(self):

        for line_ndx in range(len(self._lines) - 1):
            if not self.lines[line_ndx].pts[-1].spat_coincident_with(self.lines[line_ndx + 1].pts[0]) or \
                    not self.lines[line_ndx].pts[-1].spat_coincident_with(self.lines[line_ndx + 1].pts[-1]):
                return False
        return True

    @property
    def x_min(self):

        return min([line.x_min for line in self.lines])

    @property
    def x_max(self):

        return max([line.x_max for line in self.lines])

    @property
    def y_min(self):

        return min([line.y_min for line in self.lines])

    @property
    def y_max(self):

        return max([line.y_max for line in self.lines])

    def is_unidirectional(self):

        for line_ndx in range(len(self.lines) - 1):
            if not self.lines[line_ndx].pts[-1].spat_coincident_with(self.lines[line_ndx + 1].pts[0]):
                return False

        return True

    def as_line2dt(self):

        return CartesianLine2DT([point for line in self.lines for point in line.pts])

    def crs_project(self, srcCrs, destCrs):

        lines = []
        for line_2d in self.lines:
            lines.append(line_2d.crs_project(srcCrs, destCrs))

        return CartesianMultiLine2DT(lines)

    def densify(self, sample_distance):

        densified_multiline_2d_list = []
        for line_2d in self.lines:
            densified_multiline_2d_list.append(line_2d.densify(sample_distance))

        return CartesianMultiLine2DT(densified_multiline_2d_list)

    def remove_coincident_points(self):

        cleaned_lines = []
        for line_2d in self.lines:
            cleaned_lines.append(line_2d.remove_coincident_successive_points())

        return CartesianMultiLine2DT(cleaned_lines)


class CartesianPoint3DT(object):
    def __init__(self, x=np.nan, y=np.nan, z=np.nan, t=None):

        self._x = x
        self._y = y
        self._z = z
        self._t = t

    @property
    def p_x(self):

        return self._x

    @property
    def p_y(self):

        return self._y

    @property
    def p_z(self):

        return self._z

    @property
    def p_t(self):

        return self._t

    def clone(self):

        return CartesianPoint3DT(self.p_x, self.p_y, self.p_z, self.p_t)

    def spat_distance(self, another):
        """
        Calculate Euclidean spatial distance between two points.

        @param  another:  the CartesianPoint3DT instance for which the spatial distance should be calculated
        @type  another:  CartesianPoint3DT.
        
        @return:  spatial distance between the two points - float.
        """

        return sqrt((self.p_x - another.p_x) ** 2 + (self.p_y - another.p_y) ** 2 + (self.p_z - another.p_z) ** 2)

    def distance_2d(self, another):

        return sqrt((self.p_x - another.p_x) ** 2 + (self.p_y - another.p_y) ** 2)

    def spat_coincident_with(self, another, tolerance=MINIMUM_SEPARATION_THRESHOLD):

        if self.spat_distance(another) > tolerance:
            return False
        else:
            return True

    def translate(self, sx=0.0, sy=0.0, sz=0.0):
        """
        Create a new point shifted by given amount from the self instance.

        @param  sx:  the shift to be applied along the x axis.
        @type  sx:  float.
        @param  sy:  the shift to be applied along the y axis.
        @type  sy:  float.
        @param  sz:  the shift to be applied along the z axis.
        @type  sz:  float.
                
        @return:  a new CartesianPoint3DT instance shifted by the given amounts with respect to the original one.
        """

        return CartesianPoint3DT(self.p_x + sx, self.p_y + sy, self.p_z + sz, self.p_t)

    def translate_with_vector(self, displacement_vector):

        return CartesianPoint3DT(self.p_x + displacement_vector.x, self.p_y + displacement_vector.y,
                                 self.p_z + displacement_vector.z, self.p_t)

    def as_vector3d(self):

        return CartesianVector3D(self.p_x, self.p_y, self.p_z)

    def delta_time(self, another):

        return another.p_t - self.p_t

    def speed(self, another):

        try:
            return self.spat_distance(another) / self.delta_time(another)
        except:
            return np.nan


class CartesianSegment3DT(object):
    def __init__(self, start_point, end_point):

        self._start_pt = start_point.clone()
        self._end_pt = end_point.clone()

    @property
    def start_pt(self):

        return self._start_pt

    @property
    def end_pt(self):

        return self._end_pt

    def clone(self):

        return CartesianSegment3DT(self.start_pt, self.end_pt)

    def as_vector3d(self):

        return CartesianVector3D(self.end_pt.p_x - self.start_pt.p_x,
                                 self.end_pt.p_y - self.start_pt.p_y,
                                 self.end_pt.p_z - self.start_pt.p_z)

    @property
    def length(self):

        return self.start_pt.spat_distance(self.end_pt)

    def trend_and_plunge(self):

        as_geol_axis = self.as_vector3d().as_geolaxis()

        return as_geol_axis.trend, as_geol_axis.plunge

    def vertical_cartes_plane(self):
        """
        Creates a vertical Cartesian plane passing through the self CartesianSegment3DT
        """

        trend, _ = self.trend_and_plunge()
        dip_dir = trend + 90.0
        if dip_dir >= 360.0:
            dip_dir -= 360.0

        return GeolPlane(dip_dir, 90.0).as_cartesplane(self.start_pt)

    def densify(self, densify_distance):

        length = self.length

        assert length > 0.0

        generator_vector = self.as_vector3d().as_versor3d().scale(densify_distance)

        interpolated_line = CartesianLine3DT([self.start_pt])
        n = 0
        while True:
            n += 1
            new_pt = self.start_pt.translate_with_vector(generator_vector.scale(n))
            if self.start_pt.spat_distance(new_pt) >= length:
                break
            interpolated_line.add_pt(new_pt)
        interpolated_line.add_pt(self.end_pt)

        return interpolated_line

    def is_point_projection_in_segment(self, pt_3d):
        """ 
        return a boolean value depending on whether the 
        projection of a point lies into the considered segment.
        The determination uses scalar product between the segment, considered
        as a vector, and the point, transformed into a vector with the start point
        given by the segment start, so that:
           0 <= scalar product <= b**2
        where b is the segment length 
        """

        pt_vector = CartesianSegment3DT(self.start_pt, pt_3d).as_vector3d()
        scal_prod = self.as_vector3d().scalar_product(pt_vector)

        if 0 <= scal_prod <= self.length ** 2:
            return True
        else:
            return False


class CartesianVector3D(object):
    def __init__(self, x=np.nan, y=np.nan, z=np.nan):

        self._x = x
        self._y = y
        self._z = z

    @property
    def x(self):

        return self._x

    @property
    def y(self):

        return self._y

    @property
    def z(self):

        return self._z

    def clone(self):

        return CartesianVector3D(self.x, self.y, self.z)

    @property
    def length(self):

        return sqrt(self.x * self.x + self.y * self.y + self.z * self.z)

    @property
    def length_horiz(self):

        return sqrt(self.x * self.x + self.y * self.y)

    def scale(self, scale_factor):

        return CartesianVector3D(self.x * scale_factor,
                                 self.y * scale_factor,
                                 self.z * scale_factor)

    def as_versor3d(self):

        return self.scale(1.0 / self.length)

    def as_downvector3d(self):

        if self.z > 0.0:
            return self.scale(-1.0)
        else:
            return self.clone()

    def add(self, another):

        return CartesianVector3D(self.x + another.x,
                                 self.y + another.y,
                                 self.z + another.z)

    def slope_radians(self):

        return atan(self.z / self.length_horiz)

    def as_geolaxis(self):

        if self.length < MINIMUM_VECTOR_MAGNITUDE:
            return None

        unit_vect = self.as_versor3d()

        plunge = - degrees(asin(unit_vect.z))  # upward negative, downward positive

        trend = 90.0 - degrees(atan2(unit_vect.y, unit_vect.x))
        if trend < 0.0:
            trend += 360.0
        elif trend > 360.0:
            trend -= 360.0

        assert 0.0 <= trend < 360.0
        assert -90.0 <= plunge <= 90.0

        return GeolAxis(trend, plunge)

    def scalar_product(self, another):

        return self.x * another.x + self.y * another.y + self.z * another.z

    def vectors_cos_angle(self, another):

        try:
            return self.scalar_product(another) / (self.length * another.length)
        except ZeroDivisionError:
            return np.nan

    def angle_degr(self, another):
        """
        angle between two vectors,
        in 0 - pi range
        """

        return degrees(acos(self.vectors_cos_angle(another)))

    def vector_product(self, another):

        x = self.y * another.z - self.z * another.y
        y = self.z * another.x - self.x * another.z
        z = self.x * another.y - self.y * another.x

        return CartesianVector3D(x, y, z)

    def by_matrix(self, matrix3x3):

        vx = matrix3x3[0, 0] * self.x + matrix3x3[0, 1] * self.y + matrix3x3[0, 2] * self.z
        vy = matrix3x3[1, 0] * self.x + matrix3x3[1, 1] * self.y + matrix3x3[1, 2] * self.z
        vz = matrix3x3[2, 0] * self.x + matrix3x3[2, 1] * self.y + matrix3x3[2, 2] * self.z

        return CartesianVector3D(vx, vy, vz)

    def as_point3dt(self):

        return CartesianPoint3DT(self.x, self.y, self.z, None)


class CartesianLine3DT(object):
    # CartesianLine3DT is a list of CartesianPoint3DT objects

    def __init__(self, pts_3dt=None):

        if pts_3dt is None:
            pts_3dt = []
        self._pts = [pt_3dt.clone() for pt_3dt in pts_3dt]

    @property
    def pts(self):

        return self._pts

    @property
    def num_pts(self):

        return len(self.pts)

    def clone(self):

        return CartesianLine3DT(self.pts)

    def add_pt(self, pt):

        self.pts.append(pt)

    def add_pts(self, pt_list):

        self._pts += pt_list

    def remove_coincident_successive_points(self):

        new_line = CartesianLine3DT(self.pts[: 1])
        for ndx in range(1, self.num_pts):
            if not self.pts[ndx].spat_coincident_with(new_line.pts[-1]):
                new_line = new_line.add_point(self.pts[ndx])
        return new_line

    def join(self, another):
        """
        Joins together two lines and returns the join as a new line without point changes, 
        with possible overlapping points 
        and orientation mismatches between the two original lines
        """

        return CartesianLine3DT(self.pts + another.pts)

    @property
    def length_3d(self):

        length = 0.0
        for ndx in range(self.num_pts - 1):
            length += self.pts[ndx].spat_distance(self.pts[ndx + 1])
        return length

    @property
    def length_2d(self):

        length = 0.0
        for ndx in range(self.num_pts - 1):
            length += self.pts[ndx].distance_2d(self.pts[ndx + 1])
        return length

    def zs(self):

        return np.array(map(lambda pt: pt.p_z, self.pts))

    def zs_not_nan(self):

        return np.array(filter(lambda pt: not isnan(pt.p_z), self.pts))

    @property
    def z_min(self):

        return np.nanmin(self.zs())

    @property
    def z_max(self):

        return np.nanmax(self.zs())

    @property
    def z_mean(self):

        return np.nanmean(self.zs())

    @property
    def z_var(self):

        return np.nanvar(self.zs())

    @property
    def z_std(self):

        return np.nanstd(self.zs())

    def incremental_length_3d(self):

        incremental_length_list = []
        length = 0.0
        incremental_length_list.append(length)
        for ndx in range(self.num_pts - 1):
            length += self.pts[ndx].spat_distance(self.pts[ndx + 1])
            incremental_length_list.append(length)

        return incremental_length_list

    def incremental_length_2d(self):

        incremental_length_list = []
        length = 0.0
        incremental_length_list.append(length)
        for ndx in range(self.num_pts - 1):
            length += self.pts[ndx].distance_2d(self.pts[ndx + 1])
            incremental_length_list.append(length)

        return incremental_length_list

    def reverse_direction(self):

        new_line = self.clone()
        new_line.pts.reverse()  # in-place operation on new_line

        return new_line

    def slopes_list(self):

        slopes_list = []
        for ndx in range(self.num_pts - 1):
            vector = CartesianSegment3DT(self.pts[ndx], self.pts[ndx + 1]).as_vector3d()
            slopes_list.append(degrees(vector.slope_radians()))
        slopes_list.append(np.nan)  # slope value for last point is unknown

        return slopes_list

    def slopes_absolute_list(self):

        slopes_list = []
        for ndx in range(self.num_pts - 1):
            vector = CartesianSegment3DT(self.pts[ndx], self.pts[ndx + 1]).as_vector3d()
            slopes_list.append(abs(degrees(vector.slope_radians())))
        slopes_list.append(np.nan)  # slope value for last point is undefined

        return slopes_list


class CartesianMultiLine3DT(object):
    # CartesianMultiLine3DT is a list of CartesianLine3DT objects


    def __init__(self, lines_list):

        self._lines = lines_list

    @property
    def lines(self):

        return self._lines

    @property
    def num_parts(self):

        return len(self.lines)

    @property
    def num_points(self):

        num_points = 0
        for line in self.lines:
            num_points += line.num_pts

        return num_points

    def is_continuous(self):

        for line_ndx in range(len(self._lines) - 1):
            if not self.lines[line_ndx].pts[-1].spat_coincident_with(self.lines[line_ndx + 1].pts[0]) or \
                    not self.lines[line_ndx].pts[-1].spat_coincident_with(self.lines[line_ndx + 1].pts[-1]):
                return False

        return True

    def is_unidirectional(self):

        for line_ndx in range(len(self.lines) - 1):
            if not self.lines[line_ndx].pts[-1].spat_coincident_with(self.lines[line_ndx + 1].pts[0]):
                return False

        return True

    def to_line3dt(self):

        return CartesianLine3DT([point for line in self.lines for point in line.pts])


class CartesianParamLine(object):
    # parametric line
    # srcPt: source CartesianPoint3DT
    # l, m, n: .....


    def __init__(self, srcPt, l, m, n):

        assert -1.0 <= l <= 1.0
        assert -1.0 <= m <= 1.0
        assert -1.0 <= n <= 1.0

        self._srcPt = srcPt
        self._l = l
        self._m = m
        self._n = n

    def intersect_cartes_plane(self, cartes_plane):
        """
        Return intersection point between parametric line and Cartesian plane
        """

        # line parameters
        x1, y1, z1 = self._srcPt.p_x, self._srcPt.p_y, self._srcPt.p_z
        l, m, n = self._l, self._m, self._n

        # Cartesian plane parameters
        a, b, c, d = cartes_plane.a, cartes_plane.b, cartes_plane.c, cartes_plane.d

        try:
            k = (a * x1 + b * y1 + c * z1 + d) / (a * l + b * m + c * n)
        except ZeroDivisionError:
            return None

        return CartesianPoint3DT(x1 - l * k,
                                 y1 - m * k,
                                 z1 - n * k)


class GeolAxis(object):
    """
    Structural axis,
    defined by trend and plunge (both in degrees)
    Trend range: [0.0, 360.0[ clockwise, from 0 (North) 
    Plunge: [-90.0, 90.0], negative value: upward axis, positive values: downward axis
    """

    def __init__(self, srcTrend, srcPlunge):

        assert 0.0 <= srcTrend < 360.0
        assert -90.0 <= srcPlunge <= 90.0

        self._trend = srcTrend
        self._plunge = srcPlunge

    @property
    def trend(self):

        return self._trend

    @property
    def plunge(self):

        return self._plunge

    def versor_3d(self):

        north_coord = cos(radians(self.plunge)) * cos(radians(self.trend))
        east_coord = cos(radians(self.plunge)) * sin(radians(self.trend))
        down_coord = sin(radians(self.plunge))

        return CartesianVector3D(east_coord, north_coord, -down_coord)

    def as_downgeolaxis(self):

        trend, plunge = self.trend, self.plunge
        if plunge < 0.0:
            trend += 180.0
            if trend > 360.0:
                trend -= 360.0
            plunge = - plunge

        return GeolAxis(trend, plunge)

    def as_normalgeolplane(self):

        down_axis = self.as_downgeolaxis()

        dipdir = down_axis.trend + 180.0
        if dipdir >= 360.0:
            dipdir -= 360.0
        dipangle = 90.0 - down_axis.plunge

        return GeolPlane(dipdir, dipangle)


class CartesianPlane(object):
    """
    Cartesian plane, expressed by equation:
    ax + by + cz + d = 0

    """

    def __init__(self, a=None, b=None, c=None, d=None):
        self._a = a
        self._b = b
        self._c = c
        self._d = d

    @property
    def a(self):
        return self._a

    @property
    def b(self):
        return self._b

    @property
    def c(self):
        return self._c

    @property
    def d(self):
        return self._d

    @classmethod
    def from_points(cls, pt1, pt2, pt3):
        matr_a = np.array([[pt1.p_y, pt1.p_z, 1],
                           [pt2.p_y, pt2.p_z, 1],
                           [pt3.p_y, pt3.p_z, 1]])

        matr_b = - np.array([[pt1.p_x, pt1.p_z, 1],
                             [pt2.p_x, pt2.p_z, 1],
                             [pt3.p_x, pt3.p_z, 1]])

        matr_c = np.array([[pt1.p_x, pt1.p_y, 1],
                           [pt2.p_x, pt2.p_y, 1],
                           [pt3.p_x, pt3.p_y, 1]])

        matr_d = - np.array([[pt1.p_x, pt1.p_y, pt1.p_z],
                             [pt2.p_x, pt2.p_y, pt2.p_z],
                             [pt3.p_x, pt3.p_y, pt3.p_z]])

        return cls(np.linalg.det(matr_a),
                   np.linalg.det(matr_b),
                   np.linalg.det(matr_c),
                   np.linalg.det(matr_d))

    def normal_versor3d(self):
        """
        return the normal versor to the cartesian plane
        """

        return CartesianVector3D(self.a, self.b, self.c).as_versor3d()

    def as_geolplane_and_point_3d(self):
        """
        converts a cartesian plane into a geological plane
        and a point lying in the plane (non-unique solution)
        """

        geol_plane = self.normal_versor3d().as_geolaxis().as_normalgeolplane()
        point = CartesianPoint3DT(point_solution(np.array([[self.a, self.b, self.c]]),
                                                 np.array([-self.d])))
        return geol_plane, point

    def intersection_versor3d(self, another):
        """
        return intersection versor for two intersecting planes
        """

        return self.normal_versor3d().vector_product(another.normal_versor3d()).as_versor3d()

    def intersection_point3dt(self, another):
        """
        return point on intersection line (obviously non-unique solution)
        for two planes
        """

        # find a point lying on the intersection line (this is a non-unique solution)    
        a = np.array([[self.a, self.b, self.c], [another.a, another.b, another.c]])
        b = np.array([-self.d, -another.d])
        x, y, z = point_solution(a, b)

        return CartesianPoint3DT(x, y, z)

    def set_point_inside(self, pt):
        return self.a * pt.p_x + self.b * pt.p_y + self.c * pt.p_z + self.d

    def angle_degr(self, another):
        angle_degr = self.normal_versor3d().angle_degr(another.normal_versor3d())

        assert angle_degr > 0.0

        if angle_degr > 90.0:
            angle_degr = 180.0 - angle_degr

        return angle_degr


def rhrstrike2dipdir(rhr_strk):

    return (rhr_strk + 90.0) % 360.0


class GeolPlane(object):
    """
    Structural plane, following geological conventions:
    dip direction and dip angle.
    
    """

    def __init__(self, srcAzimuth, srcDipAngle, isRHRStrike=False):
        """
        Class constructor

        @param  srcAzimuth:  Azimuth of the plane (RHR strike or dip direction).
        @type  srcAzimuth:  number or string convertible to float.
        @param  srcDipAngle:  Dip angle of the plane (0-90�).
        @type  srcDipAngle:  number or string convertible to float.
           
        @return:  GeolPlane.
    
        """

        if isRHRStrike:
            self._dipdir = rhrstrike2dipdir(srcAzimuth)
        else:
            self._dipdir = srcAzimuth % 360.0
        self._dipangle = srcDipAngle

    @property
    def dipdir(self):
        return self._dipdir

    @property
    def dipangle(self):
        return self._dipangle

    def as_normalgeolaxis(self):
        trend = self.dipdir + 180.0
        if trend >= 360.0:
            trend -= 360.0
        plunge = 90.0 - self.dipangle

        return GeolAxis(trend, plunge)

    def plane_x_coeff(self):
        """
        Calculate the slope of a given plane along the x direction.
        The plane orientation  is expressed following the geological convention. 
               
        @return:  slope - float.    
        """
        return - sin(radians(self.dipdir)) * tan(radians(self.dipangle))

    def plane_y_coeff(self):
        """
        Calculate the slope of a given plane along the y direction.
        The plane orientation  is expressed following the geological convention. 
               
        @return:  slope - float.     
        """
        return - cos(radians(self.dipdir)) * tan(radians(self.dipangle))

    def plane_from_geo(self, or_Pt):
        """
        Closure that embodies the analytical formula for a given, non-vertical plane.
        This closure is used to calculate the z value from given horizontal coordinates (x, y).
    
        @param  or_Pt:  CartesianPoint3DT instance expressing a location point contained by the plane.
        @type  or_Pt:  CartesianPoint3DT.
        
        @return:  lambda (closure) expressing an analytical formula for deriving z given x and y values.
        """

        x0 = or_Pt.p_x
        y0 = or_Pt.p_y
        z0 = or_Pt.p_z

        # slope of the line parallel to the x axis and contained by the plane
        a = self.plane_x_coeff()

        # slope of the line parallel to the y axis and contained by the plane
        b = self.plane_y_coeff()

        return lambda x, y: a * (x - x0) + b * (y - y0) + z0

    def as_cartesplane(self, point):
        normal_versor = self.as_normalgeolaxis().as_downgeolaxis().versor_3d()
        a, b, c = normal_versor.x, normal_versor.y, normal_versor.z
        d = - (a * point.p_x + b * point.p_y + c * point.p_z)
        return CartesianPlane(a, b, c, d)


def eq_xy_pair(xy_pair_1, xy_pair_2):
    if xy_pair_1[0] == xy_pair_2[0] and xy_pair_1[1] == xy_pair_2[1]:
        return True

    return False


def remove_equal_consecutive_xypairs(xy_list):
    out_xy_list = [xy_list[0]]

    for n in range(1, len(xy_list)):
        if not eq_xy_pair(xy_list[n], out_xy_list[-1]):
            out_xy_list.append(xy_list[n])

    return out_xy_list


def xytuple_list_to_Line2D(xy_list):
    return CartesianLine2DT([CartesianPoint2DT(x, y) for (x, y) in xy_list])


def xytuple_list2_to_MultiLine2D(xytuple_list2):
    # input is a list of list of (x,y) values

    assert len(xytuple_list2) > 0
    lines_list = []
    for xy_list in xytuple_list2:
        assert len(xy_list) > 0
        lines_list.append(xytuple_list_to_Line2D(xy_list))

    return CartesianMultiLine2DT(lines_list)


def list2_to_list(list2):
    """
    input: a list of list of (x,y) tuples
    output: a list of (x,y) tuples
    """

    out_list = []
    for list1 in list2:
        for el in list1:
            out_list.append(el)

    return out_list


def list3_to_list(list3):
    """
    input: a list of list of (x,y) tuples
    output: a list of (x,y) tuples
    """

    out_list = []
    for list2 in list3:
        for list1 in list2:
            out_list += list1

    return out_list


def merge_lines(lines, progress_ids):
    """
    lines: a list of list of (x,y,z) tuples for multilines
    """

    sorted_line_list = [line for (_, line) in sorted(zip(progress_ids, lines))]

    line_list = []
    for line in sorted_line_list:

        line_type, line_geometry = line

        if line_type == 'multiline':
            path_line = xytuple_list2_to_MultiLine2D(line_geometry).as_line2dt()
        elif line_type == 'line':
            path_line = xytuple_list_to_Line2D(line_geometry)
        else:
            continue

        line_list.append(path_line)  # now a list of Lines     

    # now the list of Lines is transformed into a single CartesianLine2DT
    return CartesianMultiLine2DT(line_list).as_line2dt().remove_coincident_successive_points()


class ArrCoord(object):
    """
    2D Array coordinates.
    Manages coordinates in the raster (array) space. 
    
    """

    def __init__(self, ival=0.0, jval=0.0):
        """
        @param  ival:  the i (-y) array coordinate of the point.
        @type  ival:  number or string convertible to float.
        @param  jval:  the j (x) array coordinate of the point.
        @type  jval:  number or string convertible to float.
               
        @return:  self.                
        """
        self._i = float(ival)
        self._j = float(jval)

    def g_i(self):
        """
        Get i (row) coordinate value.
        
        @return:  the i (-y) array coordinate of the point - float.
        """
        return self._i

    def s_i(self, ival):
        """
        Set i (row) coordinate value.
        
        @param  ival:  the i (-y) array coordinate of the point.
        @type  ival:  number or string convertible to float.
        
        @return:  self.        
        """
        self._i = float(ival)

    # set property for i
    i = property(g_i, s_i)

    def g_j(self):
        """
        Get j (column) coordinate value.
        
        @return:  the j (x) array coordinate of the point - float.
        """
        return self._j

    def s_j(self, jval):
        """
        Set j (column) coordinate value.
        
        @param  jval:  the j (x) array coordinate of the point.
        @type  jval:  number or string convertible to float.
        
        @return:  self.         
        """
        self._j = jval

    j = property(g_j, s_j)

    def grid2geogcoord(self, currGeoGrid):
        currPt_geogr_y = currGeoGrid.domain.trcorner.p_y - self.i * currGeoGrid.cellsize_y
        currPt_geogr_x = currGeoGrid.domain.llcorner.p_x + self.j * currGeoGrid.cellsize_x

        return CartesianPoint3DT(currPt_geogr_x, currPt_geogr_y)


class RectangularDomain(object):
    """
    Rectangular spatial domain class.
    
    """

    def __init__(self, pt_llc=None, pt_trc=None):
        """
        Class constructor.
        
        @param  pt_llc:  lower-left corner of the domain.
        @type  pt_llc:  CartesianPoint3DT.
        @param  pt_trc:  top-right corner of the domain.
        @type  pt_trc:  CartesianPoint3DT.
                        
        @return:  RectangularDomain instance.
        """
        self._llcorner = pt_llc
        self._trcorner = pt_trc

    @property
    def llcorner(self):
        """
        Get lower-left corner of the spatial domain.
        
        @return:  lower-left corner of the spatial domain - CartesianPoint3DT.
        """
        return self._llcorner

    @property
    def trcorner(self):
        """
        Get top-right corner of the spatial domain.
        
        @return:  top-right corner of the spatial domain - CartesianPoint3DT.
        """
        return self._trcorner

    @property
    def xrange(self):
        """
        Get x range of spatial domain.
        
        @return:  x range - float.
        """
        return self.trcorner.p_x - self.llcorner.p_x

    @property
    def yrange(self):
        """
        Get y range of spatial domain.
        
        @return:  y range - float.
        """
        return self.trcorner.p_y - self.llcorner.p_y

    @property
    def zrange(self):
        """
        Get z range of spatial domain.
        
        @return:  z range - float.
        """
        return self.trcorner.p_z - self.llcorner.p_z

    @property
    def horiz_area(self):
        """
        Get horizontal area of spatial domain.
        
        @return:  area - float.
        """
        return self.xrange * self.yrange


class Grid(object):
    """
    Grid class.
    Stores and manages the most of data and processing.
    
    """

    def __init__(self, source_filename=None, grid_params=None, grid_data=None):
        """
        Grid class constructor.
        
        @param  source_filename:  name of file from which data and geo-parameters derive.
        @type  source_filename:  string.
        @param  grid_params:  the geo-parameters of the grid.
        @type  grid_params:  class GDALParameters.
        @param  grid_data:  the array storing the data.
        @type  grid_data:  2D np.array.
               
        @return:  self.
        """
        self._sourcename = source_filename

        if grid_params is not None:
            pt_llc = grid_params.llcorner
            pt_trc = grid_params.trcorner
        else:
            pt_llc = None
            pt_trc = None

        self._grid_domain = RectangularDomain(pt_llc, pt_trc)

        if grid_data is not None:
            self._grid_data = grid_data.copy()
        else:
            self._grid_data = None

    def s_domain(self, domain):
        """
        Set spatial domain.
        
        @param  domain:  Spatial domain to be attributed to the current Grid instance.
        @type  domain:  class RectangularDomain.
        
        @return: self
        """

        del self._grid_domain
        self._grid_domain = copy.deepcopy(domain)

    def g_domain(self):
        """
        Get spatial domain.
        
        @return: the spatial domain of the current Grid instance - class RectangularDomain.
        """

        return self._grid_domain

    def d_domain(self):
        """
        Delete current spatial domain of the Grid instance.
        
        @return: self
        """

        del self._grid_domain

    # set property for spatial domain
    domain = property(g_domain, s_domain, d_domain)

    def s_grid_data(self, data_array):
        """
        Set grid data array.
        
        @param data_array: numpy.array of data values.
        @param type: 2D numpy.array.
        
        @return: self.
        """

        if self._grid_data is not None:
            del self._grid_data

        self._grid_data = data_array.copy()

    def g_grid_data(self):
        """
        Get grid data array.
        
        @return: 2D numpy.array.
        """

        return self._grid_data

    def d_grid_data(self):
        """
        Delete grid data array.
        
        @return: self.
        """

        del self._grid_data

    data = property(g_grid_data, s_grid_data, d_grid_data)

    def grid_extent(self):
        """
        Return the xmin, xmax and ymin, ymax values as a dictionary
        """

        return dict(xmin=self.domain.llcorner.p_x,
                    xmax=self.domain.trcorner.p_x,
                    ymin=self.domain.llcorner.p_y,
                    ymax=self.domain.trcorner.p_y)

    @property
    def xmin(self):

        return self.grid_extent()['xmin']

    @property
    def xmax(self):

        return self.grid_extent()['xmax']

    @property
    def ymin(self):

        return self.grid_extent()['ymin']

    @property
    def ymax(self):

        return self.grid_extent()['ymax']

    @property
    def row_num(self):
        """
        Get row number of the grid domain.     
        
        @return: number of rows of data array - int. 
        """

        return np.shape(self.data)[0]

    @property
    def col_num(self):
        """
        Get column number of the grid domain.
        
        @return: number of columns of data array - int. 
        """

        return np.shape(self.data)[1]

    @property
    def cellsize_x(self):
        """
        Get the cell size of the grid in the x direction.
        
        @return: cell size in the x (j) direction - float.
        """

        return self.domain.xrange / float(self.col_num)

    @property
    def cellsize_y(self):
        """
        Get the cell size of the grid in the y direction.
        
        @return: cell size in the y (-i) direction - float.
        """

        return self.domain.yrange / float(self.row_num)

    @property
    def cellsize_h(self):
        """
        Get the mean horizontal cell size.
        
        @return: mean horizontal cell size - float.
        """

        return (self.cellsize_x + self.cellsize_y) / 2.0

    def geog2array_coord(self, curr_Pt):
        """
        Converts from geographic to raster (array) coordinates.
        
        @param curr_Pt: point whose geographical coordinates will be converted to raster (array) ones.
        @type curr_Pt: CartesianPoint3DT.
        
        @return: point coordinates in raster (array) frame - class ArrCoord.
        """
        currArrCoord_grid_i = (self.domain.trcorner.p_y - curr_Pt.p_y) / self.cellsize_y
        currArrCoord_grid_j = (curr_Pt.p_x - self.domain.llcorner.p_x) / self.cellsize_x

        return ArrCoord(currArrCoord_grid_i, currArrCoord_grid_j)

    def x(self):
        """
        Creates an array storing the geographical coordinates of the cell centers along the x axis.
        Direction is from left to right.
        
        @return: numpy.array, shape: 1 x col_num.
        """

        x_values = self.domain.llcorner.p_x + self.cellsize_x * (0.5 + np.arange(self.col_num))

        return x_values[np.newaxis, :]

    def y(self):
        """
        Creates an array storing the geographical coordinates of the cell centers along the y axis.
        Direction is from top to bottom.
        
        @return: numpy.array, shape: row_num x 1.
        """

        y_values = self.domain.trcorner.p_y - self.cellsize_y * (0.5 + np.arange(self.row_num))

        return y_values[:, np.newaxis]

    def grad_forward_y(self):
        """
        Return an array representing the forward gradient in the y direction (top-wards), with values scaled by cell size.
        
        @return: numpy.array, same shape as current Grid instance        
        """

        gf = np.zeros(np.shape(self.data)) * np.NaN
        gf[1:, :] = self.data[:-1, :] - self.data[1:, :]

        return gf / float(self.cellsize_y)

    def grad_forward_x(self):
        """
        Return an array representing the forward gradient in the x direction (right-wards), with values scaled by cell size.

        @return: numpy.array, same shape as current Grid instance
        """

        gf = np.zeros(np.shape(self.data), ) * np.NaN
        gf[:, :-1] = self.data[:, 1:] - self.data[:, :-1]

        return gf / float(self.cellsize_x)

    def interpolate_bilinear(self, curr_Pt_array_coord):
        """
        Interpolate the z value at a point, given its array coordinates.
        Interpolation method: bilinear.
        
        @param curr_Pt_array_coord: array coordinates of the point for which the interpolation will be made.
        @type curr_Pt_array_coord: class ArrCoord.
        
        @return: interpolated z value - float.
        """

        currPt_cellcenter_i = curr_Pt_array_coord.i - 0.5
        currPt_cellcenter_j = curr_Pt_array_coord.j - 0.5

        assert currPt_cellcenter_i > 0, currPt_cellcenter_j > 0

        grid_val_00 = self.data[int(floor(currPt_cellcenter_i)), int(floor(currPt_cellcenter_j))]
        grid_val_01 = self.data[int(floor(currPt_cellcenter_i)), int(ceil(currPt_cellcenter_j))]
        grid_val_10 = self.data[int(ceil(currPt_cellcenter_i)), int(floor(currPt_cellcenter_j))]
        grid_val_11 = self.data[int(ceil(currPt_cellcenter_i)), int(ceil(currPt_cellcenter_j))]

        delta_i = currPt_cellcenter_i - floor(currPt_cellcenter_i)
        delta_j = currPt_cellcenter_j - floor(currPt_cellcenter_j)

        grid_val_y0 = grid_val_00 + (grid_val_10 - grid_val_00) * delta_i
        grid_val_y1 = grid_val_01 + (grid_val_11 - grid_val_01) * delta_i

        grid_val_interp = grid_val_y0 + (grid_val_y1 - grid_val_y0) * delta_j

        return grid_val_interp

    def intersection_with_surface(self, surf_type, srcPt, srcPlaneAttitude):
        """
        Calculates the intersections (as points) between DEM (the self object) and an analytical surface.
        Currently it works only with planes.
        
        @param surf_type: type of considered surface (i.e., plane, the only case implemented at present).
        @type surf_type: String.
        @param srcPt: point, expressed in geographical coordinates, that the plane must contain.
        @type srcPt: CartesianPoint3DT.
        @param srcPlaneAttitude: orientation of the surface (currently only planes).
        @type srcPlaneAttitude: class GeolPlane.
        
        @return: tuple of four arrays
        """

        if surf_type == 'plane':

            # closures to compute the geographic coordinates (in x- and y-) of a cell center
            # the grid coordinates of the cell center are expressed by i and j 
            grid_coord_to_geogr_coord_x_closure = lambda j: self.domain.llcorner.p_x + self.cellsize_x * (0.5 + j)
            grid_coord_to_geogr_coord_y_closure = lambda i: self.domain.trcorner.p_y - self.cellsize_y * (0.5 + i)

            # arrays storing the geographical coordinates of the cell centers along the x- and y- axes
            cell_center_x_array = self.x()
            cell_center_y_array = self.y()

            ycoords_x, xcoords_y = np.broadcast_arrays(cell_center_x_array, cell_center_y_array)

            #### x-axis direction intersections

            # 2D array of DEM segment parameters                         
            x_dem_m = self.grad_forward_x()
            x_dem_q = self.data - cell_center_x_array * x_dem_m

            # closure for the planar surface that, given (x,y), will be used to derive z
            plane_z_closure = srcPlaneAttitude.plane_from_geo(srcPt)

            # 2D array of plane segment parameters
            x_plane_m = srcPlaneAttitude.plane_x_coeff()
            x_plane_q = array_from_function(self.row_num(), 1, lambda j: 0, grid_coord_to_geogr_coord_y_closure,
                                            plane_z_closure)

            # 2D array that defines denominator for intersections between local segments
            x_inters_denomin = np.where(x_dem_m != x_plane_m, x_dem_m - x_plane_m, np.NaN)

            coincident_x = np.where(x_dem_q != x_plane_q, np.NaN, ycoords_x)

            xcoords_x = np.where(x_dem_m != x_plane_m, (x_plane_q - x_dem_q) / x_inters_denomin, coincident_x)
            xcoords_x = np.where(xcoords_x < ycoords_x, np.NaN, xcoords_x)
            xcoords_x = np.where(xcoords_x >= ycoords_x + self.cellsize_x, np.NaN, xcoords_x)

            #### y-axis direction intersections

            # 2D array of DEM segment parameters  
            y_dem_m = self.grad_forward_y()
            y_dem_q = self.data - cell_center_y_array * y_dem_m

            # 2D array of plane segment parameters
            y_plane_m = srcPlaneAttitude.plane_y_coeff()
            y_plane_q = array_from_function(1, self.col_num, grid_coord_to_geogr_coord_x_closure, lambda i: 0,
                                            plane_z_closure)

            # 2D array that defines denominator for intersections between local segments
            y_inters_denomin = np.where(y_dem_m != y_plane_m, y_dem_m - y_plane_m, np.NaN)
            coincident_y = np.where(y_dem_q != y_plane_q, np.NaN, xcoords_y)

            ycoords_y = np.where(y_dem_m != y_plane_m, (y_plane_q - y_dem_q) / y_inters_denomin, coincident_y)

            # filter out cases where intersection is outside cell range
            ycoords_y = np.where(ycoords_y < xcoords_y, np.NaN, ycoords_y)
            ycoords_y = np.where(ycoords_y >= xcoords_y + self.cellsize_y, np.NaN, ycoords_y)

            for i in xrange(xcoords_x.shape[0]):
                for j in xrange(xcoords_x.shape[1]):
                    if abs(xcoords_x[i, j] - ycoords_x[i, j]) < MINIMUM_SEPARATION_THRESHOLD and abs(
                                    ycoords_y[i, j] - xcoords_y[i, j]) < MINIMUM_SEPARATION_THRESHOLD:
                        ycoords_y[i, j] = np.NaN

            return xcoords_x, xcoords_y, ycoords_x, ycoords_y


class TriangBeam(object):
    """
    represents a 'fascio', a 2D semi-infinite geometrical object,
    defined by an apex (a CartesianPoint3DT object) and two semi-infinite segments, originating
    from the apex and defined by two versors (CartesianVector3D objects, with init name 'versor_1' and 'versor_2').
    Its maximum width 
    """

    def __init__(self, apex_pt3d, vector_1, vector_2):
        """
        assert almost_zero(versor_1.length() - 1.0)
        assert almost_zero(versor_2.length() - 1.0)
        """

        self._apex = apex_pt3d
        self._versor_1 = vector_1.as_versor3d()
        self._versor_2 = vector_2.as_versor3d()

    def fangle_degr(self):
        """
        angle 'sotteso' by the 'fascio'
        in the 0 - 180 degrees range
        """

        return self._versor_1.angle_degr(self._versor_2)

    def point_fangles_degr(self, pt_3d):
        """
        angles
        """

        vector_pt = CartesianSegment3DT(self._apex, pt_3d).as_vector3d()

        angle_side_1 = self._versor_1.angle_degr(vector_pt)
        angle_side_2 = self._versor_2.angle_degr(vector_pt)

        return angle_side_1, angle_side_2

    def is_within_fascio(self, pt_3d):
        apertura = self.fangle_degr()

        assert apertura < 180.0

        ang1, ang2 = self.point_fangles_degr(pt_3d)
        angle_sum = ang1 + ang2

        return almost_zero(apertura - angle_sum)


class CartesianTriangle(object):
    def __init__(self, pt_3d_1, pt_3d_2, pt_3d_3):

        self._pt_1 = pt_3d_1
        self._pt_2 = pt_3d_2
        self._pt_3 = pt_3d_3

    def is_pt_within(self, pt_3d):

        def versor3d(pt_1, pt_2):

            return CartesianSegment3DT(pt_1, pt_2).as_vector3d().as_versor3d()

        def is_pt_in_fascio(pt_1, pt_2, pt_3):

            apex = pt_1
            versor_1 = versor3d(pt_1, pt_2)
            versor_2 = versor3d(pt_1, pt_3)

            fascio = TriangBeam(apex, versor_1, versor_2)
            if not fascio.is_within_fascio(pt_3d):
                return False
            else:
                return True

        if not (is_pt_in_fascio(self._pt_1, self._pt_2, self._pt_3) and
                    is_pt_in_fascio(self._pt_2, self._pt_1, self._pt_3)):
            return False
        else:
            return True


class AnalyticGeosurface(object):
    def __init__(self, analytical_params, geogr_params, deform_params):

        self.analytical_params = analytical_params
        self.geographical_params = geogr_params
        self.deformational_params = deform_params

        # extract array params
        self.anal_param_values = self.get_analytical_param_values()
        array_range, array_size, formula = self.anal_param_values
        a_min, a_max, b_min, b_max = array_range
        a_range, b_range = a_max - a_min, b_max - b_min

        # calculate array from formula    
        try:
            self.X, self.Y, self.Z = formula_to_grid(array_range, array_size, formula)
        except AnaliticSurfaceCalcException, msg:
            raise AnaliticSurfaceCalcException, msg

        # calculate geographic transformations to surface
        self.geographical_values = self.get_geographical_param_values()
        (geog_x_min, geog_y_min), (area_height, area_width), area_rot_ang_deg = self.geographical_values
        geog_scale_matr = geographic_scale_matrix(a_range, b_range, area_height, area_width)
        geogr_rot_matrix = geographic_rotation_matrix(area_rot_ang_deg)

        self.geographic_transformation_matrix = np.dot(geogr_rot_matrix, geog_scale_matr)

        self.geographic_offset_matrix = geographic_offset(self.geographic_transformation_matrix,
                                                          np.array([a_min, b_min, 0.0]),
                                                          np.array([geog_x_min, geog_y_min, 0.0]))

        # apply total transformations to grid points 
        self.deformations = deformation_matrices(self.deformational_params)

    def geosurface_center(self):

        array_range, _, _ = self.anal_param_values
        a_min, a_max, b_min, b_max = array_range

        x = (a_min + a_max) / 2.0
        y = (b_min + b_max) / 2.0
        z = (min(self.Z) + max(self.Z)) / 2.0

        return self.transform_loc(x, y, z)

    def geosurface_XYZ(self):

        geosurface_X = []
        geosurface_Y = []
        geosurface_Z = []

        for x, y, z in zip(self.X, self.Y, self.Z):
            pt = self.transform_loc(x, y, z)
            geosurface_X.append(pt[0])
            geosurface_Y.append(pt[1])
            geosurface_Z.append(pt[2])

        return geosurface_X, geosurface_Y, geosurface_Z

    def get_analytical_param_values(self):

        try:
            a_min = float(self.analytical_params['a min'])
            a_max = float(self.analytical_params['a max'])
            grid_cols = int(self.analytical_params['grid cols'])

            b_min = float(self.analytical_params['b min'])
            b_max = float(self.analytical_params['b max'])
            grid_rows = int(self.analytical_params['grid rows'])

            formula = str(self.analytical_params['formula'])
        except:
            raise AnaliticSurfaceIOException, "Analytical value error"

        if a_min >= a_max or b_min >= b_max:
            raise AnaliticSurfaceIOException, "Input a and b value error"

        if grid_cols <= 0 or grid_rows <= 0:
            raise AnaliticSurfaceIOException, "Grid column/row value error"

        if formula == '':
            raise AnaliticSurfaceIOException, "Input analytical formula error"

        return (a_min, a_max, b_min, b_max), (grid_rows, grid_cols), formula

    def get_geographical_param_values(self):

        try:
            geog_x_min = float(self.geographical_params['geog x min'])
            geog_y_min = float(self.geographical_params['geog y min'])
            grid_height = float(self.geographical_params['grid height'])
            grid_width = float(self.geographical_params['grid width'])
            grid_rot_angle_degr = float(self.geographical_params['grid rot angle degr'])
        except:
            raise AnaliticSurfaceIOException, "Input geographic value error"

        return (geog_x_min, geog_y_min), (grid_height, grid_width), grid_rot_angle_degr

    def transform_loc(self, x, y, z):

        pt = np.dot(self.geographic_transformation_matrix, np.array([x, y, z])) + self.geographic_offset_matrix
        for deformation in self.deformations:
            if deformation['increment'] == 'additive':
                pt = pt + deformation['matrix']
            elif deformation['increment'] == 'multiplicative':
                pt = pt - deformation['shift_pt']
                pt = np.dot(deformation['matrix'], pt)
                pt = pt + deformation['shift_pt']
        return pt


def geographic_scale_matrix(a_range, b_range, grid_height, grid_width):
    assert a_range > 0.0
    assert b_range > 0.0
    assert grid_height > 0.0
    assert grid_width > 0.0

    sx = grid_width / a_range
    sy = grid_height / b_range
    sz = 1

    return np.array([(sx, 0.0, 0.0), (0.0, sy, 0.0), (0.0, 0.0, sz)])


def geographic_rotation_matrix(grid_rot_angle_degr):
    grid_rot_angle_rad = radians(grid_rot_angle_degr)
    sin_rot_angle = sin(grid_rot_angle_rad)
    cos_rot_angle = cos(grid_rot_angle_rad)

    return np.array([(cos_rot_angle, -sin_rot_angle, 0.0),
                     (sin_rot_angle, cos_rot_angle, 0.0),
                     (0.0, 0.0, 1.0)])


def geographic_offset(transformation_matrix, llc_point_matr, llc_point_geog):
    return llc_point_geog - np.dot(transformation_matrix, llc_point_matr)


def rotation_matrix(rot_axis_trend, rot_axis_plunge, rot_angle):
    phi = radians(rot_angle)

    rotation_versor = GeolAxis(rot_axis_trend, rot_axis_plunge).versor_3d()

    l = rotation_versor.x
    m = rotation_versor.y
    n = rotation_versor.z

    cos_phi = cos(phi)
    sin_phi = sin(phi)

    a11 = cos_phi + ((l * l) * (1 - cos_phi))
    a12 = ((l * m) * (1 - cos_phi)) - (n * sin_phi)
    a13 = ((l * n) * (1 - cos_phi)) + (m * sin_phi)

    a21 = ((l * m) * (1 - cos_phi)) + (n * sin_phi)
    a22 = cos_phi + ((m * m) * (1 - cos_phi))
    a23 = ((m * n) * (1 - cos_phi)) - (l * sin_phi)

    a31 = ((l * n) * (1 - cos_phi)) - (m * sin_phi)
    a32 = ((m * n) * (1 - cos_phi)) + (l * sin_phi)
    a33 = cos_phi + ((n * n) * (1 - cos_phi))

    return np.array([(a11, a12, a13),
                     (a21, a22, a23),
                     (a31, a32, a33)])


def scaling_matrix(scale_factor_x, scale_factor_y, scale_factor_z):
    return np.array([(scale_factor_x, 0.0, 0.0),
                     (0.0, scale_factor_y, 0.0),
                     (0.0, 0.0, scale_factor_z)])


def simple_shear_horiz_matrix(phi_angle_degr, alpha_angle_degr):
    phi_angle_rad = radians(phi_angle_degr)
    alpha_angle_rad = radians(alpha_angle_degr)

    gamma = tan(phi_angle_rad)
    sin_a = sin(alpha_angle_rad)
    cos_a = cos(alpha_angle_rad)

    return np.array([(1.0 - gamma * sin_a * cos_a, gamma * cos_a * cos_a, 0.0),
                     (-gamma * sin_a * sin_a, 1.0 + gamma * sin_a * cos_a, 0.0),
                     (0.0, 0.0, 1.0)])


def simple_shear_vert_matrix(phi_angle_degr, alpha_angle_degr):
    phi_angle_rad = radians(phi_angle_degr)
    alpha_angle_rad = radians(alpha_angle_degr)

    gamma = tan(phi_angle_rad)
    sin_a = sin(alpha_angle_rad)
    cos_a = cos(alpha_angle_rad)

    return np.array([(1.0, 0.0, gamma * cos_a),
                     (0.0, 1.0, gamma * sin_a),
                     (0.0, 0.0, 1.0)])


def deformation_matrices(deform_params):
    deformation_matrices = []

    for deform_param in deform_params:
        if deform_param['type'] == 'displacement':
            displ_x = deform_param['parameters']['delta_x']
            displ_y = deform_param['parameters']['delta_y']
            displ_z = deform_param['parameters']['delta_z']
            deformation = {'increment': 'additive',
                           'matrix': np.array([displ_x, displ_y, displ_z])}
        elif deform_param['type'] == 'rotation':
            rot_matr = rotation_matrix(deform_param['parameters']['rotation axis trend'],
                                       deform_param['parameters']['rotation axis plunge'],
                                       deform_param['parameters']['rotation angle'])
            deformation = {'increment': 'multiplicative',
                           'matrix': rot_matr,
                           'shift_pt': np.array([deform_param['parameters']['center x'],
                                                 deform_param['parameters']['center y'],
                                                 deform_param['parameters']['center z']])}
        elif deform_param['type'] == 'scaling':
            scal_matr = scaling_matrix(deform_param['parameters']['x factor'],
                                       deform_param['parameters']['y factor'],
                                       deform_param['parameters']['z factor'])
            deformation = {'increment': 'multiplicative',
                           'matrix': scal_matr,
                           'shift_pt': np.array([deform_param['parameters']['center x'],
                                                 deform_param['parameters']['center y'],
                                                 deform_param['parameters']['center z']])}
        elif deform_param['type'] == 'simple shear - horizontal':
            simple_shear_horiz_matr = simple_shear_horiz_matrix(deform_param['parameters']['psi angle (degr.)'],
                                                                deform_param['parameters']['alpha angle (degr.)'])
            deformation = {'increment': 'multiplicative',
                           'matrix': simple_shear_horiz_matr,
                           'shift_pt': np.array([deform_param['parameters']['center x'],
                                                 deform_param['parameters']['center y'],
                                                 deform_param['parameters']['center z']])}
        elif deform_param['type'] == 'simple shear - vertical':
            simple_shear_vert_matr = simple_shear_vert_matrix(deform_param['parameters']['psi angle (degr.)'],
                                                              deform_param['parameters']['alpha angle (degr.)'])
            deformation = {'increment': 'multiplicative',
                           'matrix': simple_shear_vert_matr,
                           'shift_pt': np.array([deform_param['parameters']['center x'],
                                                 deform_param['parameters']['center y'],
                                                 deform_param['parameters']['center z']])}
        else:
            continue

        deformation_matrices.append(deformation)

    return deformation_matrices
