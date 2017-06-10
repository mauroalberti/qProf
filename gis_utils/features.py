
from math import sqrt, degrees
import numpy as np

from .qgs_tools import qgs_point_2d, project_qgs_point
from ..gsf.geometry import Vect, GPlane, Point

MIN_2D_SEPARATION_THRESHOLD = 1e-10


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

    def translate_with_vector(self, displacement_vector):

        return CartesianPoint2DT(self.p_x + displacement_vector.x, self.p_y + displacement_vector.y, self.p_t)

    def spat_coincident_with(self, another, tolerance=MIN_2D_SEPARATION_THRESHOLD):

        if self.spat_distance(another) > tolerance:
            return False
        else:
            return True

    def crs_project_2d(self, srcCrs, destCrs):

        qgis_pt = qgs_point_2d(self.p_x, self.p_y)
        destCrs_qgis_pt = project_qgs_point(qgis_pt, srcCrs, destCrs)

        return CartesianPoint2DT(destCrs_qgis_pt.x(), destCrs_qgis_pt.y(), self.p_t)


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

        return min([x for x in self.x_list if not np.isnan(x)])

    @property
    def x_max(self):

        return max([x for x in self.x_list if not np.isnan(x)])

    @property
    def y_min(self):

        return min([y for y in self.y_list if not np.isnan(y)])

    @property
    def y_max(self):

        return max([y for y in self.y_list if not np.isnan(y)])

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

        return Vect(self.end_pt.p_x - self.start_pt.p_x,
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

        return GPlane(dip_dir, 90.0).as_cartesplane(self.start_pt)

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


class CartesianLine3DT(object):
    # CartesianLine3DT is a list of Point objects

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

        return np.array(filter(lambda pt: not np.isnan(pt.p_z), self.pts))

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
    # srcPt: source Point
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

        return Point(x1 - l * k,
                     y1 - m * k,
                     z1 - n * k)


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
