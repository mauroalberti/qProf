
import functools
import abc
import itertools
import numbers
import random
from array import array

from pygsf.mathematics.vectors import *
from pygsf.geometries.shapes.space3d import *


class Shape2D(object, metaclass=abc.ABCMeta):

    @abc.abstractmethod
    def area(self):
        """Calculate shape area"""
        
    @abc.abstractmethod
    def length(self):
        """Calculate shape area"""

    '''
    @abc.abstractmethod
    def clone(self):
        """Create a clone of the shape"""
    '''


class Point2D(Shape2D):
    """
    Cartesian point.
    Dimensions: 2D
    """

    def __init__(self,
                 x: numbers.Real,
                 y: numbers.Real
                 ):
        """
        Construct a Point instance.

        :param x: point x coordinate.
        :type x: numbers.Real.
        :param y: point y coordinate.
        :type y: numbers.Real.
        """

        vals = [x, y]

        if any(map(lambda val: not isinstance(val, numbers.Real), vals)):
            raise Exception("Input values must be integer or float type")

        if not all(map(math.isfinite, vals)):
            raise Exception("Input values must be finite")

        self._x = float(x)
        self._y = float(y)

    def area(self):
        return 0.0

    def length(self):
        return 0.0

    @property
    def x(self) -> numbers.Real:
        """
        Return the x coordinate of the current point.

        :return: x coordinate.
        :rtype: numbers.Real

        Examples:
          >>> Point2D(4, 3).x
          4.0
          >>> Point2D(-0.39, 3).x
          -0.39
        """

        return self._x

    @property
    def y(self) -> numbers.Real:
        """
        Return the y coordinate of the current point.

        :return: y coordinate.
        :rtype: numbers.Real

        Examples:
          >>> Point2D(4, 3).y
          3.0
          >>> Point2D(-0.39, 17.42).y
          17.42
        """

        return self._y

    def __iter__(self):
        """
        Return the elements of a Point.

        :return:

        Examples;
          >>> x, y = Point2D(1,1)
          >>> x == 1
          True
          >>> y == 1
          True

        """

        return (i for i in self.a())

    def __repr__(self) -> str:

        return "Point2D({:.4f}, {:.4f})".format(self.x, self.y)

    def __eq__(self,
        another: 'Point2D'
    ) -> bool:
        """
        Return True if objects are equal.

        :param another: another point.
        :type another: Point.
        :raise: Exception.

        Example:
          >>> Point2D(1., 1.) == Point2D(1, 1)
          True
          >>> Point2D(1., 1.) == Point2D(1, 1)
          True
          >>> Point2D(1., 1.) == Point2D(1, -1)
          False
        """

        if not isinstance(another, Point2D):
            raise Exception("Another instance must be a Point")

        return all([
            self.x == another.x,
            self.y == another.y,
            ]
        )

    def __ne__(self,
        another: 'Point2D'
    ) -> bool:
        """
        Return False if objects are equal.

        Example:
          >>> Point2D(1., 1.) != Point2D(0., 0.)
          True
          >>> Point2D(1., 1.) != Point2D(1, 1)
          False
        """

        return not (self == another)

    def a(self) -> Tuple[numbers.Real, numbers.Real]:
        """
        Return the individual values of the point.

        :return: double array of x, y values

        Examples:
          >>> Point2D(4, 3).a()
          (4.0, 3.0)
        """

        return self.x, self.y

    def __add__(self, another: 'Point2D') -> 'Point2D':
        """
        Sum of two points.

        :param another: the point to add
        :type another: Point2D
        :return: the sum of the two points
        :rtype: Point2D
        :raise: Exception

        Example:
          >>> Point2D(1, 0) + Point2D(0, 1)
          Point(1.0000, 1.0000)
          >>> Point2D(1, 1) + Point2D(-1, -1)
          Point(0.0000, 0.0000)
        """

        check_type(another, "Second point", Point2D)

        x0, y0 = self
        x1, y1 = another

        return Point2D(
            x=x0+x1,
            y=y0+y1
        )

    def __sub__(self,
        another: 'Point2D'
    ) -> 'Point2D':
        """Subtract two points.

        :param another: the point to subtract
        :type another: Point2D
        :return: the difference between the two points
        :rtype: Point2D
        :raise: Exception

        Example:
          >>> Point2D(1., 1.) - Point2D(1., 1.)
          Point(0.0000, 0.0000)
          >>> Point2D(1., 1.) - Point2D(1., 1.)
          Point(0.0000, 0.0000)
        """

        check_type(another, "Second point", Point2D)

        x0, y0 = self
        x1, y1 = another

        return Point2D(
            x=x0 - x1,
            y=y0 - y1
        )

    def clone(self) -> 'Point2D':
        """
        Clone a point.

        :return: a new point.
        :rtype: Point.
        """

        return Point2D(*self.a())

    def toXY(self) -> Tuple[numbers.Real, numbers.Real]:
        """
        Returns the spatial components as a tuple of two values.

        :return: the spatial components (x, y).
        :rtype: a tuple of two floats.

        Examples:
          >>> Point2D(1, 0).toXY()
          (1.0, 0.0)
        """

        return self.x, self.y

    def toArray(self) -> np.ndarray:
        """
        Return a Numpy array representing the point values.

        :return: Numpy array

        Examples:
          >>> np.allclose(Point2D(1, 2).toArray(), np.array([ 1., 2.]))
          True
        """

        return np.asarray(self.toXY())

    def deltaX(self,
        another: 'Point2D'
    ) -> Optional[numbers.Real]:
        """
        Delta between x components of two Point Instances.

        :return: x coordinates difference value.
        :rtype: optional numbers.Real.
        :raise: Exception

        Examples:
          >>> Point2D(1, 2).deltaX(Point2D(4, 7))
          3.0
        """

        return another.x - self.x

    def deltaY(self,
        another: 'Point2D'
    ) -> Optional[numbers.Real]:
        """
        Delta between y components of two Point Instances.

        :return: y coordinates difference value.
        :rtype: optional numbers.Real.

        Examples:
          >>> Point2D(1, 2).deltaY(Point2D(4, 7))
          5.0
        """

        return another.y - self.y

    def dist2DWith(self,
                 another: 'Point2D'
                 ) -> numbers.Real:
        """
        Calculate horizontal (2D) distance between two points.
        TODO: consider case of polar CRS

        :param another: another Point instance.
        :type another: Point.
        :return: the 2D distance (when the two points have the same CRS).
        :rtype: numbers.Real.
        :raise: Exception.

        Examples:
          >>> Point2D(1., 1.).dist2DWith(Point2D(4., 5.))
          5.0
        """

        check_type(another, "Second point", Point2D)

        return math.sqrt((self.x - another.x) ** 2 + (self.y - another.y) ** 2)

    def scale(self,
        scale_factor: numbers.Real
    ) -> 'Point2D':
        """
        Create a scaled object.
        Note: it does not make sense for polar coordinates.
        TODO: manage polar coordinates cases OR deprecate and remove - after dependency check.

        Example;
          >>> Point2D(1, 0).scale(2.5)
          Point(2.5000, 0.0000)
          >>> Point2D(1, 0).scale(2.5)
          Point(2.5000, 0.0000)
        """

        x, y = self.x * scale_factor, self.y * scale_factor
        return Point2D(x, y)

    def invert(self) -> 'Point2D':
        """
        Create a new object with inverted direction.
        Note: it depends on scale method, that could be deprecated/removed.

        Examples:
          >>> Point2D(1, 1).invert()
          Point(-1.0000, -1.0000)
          >>> Point2D(2, -1).invert()
          Point(-2.0000, 1.0000)
        """

        return self.scale(-1)

    def is_coincident(self,
                      another: 'Point2D',
                      tolerance: numbers.Real = MIN_SEPARATION_THRESHOLD
                      ) -> bool:
        """
        Check spatial coincidence of two points, limiting to the horizontal (XY) plane.

        :param another: the point to compare.
        :type another: Point.
        :param tolerance: the maximum allowed distance between the two points.
        :type tolerance: numbers.Real.
        :return: whether the two points are coincident.
        :rtype: bool.
        :raise: Exception.

        Example:

        """

        check_type(another, "Second point", Point2D)

        return self.dist2DWith(another) <= tolerance

    def already_present(self,
                        pt_list: List['Point2D'],
                        tolerance: numbers.Real = MIN_SEPARATION_THRESHOLD
                        ) -> Optional[bool]:
        """
        Determines if a point is already in a given point list, using an optional distance separation,

        :param pt_list: list of points. May be empty.
        :type pt_list: List of Points.
        :param tolerance: optional maximum distance between near-coincident point pair.
        :type tolerance: numbers.Real.
        :return: True if already present, False otherwise.
        :rtype: optional boolean.
        """

        for pt in pt_list:
            if self.is_coincident(pt, tolerance=tolerance):
                return True
        return False

    def shift(self,
        sx: numbers.Real,
        sy: numbers.Real
    ) -> Optional['Point2D']:
        """
        Create a new object shifted by given amount from the self instance.

        Example:
          >>> Point2D(1, 1).shift(0.5, 1.)
          Point(1.5000, 2.0000)
          >>> Point2D(1, 2).shift(0.5, 1.)
          Point(1.5000, 3.0000)
       """

        return Point2D(self.x + sx, self.y + sy)

    @classmethod
    def random(cls,
        lower_boundary: float = -MAX_SCALAR_VALUE,
        upper_boundary: float =  MAX_SCALAR_VALUE
    ):
        """
        Creates a random point.

        :return: random point
        :rtype: Point2D
        """

        vals = [random.uniform(lower_boundary, upper_boundary) for _ in range(2)]
        return cls(*vals)


class Segment2D(Shape2D):
    """
    Segment2D is a geometric object defined by the straight line between
    two vertices.
    """

    def __init__(self, start_pt: Point2D, end_pt: Point2D):
        """
        Creates a segment instance provided the two points have the same CRS code.

        :param start_pt: the start point.
        :param end_pt: the end point.
        :return: the new segment instance if both points have the same georeferenced.
        :raises: CRSCodeException.
        """

        check_type(start_pt, "Start point", Point2D)
        check_type(end_pt, "End point", Point2D)

        if start_pt.dist2DWith(end_pt) == 0.0:
            raise Exception("Source points cannot be coincident")

        self._start_pt = start_pt.clone()
        self._end_pt = end_pt.clone()

    def __repr__(self) -> str:
        """
        Represents a Segment2D instance.

        :return: the Segment2D representation.
        :rtype: str.
        """

        return "Segment2D(start_pt={}, end_pt={})".format(
            self.start_pt,
            self.end_pt
        )

    @property
    def start_pt(self) -> Point2D:

        return self._start_pt

    @property
    def end_pt(self) -> Point2D:

        return self._end_pt

    def asPoints(self) -> List[Point2D]:
        """
        Return the segments as points.
        """

        return [self.start_pt, self.end_pt]

    def length(self) -> numbers.Real:
        """
        Returns the horizontal length of the segment.

        :return: the horizontal length of the segment.
        :rtype: numbers.Real.
        """

        return self.start_pt.dist2DWith(self.end_pt)

    def area(self):
        return 0.0

    def __iter__(self):
        """
        Return the elements of a Segment2D, i.e., start and end point.
        """

        return (i for i in [self.start_pt, self.end_pt])

    def clone(self) -> 'Segment2D':

        return Segment2D(self._start_pt, self._end_pt)

    def increasing_x(self) -> 'Segment2D':

        if self.end_pt.x < self.start_pt.x:
            return Segment2D(self.end_pt, self.start_pt)
        else:
            return self.clone()

    def x_range(self) -> Tuple[numbers.Real, numbers.Real]:

        if self.start_pt.x < self.end_pt.x:
            return self.start_pt.x, self.end_pt.x
        else:
            return self.end_pt.x, self.start_pt.x

    def y_range(self) -> Tuple[numbers.Real, numbers.Real]:

        if self.start_pt.y < self.end_pt.y:
            return self.start_pt.y, self.end_pt.y
        else:
            return self.end_pt.y, self.start_pt.y

    def delta_x(self) -> numbers.Real:

        return self.end_pt.x - self.start_pt.x

    def delta_y(self) -> numbers.Real:

        return self.end_pt.y - self.start_pt.y

    def segment_2d_m(self) -> Optional[numbers.Real]:

        denom = self.end_pt.x - self.start_pt.x

        if denom == 0.0:
            return None

        return (self.end_pt.y - self.start_pt.y) / denom

    def segment_2d_p(self) -> Optional[numbers.Real]:

        s2d_m = self.segment_2d_m()

        if s2d_m is None:
            return None

        return self.start_pt.y - s2d_m * self.start_pt.x

    def intersection_2d_pt(self,
                           another: 'Segment2D'
                           ) -> Optional[Point2D]:
        """

        :param another:
        :return:
        """

        check_type(another, "Second segment", Segment2D)

        #check_crs(self, another)

        s_len2d = self.length()
        a_len2d = another.length()

        if s_len2d == 0.0 or a_len2d == 0.0:
            return None

        if self.start_pt.x == self.end_pt.x:  # self segment parallel to y axis
            x0 = self.start_pt.x
            m1, p1 = another.segment_2d_m(), another.segment_2d_p()
            if m1 is None:
                return None
            y0 = m1 * x0 + p1
        elif another.start_pt.x == another.end_pt.x:  # another segment parallel to y axis
            x0 = another.start_pt.x
            m1, p1 = self.segment_2d_m(), self.segment_2d_p()
            if m1 is None:
                return None
            y0 = m1 * x0 + p1
        else:  # no segment parallel to y axis
            m0, p0 = self.segment_2d_m(), self.segment_2d_p()
            m1, p1 = another.segment_2d_m(), another.segment_2d_p()
            if m0 is None or m1 is None:
                return None
            x0 = (p1 - p0) / (m0 - m1)
            y0 = m0 * x0 + p0

        return Point2D(x0, y0)

    def contains_pt(self,
                    pt: Point2D
                    ) -> bool:
        """
        Checks whether a point is contained in a segment.

        :param pt: the point for which to check containement.
        :return: bool.
        :raise: Exception.

        Examples:
          >>> segment = Segment2D(Point2D(0, 0), Point2D(1, 0))
          >>> segment.contains_pt(Point2D(0, 0))
          True
          >>> segment.contains_pt(Point2D(1, 0))
          True
          >>> segment.contains_pt(Point2D(0.5, 0))
          True
          >>> segment.contains_pt(Point2D(0.5, 0.00001))
          False
          >>> segment.contains_pt(Point2D(1.00001, 0))
          False
          >>> segment.contains_pt(Point2D(0.000001, 0))
          True
          >>> segment.contains_pt(Point2D(-0.000001, 0))
          False
          >>> segment.contains_pt(Point2D(0.5, 1000))
          False
          >>> segment = Segment2D(Point2D(0, 0), Point2D(0, 1))
          >>> segment.contains_pt(Point2D(0, 0))
          True
          >>> segment.contains_pt(Point2D(0, 0.5))
          True
          >>> segment.contains_pt(Point2D(0, 1))
          True
          >>> segment.contains_pt(Point2D(0, 1.5))
          False
          >>> segment = Segment2D(Point2D(0, 0), Point2D(1, 1))
          >>> segment.contains_pt(Point2D(0.5, 0.5))
          True
          >>> segment.contains_pt(Point2D(1, 1))
          True
          >>> segment = Segment2D(Point2D(1, 2), Point2D(9, 8))
          >>> segment.contains_pt(segment.pointAt(0.745))
          True
          >>> segment.contains_pt(segment.pointAt(1.745))
          False
          >>> segment.contains_pt(segment.pointAt(-0.745))
          False
          >>> segment.contains_pt(segment.pointAt(0))
          True
        """

        check_type(pt, "Point", Point2D)

        segment_length = self.length()
        length_startpt_pt = self.start_pt.dist2DWith(pt)
        length_endpt_pt = self.end_pt.dist2DWith(pt)

        return areClose(
            a=segment_length,
            b=length_startpt_pt + length_endpt_pt
        )

    def fast_2d_contains_pt(self,
                            pt2d
                            ) -> bool:
        """
        Deprecated. Use 'contains_pt'.

        to work properly, this function requires that the pt lies on the line defined by the segment
        """

        range_x = self.x_range
        range_y = self.y_range

        if range_x()[0] <= pt2d.x <= range_x()[1] or \
                range_y()[0] <= pt2d.y <= range_y()[1]:
            return True
        else:
            return False

    def pointAt(self,
                scale_factor: numbers.Real
                ) -> Point2D:
        """
        Returns a point aligned with the segment
        and lying at given scale factor, where 1 is segment length
        ans 0 is segment start.

        :param scale_factor: the scale factor, where 1 is the segment length.
        :type scale_factor: numbers.Real
        :return: Point at scale factor
        :rtype: Point2D

        Examples:
          >>> s = Segment2D(Point2D(0,0), Point2D(1,0))
          >>> s.pointAt(0)
          Point(0.0000, 0.0000)
          >>> s.pointAt(0.5)
          Point(0.5000, 0.0000)
          >>> s.pointAt(1)
          Point(1.0000, 0.0000)
          >>> s.pointAt(-1)
          Point(-1.0000, 0.0000)
          >>> s.pointAt(-2)
          Point(-2.0000, 0.0000)
          >>> s.pointAt(2)
          Point(2.0000, 0.0000)
          >>> s = Segment2D(Point2D(0,0), Point2D(1,1))
          >>> s.pointAt(0.5)
          Point(0.5000, 0.5000)
          >>> s = Segment2D(Point2D(0,0), Point2D(4,0))
          >>> s.pointAt(7.5)
          Point(30.0000, 0.0000)
        """

        dx = self.delta_x() * scale_factor
        dy = self.delta_y() * scale_factor

        return Point2D(
            x=self.start_pt.x + dx,
            y=self.start_pt.y + dy
        )

    '''
    def pointProjection(self,
                        point: Point
                        ) -> Point:
        """
        Return the point projection on the segment.

        Examples:
          >>> s = Segment2D(start_pt=Point(0,0), end_pt=Point(1,0))
          >>> p = Point(0.5, 1)
          >>> s.pointProjection(p)
          Point(0.5000, 0.0000)
          >>> s = Segment2D(start_pt=Point(0,0), end_pt=Point(4,0))
          >>> p = Point(7.5, 19.2)
          >>> s.pointProjection(p)
          Point(7.5000, 0.0000)
        """

        check_type(point, "Input point", Point)

        #check_crs(self, point)

        other_segment = Segment2D(
            self.start_pt,
            point
        )

        scale_factor = self.vector().scalarProjection(other_segment.vector()) / self.length2D()
        return self.pointAt(scale_factor)
    '''

    '''
    def pointDistance(self,
                      point: Point
                      ) -> numbers.Real:
        """
        Returns the point distance to the segment.

        :param point: the point to calculate the distance with
        :type point: Point
        :return: the distance of the point to the segment
        :rtype: numbers.Real

        Examples:
          >>> s = Segment2D(Point(0, 0), Point(0, 0))
          >>> s.pointDistance(Point(-17.2, 0.0))
          17.2
          >>> s.pointDistance(Point(-17.2, 1.22))
          17.24321315764553
        """

        check_type(point, "Input point", Point)

        # check_crs(self, point)

        point_projection = self.pointProjection(point)

        return point.dist2DWith(point_projection)
    '''

    def pointS(self,
               point: Point2D
               ) -> Optional[numbers.Real]:
        """
        Calculates the optional distance of the point along the segment.
        A zero value is for a point coinciding with the start point.
        Returns None if the point is not contained in the segment.

        :param point: the point to calculate the optional distance in the segment.
        :type point: Point2D
        :return: the the optional distance of the point along the segment.
        """

        check_type(point, "Input point", Point2D)

        # check_crs(self, point)

        if not self.contains_pt(point):
            return None

        return self.start_pt.dist2DWith(point)

    def scale(self,
              scale_factor
              ) -> 'Segment2D':
        """
        Scale a segment by the given scale_factor.
        Start point does not change.

        :param scale_factor: the scale factor, where 1 is the segment length.
        :type scale_factor: numbers.Real
        :return: Point at scale factor
        :rtype: Point2D
        """

        end_pt = self.pointAt(scale_factor)

        return Segment2D(
            self.start_pt,
            end_pt)

    def as_vector(self) -> Vect:
        """
        Convert a segment to a vector.
        """

        return Vect(
            x=self.delta_x(),
            y=self.delta_y()
        )

    def densify2d_asLine(self, densify_distance) -> 'Line2D':
        """
        Densify a segment by adding additional points
        separated a distance equal to densify_distance.
        The result is no longer a Segment instance, instead it is a Line instance.
        :param densify_distance: float
        :return: Line
        """

        length2d = self.length()

        vect = self.as_vector()
        vers_2d = vect.versor_2d()
        generator_vector = vers_2d.scale(densify_distance)

        interpolated_line = Line2D(
            pts=[self.start_pt])

        n = 0
        while True:
            n += 1
            shift_vector = generator_vector.scale(n)
            new_pt = self.start_pt.shift(shift_vector.x, shift_vector.y)
            distance = self.start_pt.dist2DWith(new_pt)
            if distance >= length2d:
                break
            interpolated_line.add_pt(new_pt)
        interpolated_line.add_pt(self.end_pt)

        return interpolated_line

    def densify2d_asPts(self, densify_distance) -> List[Point2D]:

        return self.densify2d_asLine(densify_distance=densify_distance).pts()

    def densify2d_asSteps(self,
                          densify_distance: numbers.Real
                          ) -> array:
        """
        Defines the array storing the incremental lengths according to the provided densify distance.

        :param densify_distance: the step distance.
        :type densify_distance: numbers.Real.
        :return: array storing incremental steps, with the last step being equal to the segment length.
        :rtype: array.
        """

        if not isinstance(densify_distance, numbers.Real):
            raise Exception("Densify distance must be float or int")

        if not math.isfinite(densify_distance):
            raise Exception("Densify distance must be finite")

        if densify_distance <= 0.0:
            raise Exception("Densify distance must be positive")

        segment_length = self.length()

        s_list = []
        n = 0
        length = n * densify_distance

        while length < segment_length:
            s_list.append(length)
            n += 1
            length = n * densify_distance

        s_list.append(segment_length)

        return array('d', s_list)

    '''
    def densify2d_asPts(self,
                        densify_distance
                        ) -> List[Point]:
        """
        Densify a segment by adding additional points
        separated a distance equal to densify_distance.
        The result is no longer a Segment2D instance, instead it is a Line instance.

        :param densify_distance: the distance with which to densify the segment.
        :type densify_distance: numbers.Real.
        :return: the set of densified points.
        :rtype: List[Point].
        """

        if not isinstance(densify_distance, numbers.Real):
            raise Exception("Input densify distance must be float or integer")

        if not math.isfinite(densify_distance):
            raise Exception("Input densify distance must be finite")

        if densify_distance <= 0.0:
            raise Exception("Input densify distance must be positive")

        length2d = self.length2D()

        vect = self.vector()
        vers_2d = vect.versor2D()
        generator_vector = vers_2d.scale(densify_distance)

        pts = [self.start_pt]

        n = 0
        while True:
            n += 1
            new_pt = self.start_pt.shiftByVect(generator_vector.scale(n))
            distance = self.start_pt.dist2DWith(new_pt)
            if distance >= length2d:
                break
            pts.append(new_pt)

        pts.append(self.end_pt)

        return pts
    '''

    '''
    def densify2d_asLine(self,
                         densify_distance
                         ) -> 'Points':
        """
        Densify a segment by adding additional points
        separated a distance equal to densify_distance.
        The result is no longer a Segment2D instance, instead it is a Line instance.

        :param densify_distance: numbers.Real
        :return: Line
        """

        pts = self.densify2d_asPts(densify_distance=densify_distance)

        return Points(
            pts=pts)
    '''

    '''
    def vertical_plane(self) -> Optional[CPlane]:
        """
        Returns the vertical Cartesian plane containing the segment.

        :return: the vertical Cartesian plane containing the segment.
        :rtype: Optional[CPlane].
        """

        if self.length2D() == 0.0:
            return None

        # arbitrary point on the same vertical as end point

        section_final_pt_up = self.end_pt.shift(
            sx=0.0,
            sy=0.0,
            sz=1000.0)

        return CPlane.fromPoints(
            pt1=self.start_pt,
            pt2=self.end_pt,
            pt3=section_final_pt_up)
    '''

    def same_start(self,
                   another: 'Segment2D',
                   tol: numbers.Real = 1e-12
                   ) -> bool:
        """
        Check whether the two segments have the same start point.

        :param another: a segment to check for.
        :type another: Segment2D.
        :param tol: tolerance for distance between points.
        :type tol: numbers.Real.
        :return: whether the two segments have the same start point.
        :rtype: bool.

        Examples:
          >>> s1 = Segment2D(Point2D(0,0), Point2D(1,0))
          >>> s2 = Segment2D(Point2D(0,0), Point2D(0,1))
          >>> s1.same_start(s2)
          True
        """

        return self.start_pt.is_coincident(
            another=another.start_pt,
            tolerance=tol
        )

    def same_end(self,
                 another: 'Segment2D',
                 tol: numbers.Real = 1e-12
                 ) -> bool:
        """
        Check whether the two segments have the same end point.

        :param another: a segment to check for.
        :type another: Segment2D.
        :param tol: tolerance for distance between points.
        :type tol: numbers.Real.
        :return: whether the two segments have the same end point.
        :rtype: bool.

        Examples:
          >>> s1 = Segment2D(Point2D(0,0), Point2D(1,0))
          >>> s2 = Segment2D(Point2D(2,0), Point2D(1,0))
          >>> s1.same_end(s2)
          True
        """

        return self.end_pt.is_coincident(
            another=another.end_pt,
            tolerance=tol)

    def conn_to_other(self,
                      another: 'Segment2D',
                      tol: numbers.Real = 1e-12
                      ) -> bool:
        """
        Check whether the first segment is sequentially connected to the second one.

        :param another: a segment to check for.
        :type another: Segment2D.
        :param tol: tolerance for distance between points.
        :type tol: numbers.Real.
        :return: whether the first segment is sequentially connected to the second one.
        :rtype: bool.

        Examples:
          >>> s1 = Segment2D(Point2D(0,0), Point2D(1,0))
          >>> s2 = Segment2D(Point2D(1,0), Point2D(2,0))
          >>> s1.conn_to_other(s2)
          True
        """

        return self.end_pt.is_coincident(
            another=another.start_pt,
            tolerance=tol)

    def other_connected(self,
                        another: 'Segment2D',
                        tol: numbers.Real = 1e-12
                        ) -> bool:
        """
        Check whether the second segment is sequentially connected to the first one.

        :param another: a segment to check for.
        :type another: Segment2D.
        :param tol: tolerance for distance between points.
        :type tol: numbers.Real.
        :return: whether the second segment is sequentially connected to the first one.
        :rtype: bool.

        Examples:
          >>> s1 = Segment2D(Point2D(0,0), Point2D(1,0))
          >>> s2 = Segment2D(Point2D(-1,0), Point2D(0,0))
          >>> s1.other_connected(s2)
          True
        """

        return another.end_pt.is_coincident(
            another=self.start_pt,
            tolerance=tol)

    def segment_start_in(self,
                         another: 'Segment2D'
                         ) -> bool:
        """
        Check whether the second segment contains the first segment start point.

        :param another: a segment to check for.
        :type another: Segment2D.
        :return: whether the second segment contains the first segment start point.
        :rtype: bool.

        Examples:
          >>> s1 = Segment2D(Point2D(0,0), Point2D(1,0))
          >>> s2 = Segment2D(Point2D(-0.5,0), Point2D(0.5,0))
          >>> s1.segment_start_in(s2)
          True
          >>> s1 = Segment2D(Point2D(0,0), Point2D(1,1))
          >>> s1.segment_start_in(s2)
          True
          >>> s1 = Segment2D(Point2D(0,1), Point2D(1,1))
          >>> s1.segment_start_in(s2)
          False
          >>> s1 = Segment2D(Point2D(-1,-1), Point2D(1,1))
          >>> s1.segment_start_in(s2)
          False
        """

        return another.contains_pt(self.start_pt)

    def segment_end_in(self,
                       another: 'Segment2D'
                       ) -> bool:
        """
        Check whether the second segment contains the first segment end point.

        :param another: a segment to check for.
        :type another: Segment2D.
        :return: whether the second segment contains the first segment end point.
        :rtype: bool.

        Examples:
          >>> s1 = Segment2D(Point2D(0,0), Point2D(1,0))
          >>> s2 = Segment2D(Point2D(-0.5,0), Point2D(0.5,0))
          >>> s1.segment_end_in(s2)
          False
          >>> s1 = Segment2D(Point2D(0,0), Point2D(1,1))
          >>> s1.segment_end_in(s2)
          False
          >>> s1 = Segment2D(Point2D(0,1), Point2D(1,1))
          >>> s2 = Segment2D(Point2D(1,1), Point2D(0.5,0))
          >>> s1.segment_end_in(s2)
          True
          >>> s1 = Segment2D(Point2D(-1,-1), Point2D(1,1))
          >>> s2 = Segment2D(Point2D(0,2), Point2D(2,0))
          >>> s1.segment_end_in(s2)
          True
        """

        return another.contains_pt(self.end_pt)

    @classmethod
    def random(cls,
               lower_boundary: float = -MAX_SCALAR_VALUE,
               upper_boundary: float = MAX_SCALAR_VALUE):
        """
        Creates a random segment.

        :return: random segment
        :rtype: Segment2D
        """

        return cls(
            start_pt=Point2D.random(lower_boundary, upper_boundary),
            end_pt=Point2D.random(lower_boundary, upper_boundary)
        )

    def vertical_plane(self) -> Optional[CPlane3D]:
        """
        Returns the vertical Cartesian plane containing the segment.

        :return: the vertical Cartesian plane containing the segment.
        :rtype: Optional[CPlane].
        """

        if self.length() == 0.0:  # collapsed segment
            return None

        # arbitrary point on the same vertical as end point

        p1 = Point3D(
            x=self.start_pt.x,
            y=self.start_pt.y,
            z=0.0
        )
        p2 = Point3D(
            x=self.end_pt.x,
            y=self.end_pt.y,
            z=0.0
        )
        p3 = Point3D(
            x=self.end_pt.x,
            y=self.end_pt.y,
            z=1000.0
        )
        return CPlane3D.fromPoints(
            pt1=p1,
            pt2=p2,
            pt3=p3)

    def vector(self) -> Vect:

        return Vect(self.delta_x(),
                    self.delta_y(),
                    0
        )

class Line2D(Shape2D):
    """
    A list of Point objects.
    """

    def area(self):
        pass

    def length(self):
        pass

    def __init__(self,
        pts: Optional[List[Point2D]] = None
    ):
        """
        Creates the Line2D instance.

        :param pts: a list of points
        :return: a Line2D instance.
        """

        if pts is None:
            pts = []

        for pt in pts:
            if not isinstance(pt, Point2D):
                raise Exception("All input data must be point")

        self._x = array('d', [pt.x for pt in pts])
        self._y = array('d', [pt.y for pt in pts])

    @classmethod
    def fromArrays(cls,
        xs: array,
        ys: array
    ) -> 'Line2D':
        """
        Create a Line2D instance from a list of x and y values.

        Example:
          >>> Line2D.fromArrays(xs=array('d',[1,2,3]), ys=array('d', [3,4,5]))
          Line with 3 points: (1.0000, 3.0000) ... (3.0000, 5.0000)
          >>> Line2D.fromArrays(xs=array('d',[1,2,3]), ys=array('d', [3,4,5]))
          Line with 3 points: (1.0000, 3.0000) ... (3.0000, 5.0000)
        """

        if not isinstance(xs, array):
            raise Exception("X values have type {} instead of array".format(type(xs)))

        if not isinstance(ys, array):
            raise Exception("Y values have type {} instead of array".format(type(ys)))

        num_vals = len(xs)
        if len(ys) != num_vals:
            raise Exception("Y array has length {} while x array has length {}".format(len(ys), num_vals))

        self = cls()

        self._x = xs
        self._y = ys

        return self

    @classmethod
    def fromPointList(cls,
        pt_list: List[List[numbers.Real]]
    ) -> 'Line2D':
        """
        Create a Line2D instance from a list of x and y values.

        Example:
          >>> Line2D.fromPointList([[0, 0], [1, 0], [0, 1]])
          Line with 3 points: (0.0000, 0.0000) ... (0.0000, 1.0000)
        """

        pts = []
        for vals in pt_list:
            if len(vals) == 2:
                pt = Point2D(
                    x=vals[0],
                    y=vals[1]
                )
            else:
                raise Exception(f"Point input values should be 2 lists, got {len(vals)} lists")

            pts.append(pt)

        return cls(pts)

    def pt(self, pt_ndx: numbers.Integral) -> Point2D:
        """
        Extract the point at index pt_ndx.

        :param pt_ndx: point index.
        :type pt_ndx: numbers.Integral.
        :return: the extracted Point instance.
        :rtype: Point.

        Examples:
        """

        return Point2D(
            x=self._x[pt_ndx],
            y=self._y[pt_ndx]
        )

    def values_at(self,
        ndx: numbers.Integral
    ) -> Tuple[float, float]:
        """
        Return the values at given index.

        :param ndx: the index of the point values to extract
        :type ndx: numbers.Integral
        :return: the x and y values
        """

        return (
            self._x[ndx],
            self._y[ndx]
        )

    def pts(self):

        return [Point2D(*self.values_at(ndx)) for ndx in range(self.num_pts())]

    def segment(self,
        ndx: numbers.Integral
    ) -> Optional[Segment2D]:
        """
        Returns the optional segment at index ndx.

        :param ndx: the segment index.
        :type ndx: numbers.Integral
        :return: the optional segment
        :rtype: Optional[Segment]
        """

        start_pt = self.pt(ndx)
        end_pt = self.pt(ndx + 1)

        if start_pt.is_coincident(end_pt):
            return None
        else:
            return Segment2D(
                start_pt=self.pt(ndx),
                end_pt=self.pt(ndx + 1)
            )

    def num_pts(self):

        return len(self._x)

    def start_pt(self) -> Optional[Point2D]:
        """
        Return the first point of a Line or None when no points.

        :return: the first point or None.
        :rtype: optional Point instance.
        """

        return self.pt(0) if self.num_pts() > 0 else None

    def end_pt(self) -> Optional[Point2D]:
        """
        Return the last point of a Line or None when no points.

        :return: the last point or None.
        :rtype: optional Point instance.
        """

        return self.pt(-1) if self.num_pts() > 0 else None

    def __iter__(self):
        """
        Return each element of a Line, i.e., its segments.
        """

        return (self.segment(i) for i in range(self.num_pts()-1))

    def __repr__(self) -> str:
        """
        Represents a Line instance as a shortened text.

        :return: a textual shortened representation of a Line instance.
        :rtype: str.
        """

        num_points = self.num_pts()

        if num_points == 0:
            txt = "Empty Line2D"
        else:
            x1, y1 = self.start_pt()
            if num_points == 1:
                txt = "Line with unique point: {.4f}, {.4f}".format(x1, y1)
            else:
                x2, y2 = self.end_pt()
                txt = "Line with {} points: ({:.4f}, {:.4f}) ... ({:.4f}, {:.4f})".format(num_points, x1, y1, x2, y2)

        return txt

    def add_pt(self, pt) -> bool:
        """
        In-place transformation of the original Line2D instance
        by adding a new point at the end.

        :param pt: the point to add
        :type pt: Point.
        :return: status of addition. True when added, False otherwise.
        :rtype: bool.
        """

        self._x.append(pt.x)
        self._y.append(pt.y)
        return True

    def add_pts(self, pt_list) -> numbers.Integral:
        """
        In-place transformation of the original Line instance
        by adding a new set of points at the end.

        :param pt_list: list of Points.
        :type pt_list: List of Point instances.
        :return: number of added points
        :rtype: numbers.Integral.
        """

        num_added = 0
        for pt in pt_list:
            success = self.add_pt(pt)
            if success:
                num_added += 1

        return num_added

    def x_list(self) -> List[numbers.Real]:

        return list(self._x)

    def y_list(self) -> List[numbers.Real]:

        return list(self._y)

    def xy_lists(self) -> Tuple[List[numbers.Real], List[numbers.Real]]:

        return self.x_list(), self.y_list()

    def xy_zipped(self) -> List[Tuple[numbers.Real, numbers.Real]]:

        return [(x, y) for x, y in zip(self.x_list(), self.y_list())]

    def x_min(self) -> Optional[numbers.Real]:
        """
        Optional minimum of x values.

        :return: the optional minimum of x values.
        :rtype: Optional[numbers.Real]

        Examples:
          >>> l = Line2D.fromPointList([[0, 0], [1, 0], [0, 1]])
          >>> l.x_min()
          0.0
        """

        return min(self._x) if self.num_pts() > 0 else None

    def x_max(self) -> Optional[numbers.Real]:
        """
        Optional maximum of x values.

        :return: the optional maximum of x values.
        :rtype: Optional[numbers.Real]

        Examples:
          >>> l = Line2D.fromPointList([[0, 0], [1, 0], [0, 1]])
          >>> l.x_max()
          1.0
        """

        return max(self._x) if self.num_pts() > 0 else None

    def y_min(self) -> Optional[numbers.Real]:
        """
        Optional minimum of y values.

        :return: the optional minimum of y values.
        :rtype: Optional[numbers.Real]

        Examples:
          >>> l = Line2D.fromPointList([[0, 0], [1, 0], [0, 1]])
          >>> l.y_min()
          0.0
        """

        return min(self._y) if self.num_pts() > 0 else None

    def y_max(self) -> Optional[numbers.Real]:
        """
        Optional maximum of y values.

        :return: the optional maximum of y values.
        :rtype: Optional[numbers.Real]

        Examples:
          >>> l = Line2D.fromPointList([[0, 0], [1, 0], [0, 1]])
          >>> l.y_max()
          1.0
        """

        return max(self._y) if self.num_pts() > 0 else None

    def remove_coincident_points(self) -> Optional['Line2D']:
        """
        Remove coincident successive points

        :return: Line2D instance
        """

        if self.num_pts() == 0:
            return

        new_line = Line2D(
            pts=[self.pt(0)]
        )

        for ndx in range(1, self.num_pts()):
            if not self.pt(ndx).is_coincident(new_line.pt(-1)):
                new_line.add_pt(self.pt(ndx))

        return new_line

    def as_segments(self):
        """
        Convert to a list of segments 2d.

        :return: list of Segment2D objects
        """

        pts_pairs = zip(self.pts()[:-1], self.pts()[1:])

        segments = [Segment2D(pt_a, pt_b) for (pt_a, pt_b) in pts_pairs]

        return segments

    def densify_2d_line(self, sample_distance) -> 'Line2D':
        """
        Densify a line into a new line instance,
        using the provided sample distance.
        Returned Line instance has coincident successive points removed.

        :param sample_distance: numbers.Real
        :return: Line2D instance
        """

        if sample_distance <= 0.0:
            raise Exception(f"Sample distance must be positive. {sample_distance} received")

        segments = self.as_segments()

        densified_line_list = [segment.densify2d_asLine(sample_distance) for segment in segments]

        densifyied_multiline = MultiLine(densified_line_list)

        densifyied_line = densifyied_multiline.to_line()

        densifyied_line_wo_coinc_pts = densifyied_line.remove_coincident_points()

        return densifyied_line_wo_coinc_pts

    def join(self, another) -> 'Line2D':
        """
        Joins together two lines and returns the join as a new line without point changes,
        with possible overlapping points
        and orientation mismatches between the two original lines
        """

        return Line2D(self.pts() + another.pts())

    def length_2d(self) -> numbers.Real:

        length = 0.0
        for ndx in range(self.num_pts() - 1):
            length += self.pt(ndx).dist2DWith(self.pt(ndx + 1))
        return length

    def step_lengths_2d(self) -> List[numbers.Real]:
        """
        Returns the point-to-point 2D distances.
        It is the distance between a point and its previous one.
        The list has the same length as the source point list.

        :return: the individual 2D segment lengths.
        :rtype: list of floats.

        Examples:
        """

        step_length_list = [0.0]
        for ndx in range(1, self.num_pts()):
            length = self.pt(ndx).dist2DWith(self.pt(ndx - 1))
            step_length_list.append(length)

        return step_length_list

    def incremental_length_2d(self) -> List[numbers.Real]:
        """
        Returns the accumulated 2D segment lengths.

        :return: accumulated 2D segment lenghts
        :rtype: list of floats.
        """

        return list(itertools.accumulate(self.step_lengths_2d()))

    def reversed(self) -> 'Line2D':
        """
        Return a Line instance with reversed point list.

        :return: a new Line instance.
        :rtype: Line.
        """

        pts = [pt.clone() for pt in self.pts()]
        pts.reverse()

        return Line2D(
            pts=pts
        )

    def extremes_distance_2d(self) -> numbers.Real:
        """
        Calculate the 2D distance between start and end points.

        :return: the 2D distance between start and end points
        """

        return self.end_pt().dist2DWith(self.start_pt())

    def isClosed_2d(self,
        tolerance: numbers.Real = MIN_SEPARATION_THRESHOLD
    ) -> bool:
        """
        Determine if the line is 2D-closed.

        :param tolerance: the tolerance for considering the line closed
        :return: whether the line is to be considered 2D-closed
        """

        return self.end_pt().is_coincident(self.start_pt(), tolerance=tolerance)

    def walk_backward(self) -> 'Line2D':
        """
        Create a new line by walking the line backward from the last point up to the first and thus closing it.

        :return: a closed line with zero area
        :rtype: 'Line'
        """

        return Line2D(self.pts() + self.reversed().pts()[1:])

    def clone(self) -> 'Line2D':
        """
        Clone a line.

        :return: the cloned line
        :rtype: Line2D
        """

        return Line2D(self.pts())

    def close_2d(self) -> 'Line2D':
        """
        Return a line that is 2D-closed.

        :return: a 2D-closed line
        :rtype: Line2D
        """

        line = self.clone()

        if not line.isClosed_2d():

            line.add_pt(line.start_pt())

        return line


class Ellipse2D(Shape2D):

    def area(self):
        pass

    def length(self):
        pass


class Circle2D(Ellipse2D):

    def __init__(self,
                 x: numbers.Real,
                 y: numbers.Real,
                 r: numbers.Real
                 ):

        self._x = float(x)
        self._y = float(y)
        self._r = float(r)

    @property
    def x(self) -> numbers.Real:
        """
        Return the x coordinate of the current circle.

        :return: x coordinate.

        Examples:
          >>> Circle2D(4, 3, 2).x
          4.0
          >>> Circle2D(-0.39, 3, 7).x
          -0.39
        """

        return self._x

    @property
    def y(self) -> numbers.Real:
        """
        Return the y coordinate of the current circle.

        :return: y coordinate.

        Examples:
          >>> Point2D(4, 3).y
          3.0
          >>> Point2D(-0.39, 17.42).y
          17.42
        """

        return self._y

    @property
    def radius(self):
        return self._r

    def area(self):
        return math.pi * self._r * self._r

    def length(self):
        return 2.0 * math.pi * self._r

    def clone(self):
        return Circle2D(self._x, self._y, self._r)


class Polygon2D(Shape2D, metaclass=abc.ABCMeta):

    @abc.abstractmethod
    def num_side(self):
        """Return numer of sides"""

    def area(self):
        pass

    def length(self):
        pass


class Triangle2D(Polygon2D):

    def area(self):
        pass

    def length(self):
        pass

    @property
    def num_side(self):
        """Return numer of sides"""

        return 3


class Quadrilateral2D(Polygon2D, metaclass=abc.ABCMeta):

    def area(self):
        pass

    def length(self):
        pass

    @property
    def num_side(self):
        """Return numer of sides"""

        return 4


class Square2D(Quadrilateral2D):

    def __init__(self,
                 x: numbers.Real,
                 y: numbers.Real,
                 side: numbers.Real,
                 rotation: numbers.Real
                 ):
        """

        :param x: x coordinate of square center
        :param y: y coordinate of square center
        :param side: square side
        :param rotation: square rotation, counter-clockwise, decimal degrees
        """

        self._x = x
        self._y = y
        self._side = side
        self._cc_rotat = rotation % 360.0

    def area(self):
        return self._side * self._side

    def length(self):
        return 4.0 * self._side


@functools.singledispatch
def center(
        shape: Shape2D
) -> Point2D:
    """
    The 2D shape center as a point (2D)
    """

    raise NotImplementedError("Center method is not implemented for {type(shape)}")


@center.register(Point2D)
def _(
        shape: Point2D
) -> Point2D:

    return Point2D(
        x=shape.x,
        y=shape.y
    )


@center.register(Segment2D)
def _(
        shape: Segment2D
) -> Point2D:
    """
    Segment mean point.

    Examples:
      >>> center(Segment2D(Point2D(3, 2), Point2D(1, 1)))
      Point(2.0, 1.5)
    """

    x0, y0 = shape.start_pt.a()
    x1, y1 = shape.end_pt.a()

    return Point2D(
        x=(x0 + x1)/2.0,
        y=(y0 + y1)/2.0
    )


@center.register(Circle2D)
def _(
        shape: Circle2D
) -> Point2D:

    return Point2D(
        x=shape.x,
        y=shape.y
    )


if __name__ == "__main__":

    import doctest
    doctest.testmod()
