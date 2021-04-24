
import abc

from typing import Optional, Callable

from math import fabs

import itertools

from copy import copy

import numbers
from array import array

from .space2d import Point2D, Segment2D
from ...orientations.orientations import *
from ...mathematics.statistics import *
from ...mathematics.quaternions import *
from ...utils.types import check_type



class Shape3D(object, metaclass=abc.ABCMeta):

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


class Point3D:
    """
    Cartesian point.
    Dimensions: 3D
    """

    def __init__(
        self,
        x: numbers.Real,
        y: numbers.Real,
        z: numbers.Real = 0.0
    ):
        """
        Construct a Point instance.

        :param x: point x coordinate.
        :type x: numbers.Real.
        :param y: point y coordinate.
        :type y: numbers.Real.
        :param z: point z coordinate.
        :type z: numbers.Real.
        """

        vals = [x, y]
        if any(map(lambda val: not isinstance(val, numbers.Real), vals)):
            raise Exception("X and y input values must be integer or float type")
        if not all(map(math.isfinite, vals)):
            raise Exception("X and y input values must be finite (#03)")

        self._x = float(x)
        self._y = float(y)
        self._z = float(z)

    @classmethod
    def fromVect(cls,
                 vect: Vect3D) -> 'Point3D':
        """

        :param vect:
        :return:
        """

        return cls(
            x=vect.x,
            y=vect.y,
            z=vect.z
        )

    @property
    def x(self) -> numbers.Real:
        """
        Return the x coordinate of the current point.

        :return: x coordinate.
        :rtype: numbers.Real

        Examples:
          >>> Point3D(4, 3, 7).x
          4.0
          >>> Point3D(-0.39, 3, 7).x
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
          >>> Point3D(4, 3, 7).y
          3.0
          >>> Point3D(-0.39, 17.42, 7).y
          17.42
        """

        return self._y

    @property
    def z(self) -> numbers.Real:
        """
        Return the z coordinate of the current point.

        :return: z coordinate.
        :rtype: numbers.Real

        Examples:
          >>> Point3D(4, 3, 7).z
          7.0
          >>> Point3D(-0.39, 17.42, 8.9).z
          8.9
        """

        return self._z

    def __iter__(self):
        """
        Return the elements of a Point.

        :return:

        Examples;
          >>> x, y, z = Point3D(1,1)
          >>> x == 1
          True
          >>> y == 1
          True

        """

        return (i for i in self.a())

    def __repr__(self) -> str:

        return "Point3D({:.4f}, {:.4f}, {:.4f})".format(self.x, self.y, self.z)

    def __eq__(self,
        another: 'Point3D'
    ) -> bool:
        """
        Return True if objects are equal.

        :param another: another point.
        :type another: Point.
        :raise: Exception.

        Example:
          >>> Point3D(1., 1., 1.) == Point3D(1, 1, 1)
          True
          >>> Point3D(1., 1., 1.) == Point3D(1, 1, 1)
          True
          >>> Point3D(1., 1., 1.) == Point3D(1, 1, -1)
          False
        """

        if not isinstance(another, Point3D):
            raise Exception("Another instance must be a Point")

        return all([
            self.x == another.x,
            self.y == another.y,
            self.z == another.z
            ]
        )

    def __ne__(self,
        another: 'Point3D'
    ) -> bool:
        """
        Return False if objects are equal.

        Example:
          >>> Point3D(1., 1., 1.) != Point3D(0., 0., 0.)
          True
          >>> Point3D(1., 1., 1.) != Point3D(1, 1, 1)
          False
        """

        return not (self == another)

    def a(self) -> Tuple[numbers.Real, numbers.Real, numbers.Real]:
        """
        Return the individual values of the point.

        :return: double array of x, y, z values

        Examples:
          >>> Point3D(4, 3, 7).a()
          (4.0, 3.0, 7.0)
        """

        return self.x, self.y, self.z

    def __add__(self, another: 'Point3D') -> 'Point3D':
        """
        Sum of two points.

        :param another: the point to add
        :type another: Point3D
        :return: the sum of the two points
        :rtype: Point3D
        :raise: Exception

        Example:
          >>> Point3D(1, 0, 0) + Point3D(0, 1, 1)
          Point3D(1.0000, 1.0000, 1.0000)
          >>> Point3D(1, 1, 1) + Point3D(-1, -1, -1)
          Point3D(0.0000, 0.0000, 0.0000)
        """

        check_type(another, "Second point", Point3D)

        x0, y0, z0 = self
        x1, y1, z1 = another

        return Point3D(
            x=x0+x1,
            y=y0+y1,
            z=z0+z1
        )

    def __sub__(self,
        another: 'Point3D'
    ) -> 'Point3D':
        """Subtract two points.

        :param another: the point to subtract
        :type another: Point3D
        :return: the difference between the two points
        :rtype: Point3D
        :raise: Exception

        Example:
          >>> Point3D(1., 1., 1.) - Point3D(1., 1., 1.)
          Point3D(0.0000, 0.0000, 0.0000)
          >>> Point3D(1., 1., 3.) - Point3D(1., 1., 2.2)
          Point3D(0.0000, 0.0000, 0.8000)
        """

        check_type(another, "Second point", Point3D)

        x0, y0, z0 = self
        x1, y1, z1 = another

        return Point3D(
            x=x0 - x1,
            y=y0 - y1,
            z=z0 - z1
        )

    def clone(self) -> 'Point3D':
        """
        Clone a point.

        :return: a new point.
        :rtype: Point.
        """

        return Point3D(*self.a())

    def toXYZ(self) -> Tuple[numbers.Real, numbers.Real, numbers.Real]:
        """
        Returns the spatial components as a tuple of three values.

        :return: the spatial components (x, y, z).
        :rtype: a tuple of three floats.

        Examples:
          >>> Point3D(1, 0, 3).toXYZ()
          (1.0, 0.0, 3.0)
        """

        return self.x, self.y, self.z

    def toArray(self) -> np.ndarray:
        """
        Return a Numpy array representing the point values.

        :return: Numpy array

        Examples:
          >>> np.allclose(Point3D(1, 2, 3).toArray(), np.array([ 1., 2., 3.]))
          True
        """

        return np.asarray(self.toXYZ())

    def to2d(self) -> Point2D:
        """
        Projection on the x-y plane as a 2D point.

        Examples:
          >>> Point3D(2, 3, 4).to2d()
          Point2D(2.0000, 3.0000)
        """

        return Point2D(
            x=self.x,
            y=self.y
        )

    def pXY(self) -> 'Point3D':
        """
        Projection on the x-y plane

        :return: projected object instance

        Examples:
          >>> Point3D(2, 3, 4).pXY()
          Point3D(2.0000, 3.0000, 0.0000)
        """

        return Point3D(self.x, self.y, 0.0)

    def pXZ(self) -> 'Point3D':
        """
        Projection on the x-z plane

        :return: projected object instance

        Examples:
          >>> Point3D(2, 3, 4).pXZ()
          Point3D(2.0000, 0.0000, 4.0000)
        """

        return Point3D(self.x, 0.0, self.z)

    def pYZ(self) -> 'Point3D':
        """
        Projection on the y-z plane

        :return: projected object instance

        Examples:
          >>> Point3D(2, 3, 4).pYZ()
          Point3D(0.0000, 3.0000, 4.0000)
        """

        return Point3D(0.0, self.y, self.z)

    def deltaX(self,
        another: 'Point3D'
    ) -> Optional[numbers.Real]:
        """
        Delta between x components of two Point Instances.

        :return: x coordinates difference value.
        :rtype: optional numbers.Real.
        :raise: Exception

        Examples:
          >>> Point3D(1, 2, 3).deltaX(Point3D(4, 7, 1))
          3.0
        """

        return another.x - self.x

    def deltaY(self,
        another: 'Point3D'
    ) -> Optional[numbers.Real]:
        """
        Delta between y components of two Point Instances.

        :return: y coordinates difference value.
        :rtype: optional numbers.Real.

        Examples:
          >>> Point3D(1, 2, 3).deltaY(Point3D(4, 7, 1))
          5.0
        """

        return another.y - self.y

    def deltaZ(self,
        another: 'Point3D'
    ) -> Optional[numbers.Real]:
        """
        Delta between z components of two Point Instances.

        :return: z coordinates difference value.
        :rtype: optional numbers.Real.

        Examples:
          >>> Point3D(1, 2, 3).deltaZ(Point3D(4, 7, 1))
          -2.0
        """

        return another.z - self.z

    def distance(self,
                 another: 'Point3D'
                 ) -> numbers.Real:
        """
        Calculate Euclidean spatial distance between two points.
        TODO: consider case of polar CRS

        :param another: another Point instance.
        :type another: Point.
        :return: the distance (when the two points have the same CRS).
        :rtype: numbers.Real.
        :raise: Exception.

        Examples:
          >>> Point3D(1., 1., 1.).distance(Point3D(4., 5., 1))
          5.0
          >>> Point3D(1, 1, 1).distance(Point3D(4, 5, 1))
          5.0
          >>> Point3D(1, 1, 1).distance(Point3D(4, 5, 1))
          5.0
        """

        check_type(another, "Point", Point3D)

        return sqrt((self.x - another.x) ** 2 + (self.y - another.y) ** 2 + (self.z - another.z) ** 2)

    def dist_2d(self,
                another: 'Point3D'
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
          >>> Point3D(1., 1., 1.).dist_2d(Point3D(4., 5., 7.))
          5.0
        """

        check_type(another, "Second point", Point3D)

        return sqrt((self.x - another.x) ** 2 + (self.y - another.y) ** 2)

    def scale(self,
        scale_factor: numbers.Real
    ) -> 'Point3D':
        """
        Create a scaled object.
        Note: it does not make sense for polar coordinates.
        TODO: manage polar coordinates cases OR deprecate and remove - after dependency check.

        Example;
          >>> Point3D(1, 0, 1).scale(2.5)
          Point3D(2.5000, 0.0000, 2.5000)
          >>> Point3D(1, 0, 1).scale(2.5)
          Point3D(2.5000, 0.0000, 2.5000)
        """

        x, y, z = self.x * scale_factor, self.y * scale_factor, self.z * scale_factor
        return Point3D(x, y, z)

    def invert(self) -> 'Point3D':
        """
        Create a new object with inverted direction.
        Note: it depends on scale method, that could be deprecated/removed.

        Examples:
          >>> Point3D(1, 1, 1).invert()
          Point3D(-1.0000, -1.0000, -1.0000)
          >>> Point3D(2, -1, 4).invert()
          Point3D(-2.0000, 1.0000, -4.0000)
        """

        return self.scale(-1)

    def reflect_vertical(self) -> 'Point3D':
        """
        Reflect a point along a vertical axis.

        :return: reflected point.
        :rtype: Point3D

        Examples:
          >>> Point3D(1,1,1).reflect_vertical()
          Point3D(-1.0000, -1.0000, 1.0000)
        """

        x, y, z = self

        return Point3D(
            x=-x,
            y=-y,
            z=z
        )

    def is_coincident(self,
                      another: 'Point3D',
                      tolerance: numbers.Real = MIN_SEPARATION_THRESHOLD
                      ) -> bool:
        """
        Check spatial coincidence of two points

        :param another: the point to compare.
        :type another: Point.
        :param tolerance: the maximum allowed distance between the two points.
        :type tolerance: numbers.Real.
        :return: whether the two points are coincident.
        :rtype: bool.
        :raise: Exception.

        Example:
          >>> Point3D(1., 0., -1.).is_coincident(Point3D(1., 1.5, -1.))
          False
          >>> Point3D(1., 0., 0.).is_coincident(Point3D(1., 0., 0.))
          True
        """

        check_type(another, "Second point", Point3D)

        return self.distance(another) <= tolerance

    def already_present(self,
                        pt_list: List['Point3D'],
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
        sy: numbers.Real,
        sz: numbers.Real
    ) -> Optional['Point3D']:
        """
        Create a new object shifted by given amount from the self instance.

        Example:
          >>> Point3D(1, 1, 1).shift(0.5, 1., 1.5)
          Point3D(1.5000, 2.0000, 2.5000)
          >>> Point3D(1, 2, -1).shift(0.5, 1., 1.5)
          Point3D(1.5000, 3.0000, 0.5000)
       """

        return Point3D(self.x + sx, self.y + sy, self.z + sz)

    def shiftByVect(self,
        v: Vect3D
    ) -> 'Point3D':
        """
        Create a new point shifted from the self instance by given vector.

        :param v: the shift vector.
        :type v: Vect.
        :return: the shifted point.
        :rtype: Point.
        :raise: Exception

        Example:
          >>> Point3D(1, 1, 1).shiftByVect(Vect3D(0.5, 1., 1.5))
          Point3D(1.5000, 2.0000, 2.5000)
          >>> Point3D(1, 2, -1).shiftByVect(Vect3D(0.5, 1., 1.5))
          Point3D(1.5000, 3.0000, 0.5000)
       """

        x, y, z = self

        sx, sy, sz = v.toXYZ()

        return Point3D(x + sx, y + sy, z + sz)

    def asVect(self) -> 'Vect3D':
        """
        Create a vector based on the point coordinates

        Example:
          >>> Point3D(1, 1, 0).asVect()
          Vect3D(1.0000, 1.0000, 0.0000)
          >>> Point3D(0.2, 1, 6).asVect()
          Vect3D(0.2000, 1.0000, 6.0000)
        """

        return Vect3D(self.x, self.y, self.z)

    def rotate(self,
        rotation_axis: RotationAxis,
        center_point: 'Point3D' = None
        ) -> 'Point3D':
        """
        Rotates a point.
        :param rotation_axis:
        :param center_point:
        :return: the rotated point
        :rtype: Point3D

        Examples:
          >>> pt = Point3D(0,0,1)
          >>> rot_axis = RotationAxis(0,0,90)
          >>> center_pt = Point3D(0,0,0.5)
          >>> pt.rotate(rotation_axis=rot_axis, center_point=center_pt)
          Point3D(0.5000, 0.0000, 0.5000)
          >>> center_pt = Point3D(0,0,1)
          >>> pt.rotate(rotation_axis=rot_axis, center_point=center_pt)
          Point3D(0.0000, 0.0000, 1.0000)
          >>> center_pt = Point3D(0, 0, 2)
          >>> pt.rotate(rotation_axis=rot_axis, center_point=center_pt)
          Point3D(-1.0000, 0.0000, 2.0000)
          >>> rot_axis = RotationAxis(0,0,180)
          >>> pt.rotate(rotation_axis=rot_axis, center_point=center_pt)
          Point3D(-0.0000, 0.0000, 3.0000)
          >>> pt.rotate(rotation_axis=rot_axis)
          Point3D(0.0000, 0.0000, -1.0000)
          >>> pt = Point3D(1,1,1)
          >>> rot_axis = RotationAxis(0,90,90)
          >>> pt.rotate(rotation_axis=rot_axis)
          Point3D(1.0000, -1.0000, 1.0000)
          >>> rot_axis = RotationAxis(0,90,180)
          >>> pt.rotate(rotation_axis=rot_axis)
          Point3D(-1.0000, -1.0000, 1.0000)
          >>> center_pt = Point3D(1,1,1)
          >>> pt.rotate(rotation_axis=rot_axis, center_point=center_pt)
          Point3D(1.0000, 1.0000, 1.0000)
          >>> center_pt = Point3D(2,2,10)
          >>> pt.rotate(rotation_axis=rot_axis, center_point=center_pt)
          Point3D(3.0000, 3.0000, 1.0000)
          >>> pt = Point3D(1, 1, 2)
          >>> rot_axis = RotationAxis(135, 0, 180)
          >>> center_pt = Point3D(0,0,1)
          >>> pt.rotate(rotation_axis=rot_axis, center_point=center_pt)
          Point3D(-1.0000, -1.0000, 0.0000)
        """

        if not center_point:

            center_point = Point3D(
                x=0.0,
                y=0.0,
                z=0.0
            )

        check_type(center_point, "Center point", Point3D)

        p_diff = self - center_point

        p_vect = p_diff.asVect()

        rot_vect = rotVectByAxis(
            v=p_vect,
            rot_axis=rotation_axis
        )

        x, y, z = rot_vect

        rot_pt = Point3D(
            x=x,
            y=y,
            z=z
        )

        transl_pt = center_point + rot_pt

        return transl_pt

    @classmethod
    def random(cls,
        lower_boundary: float = -MAX_SCALAR_VALUE,
        upper_boundary: float = MAX_SCALAR_VALUE
    ):
        """
        Creates a random point.

        :return: random point
        :rtype: Point3D
        """

        vals = [random.uniform(lower_boundary, upper_boundary) for _ in range(3)]
        return cls(*vals)


def pack_to_points(
    xs: array,
    ys: array,
    zs: Optional[array] = None,
) -> List[Point3D]:
    # Side effects: None
    """
    Create a list of points given a set
    of input arrays.

    :param xs: array of x values
    :param ys: array of y values
    :param zs: optional array of z values
    :return: a list of Point3D instances
    """

    if zs is None:
        zs = [0.0] * len(xs)
    pts = []
    for x, y, z, t in zip(xs, ys, zs):
        pts.append(
            Point3D(
                x,
                y,
                z
            )
        )

    return pts


class Segment3D:
    """
    Segment is a geometric object defined by the straight line between
    two vertices.
    """

    def __init__(self,
                 start_pt: Point3D,
                 end_pt: Point3D):
        """
        Creates a segment instance provided the two points have the same CRS code.

        :param start_pt: the start point.
        :type: Point.
        :param end_pt: the end point.
        :type end_pt: Point.
        :return: the new segment instance if both points have the same georeferenced.
        :raises: CRSCodeException.
        """

        check_type(start_pt, "Start point", Point3D)

        check_type(end_pt, "End point", Point3D)

        if start_pt.distance(end_pt) == 0.0:
            raise Exception("Source points cannot be coincident")

        self._start_pt = start_pt.clone()
        self._end_pt = end_pt.clone()

    @classmethod
    def fromVector(cls,
                   point: Point3D,
                   dir_vector: Vect3D):

        check_type(point, "Input point", Point3D)
        check_type(dir_vector, "Directional vector", Vect3D)

        start_pt = point
        end_pt = start_pt.shiftByVect(dir_vector)

        return cls(
            start_pt=start_pt,
            end_pt=end_pt
        )

    @classmethod
    def from2D(cls,
        segment: Segment2D):

        check_type(segment, "Input segment", Segment2D)

        start_pt = Point3D(
            x=segment.start_pt.x,
            y=segment.start_pt.y,
            z=0.0
        )

        end_pt = Point3D(
            x=segment.end_pt.x,
            y=segment.end_pt.y,
            z=0.0
        )

        return cls(
            start_pt=start_pt,
            end_pt=end_pt
        )

    def __repr__(self) -> str:
        """
        Represents a Segment instance.

        :return: the Segment representation.
        :rtype: str.
        """

        return "Segment3D(start_pt={}, end_pt={})".format(
            self.start_pt,
            self.end_pt
        )

    @property
    def start_pt(self) -> Point3D:

        return self._start_pt

    @property
    def end_pt(self) -> Point3D:

        return self._end_pt

    def __iter__(self):
        """
        Return the elements of a Segment, i.e., start and end point.
        """

        return (i for i in [self.start_pt, self.end_pt])

    def clone(self) -> 'Segment3D':

        return Segment3D(self._start_pt, self._end_pt)

    def increasing_x(self) -> 'Segment3D':

        if self.end_pt.x < self.start_pt.x:
            return Segment3D(self.end_pt, self.start_pt)
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

    def z_range(self) -> Tuple[numbers.Real, numbers.Real]:

        if self.start_pt.z < self.end_pt.z:
            return self.start_pt.z, self.end_pt.z
        else:
            return self.end_pt.z, self.start_pt.z

    def delta_x(self) -> numbers.Real:
        """
        X delta between segment end point and start point.

        :return: the horizontal, x-parallel distance between segment end point and start point.
        """

        return self.end_pt.x - self.start_pt.x

    def delta_y(self) -> numbers.Real:
        """
        Y delta between segment end point and start point.

        :return: the horizontal, y-parallel distance between segment end point and start point.
        """

        return self.end_pt.y - self.start_pt.y

    def delta_z(self) -> numbers.Real:
        """
        Z delta between segment end point and start point.

        :return: the vertical distance between segment end point and start point.
        """

        return self.end_pt.z - self.start_pt.z

    def as_vector(self) -> Vect3D:
        """
        Convert a segment to a vector.
        """

        return Vect3D(
            x=self.delta_x(),
            y=self.delta_y(),
            z=self.delta_z()
        )

    def length_horizontal(self) -> numbers.Real:

        return self.start_pt.dist_2d(self.end_pt)

    def length(self) -> numbers.Real:

        return self.start_pt.distance(self.end_pt)

    def ratio_delta_zs(self) -> Optional[numbers.Real]:
        """
        Calculates the delta z - delta s ratio of a segment.

        :return: optional numbers.Real.
        """

        len2d = self.length_horizontal()

        if len2d == 0.0:
            return None

        return self.delta_z() / len2d

    def slope_rad(self) -> Optional[numbers.Real]:
        """
        Calculates the slope in radians of the segment.
        Positive is downward point, negative upward pointing.

        :return: optional numbers.Real.
        """

        delta_zs = self.ratio_delta_zs()

        if delta_zs is None:
            return None
        else:
            return - math.atan(delta_zs)

    def vector(self) -> Vect3D:

        return Vect3D(self.delta_x(),
                      self.delta_y(),
                      self.delta_z()
                      )

    def antivector(self) -> Vect3D:
        """
        Returns the vector pointing from the segment end to the segment start.

        :return: the vector pointing from the segment end to the segment start.
        :rtype: Vect.
        """

        return self.vector().invert()

    def contains_pt(self,
        pt: Point3D
    ) -> bool:
        """
        Checks whether a point is contained in a segment.

        :param pt: the point for which to check containement.
        :return: bool.
        :raise: Exception.

        Examples:
          >>> segment = Segment3D(Point3D(0, 0, 0), Point3D(1, 0, 0))
          >>> segment.contains_pt(Point3D(0, 0, 0))
          True
          >>> segment.contains_pt(Point3D(1, 0, 0))
          True
          >>> segment.contains_pt(Point3D(0.5, 0, 0))
          True
          >>> segment.contains_pt(Point3D(0.5, 0.00001, 0))
          False
          >>> segment.contains_pt(Point3D(0.5, 0, 0.00001))
          False
          >>> segment.contains_pt(Point3D(1.00001, 0, 0))
          False
          >>> segment.contains_pt(Point3D(0.000001, 0, 0))
          True
          >>> segment.contains_pt(Point3D(-0.000001, 0, 0))
          False
          >>> segment.contains_pt(Point3D(0.5, 1000, 1000))
          False
          >>> segment = Segment3D(Point3D(0, 0, 0), Point3D(0, 1, 0))
          >>> segment.contains_pt(Point3D(0, 0, 0))
          True
          >>> segment.contains_pt(Point3D(0, 0.5, 0))
          True
          >>> segment.contains_pt(Point3D(0, 1, 0))
          True
          >>> segment.contains_pt(Point3D(0, 1.5, 0))
          False
          >>> segment = Segment3D(Point3D(0, 0, 0), Point3D(1, 1, 1))
          >>> segment.contains_pt(Point3D(0.5, 0.5, 0.5))
          True
          >>> segment.contains_pt(Point3D(1, 1, 1))
          True
          >>> segment = Segment3D(Point3D(1,2,3), Point3D(9,8,2))
          >>> segment.contains_pt(segment.pointAt(0.745))
          True
          >>> segment.contains_pt(segment.pointAt(1.745))
          False
          >>> segment.contains_pt(segment.pointAt(-0.745))
          False
          >>> segment.contains_pt(segment.pointAt(0))
          True
        """

        check_type(pt, "Point", Point3D)

        segment_length = self.length()
        length_startpt_pt = self.start_pt.distance(pt)
        length_endpt_pt = self.end_pt.distance(pt)

        return areClose(
            a=segment_length,
            b=length_startpt_pt + length_endpt_pt
        )

    def pointAt(self,
        scale_factor: numbers.Real
    ) -> Point3D:
        """
        Returns a point aligned with the segment
        and lying at given scale factor, where 1 is segment length
        ans 0 is segment start.

        :param scale_factor: the scale factor, where 1 is the segment length.
        :type scale_factor: numbers.Real
        :return: Point at scale factor
        :rtype: Point3D

        Examples:
          >>> s = Segment3D(Point3D(0,0,0), Point3D(1,0,0))
          >>> s.pointAt(0)
          Point3D(0.0000, 0.0000, 0.0000)
          >>> s.pointAt(0.5)
          Point3D(0.5000, 0.0000, 0.0000)
          >>> s.pointAt(1)
          Point3D(1.0000, 0.0000, 0.0000)
          >>> s.pointAt(-1)
          Point3D(-1.0000, 0.0000, 0.0000)
          >>> s.pointAt(-2)
          Point3D(-2.0000, 0.0000, 0.0000)
          >>> s.pointAt(2)
          Point3D(2.0000, 0.0000, 0.0000)
          >>> s = Segment3D(Point3D(0,0,0), Point3D(0,0,1))
          >>> s.pointAt(0)
          Point3D(0.0000, 0.0000, 0.0000)
          >>> s.pointAt(0.5)
          Point3D(0.0000, 0.0000, 0.5000)
          >>> s.pointAt(1)
          Point3D(0.0000, 0.0000, 1.0000)
          >>> s.pointAt(-1)
          Point3D(0.0000, 0.0000, -1.0000)
          >>> s.pointAt(-2)
          Point3D(0.0000, 0.0000, -2.0000)
          >>> s.pointAt(2)
          Point3D(0.0000, 0.0000, 2.0000)
          >>> s = Segment3D(Point3D(0,0,0), Point3D(1,1,1))
          >>> s.pointAt(0.5)
          Point3D(0.5000, 0.5000, 0.5000)
          >>> s = Segment3D(Point3D(0,0,0), Point3D(4,0,0))
          >>> s.pointAt(7.5)
          Point3D(30.0000, 0.0000, 0.0000)
        """

        dx = self.delta_x() * scale_factor
        dy = self.delta_y() * scale_factor
        dz = self.delta_z() * scale_factor

        return Point3D(
            x=self.start_pt.x + dx,
            y=self.start_pt.y + dy,
            z=self.start_pt.z + dz
        )

    def pointProjection(self,
        point: Point3D
    ) -> Point3D:
        """
        Return the point projection on the segment.

        Examples:
          >>> s = Segment3D(start_pt=Point3D(0,0,0), end_pt=Point3D(1,0,0))
          >>> p = Point3D(0.5, 1, 4)
          >>> s.pointProjection(p)
          Point3D(0.5000, 0.0000, 0.0000)
          >>> s = Segment3D(start_pt=Point3D(0,0,0), end_pt=Point3D(4,0,0))
          >>> p = Point3D(7.5, 19.2, -14.72)
          >>> s.pointProjection(p)
          Point3D(7.5000, 0.0000, 0.0000)
        """

        check_type(point, "Input point", Point3D)

        other_segment = Segment3D(
            self.start_pt,
            point
        )

        scale_factor = self.vector().scalar_projection(other_segment.vector()) / self.length()
        return self.pointAt(scale_factor)

    def pointDistance(self,
        point: Point3D
    ) -> numbers.Real:
        """
        Returns the point distance to the segment.

        :param point: the point to calculate the distance with
        :type point: Point3D
        :return: the distance of the point to the segment
        :rtype: numbers.Real

        Examples:
          >>> s = Segment3D(Point3D(0,0,0), Point3D(0,0,4))
          >>> s.pointDistance(Point3D(-17.2, 0.0, -49))
          17.2
          >>> s.pointDistance(Point3D(-17.2, 1.22, -49))
          17.24321315764553
        """

        check_type(point, "Input point", Point3D)

        #check_crs(self, point)

        point_projection = self.pointProjection(point)

        return point.distance(point_projection)

    def point_s(self,
                point: Point3D
                ) -> Optional[numbers.Real]:
        """
        Calculates the optional distance of the point along the segment.
        A zero value is for a point coinciding with the start point.
        Returns None if the point is not contained in the segment.

        :param point: the point to calculate the optional distance in the segment.
        :type point: Point3D
        :return: the the optional distance of the point along the segment.
        """

        check_type(point, "Input point", Point3D)

        #check_crs(self, point)

        if not self.contains_pt(point):
            return None

        return self.start_pt.distance(point)

    def scale(self,
        scale_factor
    ) -> 'Segment3D':
        """
        Scale a segment by the given scale_factor.
        Start point does not change.

        :param scale_factor: the scale factor, where 1 is the segment length.
        :type scale_factor: numbers.Real
        :return: Point at scale factor
        :rtype: Point3D
        """

        end_pt = self.pointAt(scale_factor)

        return Segment3D(
            self.start_pt,
            end_pt)

    def vertical_plane(self) -> Optional['CPlane3D']:
        """
        Returns the vertical Cartesian plane containing the segment.

        :return: the vertical Cartesian plane containing the segment.
        :rtype: Optional[CPlane3D].
        """

        if self.length_horizontal() == 0.0:  # collapsed segment
            return None
        elif self.length_horizontal() == 0.0:  # vertical segment
            return None

        # arbitrary point on the same vertical as end point

        section_final_pt_up = self.end_pt.shift(
            sx=0.0,
            sy=0.0,
            sz=1000.0)

        return CPlane3D.fromPoints(
            pt1=self.start_pt,
            pt2=self.end_pt,
            pt3=section_final_pt_up)

    def same_start(self,
                   another: 'Segment3D',
                   tol: numbers.Real = 1e-12
                   ) -> bool:
        """
        Check whether the two segments have the same start point.

        :param another: a segment to check for.
        :type another: Segment.
        :param tol: tolerance for distance between points.
        :type tol: numbers.Real.
        :return: whether the two segments have the same start point.
        :rtype: bool.

        Examples:
          >>> s1 = Segment3D(Point3D(0,0,0), Point3D(1,0,0))
          >>> s2 = Segment3D(Point3D(0,0,0), Point3D(0,1,0))
          >>> s1.same_start(s2)
          True
        """

        return self.start_pt.is_coincident(
            another=another.start_pt,
            tolerance=tol
        )

    def same_end(self,
                 another: 'Segment3D',
                 tol: numbers.Real = 1e-12
                 ) -> bool:
        """
        Check whether the two segments have the same end point.

        :param another: a segment to check for.
        :type another: Segment.
        :param tol: tolerance for distance between points.
        :type tol: numbers.Real.
        :return: whether the two segments have the same end point.
        :rtype: bool.

        Examples:
          >>> s1 = Segment3D(Point3D(0,0,0), Point3D(1,0,0))
          >>> s2 = Segment3D(Point3D(2,0,0), Point3D(1,0,0))
          >>> s1.same_end(s2)
          True
        """

        return self.end_pt.is_coincident(
            another=another.end_pt,
            tolerance=tol)

    def conn_to_other(self,
                      another: 'Segment3D',
                      tol: numbers.Real = 1e-12
                      ) -> bool:
        """
        Check whether the first segment is sequentially connected to the second one.

        :param another: a segment to check for.
        :type another: Segment.
        :param tol: tolerance for distance between points.
        :type tol: numbers.Real.
        :return: whether the first segment is sequentially connected to the second one.
        :rtype: bool.

        Examples:
          >>> s1 = Segment3D(Point3D(0,0,0), Point3D(1,0,0))
          >>> s2 = Segment3D(Point3D(1,0,0), Point3D(2,0,0))
          >>> s1.conn_to_other(s2)
          True
        """

        return self.end_pt.is_coincident(
            another=another.start_pt,
            tolerance=tol)

    def other_connected(self,
                        another: 'Segment3D',
                        tol: numbers.Real = 1e-12
                        ) -> bool:
        """
        Check whether the second segment is sequentially connected to the first one.

        :param another: a segment to check for.
        :type another: Segment.
        :param tol: tolerance for distance between points.
        :type tol: numbers.Real.
        :return: whether the second segment is sequentially connected to the first one.
        :rtype: bool.

        Examples:
          >>> s1 = Segment3D(Point3D(0,0,0), Point3D(1,0,0))
          >>> s2 = Segment3D(Point3D(-1,0,0), Point3D(0,0,0))
          >>> s1.other_connected(s2)
          True
        """

        return another.end_pt.is_coincident(
            another=self.start_pt,
            tolerance=tol)

    def segment_start_in(self,
        another: 'Segment3D'
    ) -> bool:
        """
        Check whether the second segment contains the first segment start point.

        :param another: a segment to check for.
        :type another: Segment.
        :return: whether the second segment contains the first segment start point.
        :rtype: bool.

        Examples:
          >>> s1 = Segment3D(Point3D(0,0,0), Point3D(1,0,0))
          >>> s2 = Segment3D(Point3D(-0.5,0,0), Point3D(0.5,0,0))
          >>> s1.segment_start_in(s2)
          True
          >>> s1 = Segment3D(Point3D(0,0,0), Point3D(1,1,1))
          >>> s1.segment_start_in(s2)
          True
          >>> s1 = Segment3D(Point3D(0,1,0), Point3D(1,1,1))
          >>> s1.segment_start_in(s2)
          False
          >>> s1 = Segment3D(Point3D(-1,-1,-1), Point3D(1,1,1))
          >>> s1.segment_start_in(s2)
          False
        """

        return another.contains_pt(self.start_pt)

    def segment_end_in(self,
        another: 'Segment3D'
    ) -> bool:
        """
        Check whether the second segment contains the first segment end point.

        :param another: a segment to check for.
        :type another: Segment.
        :return: whether the second segment contains the first segment end point.
        :rtype: bool.

        Examples:
          >>> s1 = Segment3D(Point3D(0,0,0), Point3D(1,0,0))
          >>> s2 = Segment3D(Point3D(-0.5,0,0), Point3D(0.5,0,0))
          >>> s1.segment_end_in(s2)
          False
          >>> s1 = Segment3D(Point3D(0,0,0), Point3D(1,1,1))
          >>> s1.segment_end_in(s2)
          False
          >>> s1 = Segment3D(Point3D(0,1,0), Point3D(1,1,1))
          >>> s2 = Segment3D(Point3D(1,1,1), Point3D(0.5,0,0))
          >>> s1.segment_end_in(s2)
          True
          >>> s1 = Segment3D(Point3D(-1,-1,3), Point3D(1,1,3))
          >>> s2 = Segment3D(Point3D(0,2,3), Point3D(2,0,3))
          >>> s1.segment_end_in(s2)
          True
        """

        return another.contains_pt(self.end_pt)

    def rotate(self,
        rotation_axis: 'RotationAxis',
        center_point: 'Point3D' = None
        ) -> 'Segment3D':
        """
        Rotates a segment.
        :param rotation_axis:
        :param center_point:
        :return: the rotated segment
        :rtype: Segment3D

        Examples:
        >>> seg = Segment3D(Point3D(0,0,0), Point3D(0,0,1))
        >>> rot_ax = RotationAxis(0, 0, 90)
        >>> seg.rotate(rot_ax)
        Segment3D(start_pt=Point3D(0.0000, 0.0000, 0.0000), end_pt=Point3D(1.0000, 0.0000, 0.0000))
        >>> rot_ax = RotationAxis(0, 0, 180)
        >>> centr_pt = Point3D(0,0,0.5)
        >>> seg.rotate(rotation_axis=rot_ax, center_point=centr_pt)
        Segment3D(start_pt=Point3D(-0.0000, 0.0000, 1.0000), end_pt=Point3D(0.0000, 0.0000, 0.0000))
        >>> seg = Segment3D(Point3D(0,0,0), Point3D(1,1,0))
        >>> centr_pt = Point3D(1,0,0)
        >>> rot_ax = RotationAxis(0, 90, 90)
        >>> seg.rotate(rotation_axis=rot_ax, center_point=centr_pt)
        Segment3D(start_pt=Point3D(1.0000, 1.0000, 0.0000), end_pt=Point3D(2.0000, 0.0000, -0.0000))
        >>> seg = Segment3D(Point3D(1,1,1), Point3D(0,0,0))
        >>> rot_ax = RotationAxis(135, 0, 180)
        >>> centr_pt = Point3D(0.5,0.5,0.5)
        >>> seg.rotate(rotation_axis=rot_ax, center_point=centr_pt)
        Segment3D(start_pt=Point3D(0.0000, 0.0000, 0.0000), end_pt=Point3D(1.0000, 1.0000, 1.0000))
        """

        start_pt, end_pt = self

        rotated_start_pt = start_pt.rotate(
            rotation_axis=rotation_axis,
            center_point=center_point
        )

        rotated_end_pt = end_pt.rotate(
            rotation_axis=rotation_axis,
            center_point=center_point
        )

        return Segment3D(
            start_pt=rotated_start_pt,
            end_pt=rotated_end_pt
        )

    @classmethod
    def random(cls,
        lower_boundary: float = -MAX_SCALAR_VALUE,
        upper_boundary: float = MAX_SCALAR_VALUE):
        """
        Creates a random segment.

        :return: random segment
        :rtype: Segment3D
        """

        return cls(
            start_pt=Point3D.random(lower_boundary, upper_boundary),
            end_pt=Point3D.random(lower_boundary, upper_boundary)
        )

    def densify_as_line3d(self,
                          densify_distance
                          ) -> 'Line3D':
        """
        Densify a segment by adding additional points
        separated a distance equal to densify_distance.
        The result is no longer a Segment instance, instead it is a Line instance.

        :param densify_distance: float
        :return: a Line3D
        """

        length3d = self.length()

        segment_versor = self.as_vector().versor()
        generator_vector = segment_versor.scale(densify_distance)

        interpolated_line = Line3D(
            pts=[self.start_pt])

        n = 0
        while True:
            n += 1
            shift_vector = generator_vector.scale(n)
            new_pt = self.start_pt.shift(
                shift_vector.x,
                shift_vector.y,
                shift_vector.z
            )
            distance = self.start_pt.distance(new_pt)
            if distance >= length3d:
                break
            interpolated_line.add_pt(new_pt)
        interpolated_line.add_pt(self.end_pt)

        return interpolated_line

    def densify_as_pts3d(self,
                         densify_distance
                         ) -> List[Point3D]:

        return self.densify_as_line3d(densify_distance=densify_distance).pts()

    def densify_as_steps3d(self,
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


def point_or_segment3d(
        point1: Point3D,
        point2: Point3D,
        tol: numbers.Real = PRACTICAL_MIN_DIST
) -> Union[Point3D, Segment3D]:
    """
    Creates a point or segment based on the points distance.

    :param point1: first input point.
    :type point1: Point.
    :param point2: second input point.
    :type point2: Point.
    :param tol: distance tolerance between the two points.
    :type tol: numbers.Real.
    :return: point or segment based on their distance.
    :rtype: PointOrSegment.
    :raise: Exception.
    """

    check_type(point1, "First point", Point3D)
    check_type(point2, "Second point", Point3D)

    if point1.distance(point2) <= tol:
        return Point3D(
            x=(point1.x + point2.x) / 2,
            y=(point1.y + point2.y) / 2,
            z=(point1.z + point2.z) / 2
        )
    else:
        return Segment3D(
            start_pt=point1,
            end_pt=point2
        )


def intersect_segments3d(
    segment1: Segment3D,
    segment2: Segment3D,
    tol: numbers.Real = PRACTICAL_MIN_DIST
) -> Optional[Union[Point3D, Segment3D]]:
    """
    Determines the optional point or segment intersection between the segment pair.

    :param segment1: the first segment
    :param segment2: the second segment
    :param tol: the distance tolerance for collapsing a intersection segment into a point
    :return: the optional point or segment intersection between the segment pair.

    Examples:
      >>> s2 = Segment3D(Point3D(0,0,0), Point3D(1,0,0))
      >>> s1 = Segment3D(Point3D(0,0,0), Point3D(1,0,0))
      >>> intersect_segments3d(s1, s2)
      Segment3D(start_pt=Point3D(0.0000, 0.0000, 0.0000), end_pt=Point3D(1.0000, 0.0000, 0.0000))
      >>> s1 = Segment3D(Point3D(-2,0,0), Point3D(-1,0,0))
      >>> intersect_segments3d(s1, s2) is None
      True
      >>> s1 = Segment3D(Point3D(-2,0,0), Point3D(0,0,0))
      >>> intersect_segments3d(s1, s2)
      Point3D(0.0000, 0.0000, 0.0000)
      >>> s1 = Segment3D(Point3D(-2,0,0), Point3D(0.5,0,0))
      >>> intersect_segments3d(s1, s2)
      Segment3D(start_pt=Point3D(0.0000, 0.0000, 0.0000), end_pt=Point3D(0.5000, 0.0000, 0.0000))
      >>> s1 = Segment3D(Point3D(-2,0,0), Point3D(1,0,0))
      >>> intersect_segments3d(s1, s2)
      Segment3D(start_pt=Point3D(0.0000, 0.0000, 0.0000), end_pt=Point3D(1.0000, 0.0000, 0.0000))
      >>> s1 = Segment3D(Point3D(-2,0,0), Point3D(2,0,0))
      >>> intersect_segments3d(s1, s2)
      Segment3D(start_pt=Point3D(0.0000, 0.0000, 0.0000), end_pt=Point3D(1.0000, 0.0000, 0.0000))
      >>> s1 = Segment3D(Point3D(0,0,0), Point3D(0.5,0,0))
      >>> intersect_segments3d(s1, s2)
      Segment3D(start_pt=Point3D(0.0000, 0.0000, 0.0000), end_pt=Point3D(0.5000, 0.0000, 0.0000))
      >>> s1 = Segment3D(Point3D(0.25,0,0), Point3D(0.75,0,0))
      >>> intersect_segments3d(s1, s2)
      Segment3D(start_pt=Point3D(0.2500, 0.0000, 0.0000), end_pt=Point3D(0.7500, 0.0000, 0.0000))
      >>> s1 = Segment3D(Point3D(0.25,0,0), Point3D(1,0,0))
      >>> intersect_segments3d(s1, s2)
      Segment3D(start_pt=Point3D(0.2500, 0.0000, 0.0000), end_pt=Point3D(1.0000, 0.0000, 0.0000))
      >>> s1 = Segment3D(Point3D(0.25,0,0), Point3D(1.25,0,0))
      >>> intersect_segments3d(s1, s2)
      Segment3D(start_pt=Point3D(0.2500, 0.0000, 0.0000), end_pt=Point3D(1.0000, 0.0000, 0.0000))
      >>> s1 = Segment3D(Point3D(0,0,0), Point3D(1.25,0,0))
      >>> intersect_segments3d(s1, s2)
      Segment3D(start_pt=Point3D(0.0000, 0.0000, 0.0000), end_pt=Point3D(1.0000, 0.0000, 0.0000))
      >>> s1 = Segment3D(Point3D(1,0,0), Point3D(1.25,0,0))
      >>> intersect_segments3d(s1, s2)
      Point3D(1.0000, 0.0000, 0.0000)
      >>> s2 = Segment3D(Point3D(0,0,0), Point3D(1,1,1))
      >>> s1 = Segment3D(Point3D(0.25,0.25,0.25), Point3D(0.75,0.75,0.75))
      >>> intersect_segments3d(s1, s2)
      Segment3D(start_pt=Point3D(0.2500, 0.2500, 0.2500), end_pt=Point3D(0.7500, 0.7500, 0.7500))
      >>> s1 = Segment3D(Point3D(0.25,0.25,0.25), Point3D(1.75,1.75,1.75))
      >>> intersect_segments3d(s1, s2)
      Segment3D(start_pt=Point3D(0.2500, 0.2500, 0.2500), end_pt=Point3D(1.0000, 1.0000, 1.0000))
      >>> s1 = Segment3D(Point3D(0.25,0.25,0.25), Point3D(1.75,0,1.75))
      >>> intersect_segments3d(s1, s2)
      Point3D(0.2500, 0.2500, 0.2500)
      >>> s1 = Segment3D(Point3D(0.25,1,0.25), Point3D(0.75,0.75,0.75))
      >>> intersect_segments3d(s1, s2)
      Point3D(0.7500, 0.7500, 0.7500)
      >>> s2 = Segment3D(Point3D(-1,-1,-1), Point3D(1,1,1))
      >>> s1 = Segment3D(Point3D(-1,1,1), Point3D(1,-1,-1))
      >>> intersect_segments3d(s1, s2)
      Point3D(0.0000, 0.0000, 0.0000)
    """

    check_type(segment1, "First segment", Segment3D)
    check_type(segment2, "Second segment", Segment3D)

    #check_crs(segment1, segment2)

    s1_startpt_inside = segment1.segment_start_in(segment2)
    s2_startpt_inside = segment2.segment_start_in(segment1)

    s1_endpt_inside = segment1.segment_end_in(segment2)
    s2_endpt_inside = segment2.segment_end_in(segment1)

    elements = [s1_startpt_inside, s2_startpt_inside, s1_endpt_inside, s2_endpt_inside]

    if all(elements):
        return segment1.clone()

    if s1_startpt_inside and s1_endpt_inside:
        return segment1.clone()

    if s2_startpt_inside and s2_endpt_inside:
        return segment2.clone()

    if s1_startpt_inside and s2_startpt_inside:
        return point_or_segment3d(
            segment1.start_pt,
            segment2.start_pt,
            tol=tol
        )

    if s1_startpt_inside and s2_endpt_inside:
        return point_or_segment3d(
            segment1.start_pt,
            segment2.end_pt,
            tol=tol
        )

    if s1_endpt_inside and s2_startpt_inside:
        return point_or_segment3d(
            segment2.start_pt,
            segment1.end_pt,
            tol=tol
        )

    if s1_endpt_inside and s2_endpt_inside:
        return point_or_segment3d(
            segment1.end_pt,
            segment2.end_pt,
            tol=tol
        )

    if s1_startpt_inside:
        return segment1.start_pt.clone()

    if s1_endpt_inside:
        return segment1.end_pt.clone()

    if s2_startpt_inside:
        return segment2.start_pt.clone()

    if s2_endpt_inside:
        return segment2.end_pt.clone()

    shortest_segm_or_pt = shortest_segment_or_point3d(
        segment1,
        segment2,
        tol=tol
    )

    if not shortest_segm_or_pt:
        return None

    if not isinstance(shortest_segm_or_pt, Point3D):
        return None

    inters_pt = shortest_segm_or_pt

    if not segment1.contains_pt(inters_pt):
        return None

    if not segment2.contains_pt(inters_pt):
        return None

    return inters_pt


class PointSegmentCollection3D(list):
    """
    Collection of point or segment elements.

    """

    def __init__(
            self,
            geoms: Optional[List[Union[Point3D, Segment3D]]] = None,
            # epsg_code: Optional[numbers.Integral] = None
    ):

        if geoms is not None:

            for geom in geoms:
                check_type(geom, "Spatial element", (Point3D, Segment3D))

        """
        if epsg_code is not None:
            check_type(
                var=epsg_code,
                name="EPSG code",
                expected_types=numbers.Integral
            )

        if geoms is not None and epsg_code is not None:

            for geom in geoms:
                check_epsg(
                    spatial_element=geom,
                    epsg_code=epsg_code
                )

        elif geoms is not None and len(geoms) > 0:

            epsg_code = geoms[0].epsg_code()
        """

        if geoms is not None and len(geoms) > 0:

            super(PointSegmentCollection3D, self).__init__(geoms)

        else:

            super(PointSegmentCollection3D, self).__init__()

        # self.epsg_code = epsg_code

    def append(self,
               spatial_element: Union[Point3D, Segment3D]
               ) -> None:

        check_type(
            var=spatial_element,
            name="Spatial element",
            expected_types=(Point3D, Segment3D)
        )

        """
        if self.epsg_code is not None:

            check_epsg(
                spatial_element=spatial_element,
                epsg_code=self.epsg_code
            )

        else:

            self.epsg_code = spatial_element.epsg_code()
        """

        self.append(spatial_element)


class Line3D:
    """
    A line.
    """

    def __init__(self,
                 pts: Optional[List[Point3D]] = None):
        """

        """

        if pts is not None:

            check_type(pts, "List", list)
            for el in pts:
                check_type(el, "Point3D", Point3D)

            self._pts = pts

        else:

            self._pts = []

    def __repr__(self) -> str:
        """
        Represents a Line instance as a shortened text.

        :return: a textual shortened representation of a Line instance.
        :rtype: str.
        """

        num_points = self.num_pts()

        if num_points == 0:
            txt = "Empty Line3D"
        else:
            x1, y1, z1 = self.start_pt()
            if num_points == 1:
                txt = f"Line3D with unique point: {x1:.4f}, {y1:.4f}, {z1:.4f}"
            else:
                x2, y2, z2 = self.end_pt()
                txt = f"Line3D with {self.num_pts()} points: ({x1:.4f}, {y1:.4f}, {z1:.4f}) ... ({x2:.4f}, {y2:.4f}, {z2:.4f})"

        return txt

    def pts(self):
        return self._pts

    def pt(self,
           ndx: numbers.Integral):
        """

        """

        return self._pts[ndx]

    def start_pt(self) -> Optional[Point3D]:
        """
        Return the first point of a Line or None when no points.

        :return: the first point or None.
        """

        return self.pt(0) if self.num_pts() > 0 else None

    def end_pt(self) -> Optional[Point3D]:
        """
        Return the last point of a Line or None when no points.

        :return: the last point or None.
        """

        return self.pt(-1) if self.num_pts() > 0 else None

    def add_pt(self,
               pt: Point3D):

        self._pts.append(pt)

    def num_pts(self):
        return len(self._pts)

    def segment(self,
        ndx: numbers.Integral
    ) -> Optional[Segment3D]:
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
            return Segment3D(
                start_pt=self.pt(ndx),
                end_pt=self.pt(ndx + 1)
            )

    def __iter__(self):
        """
        Return each element of a Line, i.e., its segments.
        """

        return (self.segment(i) for i in range(self.num_pts()-1))

    def x_list(self) -> List[numbers.Real]:

        return list(map(lambda pt: pt.x, self._pts))

    def y_list(self) -> List[numbers.Real]:

        return list(map(lambda pt: pt.y, self._pts))

    def x_array(self):

        return np.asarray([pt.x for pt in self.pts()])

    def y_array(self):

        return np.asarray([pt.y for pt in self.pts()])

    def z_array(self):

        return np.asarray([pt.z for pt in self.pts()])

    def xy_arrays(self):

        return self.x_array, self.y_array

    def x_min(self):
        return np.nanmin(list(map(lambda pt: pt.x, self._pts)))

    def x_max(self):
        return np.nanmax(list(map(lambda pt: pt.x, self._pts)))

    def y_min(self):
        return np.nanmin(list(map(lambda pt: pt.y, self._pts)))

    def y_max(self):
        return np.nanmax(list(map(lambda pt: pt.y, self._pts)))

    def z_min(self):
        return np.nanmin(list(map(lambda pt: pt.z, self._pts)))

    def z_max(self):
        return np.nanmax(list(map(lambda pt: pt.z, self._pts)))

    def as_segments(self):
        """
        Convert to a list of segments.

        :return: list of Segment objects
        """

        pts_pairs = zip(self.pts()[:-1], self.pts()[1:])

        segments = [Segment3D(pt_a, pt_b) for (pt_a, pt_b) in pts_pairs]

        return segments

    '''
    def densify_2d_line(self, sample_distance) -> 'Points':
        """
        Densify a line into a new line instance,
        using the provided sample distance.
        Returned Line instance has coincident successive points removed.

        :param sample_distance: numbers.Real
        :return: Line instance
        """

        if sample_distance <= 0.0:
            raise Exception(f"Sample distance must be positive. {sample_distance} received")

        segments = self.as_segments()

        densified_line_list = [segment.densify2d_asLine(sample_distance) for segment in segments]

        densifyied_multiline = MultiLine(densified_line_list)

        densifyied_line = densifyied_multiline.to_line()

        densifyied_line_wo_coinc_pts = densifyied_line.remove_coincident_points()

        return densifyied_line_wo_coinc_pts
    '''

    def join(self, another) -> 'Line3D':
        """
        Joins together two lines and returns the join as a new line without point changes,
        with possible overlapping points
        and orientation mismatches between the two original lines
        """

        return Line3D(self.pts() + another.pts())

    def length(self) -> numbers.Real:

        length = 0.0
        for ndx in range(self.num_pts() - 1):
            length += self.pt(ndx).distance(self.pt(ndx + 1))
        return length

    def length_2d(self) -> numbers.Real:

        length = 0.0
        for ndx in range(self.num_pts() - 1):
            length += self.pt(ndx).to2d().distance(self.pt(ndx + 1).to2d())
        return length

    def step_delta_z(self) -> List[numbers.Real]:
        """
        Return the difference in elevation between consecutive points:
        z[ndx+1] - z[ndx]

        :return: a list of height differences.
        :rtype: list of floats.
        """

        delta_z = [0.0]

        for ndx in range(1, self.num_pts()):
            delta_z.append(self.pt(ndx).z - self.pt(ndx - 1).z)

        return delta_z

    def step_lengths_3d(self) -> List[numbers.Real]:
        """
        Returns the point-to-point 3D distances.
        It is the distance between a point and its previous one.
        The list has the same lenght as the source point list.

        :return: the individual 3D segment lengths.
        :rtype: list of floats.

        Examples:
        """

        step_length_list = [0.0]
        for ndx in range(1, self.num_pts()):
            length = self.pt(ndx).distance(self.pt(ndx - 1))
            step_length_list.append(length)

        return step_length_list

    '''
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
    '''

    def incremental_length_2d(self):

        lIncrementalLengths = []
        length = 0.0
        lIncrementalLengths.append(length)
        for ndx in range(self.num_pts() - 1):
            length += self.pts()[ndx].dist_2d(self.pts()[ndx + 1])
            lIncrementalLengths.append(length)

        return np.asarray(lIncrementalLengths)

    def incremental_length_3d(self) -> List[numbers.Real]:
        """
        Returns the accumulated 3D segment lengths.

        :return: accumulated 3D segment lenghts
        :rtype: list of floats.
        """

        return list(itertools.accumulate(self.step_lengths_3d()))

    '''
    def incremental_length_2d(self) -> List[numbers.Real]:
        """
        Returns the accumulated 2D segment lengths.

        :return: accumulated 2D segment lenghts
        :rtype: list of floats.
        """

        return list(itertools.accumulate(self.step_lengths_2d()))
    '''

    def reversed(self) -> 'Line3D':
        """
        Return a Line instance with reversed point list.

        :return: a new Line instance.
        :rtype: Line.
        """

        pts = [pt.clone() for pt in self.pts()]
        pts.reverse()

        return Line3D(
            pts=pts
        )

    def slopes_degr(self) -> List[Optional[numbers.Real]]:
        """
        Calculates the slopes (in degrees) of each Line segment.
        The first value is the slope of the first segment.
        The last value, always None, is the slope of the segment starting at the last point.
        The number of elements is equal to the number of points in the Line.

        :return: list of slopes (degrees).
        :rtype: List[Optional[numbers.Real]].
        """

        lSlopes = []

        segments = self.as_segments()
        for segment in segments:
            vector = segment.vector()
            lSlopes.append(-vector.slope_degr())  # minus because vector convention is positive downward

        lSlopes.append(None)  # None refers to the slope of the Segment starting with the last point

        return lSlopes

    def slopes_stats(self) -> Dict:
        """
        Returns the line directional slope statistics.

        :return: the statistics parameters: min, max, mean, var, std.
        """

        return get_statistics(self.slopes_degr())

    def abs_slopes_degr(self) -> List[Optional[numbers.Real]]:

        return [abs(val) for val in self.slopes_degr()]

    def dir_slopes(self) -> np.ndarray:

        lSlopes = []
        for ndx in range(self.num_pts() - 1):
            segment_start_pt = self.pts()[ndx]
            segment_end_pt = self.pts()[ndx + 1]
            if np.isnan(segment_start_pt.z) or np.isnan(segment_end_pt.z):
                lSlopes.append(np.nan)
            else:
                vector = Segment3D(self.pts()[ndx], self.pts()[ndx + 1]).vector()
                lSlopes.append(-vector.slope_degr())  # minus because vector convention is positive downward
        lSlopes.append(np.nan)  # slope value for last point is unknown

        return np.asarray(lSlopes)

    def absolute_slopes(self) -> np.ndarray:

        return np.asarray(list(map(abs, self.dir_slopes())))

    def abs_slopes_stats(self) -> Dict:
        """
        Returns the line absolute slopes statistics.

        :return: the statistics parameters: min, max, mean, var, std.
        :rtype: Dictionary.
        """

        return get_statistics(self.abs_slopes_degr())

    def extremes_distance_3d(self) -> numbers.Real:
        """
        Calculate the 3D distance between start and end points.

        :return: the 3D distance between start and end points
        :rtype: numbers.Real
        """

        return self.pt(-1).distance(self.pt(0))

    '''
    def extremes_distance_2d(self) -> numbers.Real:
        """
        Calculate the 2D distance between start and end points.

        :return: the 2D distance between start and end points
        """

        return self.end_pt().dist2DWith(self.start_pt())
    '''

    def is_closed(self,
                  tolerance: numbers.Real = MIN_SEPARATION_THRESHOLD
                  ) -> bool:
        """
        Determine if the line is 3D-closed.

        :param tolerance: the tolerance for considering the line closed
        :type tolerance: numbers.Real
        :return: whether the line is to be considered 3D-closed
        :rtype: bool
        """

        return self.pt(-1).is_coincident(self.pt(0), tolerance=tolerance)

    '''
    def isClosed_2d(self,
        tolerance: numbers.Real = MIN_SEPARATION_THRESHOLD
    ) -> bool:
        """
        Determine if the line is 2D-closed.

        :param tolerance: the tolerance for considering the line closed
        :return: whether the line is to be considered 2D-closed
        """

        return self.end_pt().isCoinc2D(self.start_pt(), tolerance=tolerance)
    '''

    def walk_backward(self) -> 'Line3D':
        """
        Create a new line by walking the line backward from the last point up to the first and thus closing it.

        :return: a closed line with zero area
        :rtype: 'Line'
        """

        return Line3D(self.pts() + self.reversed()[1:])

    def clone(self) -> 'Line3D':
        """
        Clone a line.

        :return: the cloned line
        :rtype: Line3D
        """

        return Line3D(self.pts())

    '''
    def close_2d(self) -> 'Points':
        """
        Return a line that is 2D-closed.

        :return: a 2D-closed line
        :rtype: Points
        """

        line = self.clone()

        if not line.isClosed_2d():

            line.add_pt(line.start_pt())

        return line
    '''

    def close_3d(self) -> 'Line3D':
        """
        Return a line that is 3D-closed.

        :return: a 3D-closed line
        :rtype: Line3D
        """

        line = self.clone()

        if not line.is_closed():

            line.add_pt(line.start_pt())

        return line

    def remove_coincident_points(self) -> Optional['Line3D']:
        """
        Remove coincident successive points

        :return: Line instance
        :rtype: Optional[Line3D]
        """

        if self.num_pts() == 0:
            return

        new_line = Line3D(
            pts=[self.pt(0)]
        )

        for ndx in range(1, self.num_pts()):
            if not self.pt(ndx).is_coincident(new_line.pt(-1)):
                new_line.add_pt(self.pt(ndx))

        return new_line

    def intersectSegment(self,
        segment: Segment3D
    ) -> Optional[PointSegmentCollection3D]:
        """
        Calculates the possible intersection between the line and a provided segment.

        :param segment: the input segment
        :return: the optional intersections, points or segments
        :raise: Exception
        """

        if self.num_pts() <= 1:
            return

        check_type(segment, "Input segment", Segment3D)

        intersections = [intersect_segments3d(curr_segment, segment) for curr_segment in self if curr_segment is not None]
        intersections = list(filter(lambda val: val is not None, intersections))
        intersections = PointSegmentCollection3D(intersections)

        return intersections


class MultiLine3D:
    """
    MultiLine is a list of Line objects
    """

    def __init__(self, lines_list=None):

        if lines_list is None:
            lines_list = []
        self._lines = lines_list

    @property
    def lines(self):

        return self._lines

    def add(self, line):

        return MultiLine3D(self.lines + [line])

    def clone(self):

        return MultiLine3D(self.lines)

    @property
    def num_parts(self):

        return len(self.lines)

    @property
    def num_points(self):

        num_points = 0
        for line in self.lines:
            num_points += line.num_pts

        return num_points

    @property
    def x_min(self):

        return np.nanmin([line.x_min for line in self.lines])

    @property
    def x_max(self):

        return np.nanmax([line.x_max for line in self.lines])

    @property
    def y_min(self):

        return np.nanmin([line.y_min for line in self.lines])

    @property
    def y_max(self):

        return np.nanmax([line.y_max for line in self.lines])

    @property
    def z_min(self):

        return np.nanmin([line.z_min for line in self.lines])

    @property
    def z_max(self):

        return np.nanmax([line.z_max for line in self.lines])

    def is_continuous(self):

        for line_ndx in range(len(self._lines) - 1):
            if not self.lines[line_ndx].pts[-1].coincident(self.lines[line_ndx + 1].pts[0]) or \
               not self.lines[line_ndx].pts[-1].coincident(self.lines[line_ndx + 1].pts[-1]):
                return False

        return True

    def is_unidirectional(self):

        for line_ndx in range(len(self.lines) - 1):
            if not self.lines[line_ndx].pts[-1].coincident(self.lines[line_ndx + 1].pts[0]):
                return False

        return True

    def to_line(self):

        return Line3D([point for line in self.lines for point in line.pts])

    '''
    def crs_project(self, srcCrs, destCrs):

        lines = []
        for line in self.lines:
            lines.append(line.crs_project(srcCrs, destCrs))

        return MultiLine4D(lines)
    '''

    '''
    def densify_2d_multiline(self, sample_distance):

        lDensifiedLines = []
        for line in self.lines:
            lDensifiedLines.append(line.densify_2d_line(sample_distance))

        return MultiLine4D(lDensifiedLines)
    '''

    def remove_coincident_points(self):

        cleaned_lines = []
        for line in self.lines:
            cleaned_lines.append(line.remove_coincident_points())

        return MultiLine3D(cleaned_lines)


def shortest_segment_or_point3d(
    first_segment: Segment3D,
    second_segment: Segment3D,
    tol: numbers.Real = PRACTICAL_MIN_DIST
) -> Optional[Union[Segment3D, Point3D]]:

    """
    Calculates the optional shortest segment - or the intersection point - between two lines represented by two segments.

    Adapted from:
        http://paulbourke.net/geometry/pointlineplane/

    C code from:
        http://paulbourke.net/geometry/pointlineplane/lineline.c
[
    typedef struct {
    double x,y,z;
    } XYZ;

    /*
    Calculate the line segment PaPb that is the shortest route between
    two lines P1P2 and P3P4. Calculate also the values of mua and mub where
      Pa = P1 + mua (P2 - P1)
      Pb = P3 + mub (P4 - P3)
    Return FALSE if no solution exists.
    */
    int LineLineIntersect(
    XYZ p1,XYZ p2,XYZ p3,XYZ p4,XYZ *pa,XYZ *pb,
    double *mua, double *mub)
    {
    XYZ p13,p43,p21;
    double d1343,d4321,d1321,d4343,d2121;
    double numer,denom;

    p13.x = p1.x - p3.x;
    p13.y = p1.y - p3.y;
    p13.z = p1.z - p3.z;
    p43.x = p4.x - p3.x;
    p43.y = p4.y - p3.y;
    p43.z = p4.z - p3.z;
    if (ABS(p43.x) < EPS && ABS(p43.y) < EPS && ABS(p43.z) < EPS)
      return(FALSE);
    p21.x = p2.x - p1.x;
    p21.y = p2.y - p1.y;
    p21.z = p2.z - p1.z;
    if (ABS(p21.x) < EPS && ABS(p21.y) < EPS && ABS(p21.z) < EPS)
      return(FALSE);

    d1343 = p13.x * p43.x + p13.y * p43.y + p13.z * p43.z;
    d4321 = p43.x * p21.x + p43.y * p21.y + p43.z * p21.z;
    d1321 = p13.x * p21.x + p13.y * p21.y + p13.z * p21.z;
    d4343 = p43.x * p43.x + p43.y * p43.y + p43.z * p43.z;
    d2121 = p21.x * p21.x + p21.y * p21.y + p21.z * p21.z;

    denom = d2121 * d4343 - d4321 * d4321;
    if (ABS(denom) < EPS)
      return(FALSE);
    numer = d1343 * d4321 - d1321 * d4343;

    *mua = numer / denom;
    *mub = (d1343 + d4321 * (*mua)) / d4343;

    pa->x = p1.x + *mua * p21.x;
    pa->y = p1.y + *mua * p21.y;
    pa->z = p1.z + *mua * p21.z;
    pb->x = p3.x + *mub * p43.x;
    pb->y = p3.y + *mub * p43.y;
    pb->z = p3.z + *mub * p43.z;

    return(TRUE);
    }

    :param first_segment: the first segment
    :param second_segment: the second segment
    :param tol: tolerance value for collapsing a segment into the midpoint.
    :return: the optional shortest segment or an intersection point.
    """

    check_type(second_segment, "Second Cartesian line", Segment3D)

    p1 = first_segment.start_pt
    p2 = first_segment.end_pt

    p3 = second_segment.start_pt
    p4 = second_segment.end_pt

    p13 = Point3D(
        x=p1.x - p3.x,
        y=p1.y - p3.y,
        z=p1.z - p3.z
    )

    p43 = Point3D(
        x=p4.x - p3.x,
        y=p4.y - p3.y,
        z=p4.z - p3.z
    )

    if p43.asVect().is_close_to_zero:
        return None

    p21 = Point3D(
        x=p2.x - p1.x,
        y=p2.y - p1.y,
        z=p2.z - p1.z,
    )

    if p21.asVect().is_close_to_zero:
        return None

    d1343 = p13.x * p43.x + p13.y * p43.y + p13.z * p43.z
    d4321 = p43.x * p21.x + p43.y * p21.y + p43.z * p21.z
    d1321 = p13.x * p21.x + p13.y * p21.y + p13.z * p21.z
    d4343 = p43.x * p43.x + p43.y * p43.y + p43.z * p43.z
    d2121 = p21.x * p21.x + p21.y * p21.y + p21.z * p21.z

    denom = d2121 * d4343 - d4321 * d4321

    if fabs(denom) < MIN_SCALAR_VALUE:
        return None

    numer = d1343 * d4321 - d1321 * d4343

    mua = numer / denom
    mub = (d1343 + d4321 * mua) / d4343

    pa = Point3D(
        x=p1.x + mua * p21.x,
        y=p1.y + mua * p21.y,
        z=p1.z + mua * p21.z
    )

    pb = Point3D(
        x=p3.x + mub * p43.x,
        y=p3.y + mub * p43.y,
        z=p3.z + mub * p43.z
    )

    intersection = point_or_segment3d(
        point1=pa,
        point2=pb,
        tol=tol
    )

    return intersection

'''
class ParamLine3D(object):
    """
    parametric line
    srcPt: source Point
    l, m, n: .....
    """

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
        x1, y1, z1 = self._srcPt.x, self._srcPt.y, self._srcPt.z
        l, m, n = self._l, self._m, self._n

        # Cartesian plane parameters
        a, b, c, d = cartes_plane.a, cartes_plane.b, cartes_plane.c, cartes_plane.d

        try:
            k = (a * x1 + b * y1 + c * z1 + d) / (a * l + b * m + c * n)
        except ZeroDivisionError:
            return None

        return Point3D(x1 - l * k,
                       y1 - m * k,
                       z1 - n * k)
'''

def eq_xy_pair(xy_pair_1, xy_pair_2):

    if xy_pair_1[0] == xy_pair_2[0] and xy_pair_1[1] == xy_pair_2[1]:
        return True

    return False

'''
def remove_equal_consecutive_xypairs(xy_list):

    out_xy_list = [xy_list[0]]

    for n in range(1, len(xy_list)):
        if not eq_xy_pair(xy_list[n], out_xy_list[-1]):
            out_xy_list.append(xy_list[n])

    return out_xy_list
'''


class CPlane3D:
    """
    Cartesian plane.
    Expressed by equation:
    ax + by + cz + d = 0

    Note: CPlane3D is locational - its position in space is defined.
    This contrast with Plane, defined just by its attitude, but with undefined position

    """

    def __init__(self,
                 a: numbers.Real,
                 b: numbers.Real,
                 c: numbers.Real,
                 d: numbers.Real
                 ):

        if not isinstance(a, numbers.Real):
            raise Exception("Input value a must be float or int but is {}".format(type(a)))
        if not isinstance(b, numbers.Real):
            raise Exception("Input value b must be float or int but is {}".format(type(b)))
        if not isinstance(c, numbers.Real):
            raise Exception("Input value c must be float or int but is {}".format(type(c)))
        if not isinstance(d, numbers.Real):
            raise Exception("Input value d must be float or int but is {}".format(type(d)))

        norm = sqrt(a*a + b*b + c*c)
        self._a = float(a) / norm
        self._b = float(b) / norm
        self._c = float(c) / norm
        self._d = float(d) / norm

    def a(self) -> numbers.Real:
        """
        Return a coefficient of a CPlane3D instance.

        Example:
          >>> CPlane3D(1, 0, 0, 2).a()
          1.0
        """

        return self._a

    def b(self) -> numbers.Real:
        """
        Return b coefficient of a CPlane3D instance.

        Example:
          >>> CPlane3D(1, 4, 0, 2).b()
          0.9701425001453319
        """

        return self._b

    def c(self) -> numbers.Real:
        """
        Return a coefficient of a CPlane3D instance.

        Example:
          >>> CPlane3D(1, 0, 5.4, 2).c()
          0.9832820049844602
        """

        return self._c

    def d(self) -> numbers.Real:
        """
        Return a coefficient of a CPlane3D instance.

        Example:
          >>> CPlane3D(1, 0, 0, 2).d()
          2.0
        """

        return self._d

    def v(self) -> Tuple[numbers.Real, numbers.Real, numbers.Real, numbers.Real]:
        """
        Return coefficients of a CPlane3D instance.

        Example:
          >>> CPlane3D(1, 1, 7, -4).v()
          (0.14002800840280097, 0.14002800840280097, 0.9801960588196068, -0.5601120336112039)
        """

        return self.a(), self.b(), self.c(), self.d()

    @classmethod
    def fromPoints(cls, pt1, pt2, pt3) -> 'CPlane3D':
        """
        Create a CPlane3D from three given Point instances.

        Example:
          >>> CPlane3D.fromPoints(Point3D(0, 0, 0), Point3D(1, 0, 0), Point3D(0, 1, 0))
          CPlane3D(0.0000, 0.0000, 1.0000, 0.0000)
          >>> CPlane3D.fromPoints(Point3D(0, 0, 0), Point3D(1, 0, 0), Point3D(0, 1, 0))
          CPlane3D(0.0000, 0.0000, 1.0000, 0.0000)
          >>> CPlane3D.fromPoints(Point3D(0, 0, 0), Point3D(0, 1, 0), Point3D(0, 0, 1))
          CPlane3D(1.0000, 0.0000, 0.0000, 0.0000)
          >>> CPlane3D.fromPoints(Point3D(1,2,3), Point3D(2,3,4), Point3D(-1,7,-2))
          CPlane3D(-0.7956, 0.2387, 0.5569, -1.3524)
        """

        if not (isinstance(pt1, Point3D)):
            raise Exception("First input point should be Point but is {}".format(type(pt1)))

        if not (isinstance(pt2, Point3D)):
            raise Exception("Second input point should be Point but is {}".format(type(pt2)))

        if not (isinstance(pt3, Point3D)):
            raise Exception("Third input point should be Point but is {}".format(type(pt3)))

        matr_a = np.array(
            [[pt1.y, pt1.z, 1],
             [pt2.y, pt2.z, 1],
             [pt3.y, pt3.z, 1]])

        matr_b = - np.array(
            [[pt1.x, pt1.z, 1],
             [pt2.x, pt2.z, 1],
             [pt3.x, pt3.z, 1]])

        matr_c = np.array(
            [[pt1.x, pt1.y, 1],
             [pt2.x, pt2.y, 1],
             [pt3.x, pt3.y, 1]])

        matr_d = - np.array(
            [[pt1.x, pt1.y, pt1.z],
             [pt2.x, pt2.y, pt2.z],
             [pt3.x, pt3.y, pt3.z]])

        return cls(
            np.linalg.det(matr_a),
            np.linalg.det(matr_b),
            np.linalg.det(matr_c),
            np.linalg.det(matr_d)
        )

    @classmethod
    def from_geological_plane(cls,
                              geol_plane: Plane,
                              pt: Point3D):
        """
          Given a Plane instance and a provided Point instance,
          calculate the corresponding Plane instance.

          Example:
            >>> CPlane3D.from_geological_plane(Plane(0, 0), Point3D(0, 0, 0))
            CPlane3D(0.0000, 0.0000, 1.0000, -0.0000)
            >>> CPlane3D.from_geological_plane(Plane(90, 45), Point3D(0, 0, 0))
            CPlane3D(0.7071, 0.0000, 0.7071, -0.0000)
            >>> CPlane3D.from_geological_plane(Plane(0, 90), Point3D(0, 0, 0))
            CPlane3D(0.0000, 1.0000, -0.0000, -0.0000)
          """

        normal_versor = geol_plane.normDirectFrwrd().as_versor()
        a, b, c = normal_versor.x, normal_versor.y, normal_versor.z
        d = - (a * pt.x + b * pt.y + c * pt.z)
        return CPlane3D(a, b, c, d)

    def __repr__(self):

        return "CPlane3D({:.4f}, {:.4f}, {:.4f}, {:.4f})".format(*self.v())

    def normVersor(self) -> Optional[Vect3D]:
        """
        Return the versor normal to the cartesian plane.

        Examples:
          >>> CPlane3D(0, 0, 5, -2).normVersor()
          Vect3D(0.0000, 0.0000, 1.0000)
          >>> CPlane3D(0, 7, 0, 5).normVersor()
          Vect3D(0.0000, 1.0000, 0.0000)
        """

        return Vect3D(self.a(), self.b(), self.c()).versor()

    def toPoint(self) -> Point3D:
        """
        Returns a point lying in the plane (non-unique solution).

        Examples:
          >>> CPlane3D(0, 0, 1, -1).toPoint()
          Point3D(0.0000, 0.0000, 1.0000)
        """

        point = Point3D(
            *pointSolution(
                np.array([[self.a(), self.b(), self.c()]]),
                np.array([-self.d()]))
        )

        return point

    """
    def gplane_point(self):
        '''
        Converts a cartesian plane into a geological plane
        and a point lying in the plane (non-unique solution).

        Examples:
          >>> gpl, pt = CPlane3D(0, 0, 1, -1).gplane_point()
          >>> gpl
          GPlane(000.00, +00.00)
          >>> pt
          Point(0.0000, 0.0000, 1.0000, nan)
        '''

        geol_plane = self.normVersor().gvect.normal_gplane
        point = Point4D(*point_solution(np.array([[self.a, self.b, self.c]]),
                                        np.array([-self.d])))
        return geol_plane, point
    """

    def intersVersor(self, another) -> Optional[Vect3D]:
        """
        Return intersection versor for two intersecting planes.
        Return None for not intersecting planes.

        :param another: another Cartesian plane.
        :type another: CPlane3D.
        :return: the intersection line as a vector.
        :rtype: Optional[Vect].
        :raise: Exception.

        Examples:
          >>> a = CPlane3D(1, 0, 0, 0)
          >>> b = CPlane3D(0, 0, 1, 0)
          >>> a.intersVersor(b)
          Vect3D(0.0000, -1.0000, 0.0000)
          >>> b = CPlane3D(-1, 0, 0, 0)  # parallel plane, no intersection
          >>> a.intersVersor(b) is None
          True
        """

        check_type(another, "Input Cartesian plane", CPlane3D)

        return self.normVersor().cross_product(another.normVersor()).versor()

    def intersPoint(self,
            another) -> Optional[Point3D]:
        """
        Return point on intersection line (non-unique solution)
        for two planes.

        :param another: the second cartesian plane
        :type another: CPlane3D
        :return: the optional instersection point
        :rtype: Optional[Point]
        :raise: Exception

        Examples:
          >>> p_a = CPlane3D(1, 0, 0, 0)
          >>> p_b = CPlane3D(0, 0, 1, 0)
          >>> p_a.intersPoint(p_b)
          Point3D(0.0000, 0.0000, 0.0000)
          >>> p_b = CPlane3D(-1, 0, 0, 0)  # parallel plane, no intersection
          >>> p_a.intersPoint(p_b) is None
        """

        check_type(another, "Second plane", CPlane3D)

        # find a point lying on the intersection line (this is a non-unique solution)

        a = np.array([[self.a(), self.b(), self.c()], [another.a(), another.b(), another.c()]])
        b = np.array([-self.d(), -another.d()])
        x, y, z = pointSolution(a, b)

        if x is not None and y is not None and z is not None:
            return Point3D(x, y, z)
        else:
            return None

    def pointDistance(self,
        pt: Point3D
    ) -> numbers.Real:
        """
        Calculate the distance between a point and the cartesian plane.
        Distance expression:
        distance = a * x1 + b * y1 + c * z1 + d
        where a, b, c and d are plane parameters of the plane equation:
         a * x + b * y + c * z + d = 0
        and x1, y1, and z1 are the point coordinates.

        :param pt: the point to calculate distance with.
        :type pt: Point.
        :return: the distance value.
        :rtype: numbers.Real.
        :raise: Exception.

        Examples:
          >>> cpl = CPlane3D(0, 0, 1, 0)
          >>> pt = Point3D(0, 0, 1)
          >>> cpl.pointDistance(pt)
          1.0
          >>> pt = Point3D(0, 0, 0.5)
          >>> cpl.pointDistance(pt)
          0.5
          >>> pt = Point3D(0, 0, -0.5)
          >>> cpl.pointDistance(pt)
          -0.5
          >>> pt = Point3D(10, 20, 0.0)
          >>> cpl.pointDistance(pt)
          0.0
        """

        check_type(pt, "Input point", Point3D)

        return self.a() * pt.x + self.b() * pt.y + self.c() * pt.z + self.d()

    def isPointInPlane(self,
        pt: Union[Point3D, Point2D]
    ) -> bool:
        """
        Check whether a point lies in the current plane.

        :param pt: the point to check.
        :return: whether the point lies in the current plane.
        :raise: Exception.

        Examples:
          >>> pl = CPlane3D(0, 0, 1, 0)
          >>> pt = Point3D(0, 1, 0)
          >>> pl.isPointInPlane(pt)
          True
          >>> pl = CPlane3D(0, 0, 1, 0)
          >>> pt = Point3D(0, 1, 0)
          >>> pl.isPointInPlane(pt)
          True
        """

        check_type(pt, "Input point", (Point2D, Point3D))

        if isinstance(pt, Point2D):
            pt = Point3D(
                pt.x,
                pt.y,
                0.0
            )
        if abs(self.pointDistance(pt)) < MIN_SEPARATION_THRESHOLD:
            return True
        else:
            return False

    def angle_as_degrees(self,
         another: 'CPlane3D'
         ) -> numbers.Real:
        """
        Calculate angle (in degrees) between two planes.

        :param another: the CPlane3D instance to calculate angle with.
        :type another: CPlane3D.
        :return: the angle (in degrees) between the two planes.
        :rtype: numbers.Real.
        :raise: Exception.

        Examples:
          >>> CPlane3D(1,0,0,0).angle_as_degrees(CPlane3D(0,1,0,0))
          90.0
          >>> CPlane3D(1,0,0,0).angle_as_degrees(CPlane3D(0,1,0,0))
          90.0
          >>> CPlane3D(1,0,0,0).angle_as_degrees(CPlane3D(1,0,1,0))
          45.0
          >>> CPlane3D(1,0,0,0).angle_as_degrees(CPlane3D(1,0,0,0))
          0.0
        """

        check_type(another, "Second Cartesian plane", CPlane3D)

        angle_degr = self.normVersor().angle_as_degrees(another.normVersor())

        if angle_degr > 90.0:
            angle_degr = 180.0 - angle_degr

        return angle_degr


class ParamLine3D(object):
    """
    parametric line
    srcPt: source Point
    l, m, n: line coefficients
    """

    def __init__(self, srcPt, l, m, n):

        for v in (l, m, n):
            if not (-1.0 <= v <= 1.0):
                raise Exception("Parametric line values must be in -1 to 1 range")

        self._srcPt = srcPt.clone()
        self._l = l
        self._m = m
        self._n = n

    '''
    def epsg(self) -> numbers.Integral:
        """
        Return the EPSG code of the parametric line.
        """

        return self._srcPt.epsg_code
    '''

    def intersect_cartes_plane(self, cartes_plane) -> Optional[Point3D]:
        """
        Return intersection point between parametric line and Cartesian plane.

        :param cartes_plane: a Cartesian plane:
        :type cartes_plane: CPlane3D.
        :return: the intersection point between parametric line and Cartesian plane.
        :rtype: Point.
        :raise: Exception.
        """

        if not isinstance(cartes_plane, CPlane3D):
            raise Exception("Method argument should be a Cartesian plane but is {}".format(type(cartes_plane)))

        '''
        if cartes_plane.epsg_code != self.epsg_code:
            raise Exception("Parametric line has EPSG {} while Cartesian plane has {}".format(self.epsg_code, cartes_plane.epsg_code))
        '''

        # line parameters
        x1, y1, z1 = self._srcPt.x, self._srcPt.y, self._srcPt.z
        l, m, n = self._l, self._m, self._n
        # Cartesian plane parameters
        a, b, c, d = cartes_plane.a(), cartes_plane.b(), cartes_plane.c(), cartes_plane.d()
        try:
            k = (a * x1 + b * y1 + c * z1 + d) / (a * l + b * m + c * n)
        except ZeroDivisionError:
            return None

        return Point3D(
            x=x1 - l * k,
            y=y1 - m * k,
            z=z1 - n * k
        )


def closure_plane_from_geo(
        plane: Plane,
        src_pt: Point3D
) -> Callable:
    """
    Closure that embodies the analytical formula for a given, non-vertical plane.
    This closure is used to calculate the z value from given horizontal coordinates (x, y).

    :param plane: the geological plane
    :param src_pt: the 3D point expressing a location point contained by the plane.


    :return: lambda (closure) expressing an analytical formula for deriving z given x and y values.
    """

    x0 = src_pt.x
    y0 = src_pt.y
    z0 = src_pt.z

    # slope of the line parallel to the x axis and contained by the plane
    a = plane.slope_x_dir()

    # slope of the line parallel to the y axis and contained by the plane
    b = plane.slope_y_dir()

    return lambda x, y: a * (x - x0) + b * (y - y0) + z0


class Points3D:
    """
    Collection of points.
    """

    def __init__(self,
                 x_array: array,
                 y_array: array,
                 z_array: array
                 ):
        """
        Construct a point list from a set of array values.

        :param x_array: the array storing the x values
        :param y_array: the array storing the y values
        :param z_array: the optional array storing the z values
        """

        check_type(
            var=x_array,
            name="X array",
            expected_types=array
        )

        check_type(
            var=y_array,
            name="Y array",
            expected_types=array
        )

        array_length = len(x_array)

        if len(y_array) != array_length:
            raise Exception(f"Y array has length {len(y_array)} while X array has length {len(x_array)}")

        check_type(
            var=z_array,
            name="Z array",
            expected_types=array
        )

        if len(z_array) != array_length:
            raise Exception(f"Z array has length {len(z_array)} while X array has length {len(x_array)}")

        self._x_array = x_array
        self._y_array = y_array
        self._z_array = z_array

    def num_pts(self
                ) -> int:
        """
        Numbers of points.
        """

        return len(self._x_array)

    @classmethod
    def fromPoints(cls,
                   points: List[Point3D]
                   ):
        """

        :param points: list of points
        """

        for ndx, point in enumerate(points):

            check_type(point, "Input point {}".format(ndx), Point3D)

        return Points3D(
            x_array=array('d', [p.x for p in points]),
            y_array=array('d', [p.y for p in points]),
            z_array=array('d', [p.z for p in points])
        )

    @property
    def xs(self
           ) -> array:
        """
        Returns a copy of the points x values.

        :return: points x values
        """

        return copy(self._x_array)

    @property
    def ys(self
           ) -> array:
        """
        Returns a copy of the points y values.

        :return: points y values
        """

        return copy(self._y_array)

    @property
    def zs(self
           ) -> array:
        """
        Returns a copy of the points z values.

        :return: points z values
        """

        return copy(self._z_array)

    def pt(self, pt_ndx: numbers.Integral) -> Point3D:
        """
        Extract the point at index pt_ndx.

        :param pt_ndx: point index.
        :type pt_ndx: numbers.Integral.
        :return: the extracted Point instance.
        :rtype: Point.

        Examples:
        """

        return Point3D(
            x=self._x_array[pt_ndx],
            y=self._y_array[pt_ndx],
            z=self._z_array[pt_ndx]
        )

    def values_at(self,
        ndx: numbers.Integral
                  ) -> Tuple[float, float, float]:
        """
        Return the values at given index.

        :param ndx: the index of the point values to extract
        :type ndx: numbers.Integral
        :return: the x, y and z values
        """

        return (
            self._x_array[ndx],
            self._y_array[ndx],
            self._z_array[ndx]
        )

    def pts(self):

        return [Point3D(*self.values_at(ndx)) for ndx in range(self.num_pts())]

    def __repr__(self) -> str:
        """
        Represents a Points instance as a shortened text.

        :return: a textual shortened representation of a Points instance.
        """

        num_points = self.num_pts()

        if num_points == 0:
            txt = "Empty Points3D"
        else:
            x1, y1, z1 = self.values_at(0)
            if num_points == 1:
                txt = "Points3D with unique point: {.4f}.{.4f},{.4f}".format(x1, y1, z1)
            else:
                x2, y2, z2 = self.values_at(self.num_pts()-1)
                txt = "Points3D with {} points: ({:.4f}, {:.4f}, {:.4f}) ... ({:.4f}, {:.4f}, {:.4f})".format(
                    num_points, x1, y1, z1, x2, y2, z2)

        return txt

    def __iter__(self):
        """
        Return each point.
        """

        return (self.pt(ndx) for ndx in range(self.num_pts()))

    def asXyzArray(self):
        """
        Convert to a Numpy x-y-z array
        """

        return np.vstack(
            (
                self.xs,
                self.ys,
                self.zs
            )
        ).transpose()

    def add_pt(self, pt) -> None:
        """
        In-place transformation of the original Points3D instance
        by adding a new point at the end.

        :param pt: the point to add
        :return: nothing
        """

        self._x_array.append(pt.x)
        self._y_array.append(pt.y)
        self._z_array.append(pt.z)

    def add_pts(self,
                pts: 'Points3D'
                ):
        """
        In-place transformation of the original Points instance
        by adding a new set of points at the end.

        :param pts: list of Points.
        """

        check_type(pts, "Points", Points3D)

        self._x_array.extend(pts.xs)
        self._y_array.extend(pts.ys)
        self._z_array.extend(pts.zs)

    def x_min(self) -> Optional[numbers.Real]:
        """
        Optional minimum of x values.

        :return: the optional minimum of x values.
        :rtype: Optional[numbers.Real]
        """

        return np.nanmin(self._x_array) if self.num_pts() > 0 else None

    def x_max(self) -> Optional[numbers.Real]:
        """
        Optional maximum x value.
        """

        return np.nanmax(self._x_array) if self.num_pts() > 0 else None

    def x_mean(self) -> Optional[numbers.Real]:
        """
        Optional mean x value.
        """

        return np.nanmean(self._x_array) if self.num_pts() > 0 else None

    def y_min(self) -> Optional[numbers.Real]:
        """
        Optional minimum y value.
        """

        return np.nanmin(self._y_array) if self.num_pts() > 0 else None

    def y_max(self) -> Optional[numbers.Real]:
        """
        Optional maximum y value.
        """

        return np.nanmax(self._y_array) if self.num_pts() > 0 else None

    def y_mean(self) -> Optional[numbers.Real]:
        """
        Optional mean y value.
        """

        return np.nanmean(self._y_array) if self.num_pts() > 0 else None

    def z_min(self) -> Optional[numbers.Real]:
        """
        Optional minimum z value.
        """

        return np.nanmin(self._z_array) if self.num_pts() > 0 else None

    def z_max(self) -> Optional[numbers.Real]:
        """
        Optional maximum z value.
        """

        return np.nanmax(self._z_array) if self.num_pts() > 0 else None

    def z_mean(self) -> Optional[numbers.Real]:
        """
        Optional mean z value.
        """

        return np.nanmean(self._z_array) if self.num_pts() > 0 else None

    def z_var(self) -> Optional[numbers.Real]:
        """
        Optional variance of z values.

        :return: the optional variance of z values.
        :rtype: Optional[numbers.Real]

        Examples:
          >>> l = Points3D.fromPoints([Point3D(0, 0, 2), Point3D(1, 0, 2), Point3D(0, 1, 2)])
          >>> l.z_var()
          0.0
        """

        return np.nanvar(self._z_array) if self.num_pts() > 0 else None

    def z_std(self) -> Optional[numbers.Real]:
        """
        Optional standard deviation of z values.

        :return: the optional standard deviation of z values.
        :rtype: Optional[numbers.Real]

        Examples:
          >>> l = Points3D.fromPoints([Point3D(0, 0, 2), Point3D(1, 0, 2), Point3D(0, 1, 2)])
          >>> l.z_std()
          0.0
        """

        return np.nanstd(self._z_array) if self.num_pts() > 0 else None

    def nanmean_point(self) -> Point3D:
        """
        Returns the nan- excluded mean point of the collection.
        It is the mean point for a collection of point in a x-y-z frame (i.e., not lat-lon).

        :return: the nan- excluded mean point of the collection.
        """

        return Point3D(
            x=np.nanmean(self._x_array),
            y=np.nanmean(self._y_array),
            z=np.nanmean(self._z_array)
        )

    def segment(self,
        ndx: int
    ) -> Optional[Segment3D]:
        """
        Returns the optional segment starting at index ndx.

        :param ndx: the segment index.
        :return: the optional segment
        """

        if ndx < 0 or ndx >= self.num_pts() - 1:
            return None

        return Segment3D(
            start_pt=self.pt(ndx),
            end_pt=self.pt(ndx + 1)
        )

    def reversed(self) -> 'Points3D':
        """
        Return a Points3D instance with reversed point list.

        :return: a new Points3D instance.
        """

        xs = self._x_array.reversed()
        ys = self._y_array.reversed()
        zs = self._z_array.reversed()

        return Points3D(
            x_array=xs,
            y_array=ys,
            z_array=zs
        )

