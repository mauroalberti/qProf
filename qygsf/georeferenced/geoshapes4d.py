
import datetime

import numbers

from .crs import *
from ..geometries.shapes.space3d import *


class GeoPoint4D:

    def __init__(self,
                 x: Union[numbers.Integral, numbers.Real],
                 y: Union[numbers.Integral, numbers.Real],
                 z: Optional[Union[numbers.Integral, numbers.Real]] = 0.0,
                 t: Optional[datetime.datetime] = None,
                 epsg_cd: Optional[numbers.Integral] = -1
                 ):
        """
        Construct a GeoPoint4D instance.

        :param x: point x coordinate
        :param y: point y coordinate
        :param z: point z coordinate
        :param t: time
        :param epsg_cd: EPSG code
        """

        check_type(
            var=epsg_cd,
            name="EPSG code",
            expected_types=numbers.Integral
        )

        vals = [x, y, z]
        if any(map(lambda val: not isinstance(val, numbers.Real), vals)):
            raise Exception("Input values must be integer or float type")
        if not all(map(math.isfinite, vals)):
            raise Exception("Input values must be finite (#04)")
        if t is not None:
            check_type(
                var=t,
                name="Time",
                expected_types=datetime.datetime
            )

        self._x = float(x)
        self._y = float(y)
        self._z = float(z)
        self._t = t
        self._epsg_code = epsg_cd

    @classmethod
    def fromVect(cls,
                 vect: Vect3D,
                 epsg_cd: Optional[numbers.Integral] = -1
                 ) -> 'GeoPoint4D':

        return cls(
            x=vect.x,
            y=vect.y,
            z=vect.z,
            t=None,
            epsg_cd=epsg_cd
        )

    @property
    def x(self) -> numbers.Real:
        """
        Return the x coordinate of the current point.

        :return: x coordinate.
        :rtype: numbers.Real

        Examples:
          >>> GeoPoint4D(4, 3, 7).x
          4.0
          >>> GeoPoint4D(-0.39, 3, 7).x
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
          >>> GeoPoint4D(4, 3, 7).y
          3.0
          >>> GeoPoint4D(-0.39, 17.42, 7).y
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
          >>> GeoPoint4D(4, 3, 7).z
          7.0
          >>> GeoPoint4D(-0.39, 17.42, 8.9).z
          8.9
        """

        return self._z

    @property
    def t(self) -> datetime.datetime:
        """
        Return the time associated to the current point.

        :return: z coordinate.
        """

        return self._t

    @property
    def crs(self) -> Crs:
        """
        The points CRS.

        :return: the points CRS
        :rtype: Crs
        """

        return Crs(self.epsg_code)

    @property
    def epsg_code(self) -> numbers.Integral:
        """
        Return the EPSG code of the current point.

        :return: EPSG code.
        """

        return self._epsg_code

    def __iter__(self):
        """
        Return the elements of a Point.

        :return:

        Examples;
          >>> x, y, z, t, epsg_code = GeoPoint4D(1, 1, 2, None, 32633)
          >>> x == 1.0
          True
          >>> y == 1.0
          True
          >>> z == 2.0
          True
          >>> t is None
          True
          >>> epsg_code == 32633
          True
        """

        return (i for i in self.a())

    def __repr__(self) -> str:

        return f"GeoPoint4D(x={self.x:.4f}, y={self.y:.4f}, z={self.z:.4f}, time={self.t}, epsg_code={self.epsg_code})"

    def __eq__(self,
        another: 'GeoPoint4D'
    ) -> bool:
        """
        Return True if objects are equal.

        :param another: another point.
        :type another: Point.
        :raise: Exception.

        Example:
          >>> GeoPoint4D(1., 1., 1.) == GeoPoint4D(1, 1, 1)
          True
          >>> GeoPoint4D(1., 1., 1.) == GeoPoint4D(1, 1, 1)
          True
          >>> GeoPoint4D(1., 1., 1.) == GeoPoint4D(1, 1, -1)
          False
        """

        if not isinstance(another, GeoPoint4D):
            raise Exception("Another instance must be a Point")

        return all(
            [
                self.x == another.x,
                self.y == another.y,
                self.z == another.z,
                self.t == another.t,
                self.epsg_code == another.epsg_code
            ]
        )

    def __ne__(self,
        another: 'GeoPoint4D'
    ) -> bool:
        """
        Return False if objects are equal.

        Example:
          >>> GeoPoint4D(1., 1., 1.) != GeoPoint4D(0., 0., 0.)
          True
          >>> GeoPoint4D(1., 1., 1.) != GeoPoint4D(1, 1, 1)
          False
        """

        return not (self == another)

    def a(self
          ) -> Tuple[numbers.Real, numbers.Real, numbers.Real, Optional[datetime.datetime], numbers.Integral]:
        """
        Return the individual values of the point.

        :return: double array of x, y, z values

        Examples:
          >>> GeoPoint4D(4, 3, 7).a()
          (4.0, 3.0, 7.0, None, -1)
        """

        return self.x, self.y, self.z, self.t, self.epsg_code

    def clone(self) -> 'GeoPoint4D':
        """
        Clone a point.

        :return: a new point.
        :rtype: Point.
        """

        return GeoPoint4D(*self.a())

    def toXYZ(self) -> Tuple[numbers.Real, numbers.Real, numbers.Real]:
        """
        Returns the spatial components as a tuple of three values.

        :return: the spatial components (x, y, z).
        :rtype: a tuple of three floats.

        Examples:
          >>> GeoPoint4D(1, 0, 3).toXYZ()
          (1.0, 0.0, 3.0)
        """

        return self.x, self.y, self.z

    def toArray(self) -> np.ndarray:
        """
        Return a Numpy array representing the point values.

        :return: Numpy array

        Examples:
          >>> np.allclose(GeoPoint4D(1, 2, 3).toArray(), np.array([ 1., 2., 3.]))
          True
        """

        return np.asarray(self.toXYZ())

    def pXY(self) -> 'GeoPoint4D':
        """
        Projection on the x-y plane

        :return: projected object instance

        Examples:
          >>> GeoPoint4D(2, 3, 4).pXY()
          GeoPoint4D(x=2.0000, y=3.0000, z=0.0000, time=None, epsg_code=-1)
        """

        return GeoPoint4D(self.x, self.y, 0.0)

    def pXZ(self) -> 'GeoPoint4D':
        """
        Projection on the x-z plane

        :return: projected object instance

        Examples:
          >>> GeoPoint4D(2, 3, 4).pXZ()
          GeoPoint4D(x=2.0000, y=0.0000, z=4.0000, time=None, epsg_code=-1)
        """

        return GeoPoint4D(self.x, 0.0, self.z)

    def pYZ(self) -> 'GeoPoint4D':
        """
        Projection on the y-z plane

        :return: projected object instance

        Examples:
          >>> GeoPoint4D(2, 3, 4).pYZ()
          GeoPoint4D(x=0.0000, y=3.0000, z=4.0000, time=None, epsg_code=-1)
        """

        return GeoPoint4D(0.0, self.y, self.z)

    def deltaX(self,
        another: 'GeoPoint4D'
    ) -> Optional[numbers.Real]:
        """
        Delta between x components of two Point Instances.

        :return: x coordinates difference value.
        :rtype: optional numbers.Real.
        :raise: Exception

        Examples:
          >>> GeoPoint4D(1, 2, 3).deltaX(GeoPoint4D(4, 7, 1))
          3.0
        """

        check_crs(self, another)

        return another.x - self.x

    def deltaY(self,
        another: 'GeoPoint4D'
    ) -> Optional[numbers.Real]:
        """
        Delta between y components of two Point Instances.

        :return: y coordinates difference value.
        :rtype: optional numbers.Real.

        Examples:
          >>> GeoPoint4D(1, 2, 3).deltaY(GeoPoint4D(4, 7, 1))
          5.0
        """

        check_crs(self, another)

        return another.y - self.y

    def deltaZ(self,
        another: 'GeoPoint4D'
    ) -> Optional[numbers.Real]:
        """
        Delta between z components of two Point Instances.

        :return: z coordinates difference value.
        :rtype: optional numbers.Real.

        Examples:
          >>> GeoPoint4D(1, 2, 3).deltaZ(GeoPoint4D(4, 7, 1))
          -2.0
        """

        check_crs(self, another)

        return another.z - self.z

    def deltaT(self,
        another: 'GeoPoint4D'
    ) -> Optional[datetime.timedelta]:
        """
        Delta between t components of two Point Instances.

        :return: t coordinates difference value.
        """

        #check_crs(self, another)

        return another.t - self.t

    def distance(self,
                 another: 'GeoPoint4D'
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
          >>> GeoPoint4D(1., 1., 1.).distance(GeoPoint4D(4., 5., 1))
          5.0
        """

        check_type(another, "GeoPoint", GeoPoint4D)

        return sqrt((self.x - another.x) ** 2 + (self.y - another.y) ** 2 + (self.z - another.z) ** 2)

    def horizontal_distance(self,
                            another: 'GeoPoint4D'
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
          >>> GeoPoint4D(1., 1., 1.).horizontal_distance(GeoPoint4D(4., 5., 7.))
          5.0
        """

        check_type(another, "Second point", GeoPoint4D)

        return sqrt((self.x - another.x) ** 2 + (self.y - another.y) ** 2)

    def scale(self,
        scale_factor: numbers.Real
    ) -> 'GeoPoint4D':
        """
        Create a scaled object.
        Note: it does not make sense for polar coordinates.
        TODO: manage polar coordinates cases OR deprecate and remove - after dependency check.

        Example;
          >>> GeoPoint4D(1, 0, 1).scale(2.5)
          GeoPoint4D(x=2.5000, y=0.0000, z=2.5000, time=None, epsg_code=-1)
        """

        x, y, z = self.x * scale_factor, self.y * scale_factor, self.z * scale_factor
        return GeoPoint4D(
            x=x,
            y=y,
            z=z,
            t=self.t,
            epsg_cd=self.epsg_code
        )

    def invert(self) -> 'GeoPoint4D':
        """
        Create a new object with inverted direction.
        Note: it depends on scale method, that could be deprecated/removed.

        Examples:
          >>> GeoPoint4D(1, 1, 1).invert()
          GeoPoint4D(x=-1.0000, y=-1.0000, z=-1.0000, time=None, epsg_code=-1)
          >>> GeoPoint4D(2, -1, 4).invert()
          GeoPoint4D(x=-2.0000, y=1.0000, z=-4.0000, time=None, epsg_code=-1)
        """

        return self.scale(-1)

    def reflect_vertical(self) -> 'GeoPoint4D':
        """
        Reflect a point along a vertical axis.

        :return: reflected point.
        :rtype: GeoPoint4D

        Examples:
          >>> GeoPoint4D(1, 1, 1).reflect_vertical()
          GeoPoint4D(x=-1.0000, y=-1.0000, z=1.0000, time=None, epsg_code=-1)
        """

        x, y, z, t, epsg_cd = self

        return GeoPoint4D(
            x=-x,
            y=-y,
            z=z,
            t=t,
            epsg_cd=epsg_cd
        )

    def is_coincident(self,
                      another: 'GeoPoint4D',
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
          >>> GeoPoint4D(1., 0., -1.).is_coincident(GeoPoint4D(1., 1.5, -1.))
          False
          >>> GeoPoint4D(1., 0., 0.).is_coincident(GeoPoint4D(1., 0., 0.))
          True
        """

        check_type(another, "Second point", GeoPoint4D)

        return self.distance(another) <= tolerance

    def already_present(self,
                        pt_list: List['GeoPoint4D'],
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
        sz: numbers.Real = 0.0,
        st: Optional[datetime.timedelta] = None
    ) -> Optional['GeoPoint4D']:
        """
        Create a new object shifted by given amount from the self instance.

        Example:
          >>> GeoPoint4D(1, 1, 1).shift(0.5, 1., 1.5)
          GeoPoint4D(x=1.5000, y=2.0000, z=2.5000, time=None, epsg_code=-1)
          >>> GeoPoint4D(1, 2, -1).shift(0.5, 1., 1.5)
          GeoPoint4D(x=1.5000, y=3.0000, z=0.5000, time=None, epsg_code=-1)
       """

        if self.t is None:
            time = st
        elif st is None:
            time = self.t
        else:
            time = self.t + st
        return GeoPoint4D(
            x=self.x + sx,
            y=self.y + sy,
            z=self.z + sz,
            t=time,
            epsg_cd=self.epsg_code
        )

    def shiftByVect(self,
        v: Vect3D
    ) -> 'GeoPoint4D':
        """
        Create a new point shifted from the self instance by given vector.

        :param v: the shift vector.
        :type v: Vect.
        :return: the shifted point.
        :rtype: Point.
        :raise: Exception

        Example:
          >>> GeoPoint4D(1, 1, 1).shiftByVect(Vect3D(0.5, 1., 1.5))
          GeoPoint4D(x=1.5000, y=2.0000, z=2.5000, time=None, epsg_code=-1)
          >>> GeoPoint4D(1, 2, -1).shiftByVect(Vect3D(0.5, 1., 1.5))
          GeoPoint4D(x=1.5000, y=3.0000, z=0.5000, time=None, epsg_code=-1)
       """

        x, y, z, t, epsg_cd = self

        sx, sy, sz = v.toXYZ()

        return GeoPoint4D(
            x=x + sx,
            y=y + sy,
            z=z + sz,
            t=t,
            epsg_cd=epsg_cd
        )

    def asVect(self) -> 'Vect3D':
        """
        Create a vector based on the point coordinates

        Example:
          >>> GeoPoint4D(1, 1, 0).asVect()
          Vect3D(1.0000, 1.0000, 0.0000)
          >>> GeoPoint4D(0.2, 1, 6).asVect()
          Vect3D(0.2000, 1.0000, 6.0000)
        """

        return Vect3D(
            x=self.x,
            y=self.y,
            z=self.z
        )

    def rotate(self,
        rotation_axis: RotationAxis,
        center_point: 'GeoPoint4D' = None
        ) -> 'GeoPoint4D':
        """
        Rotates a point.
        :param rotation_axis:
        :param center_point:
        :return: the rotated point
        :rtype: GeoPoint4D

        Examples:
          >>> pt = GeoPoint4D(0,0,1)
          >>> rot_axis = RotationAxis(0,0,90)
          >>> center_pt = GeoPoint4D(0,0,0.5)
          >>> pt.rotate(rotation_axis=rot_axis, center_point=center_pt)
          GeoPoint4D(x=0.5000, y=0.0000, z=0.5000, time=None, epsg_code=-1)
          >>> center_pt = GeoPoint4D(0,0,1)
          >>> pt.rotate(rotation_axis=rot_axis, center_point=center_pt)
          GeoPoint4D(x=0.0000, y=0.0000, z=1.0000, time=None, epsg_code=-1)
          >>> center_pt = GeoPoint4D(0, 0, 2)
          >>> pt.rotate(rotation_axis=rot_axis, center_point=center_pt)
          GeoPoint4D(x=-1.0000, y=0.0000, z=2.0000, time=None, epsg_code=-1)
          >>> rot_axis = RotationAxis(0,0,180)
          >>> pt.rotate(rotation_axis=rot_axis, center_point=center_pt)
          GeoPoint4D(x=-0.0000, y=0.0000, z=3.0000, time=None, epsg_code=-1)
          >>> pt.rotate(rotation_axis=rot_axis)
          GeoPoint4D(x=0.0000, y=0.0000, z=-1.0000, time=None, epsg_code=-1)
          >>> pt = GeoPoint4D(1, 1, 1)
          >>> rot_axis = RotationAxis(0,90,90)
          >>> pt.rotate(rotation_axis=rot_axis)
          GeoPoint4D(x=1.0000, y=-1.0000, z=1.0000, time=None, epsg_code=-1)
          >>> rot_axis = RotationAxis(0,90,180)
          >>> pt.rotate(rotation_axis=rot_axis)
          GeoPoint4D(x=-1.0000, y=-1.0000, z=1.0000, time=None, epsg_code=-1)
          >>> center_pt = GeoPoint4D(1,1,1)
          >>> pt.rotate(rotation_axis=rot_axis, center_point=center_pt)
          GeoPoint4D(x=1.0000, y=1.0000, z=1.0000, time=None, epsg_code=-1)
          >>> center_pt = GeoPoint4D(2,2,10)
          >>> pt.rotate(rotation_axis=rot_axis, center_point=center_pt)
          GeoPoint4D(x=3.0000, y=3.0000, z=1.0000, time=None, epsg_code=-1)
          >>> pt = GeoPoint4D(1, 1, 2)
          >>> rot_axis = RotationAxis(135, 0, 180)
          >>> center_pt = GeoPoint4D(0,0,1)
          >>> pt.rotate(rotation_axis=rot_axis, center_point=center_pt)
          GeoPoint4D(x=-1.0000, y=-1.0000, z=0.0000, time=None, epsg_code=-1)
        """

        if not center_point:

            center_point = GeoPoint4D(
                x=0.0,
                y=0.0,
                z=0.0,
                t=self.t,
                epsg_cd=self.epsg_code
            )

        check_type(center_point, "Center point", GeoPoint4D)

        if self.distance(center_point) == 0:
            return self.clone()

        p_vect = GeoSegment4D(start_pt=center_point, end_pt=self).as_vector()

        rot_vect = rotVectByAxis(
            v=p_vect,
            rot_axis=rotation_axis
        )

        x, y, z = rot_vect

        transl_pt = center_point.shift(
            sx=x,
            sy=y,
            sz=z
        )

        return transl_pt

    @classmethod
    def random(cls,
        lower_boundary: float = -MAX_SCALAR_VALUE,
        upper_boundary: float = MAX_SCALAR_VALUE
    ):
        """
        Creates a random point.

        :return: random point
        :rtype: GeoPoint4D
        """

        vals = [random.uniform(lower_boundary, upper_boundary) for _ in range(3)]
        return cls(*vals)


class GeoSegment4D:
    """
    Segment is a geometric object defined by the straight line between
    two vertices.
    """

    def __init__(self,
                 start_pt: GeoPoint4D,
                 end_pt: GeoPoint4D):
        """
        Creates a segment instance provided the two points have the same CRS code.

        :param start_pt: the start point.
        :type: Point.
        :param end_pt: the end point.
        :type end_pt: Point.
        :return: the new segment instance if both points have the same georeferenced.
        :raises: CRSCodeException.
        """

        check_type(
            var=start_pt,
            name="Start point",
            expected_types=GeoPoint4D
        )

        check_type(
            var=end_pt,
            name="End point",
            expected_types=GeoPoint4D
        )

        check_crs(
            template_element=start_pt,
            checked_element=end_pt
        )

        if start_pt.distance(end_pt) == 0.0:
            raise Exception("Source points cannot be coincident")

        self._start_pt = start_pt.clone()
        self._end_pt = end_pt.clone()

    @classmethod
    def fromVector(cls,
                   point: GeoPoint4D,
                   dir_vector: Vect3D):

        check_type(point, "Input point", GeoPoint4D)
        check_type(dir_vector, "Directional vector", Vect3D)

        start_pt = point
        end_pt = start_pt.shiftByVect(dir_vector)

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

        return f"GeoSegment4D(start_pt={self.start_pt}, end_pt={self.end_pt})"

    @property
    def start_pt(self) -> GeoPoint4D:

        return self._start_pt

    @property
    def end_pt(self) -> GeoPoint4D:

        return self._end_pt

    def __iter__(self):
        """
        Return the elements of a Segment, i.e., start and end point.
        """

        return (i for i in [self.start_pt, self.end_pt])

    def clone(self) -> 'GeoSegment4D':

        return GeoSegment4D(
            start_pt=self._start_pt,
            end_pt=self._end_pt
        )

    @property
    def crs(self) -> Crs:
        """
        The points CRS.

        :return: the points CRS
        :rtype: Crs
        """

        return Crs(self.epsg_code)

    @property
    def epsg_code(self) -> numbers.Integral:
        """
        Return the EPSG code of the current point.

        :return: EPSG code.
        """

        return self.start_pt.epsg_code

    def increasing_x(self) -> 'GeoSegment4D':

        if self.end_pt.x < self.start_pt.x:
            return GeoSegment4D(
                start_pt=self.end_pt,
                end_pt=self.start_pt
            )
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

        return self.start_pt.horizontal_distance(self.end_pt)

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
                    pt: GeoPoint4D
                    ) -> bool:
        """
        Checks whether a point is contained in a segment.

        :param pt: the point for which to check containement.
        :return: bool.
        :raise: Exception.

        Examples:
          >>> segment = GeoSegment4D(GeoPoint4D(0, 0, 0), GeoPoint4D(1, 0, 0))
          >>> segment.contains_pt(GeoPoint4D(0, 0, 0))
          True
          >>> segment.contains_pt(GeoPoint4D(1, 0, 0))
          True
          >>> segment.contains_pt(GeoPoint4D(0.5, 0, 0))
          True
          >>> segment.contains_pt(GeoPoint4D(0.5, 0.00001, 0))
          False
          >>> segment.contains_pt(GeoPoint4D(0.5, 0, 0.00001))
          False
          >>> segment.contains_pt(GeoPoint4D(1.00001, 0, 0))
          False
          >>> segment.contains_pt(GeoPoint4D(0.000001, 0, 0))
          True
          >>> segment.contains_pt(GeoPoint4D(-0.000001, 0, 0))
          False
          >>> segment.contains_pt(GeoPoint4D(0.5, 1000, 1000))
          False
          >>> segment = GeoSegment4D(GeoPoint4D(0, 0, 0), GeoPoint4D(0, 1, 0))
          >>> segment.contains_pt(GeoPoint4D(0, 0, 0))
          True
          >>> segment.contains_pt(GeoPoint4D(0, 0.5, 0))
          True
          >>> segment.contains_pt(GeoPoint4D(0, 1, 0))
          True
          >>> segment.contains_pt(GeoPoint4D(0, 1.5, 0))
          False
          >>> segment = GeoSegment4D(GeoPoint4D(0, 0, 0), GeoPoint4D(1, 1, 1))
          >>> segment.contains_pt(GeoPoint4D(0.5, 0.5, 0.5))
          True
          >>> segment.contains_pt(GeoPoint4D(1, 1, 1))
          True
          >>> segment = GeoSegment4D(GeoPoint4D(1,2,3), GeoPoint4D(9,8,2))
          >>> segment.contains_pt(segment.pointAt(0.745))
          True
          >>> segment.contains_pt(segment.pointAt(1.745))
          False
          >>> segment.contains_pt(segment.pointAt(-0.745))
          False
          >>> segment.contains_pt(segment.pointAt(0))
          True
        """

        check_type(pt, "Point", GeoPoint4D)

        segment_length = self.length()
        length_startpt_pt = self.start_pt.distance(pt)
        length_endpt_pt = self.end_pt.distance(pt)

        return areClose(
            a=segment_length,
            b=length_startpt_pt + length_endpt_pt
        )

    def pointAt(self,
                scale_factor: numbers.Real
                ) -> GeoPoint4D:
        """
        Returns a point aligned with the segment
        and lying at given scale factor, where 1 is segment length
        ans 0 is segment start.

        :param scale_factor: the scale factor, where 1 is the segment length.
        :return: Point at scale factor

        Examples:
          >>> s = GeoSegment4D(GeoPoint4D(0,0,0), GeoPoint4D(1,0,0))
          >>> s.pointAt(0)
          GeoPoint4D(x=0.0000, y=0.0000, z=0.0000, time=None, epsg_code=-1)
          >>> s.pointAt(0.5)
          GeoPoint4D(x=0.5000, y=0.0000, z=0.0000, time=None, epsg_code=-1)
          >>> s.pointAt(1)
          GeoPoint4D(x=1.0000, y=0.0000, z=0.0000, time=None, epsg_code=-1)
          >>> s.pointAt(-1)
          GeoPoint4D(x=-1.0000, y=0.0000, z=0.0000, time=None, epsg_code=-1)
          >>> s.pointAt(-2)
          GeoPoint4D(x=-2.0000, y=0.0000, z=0.0000, time=None, epsg_code=-1)
          >>> s.pointAt(2)
          GeoPoint4D(x=2.0000, y=0.0000, z=0.0000, time=None, epsg_code=-1)
          >>> s = GeoSegment4D(GeoPoint4D(0, 0, 0), GeoPoint4D(0, 0, 1))
          >>> s.pointAt(0)
          GeoPoint4D(x=0.0000, y=0.0000, z=0.0000, time=None, epsg_code=-1)
          >>> s.pointAt(0.5)
          GeoPoint4D(x=0.0000, y=0.0000, z=0.5000, time=None, epsg_code=-1)
          >>> s.pointAt(1)
          GeoPoint4D(x=0.0000, y=0.0000, z=1.0000, time=None, epsg_code=-1)
          >>> s.pointAt(-1)
          GeoPoint4D(x=0.0000, y=0.0000, z=-1.0000, time=None, epsg_code=-1)
          >>> s.pointAt(-2)
          GeoPoint4D(x=0.0000, y=0.0000, z=-2.0000, time=None, epsg_code=-1)
          >>> s.pointAt(2)
          GeoPoint4D(x=0.0000, y=0.0000, z=2.0000, time=None, epsg_code=-1)
          >>> s = GeoSegment4D(GeoPoint4D(0, 0, 0), GeoPoint4D(1, 1, 1))
          >>> s.pointAt(0.5)
          GeoPoint4D(x=0.5000, y=0.5000, z=0.5000, time=None, epsg_code=-1)
          >>> s = GeoSegment4D(GeoPoint4D(0, 0, 0), GeoPoint4D(4, 0, 0))
          >>> s.pointAt(7.5)
          GeoPoint4D(x=30.0000, y=0.0000, z=0.0000, time=None, epsg_code=-1)
        """

        dx = self.delta_x() * scale_factor
        dy = self.delta_y() * scale_factor
        dz = self.delta_z() * scale_factor

        return GeoPoint4D(
            x=self.start_pt.x + dx,
            y=self.start_pt.y + dy,
            z=self.start_pt.z + dz
        )

    def pointProjection(self,
                        point: GeoPoint4D
                        ) -> GeoPoint4D:
        """
        Return the point projection on the segment.

        Examples:
          >>> s = GeoSegment4D(start_pt=GeoPoint4D(0,0,0), end_pt=GeoPoint4D(1,0,0))
          >>> p = GeoPoint4D(0.5, 1, 4)
          >>> s.pointProjection(p)
          GeoPoint4D(x=0.5000, y=0.0000, z=0.0000, time=None, epsg_code=-1)
          >>> s = GeoSegment4D(start_pt=GeoPoint4D(0,0,0), end_pt=GeoPoint4D(4,0,0))
          >>> p = GeoPoint4D(7.5, 19.2, -14.72)
          >>> s.pointProjection(p)
          GeoPoint4D(x=7.5000, y=0.0000, z=0.0000, time=None, epsg_code=-1)
        """

        check_type(point, "Input point", GeoPoint4D)

        check_crs(self, point)

        other_segment = GeoSegment4D(
            self.start_pt,
            point
        )

        scale_factor = self.vector().scalar_projection(other_segment.vector()) / self.length()
        return self.pointAt(scale_factor)

    def pointDistance(self,
                      point: GeoPoint4D
                      ) -> numbers.Real:
        """
        Returns the point distance to the segment.

        :param point: the point to calculate the distance with
        :type point: GeoPoint4D
        :return: the distance of the point to the segment
        :rtype: numbers.Real

        Examples:
          >>> s = GeoSegment4D(GeoPoint4D(0, 0, 0), GeoPoint4D(0, 0, 4))
          >>> s.pointDistance(GeoPoint4D(-17.2, 0.0, -49))
          17.2
          >>> s.pointDistance(GeoPoint4D(-17.2, 1.22, -49))
          17.24321315764553
        """

        check_type(point, "Input point", GeoPoint4D)

        # check_crs(self, point)

        point_projection = self.pointProjection(point)

        return point.distance(point_projection)

    def point_s(self,
                point: GeoPoint4D
                ) -> Optional[numbers.Real]:
        """
        Calculates the optional distance of the point along the segment.
        A zero value is for a point coinciding with the start point.
        Returns None if the point is not contained in the segment.

        :param point: the point to calculate the optional distance in the segment.
        :type point: GeoPoint4D
        :return: the the optional distance of the point along the segment.
        """

        check_type(point, "Input point", GeoPoint4D)

        # check_crs(self, point)

        if not self.contains_pt(point):
            return None

        return self.start_pt.distance(point)

    def scale(self,
              scale_factor
              ) -> 'GeoSegment4D':
        """
        Scale a segment by the given scale_factor.
        Start point does not change.

        :param scale_factor: the scale factor, where 1 is the segment length.
        :type scale_factor: numbers.Real
        :return: Point at scale factor
        :rtype: GeoPoint4D
        """

        end_pt = self.pointAt(scale_factor)

        return GeoSegment4D(
            self.start_pt,
            end_pt)

    def vertical_plane(self) -> Optional['CPlane3D']:
        """
        Returns the vertical Cartesian plane containing the segment.

        :return: the vertical Cartesian plane containing the segment.
        :rtype: Optional[CPlane].
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
                   another: 'GeoSegment4D',
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
          >>> s1 = GeoSegment4D(GeoPoint4D(0,0,0), GeoPoint4D(1,0,0))
          >>> s2 = GeoSegment4D(GeoPoint4D(0,0,0), GeoPoint4D(0,1,0))
          >>> s1.same_start(s2)
          True
        """

        return self.start_pt.is_coincident(
            another=another.start_pt,
            tolerance=tol
        )

    def same_end(self,
                 another: 'GeoSegment4D',
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
          >>> s1 = GeoSegment4D(GeoPoint4D(0,0,0), GeoPoint4D(1,0,0))
          >>> s2 = GeoSegment4D(GeoPoint4D(2,0,0), GeoPoint4D(1,0,0))
          >>> s1.same_end(s2)
          True
        """

        return self.end_pt.is_coincident(
            another=another.end_pt,
            tolerance=tol)

    def conn_to_other(self,
                      another: 'GeoSegment4D',
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
          >>> s1 = GeoSegment4D(GeoPoint4D(0,0,0), GeoPoint4D(1,0,0))
          >>> s2 = GeoSegment4D(GeoPoint4D(1,0,0), GeoPoint4D(2,0,0))
          >>> s1.conn_to_other(s2)
          True
        """

        return self.end_pt.is_coincident(
            another=another.start_pt,
            tolerance=tol)

    def other_connected(self,
                        another: 'GeoSegment4D',
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
          >>> s1 = GeoSegment4D(GeoPoint4D(0,0,0), GeoPoint4D(1,0,0))
          >>> s2 = GeoSegment4D(GeoPoint4D(-1,0,0), GeoPoint4D(0,0,0))
          >>> s1.other_connected(s2)
          True
        """

        return another.end_pt.is_coincident(
            another=self.start_pt,
            tolerance=tol)

    def segment_start_in(self,
                         another: 'GeoSegment4D'
                         ) -> bool:
        """
        Check whether the second segment contains the first segment start point.

        :param another: a segment to check for.
        :type another: Segment.
        :return: whether the second segment contains the first segment start point.
        :rtype: bool.

        Examples:
          >>> s1 = GeoSegment4D(GeoPoint4D(0,0,0), GeoPoint4D(1,0,0))
          >>> s2 = GeoSegment4D(GeoPoint4D(-0.5,0,0), GeoPoint4D(0.5,0,0))
          >>> s1.segment_start_in(s2)
          True
          >>> s1 = GeoSegment4D(GeoPoint4D(0,0,0), GeoPoint4D(1,1,1))
          >>> s1.segment_start_in(s2)
          True
          >>> s1 = GeoSegment4D(GeoPoint4D(0,1,0), GeoPoint4D(1,1,1))
          >>> s1.segment_start_in(s2)
          False
          >>> s1 = GeoSegment4D(GeoPoint4D(-1,-1,-1), GeoPoint4D(1,1,1))
          >>> s1.segment_start_in(s2)
          False
        """

        return another.contains_pt(self.start_pt)

    def segment_end_in(self,
                       another: 'GeoSegment4D'
                       ) -> bool:
        """
        Check whether the second segment contains the first segment end point.

        :param another: a segment to check for.
        :type another: Segment.
        :return: whether the second segment contains the first segment end point.
        :rtype: bool.

        Examples:
          >>> s1 = GeoSegment4D(GeoPoint4D(0,0,0), GeoPoint4D(1,0,0))
          >>> s2 = GeoSegment4D(GeoPoint4D(-0.5,0,0), GeoPoint4D(0.5,0,0))
          >>> s1.segment_end_in(s2)
          False
          >>> s1 = GeoSegment4D(GeoPoint4D(0,0,0), GeoPoint4D(1,1,1))
          >>> s1.segment_end_in(s2)
          False
          >>> s1 = GeoSegment4D(GeoPoint4D(0,1,0), GeoPoint4D(1,1,1))
          >>> s2 = GeoSegment4D(GeoPoint4D(1,1,1), GeoPoint4D(0.5,0,0))
          >>> s1.segment_end_in(s2)
          True
          >>> s1 = GeoSegment4D(GeoPoint4D(-1,-1,3), GeoPoint4D(1,1,3))
          >>> s2 = GeoSegment4D(GeoPoint4D(0,2,3), GeoPoint4D(2,0,3))
          >>> s1.segment_end_in(s2)
          True
        """

        return another.contains_pt(self.end_pt)

    def rotate(self,
               rotation_axis: 'RotationAxis',
               center_point: 'GeoPoint4D' = None
               ) -> 'GeoSegment4D':
        """
        Rotates a segment.
        :param rotation_axis:
        :param center_point:
        :return: the rotated segment
        :rtype: GeoSegment4D

        Examples:
        >>> seg = GeoSegment4D(GeoPoint4D(0,0,0), GeoPoint4D(0,0,1))
        >>> rot_ax = RotationAxis(0, 0, 90)
        >>> seg.rotate(rot_ax)
        GeoSegment4D(start_pt=GeoPoint4D(x=0.0000, y=0.0000, z=0.0000, time=None, epsg_code=-1), end_pt=GeoPoint4D(x=1.0000, y=0.0000, z=0.0000, time=None, epsg_code=-1))
        >>> rot_ax = RotationAxis(0, 0, 180)
        >>> seg.rotate(rot_ax)
        GeoSegment4D(start_pt=GeoPoint4D(x=0.0000, y=0.0000, z=0.0000, time=None, epsg_code=-1), end_pt=GeoPoint4D(x=0.0000, y=0.0000, z=-1.0000, time=None, epsg_code=-1))
        >>> centr_pt = GeoPoint4D(0,0,0.5)
        >>> seg.rotate(rotation_axis=rot_ax, center_point=centr_pt)
        GeoSegment4D(start_pt=GeoPoint4D(x=-0.0000, y=0.0000, z=1.0000, time=None, epsg_code=-1), end_pt=GeoPoint4D(x=0.0000, y=0.0000, z=0.0000, time=None, epsg_code=-1))
        >>> seg = GeoSegment4D(GeoPoint4D(0,0,0), GeoPoint4D(1,1,0))
        >>> centr_pt = GeoPoint4D(1,0,0)
        >>> rot_ax = RotationAxis(0, 90, 90)
        >>> seg.rotate(rotation_axis=rot_ax, center_point=centr_pt)
        GeoSegment4D(start_pt=GeoPoint4D(x=1.0000, y=1.0000, z=0.0000, time=None, epsg_code=-1), end_pt=GeoPoint4D(x=2.0000, y=0.0000, z=-0.0000, time=None, epsg_code=-1))
        >>> seg = GeoSegment4D(GeoPoint4D(1,1,1), GeoPoint4D(0,0,0))
        >>> rot_ax = RotationAxis(135, 0, 180)
        >>> centr_pt = GeoPoint4D(0.5,0.5,0.5)
        >>> seg.rotate(rotation_axis=rot_ax, center_point=centr_pt)
        GeoSegment4D(start_pt=GeoPoint4D(x=0.0000, y=0.0000, z=0.0000, time=None, epsg_code=-1), end_pt=GeoPoint4D(x=1.0000, y=1.0000, z=1.0000, time=None, epsg_code=-1))
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

        return GeoSegment4D(
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
        :rtype: GeoSegment4D
        """

        return cls(
            start_pt=GeoPoint4D.random(lower_boundary, upper_boundary),
            end_pt=GeoPoint4D.random(lower_boundary, upper_boundary)
        )

    '''
    def densify_as_line4d(self,
                          densify_distance
                          ) -> Line4D:
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
    '''

    '''
    def densify_as_pts3d(self,
                         densify_distance
                         ) -> List[GeoPoint4D]:

        return self.densify_as_line3d(densify_distance=densify_distance).pts()
    '''

    '''
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
    '''