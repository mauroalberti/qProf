# -*- coding: utf-8 -*-

from __future__ import division

from builtins import object
from math import sqrt, sin, cos, radians, acos, atan, atan2, degrees
import numpy as np

from .array_utils import point_solution
from .errors import SubparallelLineationException


MIN_SEPARATION_THRESHOLD = 1e-10
MIN_VECTOR_MAGNITUDE = 1e-10
MIN_SCALAR_VALUE = 1e-15
MIN_ANGLE_DEGR_VALUE = 1e-10
MIN_VECTOR_MAGN_DIFF = MIN_SCALAR_VALUE
MIN_ANGLE_DEGR_DISORIENTATION = 5.


class Point(object):
    """
    Cartesian point.
    Dimensions: 3D + time
    """

    def __init__(self, x=np.nan, y=np.nan, z=np.nan, t=np.nan):
        """
        Construct a Point instance given 3 or 4 float values.
        """

        self._p = np.array([x, y, z, t], dtype=np.float64)

    def __repr__(self):

        return "Point({:.4f}, {:.4f}, {:.4f}, {:.4f})".format(self.x, self.y, self.z, self.t)

    @classmethod
    def from_array(cls, a):
        """
        Class method to construct a point from a numpy 1x4 array.

        Example:
          >>> Point.from_array(np.array([1, 0, 1]))
          Point(1.0000, 0.0000, 1.0000, nan)
        """

        obj = cls()

        assert 3 <= a.size <= 4
        b = a.astype(np.float64)
        if b.size == 3:
            c = np.append(b, [np.nan])
        else:
            c = b
        obj._p = c
        return obj

    @property
    def v(self):
        """
        Return values as array

        Example:
          >>> Point(1, 0, 0).v
          array([  1.,   0.,   0.,  nan])
        """

        return self._p

    @property
    def x(self):
        """
        Return x value

        Example:
          >>> Point(1.5, 1, 1).x
          1.5
        """

        return self.v[0]

    @property
    def y(self):
        """
        Return y value

        Example:
          >>> Point(1.5, 3.0, 1).y
          3.0
        """
        return self.v[1]

    @property
    def z(self):
        """
        Return z value

        Example:
          >>> Point(1.5, 3.2, 41.).z
          41.0
        """
        return self.v[2]

    @property
    def t(self):
        """
        Return time value

        Example:
          >>> Point(1.5, 3.2, 41., 22.).t
          22.0
        """
        return self.v[3]

    def clone(self):
        """
        Clone the point.

        Example:
          >>> Point(1, 1, 1).clone()
          Point(1.0000, 1.0000, 1.0000, nan)
        """

        return Point.from_array(self.v)

    def __sub__(self, another):
        """Return point difference

        Example:
          >>> Point(1., 1., 1.) - Point(1., 1., 1.)
          Point(0.0000, 0.0000, 0.0000, nan)
        """

        return Point.from_array(self.v - another.v)

    def __abs__(self):
        """
        Point distance from frame origin.
        todo: make sure it works as intended with nan values

        Example:
          >>> abs(Point(3, 4, 0))
          5.0
          >>> abs(Point(0, 12, 5))
          13.0
        """

        return sqrt(self.x * self.x + self.y * self.y + self.z * self.z)

    def dist_3d(self, another):
        """
        Calculate Euclidean spatial distance between two points.
        todo: make sure it works as intended with nan values
        
        Examples:
          >>> Point(1., 1., 1.).dist_3d(Point(4., 5., 1,))
          5.0
          >>> Point(1, 1, 1, 4).dist_3d(Point(4, 5, 1, 14))
          5.0
        """

        return abs(self - another)

    def dist_2d(self, another):
        """
        Calculate horizontal (2D) distance between two points.
        todo: make sure it works as intended with nan values
        
        Examples:
          >>> Point(1.,1.,1.).dist_2d(Point(4.,5.,7.))
          5.0
        """

        return sqrt((self.x - another.x) ** 2 + (self.y - another.y) ** 2)

    def coincident(self, another, tolerance=MIN_SEPARATION_THRESHOLD):
        """
        Check spatial coincidence of two points
        todo: make sure it works as intended with Points with nan values
        
        Example:
          >>> Point(1., 0., -1.).coincident(Point(1., 1.5, -1.))
          False
          >>> Point(1., 0., 0.).coincident(Point(1., 0., 0.))
          True
        """

        if self.dist_2d(another) > tolerance:
            return False
        elif self.dist_3d(another) > tolerance:
            return False
        else:
            return True

    def translate(self, sx=0.0, sy=0.0, sz=0.0, st=0.0):
        """
        Create a new point shifted by given amount from the self instance.

        Example:
          >>> Point(1, 1, 1).translate(0.5, 1., 1.5)
          Point(1.5000, 2.0000, 2.5000, nan)
       """

        return Point(self.x + sx, self.y + sy, self.z + sz, self.t + st)

    def vect_offset(self, displ_vect):
        """
        Create a new point from the self, with offsets defined by a vector.

        Example:
          >>> Point(1, 2, 0).vect_offset(Vect(10, 5, 0))
          Point(11.0000, 7.0000, 0.0000, nan)
        """

        return Point(self.x + displ_vect.x,
                     self.y + displ_vect.y,
                     self.z + displ_vect.z,
                     self.t)

    @property
    def vector(self):
        """
        Create a vector based on the point coordinates

        Example:
          >>> Point(1, 1, 0, 5).vector
          Vect(1.0000, 1.0000, 0.0000)
        """

        return Vect(self.x, self.y, self.z)

    def delta_time(self, another):
        """
        Calculate the time difference between two points

        Example:
          >>> Point(1,1,1,4).delta_time(Point(1,1,2,5))
          1.0
        """

        return another.t - self.t

    def speed(self, another):
        """
        Calculate the speed required to displace self to another.

        Example:
          >>> Point(1, 1, 1, 4).speed(Point(4, 5, 1, 14))
          0.5
        """

        try:
            return self.dist_3d(another) / self.delta_time(another)
        except:
            return np.Infinity


class Vect(object):
    """
    Cartesian vector, 3D
    Right-handed rectangular Cartesian coordinate system (ENU):
    x axis -> East
    y axis -> North
    z axis -> Up
    """

    def __init__(self, x=np.nan, y=np.nan, z=np.nan):
        """
        Vect constructor
        """

        self._v = np.array([x, y, z], dtype=np.float64)

    @classmethod
    def from_array(cls, a):
        """
        Class method to construct a vector from a numpy 1x3 array.

        Example:
          >>> Vect.from_array(np.array([1,0,1]))
          Vect(1.0000, 0.0000, 1.0000)
        """

        obj = cls()

        assert a.size == 3
        b = a.astype(np.float64)
        obj._v = b
        return obj

    @property
    def v(self):
        """
        Return the vector values as array

        Example:
          >>> Vect(1, 1, 0).v
          array([ 1.,  1.,  0.])
        """

        return self._v

    @property
    def x(self):
        """
        Return x value

        Example:
          >>> Vect(1, 2, 0).x
          1.0
        """

        return self.v[0]

    @property
    def y(self):
        """
        Return y value

        Example:
          >>> Vect(1, 2, 0).y
          2.0
        """

        return self.v[1]

    @property
    def z(self):
        """
        Return z value

        Example:
          >>> Vect(1, 2, 0).z
          0.0
        """

        return self.v[2]

    def __sub__(self, another):
        """
        Return vector difference.

        Example:
          >>> Vect(1., 1., 1.) - Vect(1., 1., 1.)
          Vect(0.0000, 0.0000, 0.0000)
          >>> Vect(0., 1., 4.) - Vect(7., 3., 1.)
          Vect(-7.0000, -2.0000, 3.0000)
        """

        return Vect.from_array(self.v - another.v)

    def __eq__(self, another):
        """
        Return True if vectors are equal.

        Example:
          >>> Vect(1., 1., 1.) == Vect(1., 1., 1.)
          True
        """
        return bool(abs(self - another) < MIN_VECTOR_MAGN_DIFF)

    def __ne__(self, another):
        """
        Return False if vectors are equal.

        Example:
          >>> Vect(1., 1., 1.) != Vect(0., 0., 0.)
          True
        """

        return not self == another

    def __repr__(self):

        return "Vect({:.4f}, {:.4f}, {:.4f})".format(self.x, self.y, self.z)

    def __add__(self, another):
        """
        Sum of two vectors.

        Example:
          >>> Vect(1, 0, 0) + Vect(0, 1, 1)
          Vect(1.0000, 1.0000, 1.0000)
          >>> Vect(1, 1, 1) + Vect(-1, -1, -1)
          Vect(0.0000, 0.0000, 0.0000)
        """

        return Vect.from_array(self.v + another.v)

    def clone(self):
        """
        Clone the vector.

        Example:
          >>> Vect(1, 1, 1).clone()
          Vect(1.0000, 1.0000, 1.0000)
        """
        return Vect.from_array(self.v)

    def __abs__(self):
        """
        Vector magnitude.
        """

        return sqrt(self.x * self.x + self.y * self.y + self.z * self.z)

    @property
    def len_2d(self):
        """
        Vector length projected on the horizontal (xy) plane.

        Example:
          >>> Vect(3, 4, 7).len_2d
          5.0
        """

        return sqrt(self.x * self.x + self.y * self.y)

    @property
    def len_3d(self):
        """
        Vector length projected on the xyz space.

        """

        return sqrt(self.x * self.x + self.y * self.y + self.z * self.z)

    def scale(self, scale_factor):
        """
        Create a scaled vector.

        Example;
          >>> Vect(1,0,1).scale(2.5)
          Vect(2.5000, 0.0000, 2.5000)
        """

        return Vect.from_array(self.v * scale_factor)

    @property
    def versor_full(self):
        """
        Calculate versor.

        Example:
          >>> Vect(5, 0, 0).versor_full
          Vect(1.0000, 0.0000, 0.0000)
          >>> Vect(0, 0, -1).versor_full
          Vect(0.0000, 0.0000, -1.0000)
        """

        return self.scale(1.0 / abs(self))

    @property
    def versor_3d(self):
        """
        Calculate versor in xyz space.

        """

        return self.scale(1.0 / self.len_3d)

    @property
    def versor_2d(self):
        """
        Calculate versor in xy space.

        """

        return self.scale(1.0 / self.len_2d)

    def invert(self):
        """
        Create a new vector with inverted direction.

        Examples:
          >>> Vect(1, 1, 1).invert()
          Vect(-1.0000, -1.0000, -1.0000)
          >>> Vect(2, -1, 4).invert()
          Vect(-2.0000, 1.0000, -4.0000)
        """

        return self.scale(-1)

    @property
    def is_upward(self):
        """
        Check that a vector is upward-directed.
         
        Example:
          >>> Vect(0,0,1).is_upward
          True
          >>> Vect(0,0,-0.5).is_upward
          False
        """

        if self.z > 0.0:
            return True
        else:
            return False

    @property
    def is_downward(self):
        """
        Check that a vector is downward-directed.

        Example:
          >>> Vect(0,0,1).is_downward
          False
          >>> Vect(0,0,-0.5).is_downward
          True
        """

        if self.z < 0.0:
            return True
        else:
            return False

    @property
    def upward(self):
        """
        Calculate a new vector upward-pointing.

        Example:
          >>> Vect(1, 1, 1).upward
          Vect(1.0000, 1.0000, 1.0000)
          >>> Vect(-1, -1, -1).upward
          Vect(1.0000, 1.0000, 1.0000)
        """

        if self.z < 0.0:
            return self.scale(-1.0)
        else:
            return self.clone()

    @property
    def downward(self):
        """
        Calculate a new vector downward-pointing.

        Example:
          >>> Vect(1, 1, 1).downward
          Vect(-1.0000, -1.0000, -1.0000)
          >>> Vect(-1, -1, -1).downward
          Vect(-1.0000, -1.0000, -1.0000)
        """

        if self.z > 0.0:
            return self.scale(-1.0)
        else:
            return self.clone()

    @property
    def slope(self):
        """
        Slope of a vector expressed as degrees.
        Positive when vector is downward pointing or horizontal,
        negative when upward pointing.

        Example:
          >>> Vect(1, 0, -1).slope
          45.0
          >>> Vect(1, 0, 1).slope
          -45.0
          >>> Vect(0, 1, 0).slope
          0.0
          >>> Vect(0, 0, 1).slope
          -90.0
          >>> Vect(0, 0, -1).slope
          90.0
        """

        hlen = self.len_2d
        if hlen == 0.0:
            if self.z > 0.:
                return -90.
            elif self.z < 0.:
                return 90.
            else:
                raise Exception("Zero-valued vector")
        else:
            slope = - degrees(atan(self.z / self.len_2d))
            if abs(slope) > MIN_SCALAR_VALUE:
                return slope
            else:
                return 0.

    @property
    def gvect(self):
        """
        Calculate the geological vector parallel to the Vect instance.
        Trend range: [0°, 360°[
        Plunge range: [-90°, 90°], with negative values for upward-pointing
        geological axes and positive values for downward-pointing axes.

        Examples:
          >>> Vect(0, 1, 1).gvect
          GVect(000.00, -45.00)
          >>> Vect(1, 0, 1).gvect
          GVect(090.00, -45.00)
          >>> Vect(0, 0, 1).gvect
          GVect(000.00, -90.00)
          >>> Vect(0, 0, -1).gvect
          GVect(000.00, +90.00)
          >>> Vect(-1, 0, 0).gvect
          GVect(270.00, +00.00)
          >>> Vect(0, -1, 0).gvect
          GVect(180.00, +00.00)
          >>> Vect(-1, -1, 0).gvect
          GVect(225.00, +00.00)
        """

        if abs(self) < MIN_VECTOR_MAGNITUDE:
            raise Exception("Provided vector has near-zero magnitude")

        plunge = self.slope  # upward pointing -> negative value, downward -> positive

        unit_vect = self.versor_full
        if unit_vect.y == 0. and unit_vect.x == 0:
            trend = 0.
        else:
            trend = (90. - degrees(atan2(unit_vect.y, unit_vect.x))) % 360.

        return GVect(trend, plunge)


    @property
    def gaxis(self):
        """
        Calculate the geological axis parallel to the Vect instance.
        Trend range: [0°, 360°[
        Plunge range: [-90°, 90°], with negative values for upward-pointing
        geological axes and positive values for downward-pointing axes.

        Examples:
          >>> Vect(0, 1, 1).gaxis
          GAxis(000.00, -45.00)
          >>> Vect(1, 0, 1).gaxis
          GAxis(090.00, -45.00)
          >>> Vect(0, 0, 1).gaxis
          GAxis(000.00, -90.00)
          >>> Vect(0, 0, -1).gaxis
          GAxis(000.00, +90.00)
          >>> Vect(-1, 0, 0).gaxis
          GAxis(270.00, +00.00)
          >>> Vect(0, -1, 0).gaxis
          GAxis(180.00, +00.00)
          >>> Vect(-1, -1, 0).gaxis
          GAxis(225.00, +00.00)
        """

        return self.gvect.as_axis()

    def sp(self, another):
        """
        Vector scalar product.

        Examples:
          >>> Vect(1, 0, 0).sp(Vect(1, 0, 0))
          1.0
          >>> Vect(1, 0, 0).sp(Vect(0, 1, 0))
          0.0
          >>> Vect(1, 0, 0).sp(Vect(-1, 0, 0))
          -1.0
        """

        return self.x * another.x + self.y * another.y + self.z * another.z

    def cos_angle(self, another):
        """
        Return the cosine of the angle between two vectors.

        Examples:
          >>> Vect(1,0,0).cos_angle(Vect(0,0,1))
          0.0
          >>> Vect(1,0,0).cos_angle(Vect(-1,0,0))
          -1.0
          >>> Vect(1,0,0).cos_angle(Vect(1,0,0))
          1.0
        """

        try:
            val = self.sp(another) / (abs(self) * abs(another))
            if val > 1.0:
                return 1.0
            elif val < -1.0:
                return -1.0
            else:
                return val
        except ZeroDivisionError:
            return np.nan

    def angle(self, another):
        """
        Calculate angle between two vectors, as degrees
        in 0° - 180° range.

        Example:
          >>> Vect(1, 0, 0).angle(Vect(0, 0, 1))
          90.0
          >>> Vect(1, 0, 0).angle(Vect(-1, 0, 0))
          180.0
          >>> Vect(0, 0, 1).angle(Vect(0, 0, -1))
          180.0
          >>> Vect(1, 1, 1).angle(Vect(1, 1,1 ))
          0.0
        """

        return degrees(acos(self.cos_angle(another)))

    def vp(self, another):
        """
        Vector product.

        Examples:
          >>> Vect(1, 0, 0).vp(Vect(0, 1, 0))
          Vect(0.0000, 0.0000, 1.0000)
          >>> Vect(1, 0, 0).vp(Vect(1, 0, 0))
          Vect(0.0000, 0.0000, 0.0000)
          >>> Vect(1, 0, 0).vp(Vect(-1, 0, 0))
          Vect(0.0000, 0.0000, 0.0000)
        """

        return Vect.from_array(np.cross(self.v, another.v))

    def by_matrix(self, array3x3):
        """
        Matrix multiplication of a vector.

        """

        return Vect.from_array(array3x3.dot(self.v))


class GVect(object):
    """
    Geological vector.
    Defined by trend and plunge (both in degrees):
     - trend: [0.0, 360.0[ clockwise, from 0 (North):
     - plunge: [-90.0, 90.0].
    """

    def __init__(self, srcTrend, srcPlunge):
        """
        Geological vector constructor.
        srcTrend: Trend range: [0.0, 360.0[ clockwise, from 0 (North)
        srcPlunge: Plunge: [-90.0, 90.0], 
        negative value: upward pointing axis, positive values: downward axis;
            
        Example:
          >>> a = GVect(120, -27)
          >>> b = GVect(54, -320)
          Traceback (most recent call last):
          ...
          AssertionError: plunge must be between -90° and +90° (comprised)
        """

        assert -90.0 <= float(srcPlunge) <= 90.0, "plunge must be between -90° and +90° (comprised)"
        self._trend = float(srcTrend) % 360.0
        self._plunge = float(srcPlunge)

    @property
    def tr(self):
        """
        Return trend of the geological direction.
        Range is [0, 360[

        Example:
          >>> GVect(420, -17).tr
          60.0
          >>> GVect(-20, 49).tr
          340.0
        """

        return self._trend

    @property
    def pl(self):
        """
        Return plunge of the geological direction.
        Range is [-90, 90]

        Example:
          >>> GVect(420, -17).pl
          -17.0
        """

        return self._plunge

    @property
    def tp(self):
        """
        Return trend and plunge of the geological direction.

        Example:
          >>> GVect(-90, -45).tp
          (270.0, -45.0)
        """

        return self.tr, self.pl

    def __repr__(self):

        return "GVect({:06.2f}, {:+06.2f})".format(*self.tp)

    def copy(self):
        """
        Return a copy of the GVect instance.
        
        Example:
          >>> GVect(10, 20).copy()
          GVect(010.00, +20.00)
        """

        return GVect(*(self.tp))

    def opposite(self):
        """
        Return the opposite GVect.
        
        Example:
          >>> GVect(0, 30).opposite()
          GVect(180.00, -30.00)
          >>> GVect(315, 10).opposite()
          GVect(135.00, -10.00)
        """

        trend = (self.tr + 180.) % 360.
        plunge = -self.pl

        return GVect(trend, plunge)

    @property
    def versor(self):
        """
        Return the Vect corresponding to the geological vector.

        Examples:
          >>> GVect(0, 90).versor
          Vect(0.0000, 0.0000, -1.0000)
          >>> GVect(0, -90).versor
          Vect(0.0000, 0.0000, 1.0000)
          >>> GVect(90, 90).versor
          Vect(0.0000, 0.0000, -1.0000)
          
        """

        north_coord = cos(radians(self.pl)) * cos(radians(self.tr))
        east_coord = cos(radians(self.pl)) * sin(radians(self.tr))
        down_coord = sin(radians(self.pl))

        return Vect(east_coord, north_coord, -down_coord)

    @property
    def is_upward(self):
        """
        Check whether the instance is pointing upward or horizontal.

        Examples:
          >>> GVect(10, 15).is_upward
          False
          >>> GVect(257.4, 0.0).is_upward
          False
          >>> GVect(90, -45).is_upward
          True
        """

        if self.versor.z > 0.0:
            return True
        else:
            return False

    @property
    def is_downward(self):
        """
        Check whether the instance is pointing downward or horizontal.

        Examples:
          >>> GVect(10, 15).is_downward
          True
          >>> GVect(257.4, 0.0).is_downward
          False
          >>> GVect(90, -45).is_downward
          False
        """

        if self.versor.z < 0.0:
            return True
        else:
            return False

    @property
    def upward(self):
        """
        Return upward-point geological vector.

        Examples:
          >>> GVect(90, -45).upward
          GVect(090.00, -45.00)
          >>> GVect(180, 45).upward
          GVect(000.00, -45.00)
          >>> GVect(0, 0).upward
          GVect(180.00, -00.00)
          >>> GVect(0, 90).upward
          GVect(180.00, -90.00)
        """

        if self.is_upward:
            return self.copy()
        else:
            return self.opposite()

    @property
    def downward(self):
        """
        Return downward-pointing geological vector.

        Examples:
          >>> GVect(90, -45).downward
          GVect(270.00, +45.00)
          >>> GVect(180, 45).downward
          GVect(180.00, +45.00)
          >>> GVect(0, 0).downward
          GVect(180.00, -00.00)
          >>> GVect(0, 90).downward
          GVect(000.00, +90.00)
        """

        if self.is_downward:
            return self.copy()
        else:
            return self.opposite()

    def angle(self, another):
        """
        Calculate angle (in degrees) between the two
        GVect instances.

        Examples:
          >>> GVect(0, 90).angle(GVect(90, 0)) # doctest: +NUMBER
          90.0000000
          >>> GVect(0, 0).angle(GVect(270, 0)) # doctest: +NUMBER
          90.0000000
          >>> GVect(0, 0).angle(GVect(0, 0)) # doctest: +NUMBER 
          0.0000000
          >>> GVect(0, 0).angle(GVect(180, 0)) # doctest: +NUMBER 
          180.0000000
        """

        return self.versor.angle(another.versor_full)

    @property
    def normal_gplane(self):
        """
        Return the geological plane that is normal to the geological vector.

        Examples:
          >>> GVect(0, 45).normal_gplane
          GPlane(180.00, +45.00)
          >>> GVect(0, -45).normal_gplane
          GPlane(000.00, +45.00)
          >>> GVect(0, 90).normal_gplane
          GPlane(180.00, +00.00)
        """

        down_axis = self.downward
        dipdir = (down_axis.tr + 180.0) % 360.0
        dipangle = 90.0 - down_axis.pl

        return GPlane(dipdir, dipangle)

    def common_plane(self, another):
        """
        Calculate GPlane instance defined by the two GVect instances.

        Examples:
          >>> GVect(0, 0).common_plane(GVect(90, 0))
          GPlane(180.00, +00.00)
          >>> GVect(0, 0).common_plane(GVect(90, 90))
          GPlane(090.00, +90.00)
          >>> GVect(45, 0).common_plane(GVect(135, 45))
          GPlane(135.00, +45.00)
          >>> GVect(315, 45).common_plane(GVect(135, 45))
          GPlane(225.00, +90.00)
        """

        normal = self.versor.vp(another.versor_full)
        return normal.gvect.normal_gplane

    def as_axis(self):
        """
        Create GAxis instance with the same attitude as the self instance.

        Example:
          >>> GVect(220, 32).as_axis()
          GAxis(220.00, +32.00)
        """

        return GAxis(*self.tp)

    def vp(self, another):
        """
        Calculate the GVect instance that is normal to the two provided sources.
        Angle between sources must be larger than MIN_ANGLE_DEGR_DISORIENTATION,
        otherwise a SubparallelLineationException will be raised.
        
        Example:
          >>> GVect(0, 0).vp(GVect(4, 0))
          Traceback (most recent call last):
          ...
          SubparallelLineationException: Sources must not be sub- or anti-parallel
          >>> GVect(0, 0).vp(GVect(179, 0))
          Traceback (most recent call last):
          ...
          SubparallelLineationException: Sources must not be sub- or anti-parallel
          >>> GVect(0, 0).vp(GVect(5.1, 0))
          GVect(000.00, +90.00)
          >>> GVect(90, 45).vp(GVect(90, 0))
          GVect(180.00, +00.00)
        """

        if not MIN_ANGLE_DEGR_DISORIENTATION <= self.angle(another) <= 180. - MIN_ANGLE_DEGR_DISORIENTATION:
            raise SubparallelLineationException("Sources must not be sub- or anti-parallel")

        return self.versor.vp(another.versor_full).gvect


class GAxis(GVect):

    def __init__(self, srcTrend, srcPlunge):

        super(GAxis, self).__init__(srcTrend, srcPlunge)

    def __repr__(self):

        return "GAxis({:06.2f}, {:+06.2f})".format(*self.tp)

    def as_vect(self):
        """
        Create GVect instance with the same attitude as the self instance.
        
        Example:
          >>> GAxis(220, 32).as_vect()
          GVect(220.00, +32.00)
        """

        return GVect(*self.tp)

    def angle(self, another):
        """
        Calculate angle (in degrees) between the two
        GAxis instances.

        Examples:
          >>> GAxis(0, 90).angle(GAxis(90, 0)) # doctest: +NUMBER
          90.0000000
          >>> GAxis(0, 0).angle(GAxis(270, 0)) # doctest: +NUMBER
          90.0000000
          >>> GAxis(0, 0).angle(GAxis(0, 0)) # doctest: +NUMBER 
          0.0000000
          >>> GAxis(0, 0).angle(GAxis(180, 0)) # doctest: +NUMBER 
          0.0000000
        """

        angle_vers = self.versor.angle(another.versor_full)
        return min(angle_vers, 180. - angle_vers)

    @property
    def normal_gplane(self):
        """
        Calculate the geological plane that is normal to the given GAxis instance.
        
        Example:
          >>> GAxis(0, 90).normal_gplane
          GPlane(180.00, +00.00)
          >>> GAxis(0, 0).normal_gplane
          GPlane(000.00, +90.00)
          >>> GAxis(45, 45).normal_gplane
          GPlane(225.00, +45.00)
        """

        return self.as_vect().normal_gplane

    @property
    def is_upward(self):
        """
        Check whether the instance is pointing upward or horizontal.

        Examples:
          >>> GAxis(10, 15).is_upward
          False
          >>> GAxis(257.4, 0.0).is_upward
          False
          >>> GAxis(90, -45).is_upward
          True
        """

        return self.as_vect().is_upward

    @property
    def is_downward(self):
        """
        Check whether the instance is pointing downward or horizontal.

        Examples:
          >>> GAxis(10, 15).is_downward
          True
          >>> GAxis(257.4, 0.0).is_downward
          False
          >>> GAxis(90, -45).is_downward
          False
        """

        return self.as_vect().is_downward

    @property
    def upward(self):
        """
        Return upward-point geological axis.

        Examples:
          >>> GAxis(90, -45).upward
          GAxis(090.00, -45.00)
          >>> GAxis(180, 45).upward
          GAxis(000.00, -45.00)
          >>> GAxis(0, 0).upward
          GAxis(180.00, -00.00)
          >>> GAxis(0, 90).upward
          GAxis(180.00, -90.00)
        """

        return self.as_vect().upward.as_axis()

    @property
    def downward(self):
        """
        Return downward-pointing geological vector.

        Examples:
          >>> GAxis(90, -45).downward
          GAxis(270.00, +45.00)
          >>> GAxis(180, 45).downward
          GAxis(180.00, +45.00)
          >>> GAxis(0, 0).downward
          GAxis(180.00, -00.00)
          >>> GAxis(0, 90).downward
          GAxis(000.00, +90.00)
        """

        return self.as_vect().downward.as_axis()

    def common_plane(self, another):
        """
        Calculate GPlane instance defined by the two GAxis instances.

        Examples:
          >>> GAxis(0, 0).common_plane(GAxis(90, 0))
          GPlane(180.00, +00.00)
          >>> GAxis(0, 0).common_plane(GAxis(90, 90))
          GPlane(090.00, +90.00)
          >>> GAxis(45, 0).common_plane(GAxis(135, 45))
          GPlane(135.00, +45.00)
          >>> GAxis(315, 45).common_plane(GAxis(135, 45))
          GPlane(225.00, +90.00)
        """

        return self.as_vect().common_plane(another.as_vect())

    def vp(self, another):
        """
        Calculate the GAxis instance that is perpendicular to the two provided.
        The two source GAxis must not be subparallel (threshold is MIN_ANGLE_DEGR_DISORIENTATION),
        otherwise a SubparallelLineationException will be raised.
        
        Example:
          >>> GAxis(0, 0).vp(GAxis(4, 0))
          Traceback (most recent call last):
          ...
          SubparallelLineationException: Sources must not be sub- or anti-parallel
          >>> GAxis(0, 0).vp(GAxis(180, 0))
          Traceback (most recent call last):
          ...
          SubparallelLineationException: Sources must not be sub- or anti-parallel
          >>> GAxis(90, 0).vp(GAxis(180, 0))
          GAxis(000.00, +90.00)
          >>> GAxis(90, 45).vp(GAxis(180, 0))
          GAxis(270.00, +45.00)
          >>> GAxis(270, 45).vp(GAxis(180, 90))
          GAxis(180.00, -00.00)
        """

        return self.as_vect().vp(another.as_vect()).as_axis()


class Plane(object):
    """
    Cartesian plane.
    Expressed by equation:
    ax + by + cz + d = 0

    Note: Plane is locational - its position in space is defined.
    This contrast with GPlane, defined just by its attitude, but with undefined position

    """

    def __init__(self, a, b, c, d):

        self._a = float(a)
        self._b = float(b)
        self._c = float(c)
        self._d = float(d)

    @property
    def a(self):
        """
        Return a coefficient of a Plane instance.

        Example:
          >>> Plane(1, 0, 0, 2).a
          1.0
        """

        return self._a

    @property
    def b(self):
        """
        Return b coefficient of a Plane instance.

        Example:
          >>> Plane(1, 4, 0, 2).b
          4.0
        """

        return self._b

    @property
    def c(self):
        """
        Return a coefficient of a Plane instance.

        Example:
          >>> Plane(1, 0, 5.4, 2).c
          5.4
        """

        return self._c

    @property
    def d(self):
        """
        Return a coefficient of a Plane instance.

        Example:
          >>> Plane(1, 0, 0, 2).d
          2.0
        """

        return self._d

    @property
    def v(self):
        """
        Return coefficients of a Plane instance.

        Example:
          >>> Plane(1, 1, 7, -4).v
          (1.0, 1.0, 7.0, -4.0)
        """
        return self.a, self.b, self.c, self.d

    @classmethod
    def from_points(cls, pt1, pt2, pt3):
        """
        Create a Plane from three given Point instances.

        Example:
          >>> Plane.from_points(Point(0, 0, 0), Point(1, 0, 0), Point(0, 1, 0))
          Plane(0.0000, 0.0000, 1.0000, 0.0000)
          >>> Plane.from_points(Point(0, 0, 0), Point(0, 1, 0), Point(0, 0, 1))
          Plane(1.0000, 0.0000, 0.0000, 0.0000)
        """

        matr_a = np.array([[pt1.y, pt1.z, 1],
                           [pt2.y, pt2.z, 1],
                           [pt3.y, pt3.z, 1]])

        matr_b = - np.array([[pt1.x, pt1.z, 1],
                             [pt2.x, pt2.z, 1],
                             [pt3.x, pt3.z, 1]])

        matr_c = np.array([[pt1.x, pt1.y, 1],
                           [pt2.x, pt2.y, 1],
                           [pt3.x, pt3.y, 1]])

        matr_d = - np.array([[pt1.x, pt1.y, pt1.z],
                             [pt2.x, pt2.y, pt2.z],
                             [pt3.x, pt3.y, pt3.z]])

        return cls(np.linalg.det(matr_a),
                   np.linalg.det(matr_b),
                   np.linalg.det(matr_c),
                   np.linalg.det(matr_d))

    def __repr__(self):

        return "Plane({:.4f}, {:.4f}, {:.4f}, {:.4f})".format(*self.v)

    @property
    def nversor(self):
        """
        Return the versor normal to the cartesian plane.

        Examples:
          >>> Plane(0, 0, 5, -2).nversor
          Vect(0.0000, 0.0000, 1.0000)
          >>> Plane(0, 7, 0, 5).nversor
          Vect(0.0000, 1.0000, 0.0000)
        """

        return Vect(self.a, self.b, self.c).versor_full

    def gplane_point(self):
        """
        Converts a cartesian plane into a geological plane
        and a point lying in the plane (non-unique solution).

        Examples:
          >>> gpl, pt = Plane(0, 0, 1, -1).gplane_point()
          >>> gpl
          GPlane(000.00, +00.00)
          >>> pt
          Point(0.0000, 0.0000, 1.0000, nan)
        """

        geol_plane = self.nversor.gvect.normal_gplane
        point = Point(*point_solution(np.array([[self.a, self.b, self.c]]),
                                     np.array([-self.d])))
        return geol_plane, point

    def inters_versor(self, another):
        """
        Return intersection versor for two intersecting planes.

        >>> a = Plane(1, 0, 0, 0)
        >>> b = Plane(0, 0, 1, 0)
        >>> a.inters_versor(b)
        Vect(0.0000, -1.0000, 0.0000)
        """

        return self.nversor.vp(another.nversor).versor_full

    def inters_point(self, another):
        """
        Return point on intersection line (obviously non-unique solution)
        for two planes.

        >>> a = Plane(1, 0, 0, 0)
        >>> b = Plane(0, 0, 1, 0)
        >>> a.inters_point(b)
        Point(0.0000, 0.0000, 0.0000, nan)
        """

        # find a point lying on the intersection line (this is a non-unique solution)
        a = np.array([[self.a, self.b, self.c], [another.a, another.b, another.c]])
        b = np.array([-self.d, -another.d])
        x, y, z = point_solution(a, b)

        return Point(x, y, z)

    def is_point_inplane(self, pt):
        """
          Check whether a point lie in a plane.

          >>> pl = Plane(0, 0, 1, 0)
          >>> pt = Point(0, 1, 0)
          >>> pl.is_point_inplane(pt)
          True
        """

        if abs(self.a * pt.x + self.b * pt.y + self.c * pt.z + self.d) < MIN_SCALAR_VALUE:
            return True
        else:
            return False

    def angle(self, another):
        """
        Calculate angle (in degrees) between two planes.

        Examples:
          >>> Plane(1,0,0,0).angle(Plane(0,1,0,0))
          90.0
          >>> Plane(1,0,0,0).angle(Plane(1,0,1,0))
          45.0
          >>> Plane(1,0,0,0).angle(Plane(1,0,0,0))
          0.0
        """

        angle_degr = self.nversor.angle(another.nversor)
        if abs(angle_degr) < MIN_ANGLE_DEGR_VALUE:
            angle_degr = 0.0
        elif angle_degr > 90.0:
            angle_degr = 180.0 - angle_degr

        return angle_degr


class GPlane(object):
    """
    Geological plane.
    Defined by dip direction and dip angle (both in degrees):
     - dip direction: [0.0, 360.0[ clockwise, from 0 (North);
     - dip angle: [0, 90.0]: downward-pointing.
    """

    def __init__(self, srcAzimuth, srcDipAngle, isRHRStrike=False):
        """
        Geological plane constructor.

        @param  srcAzimuth:  Azimuth of the plane (RHR strike or dip direction).
        @type  srcAzimuth:  number or string convertible to float.
        @param  srcDipAngle:  Dip angle of the plane (0-90°).
        @type  srcDipAngle:  number or string convertible to float.
           
        @return:  GPlane.

        Example:
          >>> GPlane(0, 90)
          GPlane(000.00, +90.00)
          >>> GPlane(0, 90, isRHRStrike=True)
          GPlane(090.00, +90.00)
          >>> GPlane(0, 90, True)
          GPlane(090.00, +90.00)
        """

        def rhrstrk2dd(rhr_strk):
            """Converts RHR strike value to dip direction value.

            Example:
                >>> rhrstrk2dd(285.5)
                15.5
            """

            return (rhr_strk + 90.0) % 360.0

        if isRHRStrike:
            self._dipdir = rhrstrk2dd(srcAzimuth)
        else:
            self._dipdir = srcAzimuth % 360.0
        self._dipangle = float(srcDipAngle)

    @property
    def dd(self):
        """
        Return the dip direction of the geological plane.

        Example:
          >>> GPlane(34.2, 89.7).dd
          34.2
        """

        return self._dipdir

    @property
    def da(self):
        """
        Return the dip angle of the geological plane.

        Example:
          >>> GPlane(183, 77).da
          77.0

        """

        return self._dipangle

    @property
    def dda(self):
        """
        Return a tuple storing the dip direction and dip angle values of a geological plane.

        Example:
          >>> gp = GPlane(89.4, 17.2)
          >>> gp.dda
          (89.4, 17.2)
        """
        
        return self.dd, self.da

    @property
    def strike_rhr(self):
        """
        Return the strike according to the right-hand-rule.

        Examples:
          >>> GPlane(90, 45).strike_rhr
          0.0
          >>> GPlane(45, 89).strike_rhr
          315.0
          >>> GPlane(275, 38).strike_rhr
          185.0
          >>> GPlane(0, 38).strike_rhr
          270.0
        """

        return (self.dd - 90.0) % 360.0

    @property
    def strike_lhr(self):
        """
        Return the strike according to the left-hand-rule.

        Examples:
          >>> GPlane(90, 45).strike_lhr
          180.0
          >>> GPlane(45, 89).strike_lhr
          135.0
          >>> GPlane(275, 38).strike_lhr
          5.0
          >>> GPlane(0, 38).strike_lhr
          90.0
        """

        return (self.dd + 90.0) % 360.0

    def __repr__(self):

        return "GPlane({:06.2f}, {:+06.2f})".format(*self.dda)

    @property
    def normal(self):
        """
        Return the geological vector normal to the geological plane,
        pointing in the same direction as the geological plane.

        Example:
            >>> GPlane(90, 55).normal
            GVect(090.00, -35.00)
        """
        
        trend = self.dd % 360.0
        plunge = self.da - 90.0

        return GVect(trend, plunge)

    def plane(self, point):
        """
        Given a GPlane instance and a provided Point instance,
        calculate the corresponding Plane instance.

        Example:
          >>> GPlane(0, 0).plane(Point(0, 0, 0))
          Plane(0.0000, 0.0000, 1.0000, -0.0000)
          >>> GPlane(90, 45).plane(Point(0, 0, 0))
          Plane(0.7071, 0.0000, 0.7071, -0.0000)
          >>> GPlane(0, 90).plane(Point(0, 0, 0))
          Plane(0.0000, 1.0000, -0.0000, -0.0000)
        """

        normal_versor = self.normal.versor
        a, b, c = normal_versor.x, normal_versor.y, normal_versor.z
        d = - (a * point.x + b * point.y + c * point.z)
        return Plane(a, b, c, d)

    def angle(self, another):
        """
        Calculate angle (in degrees) between two geoplanes.

        >>> p1 = GPlane(100.0, 50.0)
        >>> p1.angle(p1)
        0.0

        >>> p2 = GPlane(300.0, 10.0)
        >>> p3 = GPlane(300.0, 90.0)
        >>> p2.angle(p3)
        80.0
        """

        return self.normal.as_axis().angle(another.normal.as_axis())

    def rake_to_gv(self, rake):
        """
        Calculate GVect given a GPlane instance and a rake value.
        The rake is defined according to the Aki and Richards, 1980 conventions:
        rake = 0° -> left-lateral
        rake = 90° -> reverse
        rake = +/- 180° -> right-lateral
        rake = -90° -> normal

        Examples:
          >>> GPlane(180, 45).rake_to_gv(0.0)
          GVect(090.00, +00.00)
          >>> GPlane(180, 45).rake_to_gv(90.0)
          GVect(000.00, -45.00)
          >>> GPlane(180, 45).rake_to_gv(-90.0)
          GVect(180.00, +45.00)
          >>> GPlane(180, 45).rake_to_gv(180.0)
          GVect(270.00, -00.00)
          >>> GPlane(180, 45).rake_to_gv(-180.0)
          GVect(270.00, +00.00)
        """

        rk = radians(rake)
        strk = radians(self.strike_rhr)
        dip = radians(self.da)

        x = cos(rk)*sin(strk)-sin(rk)*cos(dip)*cos(strk)
        y = cos(rk)*cos(strk)+sin(rk)*cos(dip)*sin(strk)
        z = sin(rk) * sin(dip)

        return Vect(x, y, z).gvect


if __name__ == "__main__":

    import doctest
    import numtest  # external module, used in doctest float checks
    doctest.testmod()
