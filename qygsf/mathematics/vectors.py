from math import sqrt, degrees, atan, atan2, acos

import numpy
from numpy.core._multiarray_umath import absolute

from qygsf import Direct
from qygsf.mathematics.defaults import MIN_VECTOR_MAGNITUDE, MIN_SCALAR_VALUE, MIN_VECTOR_MAGN_DIFF


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

    def versor(self):
        """
        Calculate versor.

        Example:
          >>> Vect(5, 0, 0).versor()
          Vect(1.0000, 0.0000, 0.0000)
          >>> Vect(0, 0, -1).versor()
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

        unit_vect = self.versor()
        if unit_vect.y == 0. and unit_vect.x == 0:
            trend = 0.
        else:
            trend = (90. - degrees(atan2(unit_vect.y, unit_vect.x))) % 360.

        return Direct(trend, plunge)


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