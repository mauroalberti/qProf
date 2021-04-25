
import numbers

from ..utils.types import *
from .arrays import *
from .defaults import *


class Vect3D(object):
    """
    Cartesian 3D vector.
    Right-handed rectangular Cartesian coordinate system (ENU):
    x axis -> East
    y axis -> North
    z axis -> Up
    """

    def __init__(self,
                 x: numbers.Real,
                 y: numbers.Real,
                 z: numbers.Real = 0.0
        ):
        """
        Vect constructor.

        Example;
          >>> Vect3D(1, 0, 1)
          Vect3D(1.0000, 0.0000, 1.0000)
          >>> Vect3D(1, np.nan, 1)
          Traceback (most recent call last):
          ...
          Exception: Input values must be finite (#06)
          >>> Vect3D(1, 0, np.inf)
          Traceback (most recent call last):
          ...
          Exception: Input values must be finite (#06)
          >>> Vect3D(0, 0, 0)
          Vect3D(0.0000, 0.0000, 0.0000)
          >>> Vect3D(2.2, -19.7)
          Vect3D(2.2000, -19.7000, 0.0000)
        """

        vals = [x, y, z]

        if any(map(lambda val: not isinstance(val, numbers.Real), vals)):
            raise Exception("Input values must be integer of float")

        if not all(map(math.isfinite, vals)):
            raise Exception("Input values must be finite (#06)")

        self._a = np.array(vals, dtype=np.float64)

    def __abs__(self):
        """
        The abs of a vector.

        :return: numbers.Real
        """

        return self.length

    def __eq__(self,
               another: 'Vect3D'
        ) -> bool:
        """
        Return True if objects are equal.

        Example:
          >>> Vect3D(1., 1., 1.) == Vect3D(1, 1, 1)
          True
          >>> Vect3D(1., 1., 1.) == Vect3D(1, 1, -1)
          False
        """

        if not isinstance(another, Vect3D):
            raise Exception("Instances must be of the same type")
        else:
            return all(
                [
                    self.x == another.x,
                    self.y == another.y,
                    self.z == another.z
                ]
            )

    def __ne__(self,
               another: 'Vect3D'
        ) -> bool:
        """
        Return False if objects are equal.

        Example:
          >>> Vect3D(1., 1., 1.) != Vect3D(0., 0., 0.)
          True
        """

        if not isinstance(another, Vect3D):
            raise Exception("Instances must be of the same type")
        else:
            return not (self == another)

    @property
    def a(self) -> np.ndarray:
        """
        Return a copy of the object inner array.

        :return: double array of x, y, z values

        Examples:
          >>> np.allclose(Vect3D(4, 3, 7).a, np.array([ 4.,  3.,  7.]))
          True
        """

        return np.copy(self._a)

    @property
    def x(self) -> numbers.Real:
        """
        Return x value

        Example:
          >>> Vect3D(1.5, 1, 1).x
          1.5
        """

        return self.a[0]

    @property
    def y(self) -> numbers.Real:
        """
        Return y value

        Example:
          >>> Vect3D(1.5, 3.0, 1).y
          3.0
        """
        return self.a[1]

    @property
    def z(self) -> numbers.Real:
        """
        Return z value

        Example:
          >>> Vect3D(1.5, 3.2, 41.).z
          41.0
        """
        return self.a[2]

    def __iter__(self):
        """
        Return the elements of a Vector.

        :return:

        """

        return (i for i in [self.x, self.y, self.z])

    def toXYZ(self
        ) -> Tuple[numbers.Real, numbers.Real, numbers.Real]:
        """
        Returns the spatial components as a tuple of three values.

        :return: the spatial components (x, y, z).
        :rtype: a tuple of three floats.

        Examples:
          >>> Vect3D(1, 0, 3).toXYZ()
          (1.0, 0.0, 3.0)
        """

        return self.x, self.y, self.z

    def toArray(self) -> np.ndarray:
        """
        Return a double Numpy array representing the point values.

        :return: Numpy array

        Examples:
          >>> np.allclose(Vect3D(1, 2, 3).toArray(), np.array([ 1., 2., 3.]))
          True
        """

        return self.a

    def pXY(self) -> 'Vect3D':
        """
        Projection on the x-y plane

        :return: projected object instance

        Examples:
          >>> Vect3D(2, 3, 4).pXY()
          Vect3D(2.0000, 3.0000, 0.0000)
        """

        return self.__class__(self.x, self.y, 0.0)

    def pXZ(self) -> 'Vect3D':
        """
        Projection on the x-z plane

        :return: projected object instance

        Examples:
          >>> Vect3D(2, 3, 4).pXZ()
          Vect3D(2.0000, 0.0000, 4.0000)
        """

        return self.__class__(self.x, 0.0, self.z)

    def pYZ(self) -> 'Vect3D':
        """
        Projection on the y-z plane

        :return: projected object instance

        Examples:
          >>> Vect3D(2, 3, 4).pYZ()
          Vect3D(0.0000, 3.0000, 4.0000)
        """

        return self.__class__(0.0, self.y, self.z)

    @property
    def length(self) -> numbers.Real:
        """
        Spatial distance of the point from the axis origin.

        :return: distance
        :rtype: numbers.Real

        Examples:
          >>> Vect3D(4.0, 3.0, 0.0).length
          5.0
        """

        return math.sqrt(self.x * self.x + self.y * self.y + self.z * self.z)

    @property
    def length_xy(self) -> numbers.Real:
        """
        2D distance of the point from the axis origin.

        Example:
          >>> Vect3D(3, 4, 0).length_xy
          5.0
          >>> Vect3D(12, 5, 3).length_xy
          13.0
        """

        return math.sqrt(self.x * self.x + self.y * self.y)

    def delta_x(self,
                another: 'Vect3D'
        ) -> Optional[numbers.Real]:
        """
        Delta between x components of two Vect Instances.

        :return: x difference value.
        :rtype: Optional[numbers.Real].

        Examples:
          >>> Vect3D(1, 2, 3).delta_x(Vect3D(4, 7, 1))
          3.0
        """

        if not isinstance(another, Vect3D):
            return None

        return another.x - self.x

    def delta_y(self,
                another: 'Vect3D'
        ) -> Optional[numbers.Real]:
        """
        Delta between y components of two Vect Instances.

        :return: y difference value.
        :rtype: Optional[numbers.Real].

        Examples:
          >>> Vect3D(1, 2, 3).delta_y(Vect3D(4, 7, 1))
          5.0
        """

        if not isinstance(another, Vect3D):
            return None

        return another.y - self.y

    def delta_z(self,
                another: 'Vect3D'
        ) -> Optional[numbers.Real]:
        """
        Delta between x components of two Vect Instances.

        :return: z difference value.
        :rtype: Optional[numbers.Real].

        Examples:
          >>> Vect3D(1, 2, 3).delta_z(Vect3D(4, 7, 1))
          -2.0
        """

        if not isinstance(another, Vect3D):
            return None

        return another.z - self.z

    def scale(self,
              scale_factor: numbers.Real
        ) -> Optional['Vect3D']:
        """
        Create a scaled object.

        Example;
          >>> Vect3D(1, 0, 1).scale(2.5)
          Vect3D(2.5000, 0.0000, 2.5000)
          >>> Vect3D(1, 0, 1).scale(2.5)
          Vect3D(2.5000, 0.0000, 2.5000)
          >>> Vect3D(1, 0, 1).scale(np.nan) is None
          True
          >>> Vect3D(1, 0, 1).scale(np.inf) is None
          True
        """

        if not isinstance(scale_factor, numbers.Real):
            return None

        if not math.isfinite(scale_factor):
            return None

        x, y, z = arrToTuple(self.a * scale_factor)
        return self.__class__(x, y, z)

    def invert(self) -> 'Vect3D':
        """
        Create a new object with inverted direction.

        Examples:
          >>> Vect3D(1, 1, 1).invert()
          Vect3D(-1.0000, -1.0000, -1.0000)
          >>> Vect3D(2, -1, 4).invert()
          Vect3D(-2.0000, 1.0000, -4.0000)
        """

        return self.scale(-1)

    def __repr__(self) -> str:

        return f"Vect3D({self.x:.4f}, {self.y:.4f}, {self.z:.4f})"

    def __add__(self,
                another: 'Vect3D'
        ) -> 'Vect3D':
        """
        Sum of two vectors.

        :param another: the vector to add
        :type another: Vect3D
        :return: the sum of the two vectors
        :rtype: Vect3D
        :raise: Exception

        Example:
          >>> Vect3D(1, 0, 0) + Vect3D(0, 1, 1)
          Vect3D(1.0000, 1.0000, 1.0000)
          >>> Vect3D(1, 1, 1) + Vect3D(-1, -1, -1)
          Vect3D(0.0000, 0.0000, 0.0000)
        """

        check_type(another, "Second vector", Vect3D)

        x, y, z = arrToTuple(self.a + another.a)
        return self.__class__(x, y, z)

    def __sub__(self,
                another: 'Vect3D'
        ) -> 'Vect3D':
        """Subtract two vectors.

        :param another: the vector to subtract
        :type another: Vect3D
        :return: the difference between the two vectors
        :rtype: Vect3D
        :raise: Exception

        Example:
          >>> Vect3D(1., 1., 1.) - Vect3D(1., 1., 1.)
          Vect3D(0.0000, 0.0000, 0.0000)
          >>> Vect3D(1., 1., 3.) - Vect3D(1., 1., 2.2)
          Vect3D(0.0000, 0.0000, 0.8000)
        """

        check_type(another, "Second vector", Vect3D)

        x, y, z = arrToTuple(self.a - another.a)
        return self.__class__(x, y, z)

    @property
    def is_close_to_zero(self) -> bool:
        """
        Check if the Vect instance length is near zero.

        :return: Boolean

        Example:
          >>> Vect3D(1, 2, 0).is_close_to_zero
          False
          >>> Vect3D(0.0, 0.0, 0.0).is_close_to_zero
          True
        """

        return areClose(self.length, 0)

    @property
    def is_close_to_1(self) -> bool:
        """
        Check if the Vect instance length is near unit.

        :return: Boolean

        Example:
          >>> Vect3D(1, 2, 0).is_close_to_1
          False
          >>> Vect3D(0.0, 1.0, 0.0).is_close_to_1
          True
        """

        return areClose(self.length, 1)

    @property
    def is_valid(self) -> bool:
        """
        Check if the Vect instance components are not all valid and the xyz not all zero-valued.

        :return: Boolean

        Example:
          >>> Vect3D(1, 2, 0).is_valid
          True
          >>> Vect3D(0.0, 0.0, 0.0).is_valid
          False
        """

        return not self.is_close_to_zero

    def versor(self) -> Optional['Vect3D']:
        """
        Calculate versor in xyz space.

        Example:
          >>> Vect3D(5, 0, 0).versor()
          Vect3D(1.0000, 0.0000, 0.0000)
          >>> Vect3D(0, 0, -1).versor()
          Vect3D(0.0000, 0.0000, -1.0000)
          >>> Vect3D(0, 0, 0).versor() is None
          True
        """

        if not self.is_valid:
            return None
        else:
            return self.scale(1.0 / self.length)

    def versor_2d(self) -> Optional['Vect3D']:
        """
        Create 2D versor version of the current vector

        :return: unit vector

        Example:
          >>> Vect3D(7, 0, 10).versor_2d()
          Vect3D(1.0000, 0.0000, 0.0000)
          >>> Vect3D(0, 0, 10).versor_2d() is None
          True
        """

        vXY = self.pXY()
        if vXY.is_valid:
            return self.pXY().versor()
        else:
            return None

    @property
    def is_upward(self) -> Optional[bool]:
        """
        Check that a vector is upward-directed.

        :return: boolean

        Example:
          >>> Vect3D(0, 0, 1).is_upward
          True
          >>> Vect3D(0, 0, -0.5).is_upward
          False
          >>> Vect3D(1, 3, 0).is_upward
          False
          >>> Vect3D(0, 0, 0).is_upward is None
          True
        """

        if not self.is_valid:
            return None
        else:
            return self.z > 0.0

    @property
    def is_downward(self) -> Optional[bool]:
        """
        Check that a vector is downward-directed.

        :return: boolean

        Example:
          >>> Vect3D(0, 0, 1).is_downward
          False
          >>> Vect3D(0, 0, -0.5).is_downward
          True
          >>> Vect3D(1, 3, 0).is_downward
          False
          >>> Vect3D(0, 0, 0).is_downward is None
          True
        """

        if not self.is_valid:
            return None
        else:
            return self.z < 0.0

    def upward(self) -> Optional['Vect3D']:
        """
        Calculate a new upward-pointing vector.

        Example:
          >>> Vect3D(1, 1, 1).upward()
          Vect3D(1.0000, 1.0000, 1.0000)
          >>> Vect3D(-1, -1, -1).upward()
          Vect3D(1.0000, 1.0000, 1.0000)
          >>> Vect3D(0, 0, 0).upward() is None
          True
        """

        if not self.is_valid:
            return None
        elif self.z < 0.0:
            return self.scale(-1.0)
        else:
            return self.scale(1.0)

    def downward(self) -> Optional['Vect3D']:
        """
        Calculate a new vector downward-pointing.

        Example:
          >>> Vect3D(1, 1, 1).downward()
          Vect3D(-1.0000, -1.0000, -1.0000)
          >>> Vect3D(-1, -1, -1).downward()
          Vect3D(-1.0000, -1.0000, -1.0000)
          >>> Vect3D(0, 0, 0).downward() is None
          True
        """

        if not self.is_valid:
            return None
        elif self.z > 0.0:
            return self.scale(-1.0)
        else:
            return self.scale(1.0)

    def slope_degr(self) -> numbers.Real:
        """
        Slope of a vector expressed as degrees.
        Positive when vector is downward pointing,
        negative when upward pointing.

        Example:
          >>> Vect3D(1, 0, -1).slope_degr()
          45.0
          >>> Vect3D(1, 0, 1).slope_degr()
          -45.0
          >>> Vect3D(0, 1, 0).slope_degr()
          0.0
          >>> Vect3D(0, 0, 1).slope_degr()
          -90.0
          >>> Vect3D(0, 0, -1).slope_degr()
          90.0
        """

        hlen = self.length_xy
        if hlen == 0.0:
            if self.z > 0.:
                return -90.
            elif self.z < 0.:
                return 90.
            else:
                raise Exception("Zero-valued vector")
        else:
            slope = - math.degrees(math.atan(self.z / self.length_xy))
            if abs(slope) > MIN_SCALAR_VALUE:
                return slope
            else:
                return 0.

    def dot_product(self,
                    another: 'Vect3D'
        ) -> numbers.Real:
        """
        Vector scalar multiplication.

        Examples:
          >>> Vect3D(1, 0, 0).dot_product(Vect3D(1, 0, 0))
          1.0
          >>> Vect3D(1, 0, 0).dot_product(Vect3D(0, 1, 0))
          0.0
          >>> Vect3D(1, 0, 0).dot_product(Vect3D(-1, 0, 0))
          -1.0
        """

        return self.x * another.x + self.y * another.y + self.z * another.z

    def cosine_of_angle(self,
                        another: 'Vect3D'
        ) -> Optional[numbers.Real]:
        """
        Return the cosine of the angle between two vectors.

        Examples:
          >>> Vect3D(1,0,0).cosine_of_angle(Vect3D(0,0,1))
          0.0
          >>> Vect3D(1,0,0).cosine_of_angle(Vect3D(-1,0,0))
          -1.0
          >>> Vect3D(1,0,0).cosine_of_angle(Vect3D(1,0,0))
          1.0
          >>> Vect3D(0, 0, 0).cosine_of_angle(Vect3D(1,0,0)) is None
          True
          >>> Vect3D(1, 0, 0).cosine_of_angle(Vect3D(0,0,0)) is None
          True
        """

        if not isinstance(another, Vect3D):
            return None

        if not (self.is_valid and another.is_valid):
            return None

        val = self.dot_product(another) / (self.length * another.length)
        if val > 1.0:
            return 1.0
        elif val < -1.0:
            return -1.0
        else:
            return val

    def scalar_projection(self,
                          another: 'Vect3D'
        ) -> Optional[numbers.Real]:
        """
        Return the scalar projection of the second vector on the first vector.

        Examples:
          >>> Vect3D(1,0,0).scalar_projection(Vect3D(0,0,1))
          0.0
          >>> Vect3D(2,0,0).scalar_projection(Vect3D(1,5,0))
          1.0
          >>> Vect3D(2,0,0).scalar_projection(Vect3D(-1,5,0))
          -1.0
          >>> Vect3D(4,0,0).scalar_projection(Vect3D(7.5, 19.2, -14.72))
          7.5
        """

        check_type(another, "Second vector", Vect3D)

        return self.cosine_of_angle(another) * another.length

    def fractional_projection(self,
                              another: 'Vect3D'
        ) -> Optional[numbers.Real]:
        """
        Return the fractional projection of the second vector on the first vector.

        Examples:
          >>> Vect3D(1,0,0).fractional_projection(Vect3D(0,0,1))
          0.0
          >>> Vect3D(2,0,0).fractional_projection(Vect3D(1,5,0))
          0.5
          >>> Vect3D(2,0,0).fractional_projection(Vect3D(-1,5,0))
          -0.5
        """

        check_type(another, "Second vector", Vect3D)

        return self.scalar_projection(another) / self.length

    def angle_as_degrees(self,
                         another: 'Vect3D',
                         unit='d'
                         ) -> Optional[numbers.Real]:
        """
        Calculate angle between two vectors,
        in 0° - 180° range (as degrees).

        Example:
          >>> Vect3D(1, 0, 0).angle_as_degrees(Vect3D(0, 0, 1))
          90.0
          >>> Vect3D(1, 0, 0).angle_as_degrees(Vect3D(-1, 0, 0))
          180.0
          >>> Vect3D(0, 0, 1).angle_as_degrees(Vect3D(0, 0, -1))
          180.0
          >>> Vect3D(1, 1, 1).angle_as_degrees(Vect3D(1, 1,1))
          0.0
          >>> Vect3D(0, 0, 0).angle_as_degrees(Vect3D(1,0,0)) is None
          True
          >>> Vect3D(1, 0, 0).angle_as_degrees(Vect3D(0,0,0)) is None
          True
        """

        if not isinstance(another, Vect3D):
            return None

        if not (self.is_valid and another.is_valid):
            return None

        angle_rad = math.acos(self.cosine_of_angle(another))
        if unit == 'd':
            return math.degrees(angle_rad)
        elif unit == 'r':
            return angle_rad
        else:
            return None

    def cross_product(self,
                      another: 'Vect3D'
        ) -> Optional['Vect3D']:
        """
        Vector product (cross product).

        Examples:
          >>> Vect3D(1, 0, 0).cross_product(Vect3D(0, 1, 0))
          Vect3D(0.0000, 0.0000, 1.0000)
          >>> Vect3D(1, 0, 0).cross_product(Vect3D(1, 0, 0))
          Vect3D(0.0000, 0.0000, 0.0000)
          >>> (Vect3D(1, 0, 0).cross_product(Vect3D(-1, 0, 0))).is_close_to_zero
          True
        """

        if not isinstance(another, Vect3D):
            raise Exception("Another instance should be Vect but is {}".format(type(another)))

        x, y, z = arrToTuple(np.cross(self.a[:3], another.a[:3]))
        return Vect3D(x, y, z)

    def multiply_by_matrix(self,
                          array3x3: np.ndarray
        ) -> 'Vect3D':
        """
        Matrix multiplication of a vector.

        """

        x, y, z = arrToTuple(array3x3.dot(self.a))
        return Vect3D(x, y, z)

    @property
    def azimuth_degr(self
        ) -> Optional[numbers.Real]:
        """
        The azimuth between the Y axis and the vector, calculated clockwise.

        :return: angle in degrees.
        :rtype: optional numbers.Real.

        Examples:
          >>> Vect3D(0, 1, 0).azimuth_degr
          0.0
          >>> Vect3D(1, 1, 0).azimuth_degr
          45.0
          >>> Vect3D(1, 0, 0).azimuth_degr
          90.0
          >>> Vect3D(1, -1, 0).azimuth_degr
          135.0
          >>> Vect3D(0, -1, 0).azimuth_degr
          180.0
          >>> Vect3D(-1, -1, 0).azimuth_degr
          225.0
          >>> Vect3D(-1, 0, 0).azimuth_degr
          270.0
          >>> Vect3D(-1, 1, 0).azimuth_degr
          315.0
          >>> Vect3D(0, 0, 1).azimuth_degr is None
          True
          >>> Vect3D(0, 0, -1).azimuth_degr is None
          True
        """

        y_axis = Vect3D(0, 1, 0)
        vector_2d = self.versor_2d()

        if not vector_2d:
            return None

        angle = vector_2d.angle_as_degrees(y_axis)

        z_comp = y_axis.cross_product(vector_2d).z

        if z_comp <= 0.0:
            return angle
        else:
            return 360.0 - angle


if __name__ == "__main__":

    import doctest
    doctest.testmod()
