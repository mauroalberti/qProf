
import numbers

from ..utils.types import *
from .arrays import *
from .defaults import *


class Vect2D(object):
    """
    Cartesian 2D vector.
    Right-handed rectangular Cartesian coordinate system (ENU):
    x axis -> East
    y axis -> North
    """

    def __init__(self,
                 x: numbers.Real,
                 y: numbers.Real
        ):
        """
        Vect constructor.

        Example;
          >>> Vect2D(1, 0)
          Vect2D(1.0000, 0.0000)
          >>> Vect2D(0, 0)
          Vect2D(0.0000, 0.0000)
          >>> Vect2D(2.2, -19.7)
          Vect2D(2.2000, -19.7000)
        """

        vals = [x, y]

        if any(map(lambda val: not isinstance(val, numbers.Real), vals)):
            raise Exception("Input values must be integer of float")

        if not all(map(math.isfinite, vals)):
            raise Exception("Input values must be finite (#05)")

        self._a = np.array(vals, dtype=np.float64)

    def __abs__(self):
        """
        The abs of a vector.

        :return: numbers.Real
        """

        return self.length

    def __eq__(self,
               another: 'Vect2D'
        ) -> bool:
        """
        Return True if objects are equal.

        Example:
          >>> Vect2D(1., 1.) == Vect2D(1, 1)
          True
          >>> Vect2D(1., -1.) == Vect2D(1, 1)
          False
        """

        if not isinstance(another, Vect2D):
            raise Exception("Instances must be of the same type")
        else:
            return all(
                [
                    self.x == another.x,
                    self.y == another.y
                ]
            )

    def __ne__(self,
               another: 'Vect2D'
        ) -> bool:
        """
        Return False if objects are equal.

        Example:
          >>> Vect2D(1., 1.) != Vect2D(0., 0.)
          True
        """

        if not isinstance(another, Vect2D):
            raise Exception("Instances must be of the same type")
        else:
            return not (self == another)

    @property
    def a(self) -> np.ndarray:
        """
        Return a copy of the object inner array.

        :return: double array of x, y, z values

        Examples:
          >>> np.allclose(Vect2D(4, 3).a, np.array([ 4.,  3.]))
          True
        """

        return np.copy(self._a)

    @property
    def x(self) -> numbers.Real:
        """
        Return x value

        Example:
          >>> Vect2D(1.5, 1).x
          1.5
        """

        return self.a[0]

    @property
    def y(self) -> numbers.Real:
        """
        Return y value

        Example:
          >>> Vect2D(1.5, 3.0).y
          3.0
        """
        return self.a[1]

    def __iter__(self):
        """
        Return the elements of a Vector.

        :return:

        """

        return (i for i in [self.x, self.y])

    def toXY(self
        ) -> Tuple[numbers.Real, numbers.Real]:
        """
        Returns the spatial components as a tuple of two values.

        :return: the spatial components (x, y).
        :rtype: a tuple of two floats.

        Examples:
          >>> Vect2D(1, 0).toXY()
          (1.0, 0.0)
        """

        return self.x, self.y

    def toArray(self) -> np.ndarray:
        """
        Return a double Numpy array representing the point values.

        :return: Numpy array

        Examples:
          >>> np.allclose(Vect2D(1, 2).toArray(), np.array([ 1., 2.]))
          True
        """

        return self.a

    @property
    def length(self) -> numbers.Real:
        """
        Spatial distance of the point from the axis origin.

        :return: distance

        Examples:
          >>> Vect2D(4.0, 3.0).length
          5.0
        """

        return math.sqrt(self.x * self.x + self.y * self.y)

    def delta_x(self,
                another: 'Vect2D'
        ) -> Optional[numbers.Real]:
        """
        Delta between x components of two Vect Instances.

        :return: x difference value.
        :rtype: Optional[numbers.Real].

        Examples:
          >>> Vect2D(1, 2).delta_x(Vect2D(4, 7))
          3.0
        """

        if not isinstance(another, Vect2D):
            return None

        return another.x - self.x

    def delta_y(self,
                another: 'Vect2D'
        ) -> Optional[numbers.Real]:
        """
        Delta between y components of two Vect Instances.

        :return: y difference value.
        :rtype: Optional[numbers.Real].

        Examples:
          >>> Vect2D(1, 2).delta_y(Vect2D(4, 7))
          5.0
        """

        if not isinstance(another, Vect2D):
            return None

        return another.y - self.y

    def scale(self,
              scale_factor: numbers.Real
        ) -> Optional['Vect2D']:
        """
        Create a scaled object.

        Example;
          >>> Vect2D(1, 0).scale(2.5)
          Vect2D(2.5000, 0.0000)
          >>> Vect2D(1, 0).scale(2.5)
          Vect2D(2.5000, 0.0000)
          >>> Vect2D(1, 0).scale(np.nan) is None
          True
          >>> Vect2D(1, 0).scale(np.inf) is None
          True
        """

        if not isinstance(scale_factor, numbers.Real):
            return None

        if not math.isfinite(scale_factor):
            return None

        x, y = arrToTuple(self.a * scale_factor)
        return self.__class__(x, y)

    def invert(self) -> 'Vect2D':
        """
        Create a new object with inverted direction.

        Examples:
          >>> Vect2D(1, 1).invert()
          Vect2D(-1.0000, -1.0000)
          >>> Vect2D(2, -1).invert()
          Vect2D(-2.0000, 1.0000)
        """

        return self.scale(-1)

    def __repr__(self) -> str:

        return f"Vect2D({self.x:.4f}, {self.y:.4f})"

    def __add__(self,
                another: 'Vect2D'
        ) -> 'Vect2D':
        """
        Sum of two vectors.

        :param another: the vector to add
        :return: the sum of the two vectors
        :raise: Exception

        Example:
          >>> Vect2D(1, 0) + Vect2D(0, 1)
          Vect2D(1.0000, 1.0000)
          >>> Vect2D(1, 1) + Vect2D(-1, -1)
          Vect2D(0.0000, 0.0000)
        """

        check_type(another, "Second vector", Vect2D)

        x, y = arrToTuple(self.a + another.a)
        return self.__class__(x, y)

    def __sub__(self,
                another: 'Vect2D'
        ) -> 'Vect2D':
        """Subtract two vectors.

        :param another: the vector to subtract
        :return: the difference between the two vectors
        :raise: Exception

        Example:
          >>> Vect2D(1., 1.) - Vect2D(1., 1.)
          Vect2D(0.0000, 0.0000)
          >>> Vect2D(1., 1.) - Vect2D(1., 1.)
          Vect2D(0.0000, 0.0000)
        """

        check_type(another, "Second vector", Vect2D)

        x, y = arrToTuple(self.a - another.a)
        return self.__class__(x, y)

    @property
    def is_close_to_zero(self) -> bool:
        """
        Check if the Vect instance length is near zero.

        :return: Boolean

        Example:
          >>> Vect2D(1, 2).is_close_to_zero
          False
          >>> Vect2D(0.0, 0.0).is_close_to_zero
          True
        """

        return areClose(self.length, 0)

    @property
    def is_close_to_1(self) -> bool:
        """
        Check if the Vect instance length is near unit.

        :return: Boolean

        Example:
          >>> Vect2D(1, 2).is_close_to_1
          False
          >>> Vect2D(0.0, 1.0).is_close_to_1
          True
        """

        return areClose(self.length, 1)

    @property
    def is_valid(self) -> bool:
        """
        Check if the Vect instance components are not all valid and the xy not all zero-valued.

        :return: Boolean

        Example:
          >>> Vect2D(1, 2).is_valid
          True
          >>> Vect2D(0.0, 0.0).is_valid
          False
        """

        return not self.is_close_to_zero

    def versor(self) -> Optional['Vect2D']:
        """
        Calculate versor in xyz space.

        Example:
          >>> Vect2D(5, 0).versor()
          Vect2D(1.0000, 0.0000)
          >>> Vect2D(0, 0).versor() is None
          True
        """

        if not self.is_valid:
            return None
        else:
            return self.scale(1.0 / self.length)

    def dot_product(self,
                    another: 'Vect2D'
        ) -> numbers.Real:
        """
        Vector scalar multiplication.

        Examples:
          >>> Vect2D(1, 0).dot_product(Vect2D(1, 0))
          1.0
          >>> Vect2D(1, 0).dot_product(Vect2D(0, 1))
          0.0
          >>> Vect2D(1, 0).dot_product(Vect2D(-1, 0))
          -1.0
        """

        return self.x * another.x + self.y * another.y

    def cosine_of_angle(self,
                        another: 'Vect2D'
        ) -> Optional[numbers.Real]:
        """
        Return the cosine of the angle between two vectors.

        Examples:
          >>> Vect2D(1,0).cosine_of_angle(Vect2D(0,1))
          0.0
          >>> Vect2D(1,0).cosine_of_angle(Vect2D(-1,0))
          -1.0
          >>> Vect2D(1,0).cosine_of_angle(Vect2D(1,0))
          1.0
        """

        if not isinstance(another, Vect2D):
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
                          another: 'Vect2D'
        ) -> Optional[numbers.Real]:
        """
        Return the scalar projection of the second vector on the first vector.

        Examples:
          >>> Vect2D(2,0).scalar_projection(Vect2D(1,5))
          1.0
          >>> Vect2D(2,0).scalar_projection(Vect2D(-1,5))
          -1.0
        """

        check_type(another, "Second vector", Vect2D)

        return self.cosine_of_angle(another) * another.length

    def fractional_projection(self,
                              another: 'Vect2D'
        ) -> Optional[numbers.Real]:
        """
        Return the fractional projection of the second vector on the first vector.

        Examples:
          >>> Vect2D(2,0).fractional_projection(Vect2D(1,5))
          0.5
          >>> Vect2D(2,0).fractional_projection(Vect2D(-1,5))
          -0.5
        """

        check_type(another, "Second vector", Vect2D)

        return self.scalar_projection(another) / self.length

    '''
    def angle_as_degrees(self,
                         another: 'Vect2D',
                         unit='d'
                         ) -> Optional[numbers.Real]:
        """
        Calculate angle between two vectors,
        in 0° - 180° range (as degrees).

        Example:
          >>> Vect2D(1, 0).angle_as_degrees(Vect2D(-1, 0))
          180.0
          >>> Vect2D(1, 1).angle_as_degrees(Vect2D(1, 1))
          0.0
        """

        if not isinstance(another, Vect2D):
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
    '''

    '''
    def cross_product(self,
                      another: 'Vect2D'
        ) -> Optional['Vect2D']:
        """
        Vector product (cross product).

        Examples:
          >>> Vect2D(1, 0).cross_product(Vect2D(0, 1))
          Vect2D(0.0000, 0.0000)
          >>> Vect2D(1, 0).cross_product(Vect2D(1, 0))
          Vect2D(0.0000, 0.0000)
          >>> (Vect2D(1, 0).cross_product(Vect2D(-1, 0))).is_close_to_zero
          True
        """

        if not isinstance(another, Vect2D):
            raise Exception("Another instance should be Vect but is {}".format(type(another)))

        x, y = arrToTuple(np.cross(self.a, another.a))
        return Vect2D(x, y)
    '''

    '''
    @property
    def azimuth_degr(self
        ) -> Optional[numbers.Real]:
        """
        The azimuth between the Y axis and the vector, calculated clockwise.

        :return: angle in degrees.
        :rtype: optional numbers.Real.

        Examples:
          >>> Vect2D(0, 1).azimuth_degr
          0.0
          >>> Vect2D(1, 1).azimuth_degr
          45.0
          >>> Vect2D(1, 0).azimuth_degr
          90.0
          >>> Vect2D(1, -1).azimuth_degr
          135.0
          >>> Vect2D(0, -1).azimuth_degr
          180.0
          >>> Vect2D(-1, -1).azimuth_degr
          225.0
          >>> Vect2D(-1, 0).azimuth_degr
          270.0
          >>> Vect2D(-1, 1).azimuth_degr
          315.0
        """

        y_axis = Vect2D(0, 1)
        vector_2d = self.versor()

        if not vector_2d:
            return None

        angle = vector_2d.angle_as_degrees(y_axis)

        z_comp = y_axis.cross_product(vector_2d).z

        if z_comp <= 0.0:
            return angle
        else:
            return 360.0 - angle
    '''

if __name__ == "__main__":

    import doctest
    doctest.testmod()
