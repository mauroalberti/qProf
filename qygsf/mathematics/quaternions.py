
from typing import Union

from math import sqrt, degrees, acos

from pygsf.mathematics.vectors import *


QUAT_NORMALIZ_TOL = 1.0e-6
QUAT_DIVISION_TOL = 1.0e-10
QUAT_MAGN_THRESH = 1.0e-6


class Quaternion(object):
    """
    Quaternion class.
    """

    def __init__(self,
                 w=np.nan,
                 x=np.nan,
                 y=np.nan,
                 z=np.nan
                 ):
        """
        Construct a Quaternion instance.

        Examples:
          >>> Quaternion(1, 0, 1, 0)
          Quaternion(1.00000, 0.00000, 1.00000, 0.00000)
          >>> Quaternion()
          Quaternion(nan, nan, nan, nan)
        """

        self.q = np.array([w, x, y, z], dtype=np.float64)

    def __repr__(self) -> str:
        """
        Instance representation.

        :return:
        :rtype: string
        """

        return "Quaternion({:.5f}, {:.5f}, {:.5f}, {:.5f})".format(self.q[0], self.q[1], self.q[2], self.q[3])

    def components(self
                   ) -> Tuple[numbers.Real, numbers.Real, numbers.Real, numbers.Real]:
        """
        Returns the quaternion xyz as a float tuple.

        :return: tuple of 4 float values
        :rtype: tuple of four float values

        Examples:
          >>> Quaternion(0, 1, 0, 0).components()
          (0.0, 1.0, 0.0, 0.0)
        """

        return self.q[0], self.q[1], self.q[2], self.q[3]

    @property
    def scalar(self) -> numbers.Real:
        """
        Return the scalar component of a quaternion.

        :return: scalar component
        :rtype: numbers.Real.

        Examples:
          >>> Quaternion(1, 2, 0, 3).scalar
          1.0
          >>> Quaternion(0.0, 4.6, 6.2, 3.1).scalar
          0.0
        """

        return self.q[0]

    def vector(self,
        ) -> Vect:
        """
        Return the vector component of the quaternion.

        :return:
        :rtype: Vect

        Examples:
          >>> Quaternion(0.1, 1.2, 3.1, 0.9).vector()
          Vect(1.2000, 3.1000, 0.9000)
          >>> Quaternion(6.1, 4.9, 1.03, 5.12).vector()
          Vect(4.9000, 1.0300, 5.1200)
        """

        return Vect(*self.components()[1:])

    @classmethod
    def fromArray(cls,
                  a
                  ) -> 'Quaternion':
        """
        Class method to construct a quaternion from a numpy 1x4 array.

        Examples:
          >>> Quaternion.fromArray(np.array([1, 0, 1, 0]))
          Quaternion(1.00000, 0.00000, 1.00000, 0.00000)
          >>> Quaternion.fromArray(np.array([7.65, -12.34, -1.0, 2.234]))
          Quaternion(7.65000, -12.34000, -1.00000, 2.23400)
          >>> Quaternion.fromArray(np.array([7.65, -12.34, -1.0]))
          Traceback (most recent call last):
          ...
          Exception: Input array for quaternion must have size of 4
        """

        if a.size != 4:
            raise Exception("Input array for quaternion must have size of 4")

        obj = cls()
        obj.q = a.astype(np.float64)

        return obj

    @classmethod
    def fromVect(cls,
                 vect
                 ) -> 'Quaternion':
        """
        Class method to construct a quaternion from a Vect.

        :param vect: Vector instance
        :return: Quaternion instance

        Examples:
          >>> Quaternion.fromVect(Vect(1, 0, 3))
          Quaternion(0.00000, 1.00000, 0.00000, 3.00000)
        """

        w, x, y, z = 0, vect.x, vect.y, vect.z

        return cls(w, x, y, z)

    @classmethod
    def fromRotMatr(cls,
                    matr
                    ) -> 'Quaternion':
        """
        Class method to construct a quaternion from a 3x3 rotation matrix.
        """

        q0 = sqrt(1 + matr[0, 0] + matr[1, 1] + matr[2, 2]) / 2.0
        q1 = sqrt(1 + matr[0, 0] - matr[1, 1] - matr[2, 2]) / 2.0
        q2 = sqrt(1 - matr[0, 0] + matr[1, 1] - matr[2, 2]) / 2.0
        q3 = sqrt(1 - matr[0, 0] - matr[1, 1] + matr[2, 2]) / 2.0

        q0q1 = (matr[2, 1] - matr[1, 2]) / 4.0
        q0q2 = (matr[0, 2] - matr[2, 0]) / 4.0
        q0q3 = (matr[1, 0] - matr[0, 1]) / 4.0
        q1q2 = (matr[0, 1] + matr[1, 0]) / 4.0
        q1q3 = (matr[0, 2] + matr[2, 0]) / 4.0
        q2q3 = (matr[1, 2] + matr[2, 1]) / 4.0

        if (3 * q0) > (q1 + q2 + q3):

            q1 = q0q1 / q0
            q2 = q0q2 / q0
            q3 = q0q3 / q0

        elif (3 * q1) > (q0 + q2 + q3):

            q0 = q0q1 / q1
            q2 = q1q2 / q1
            q3 = q1q3 / q1

        elif (3 * q2) > (q0 + q1 + q3):

            q0 = q0q2 / q2
            q1 = q1q2 / q2
            q3 = q2q3 / q2

        else:

            q0 = q0q3 / q3
            q1 = q1q3 / q3
            q2 = q2q3 / q3

        w, x, y, z = q0, q1, q2, q3

        return cls(w, x, y, z)

    @classmethod
    def zero(cls) -> 'Quaternion':
        """
        Class method to construct a zero quaternion.

        Examples:
          >>> Quaternion.zero()
          Quaternion(0.00000, 0.00000, 0.00000, 0.00000)
        """

        w, x, y, z = 0, 0, 0, 0

        return Quaternion(w, x, y, z)

    @classmethod
    def identity(cls) -> 'Quaternion':
        """
        Class method to construct an identity quaternion (i.e., zero-rotation).

        Examples:
          >>> Quaternion.identity()
          Quaternion(1.00000, 0.00000, 0.00000, 0.00000)
        """

        w, x, y, z = 1, 0, 0, 0

        return Quaternion(w, x, y, z)

    @classmethod
    def i(cls) -> 'Quaternion':
        """
        Class method to construct the i elementary quaternion.

        Examples:
          >>> Quaternion.i()
          Quaternion(0.00000, 1.00000, 0.00000, 0.00000)
        """

        w, x, y, z = 0, 1, 0, 0

        return Quaternion(w, x, y, z)

    @classmethod
    def j(cls) -> 'Quaternion':
        """
        Class method to construct the j elementary quaternion.

        Examples:
          >>> Quaternion.j()
          Quaternion(0.00000, 0.00000, 1.00000, 0.00000)
        """

        w, x, y, z = 0, 0, 1, 0

        return Quaternion(w, x, y, z)

    @classmethod
    def k(cls) -> 'Quaternion':
        """
        Class method to construct the k elementary quaternion.

        Examples:
          >>> Quaternion.k()
          Quaternion(0.00000, 0.00000, 0.00000, 1.00000)
        """

        w, x, y, z = 0, 0, 0, 1

        return Quaternion(w, x, y, z)

    def __eq__(self,
               another: 'Quaternion'
               ) -> bool:
        """
        Quaternion equality.

        :param another: a Quaternion instance
        :return: Boolean

        Examples:
          >>> Quaternion(1, 1, 3, 0) == Quaternion(0, 7, -2, 4)
          False
          >>> Quaternion(1.0, 1.0, 3.0, 0.0) == Quaternion(1.0, 1.0, 3.0, 0.0)
          True
          >>> Quaternion(1.0, 1.0, 3.0, np.nan) == Quaternion(1.0, 1.0, 3.0, np.nan)
          True
          >>> Quaternion(1.0, 1.0, 3.0, 0.0) == Quaternion(1.0, 1.0, 3.0, -1.0e-20)
          False
        """

        if not isinstance(another, Quaternion):
            raise Exception("Compared instance must be of Quaternion type")

        return ((self.q == another.q) | (np.isnan(self.q) & np.isnan(another.q))).all()

    def __ne__(self,
               another: 'Quaternion'
               ) -> bool:
        """
        Quaternion inequality.

        :param another: a Quaternion instance
        :return: Boolean

        Examples:
          >>> Quaternion(1, 1, 3, 0) != Quaternion(0, 7, -2, 4)
          True
          >>> Quaternion(1.0, 1.0, 3.0, np.nan) != Quaternion(1.0, 1.0, 3.0, np.nan)
          False
        """

        if not isinstance(another, Quaternion):
            raise Exception("Compared instance must be of Quaternion type")

        return not (self == another)

    def __add__(self,
                another: 'Quaternion'
                ) -> 'Quaternion':
        """
        Quaternion sum.

        :param another: Quaternion instance.
        :return: Quaternion instance.

        Examples:
          >>> Quaternion(1, 1, 3, 0) + Quaternion(0, 7, -2, 4)
          Quaternion(1.00000, 8.00000, 1.00000, 4.00000)
          >>> Quaternion(2, 1, np.nan, 3) + Quaternion(3, 2, -2, 1)
          Quaternion(5.00000, 3.00000, nan, 4.00000)
        """

        if not isinstance(another, Quaternion):
            raise Exception("Added instance must be of Quaternion type")

        return Quaternion.fromArray(self.q + another.q)

    def __sub__(self,
                another: 'Quaternion'
                ) -> 'Quaternion':
        """
        Quaternion difference.

        :param another: Quaternion instance.
        :return: Quaternion instance.

        Examples:
          >>> Quaternion(1, 1, 3, 0) - Quaternion(0, 7, -2, 4)
          Quaternion(1.00000, -6.00000, 5.00000, -4.00000)
          >>> Quaternion(np.inf, 1, 3, np.inf) - Quaternion(np.nan, np.nan, -1, 4)
          Quaternion(nan, nan, 4.00000, inf)
        """

        if not isinstance(another, Quaternion):
            raise Exception("Subtracted instance must be of Quaternion type")

        return Quaternion.fromArray(self.q - another.q)

    def multByScalar(self,
                     val: numbers.Real
                     ) -> 'Quaternion':
        """
        Multiplication of a quaternion by a scalar value.

        :param val: scalar multiplicand
        :type val: numbers.Real.
        :return: Quaternion instance

        Examples:
          >>> Quaternion(1, 1, 3, 0).multByScalar(4)
          Quaternion(4.00000, 4.00000, 12.00000, 0.00000)
          >>> Quaternion(1.9, -1.2, 3.6, 4.1).multByScalar(2)
          Quaternion(3.80000, -2.40000, 7.20000, 8.20000)
        """

        if not isinstance(val, numbers.Real):
            raise Exception("Multiplier must be int or float")
        
        return Quaternion.fromArray(self.q * val)

    def __neg__(self) -> 'Quaternion':
        """
        Negative of quaternion.

        :return: Quaternion instance.

        Examples:
          >>> - Quaternion(1, 1, 3, 0)
          Quaternion(-1.00000, -1.00000, -3.00000, -0.00000)
          >>> - Quaternion(1.9, -1.2, 3.6, 4.1)
          Quaternion(-1.90000, 1.20000, -3.60000, -4.10000)
        """

        return self.multByScalar(-1)

    def multByQuater(self,
                     another: 'Quaternion'
                     ) -> 'Quaternion':
        """
        Quaternion multiplication.
        Examples are taken from Kuipers, 2002, chp. 5.

        :param another: quaternion multiplier.
        :type another: Quaternion
        :return: multiplied quaternion.
        :type: Quaternion.

        Examples:
          >>> Quaternion(3, 1, -2, 1).multByQuater(Quaternion(2, -1, 2, 3))
          Quaternion(8.00000, -9.00000, -2.00000, 11.00000)
        """
        
        if not isinstance(another, Quaternion):
            raise Exception("Multiplier must be of Quaternion type")
        
        a = + (self.q[0] * another.q[0]) \
            - (self.q[1] * another.q[1]) \
            - (self.q[2] * another.q[2]) \
            - (self.q[3] * another.q[3])
        
        b = + (self.q[0] * another.q[1]) \
            + (self.q[1] * another.q[0]) \
            + (self.q[2] * another.q[3]) \
            - (self.q[3] * another.q[2])

        c = + (self.q[0] * another.q[2]) \
            - (self.q[1] * another.q[3]) \
            + (self.q[2] * another.q[0]) \
            + (self.q[3] * another.q[1])
        
        d = + (self.q[0] * another.q[3]) \
            + (self.q[1] * another.q[2]) \
            - (self.q[2] * another.q[1]) \
            + (self.q[3] * another.q[0])
                
        return Quaternion(a, b, c, d)

    def multByVect(self,
                   vect: Vect
                   ) -> 'Quaternion':
        """
        Quaternion multiplication by a Vect.

        :param vect: vector multiplier.
        :type: Vect.
        :return: quaternion multiplied by vector
        :rtype: Quaternion instance.
        
        Examples:
        """

        if not isinstance(vect, Vect):
            raise Exception("Multiplier must be of Vect type")

        return self.multByQuater(Quaternion.fromVect(vect))

    def __mul__(self,
                another: [numbers.Real, Vect, 'Quaternion']
                ) -> 'Quaternion':
        """
        Wrapper for quaternion multiplication.
        Some examples are taken from Kuipers, 2002, chp. 5.

        :param another: multiplier.
        :type another: numbers.Real, Vect or Quaternion.
        :return: multiplied quaternion.
        :rtype: Quaternion.

        Examples:
          >>> Quaternion(1, 1, 3, 0) * 3
          Quaternion(3.00000, 3.00000, 9.00000, 0.00000)
          >>> Quaternion(3, 1, -2, 1) * Quaternion(2, -1, 2, 3)
          Quaternion(8.00000, -9.00000, -2.00000, 11.00000)
          >>> Quaternion(1, 1, 3, 0) * Quaternion(1, 0, 0, 0)
          Quaternion(1.00000, 1.00000, 3.00000, 0.00000)
          >>> Quaternion.identity() * Vect(1, 3, 2)
          Quaternion(0.00000, 1.00000, 3.00000, 2.00000)
        """

        if isinstance(another, numbers.Real):
            return self.multByScalar(another)
        elif isinstance(another, Vect):
            return self.multByVect(another)
        elif isinstance(another, Quaternion):
            return self.multByQuater(another)
        else:
            raise Exception("Multiplicand is not number or quaternion")

    @property
    def conjugate(self) -> 'Quaternion':
        """
        Quaternion conjugate.

        :return: conjugate quaternion.
        :rtype: Quaternion.

        Examples:
          >>> Quaternion(1, 1, 3, 0).conjugate
          Quaternion(1.00000, -1.00000, -3.00000, -0.00000)
          >>> Quaternion(2.0, 0.0, -3.3, 17.09).conjugate
          Quaternion(2.00000, -0.00000, 3.30000, -17.09000)
          >>> Quaternion(2.0, 0.0, np.nan, 17.09).conjugate
          Quaternion(2.00000, -0.00000, nan, -17.09000)
        """
        
        a = + self.q[0]
        b = - self.q[1]
        c = - self.q[2]
        d = - self.q[3]

        return Quaternion(a, b, c, d)

    def sqrdNorm(self) -> numbers.Real:
        """
        Squared norm of a quaternion.

        :return: quaternion squared norm.
        :rtype: numbers.Real.

        Examples:
          >>> Quaternion(1, 0, 0, 0).sqrdNorm()
          1.0
          >>> Quaternion(1, 1, 0, 2).sqrdNorm()
          6.0
          >>> Quaternion(2, -1, 2, 3).sqrdNorm()
          18.0
          >>> Quaternion(2, np.nan, 2, 3).sqrdNorm()
          nan
        """

        return self.q[0]**2 + self.q[1]**2 + self.q[2]**2 + self.q[3]**2

    def __abs__(self) -> numbers.Real:
        """
        Quaternion absolute value.

        :return: absolute value (magnitude) of the quaternion.
        :rtype: numbers.Real.

        Examples:
          >>> abs(Quaternion(1, 0, 0, 0))
          1.0
          >>> areClose(abs(Quaternion(2, -1, 2, 3)), sqrt(18.0))
          True
        """

        return sqrt(self.sqrdNorm())

    @property
    def norm(self) -> numbers.Real:
        """
        The norm of the quaternion.
        Equivalent to its absolute value.

        :return: absolute value (magnitude) of the quaternion.
        :rtype: numbers.Real.

        Examples:
          >>> Quaternion(1, 0, 0, 0).norm
          1.0
        """

        return abs(self)

    @property
    def inverse(self) -> 'Quaternion':
        """
        Quaternion inverse.

        :return: quaternion inverse.
        :rtype: Quaternion.

        Examples:
          >>> Quaternion(0, 1, 0, 0).inverse
          Quaternion(0.00000, -1.00000, -0.00000, -0.00000)
          >>> Quaternion(3.2, 2.4, 7.18, 4.3).inverse * Quaternion(3.2, 2.4, 7.18, 4.3)
          Quaternion(1.00000, 0.00000, 0.00000, 0.00000)
        """

        return self.conjugate / self.sqrdNorm()

    def isNormalized(self) -> bool:
        """
        Check if a quaternion is unitary.

        :return: True or False.
        :rtype: bool

        Examples:
          >>> Quaternion(0, 1, 0, 0).isNormalized()
          True
          >>> Quaternion(1, 4, 0, -4).isNormalized()
          False
        """

        return abs(1.0 - sqrt(self.sqrdNorm())) < QUAT_NORMALIZ_TOL

    def divByScalar(self,
                    denominator: numbers.Real
                    ) -> 'Quaternion':
        """
        Division of a quaternion by a scalar.

        :param denominator: divisor.
        :type denominator: numbers.Real.
        :return: division result.
        :rtype: Quaternion.

        Examples:
          >>> Quaternion(1, 1, 3, 0).divByScalar(3)
          Quaternion(0.33333, 0.33333, 1.00000, 0.00000)
          >>> Quaternion(1, 1, 3, 0).divByScalar(1e-11)
          Traceback (most recent call last):
          ...
          Exception: Quaternion division by almost zero value
        """

        if not isinstance(denominator, numbers.Real):
            raise Exception("Quaternion divisor must be integer or float")
        elif abs(denominator) < QUAT_DIVISION_TOL:
            raise Exception("Quaternion division by almost zero value")
        else:
            return Quaternion.fromArray(self.q / denominator)

    def divByQuater(self,
                    another: 'Quaternion'
                    ) -> 'Quaternion':
        """
        Quaternion division by another quaternion.

        :param another: divisor
        :type: Quaternion
        :return: division result
        :rtype: Quaternion

        Examples:
        """

        if not isinstance(another, Quaternion):
            raise Exception("Multiplier must be of Quaternion type")

        return self * (another.conjugate.divByScalar(another.sqrdNorm()))

    def __truediv__(self,
                    another: Union[numbers.Real, 'Quaternion']
                    ) -> 'Quaternion':
        """
        Wrapper for quaternion division.
        This is only compatible with Python 3.

        :param another: divisor.
        :type another: Quaternion
        :return: division result.
        :rtype: Quaternion

        Examples:
          >>> Quaternion(1, 1, 3, 0) / 3
          Quaternion(0.33333, 0.33333, 1.00000, 0.00000)
          >>> Quaternion(1, 1, 3, 0) / Quaternion(1, 1, 3, 0)
          Quaternion(1.00000, 0.00000, 0.00000, 0.00000)
        """

        if isinstance(another, numbers.Real):
            return self.divByScalar(another)
        elif isinstance(another, Quaternion):
            return self.divByQuater(another)
        else:
            raise Exception("Denominator is not number or quaternion")

    def normalize(self) -> Optional['Quaternion']:
        """
        Normalize a quaternion.

        :return: normalized quaternion.
        :rtype: Quaternion.

        Examples:
          >>> Quaternion(0, 4, 0, 0).normalize()
          Quaternion(0.00000, 1.00000, 0.00000, 0.00000)
          >>> Quaternion(0, 4, 0, 8).normalize()
          Quaternion(0.00000, 0.44721, 0.00000, 0.89443)
          >>> areClose(abs(Quaternion(0.2, 17.9, -2.7, 4.3).normalize()), 1.0)
          True
          >>> Quaternion(0.696, 0.322, -0.152, 0.624).normalize()
          Quaternion(0.69580, 0.32191, -0.15196, 0.62382)
        """

        if areClose(self.norm, 0.0):
            raise Exception("Quaternion is null or near null")
        else:
            return self / sqrt(self.sqrdNorm())

    def isCloseTo(self,
            another: 'Quaternion',
            rtol: numbers.Real = 1e-012,
            atol:numbers.Real = 1e-12,
            equal_nan: bool = False,
            equal_inf: bool = False
        ) -> bool:
        """
        Check for quaternion equivalence.

        :param another: Quaternion instance.
        :type another: Quaternion.
        :param rtol: relative tolerance
        :type rtol: numbers.Real.
        :param atol: absolute tolerance
        :type atol: numbers.Real.
        :param equal_nan: nan values are considered equal to themselves.
        :type equal_nan: bool.
        :param equal_inf: inf values are considered equal to themselves
        :type equal_inf: bool.
        :return: True if the two quaternions are close enough, false otherwise.
        :rtype: bool.

        Examples:
          >>> Quaternion(1, 2, 3, 4).isCloseTo(Quaternion(1, 2, 3, 4))
          True
          >>> Quaternion(1, 2, 3, 4).isCloseTo(Quaternion(1, 2.01, 3, 4))
          False
          >>> Quaternion(1, 2, 3, 4).isCloseTo(Quaternion(1, 2.01, 3, 4), atol=1e-1)
          True
          >>> Quaternion(1, 2, 3, np.nan).isCloseTo(Quaternion(1, 2, 3, np.nan), equal_nan=True)
          True
        """

        return arraysAreClose(self.q, another.q, rtol, atol, equal_nan, equal_inf)

    def rotAngle(self) -> numbers.Real:
        """
        Calculate the rotation angle associated with a normalized quaternion.
        Formula from p. 710 in Kagan, Y. Y., 1991. 3-D rotation of double-couple earthquake sources.

        :return: Float

        Quaternion case for Kagan, 1991, p.712:
          >>> areClose(Quaternion(0.696, 0.322, -0.152, 0.624).rotAngle(), 91.8182771683)
          True
          >>> areClose(Quaternion(0.62471, 0.32267, 0.69465, 0.15195).rotAngle(), 102.67846140868497)
          True
        """

        return 2 * degrees(acos(self.normalize().scalar))

    def toRotMatrix(self) -> np.ndarray:
        """
        Computes the rotation matrix from the quaternion xyz.
        Formula as in:
        - Eq. 3.5 in Salamin, E., 1979. Application of quaternions to computation with rotations.
        - Eq. 10 in Kagan, Y. Y., 1991. 3-D rotation of double-couple earthquake sources.
        - Eq. 32 in Kagan, Y. Y., 2008. On geometric complexity of earthquake focal zone and fault system.

        :return: 3x3 numpy array
        """

        q0, q1, q2, q3 = self.normalize().components()

        q0q0 = q0 * q0
        q0q1 = q0 * q1
        q0q2 = q0 * q2
        q0q3 = q0 * q3

        q1q1 = q1 * q1
        q1q2 = q1 * q2
        q1q3 = q1 * q3

        q2q2 = q2 * q2
        q2q3 = q2 * q3

        q3q3 = q3 * q3

        a11 = q0q0 + q1q1 - q2q2 - q3q3
        a12 = 2*(q1q2 - q0q3)
        a13 = 2*(q0q2 + q1q3)

        a21 = 2*(q0q3 + q1q2)
        a22 = q0q0 - q1q1 + q2q2 - q3q3
        a23 = 2*(q2q3 - q0q1)

        a31 = 2*(q1q3 - q0q2)
        a32 = 2*(q0q1 + q2q3)
        a33 = q0q0 - q1q1 - q2q2 + q3q3

        return np.array([(a11, a12, a13),
                     (a21, a22, a23),
                     (a31, a32, a33)])


if __name__ == "__main__":

    import doctest
    doctest.testmod()
