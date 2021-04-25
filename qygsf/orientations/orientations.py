
import random
from math import *


from ..defaults import *

from ..mathematics.quaternions import *
from ..mathematics.utils import *

from .direct_utils import *


class Azim(object):
    """
    Azim class
    """

    def __init__(self,
        val: numbers.Real,
        unit: str = 'd'
    ):
        """
        Creates an azimuth instance.

        :param val: azimuth value
        :param unit: angle measurement unit, 'd' (default, stands for decimal degrees) or 'r' (stands for radians)

        Examples:
          >>> Azim(10)
          Azim(10.00°)
          >>> Azim(370)
          Azim(10.00°)
          >>> Azim(pi/2, unit='r')
          Azim(90.00°)
        """

        # unit check
        if unit not in ("d", "r"):
            raise Exception(f"Unit input must be 'd' or 'r' but {unit} got")

        if not (isinstance(val, numbers.Real)):
            raise Exception(f"Input azimuth value must be int/float but type {type(val)} got")
        elif not isfinite(val):
            raise Exception(f"Input azimuth value must be finite but {val} got")

        if unit == 'd':
            val = radians(val)

        self.a = val % (2*pi)

    @property
    def d(self) -> numbers.Real:
        """
        Returns the angle in decimal degrees.

        :return: angle in decimal degrees

        Example:
          >>> Azim(10).d
          10.0
          >>> Azim(pi/2, unit='r').d
          90.0
        """

        return degrees(self.a)

    @property
    def r(self) -> numbers.Real:
        """
        Returns the angle in radians.

        :return: angle in radians

        Example:
          >>> Azim(180).r
          3.141592653589793
        """

        return self.a

    @classmethod
    def fromXY(cls,
        x: numbers.Real,
        y: numbers.Real
    ) -> 'Azim':
        """
        Calculates azimuth given cartesian components.

        :param cls: class
        :param x: x component
        :param y: y component
        :return: Azim instance

        Examples:
          >>> Azim.fromXY(1, 1)
          Azim(45.00°)
          >>> Azim.fromXY(1, -1)
          Azim(135.00°)
          >>> Azim.fromXY(-1, -1)
          Azim(225.00°)
          >>> Azim.fromXY(-1, 1)
          Azim(315.00°)
          >>> Azim.fromXY(0, 0)
          Azim(0.00°)
          >>> Azim.fromXY(0, np.nan)
          Traceback (most recent call last):
          ...
          Exception: Input x and y values must be finite
        """

        # input vals checks
        vals = [x, y]
        if not all(map(lambda val: isinstance(val, numbers.Real), vals)):
            raise Exception("Input x and y values must be integer or float")
        elif not all(map(isfinite, vals)):
            raise Exception("Input x and y values must be finite")

        angle = atan2(x, y)
        return cls(angle, unit='r')

    def __repr__(self) -> str:

        return "Azim({:.2f}°)".format(self.d)

    def toXY(self
    ) -> Tuple[numbers.Real, numbers.Real]:
        """
        Converts an azimuth to x-y components.

        :return: a tuple storing x and y values:
        :type: tuple of two floats

        Examples:
          >>> apprFTuple(Azim(0).toXY())
          (0.0, 1.0)
          >>> apprFTuple(Azim(90).toXY())
          (1.0, 0.0)
          >>> apprFTuple(Azim(180).toXY())
          (0.0, -1.0)
          >>> apprFTuple(Azim(270).toXY())
          (-1.0, 0.0)
          >>> apprFTuple(Azim(360).toXY())
          (0.0, 1.0)
        """

        return sin(self.a), cos(self.a)


class Plunge(object):
    """
    Class representing a plunge
    """

    def __init__(self,
        val: numbers.Real,
        unit: str='d'
    ):
        """
        Creates a Plunge instance.

        :param val: plunge value
        :param unit: angle measurement unit, decimal degrees ('d') or radians ('r')

        Examples:
          >>> Plunge(10)
          Plunge(10.00°)
          >>> Plunge(np.nan)
          Traceback (most recent call last):
          ...
          Exception: Input plunge value must be finite
          >>> Plunge(-100)
          Traceback (most recent call last):
          ...
          Exception: Input value in degrees must be between -90° and 90°
         """

        # unit check
        if unit not in ('d', 'r'):
            raise Exception("Unit input must be 'd' (for degrees) or 'r' (for radians)")

        # val check
        if not (isinstance(val, numbers.Real)):
            raise Exception("Input plunge value must be int/float")
        elif not isfinite(val):
            raise Exception("Input plunge value must be finite")
        if unit == 'd' and not (-90.0 <= val <= 90.0):
            raise Exception("Input value in degrees must be between -90° and 90°")
        elif unit == 'r' and not (-pi/2 <= val <= pi/2):
            raise Exception("Input value in radians must be between -pi/2 and pi/2")

        if unit == 'd':
            val = radians(val)

        self.p = val

    @property
    def d(self):
        """
        Returns the angle in decimal degrees.

        :return: angle in decimal degrees

        Example:
          >>> Plunge(10).d
          10.0
          >>> Plunge(-pi/2, unit='r').d
          -90.0
        """

        return degrees(self.p)

    @property
    def r(self):
        """
        Returns the angle in radians.

        :return: angle in radians

        Example:
          >>> Plunge(90).r
          1.5707963267948966
          >>> Plunge(45).r
          0.7853981633974483
        """

        return self.p

    @classmethod
    def fromHZ(cls,
               h: numbers.Real,
               z: numbers.Real
        ) -> 'Plunge':
        """
        Calculates plunge from h and z components.

        :param cls: class
        :param h: horizontal component (always positive)
        :param z: vertical component (positive upward)
        :return: Plunge instance

        Examples:
          >>> Plunge.fromHZ(1, 1)
          Plunge(-45.00°)
          >>> Plunge.fromHZ(1, -1)
          Plunge(45.00°)
          >>> Plunge.fromHZ(0, 1)
          Plunge(-90.00°)
          >>> Plunge.fromHZ(0, -1)
          Plunge(90.00°)
          >>> Plunge.fromHZ(-1, 0)
          Traceback (most recent call last):
          ...
          Exception: Horizontal component cannot be negative
          >>> Plunge.fromHZ(0, 0)
          Traceback (most recent call last):
          ...
          Exception: Input h and z values cannot be both zero
        """

        # input vals check

        vals = [h, z]
        if not all(map(lambda val: isinstance(val, numbers.Real), vals)):
            raise Exception("Input h and z values must be integer or float")
        elif not all(map(isfinite, vals)):
            raise Exception("Input h and z values must be finite")

        if h == 0.0 and z == 0.0:
            raise Exception("Input h and z values cannot be both zero")
        elif h < 0.0:
            raise Exception("Horizontal component cannot be negative")

        angle = atan2(-z, h)

        return cls(angle, unit='r')

    def __repr__(self) -> str:

        return "Plunge({:.2f}°)".format(self.d)

    def toHZ(self):

        """
        Converts an azimuth to h-z components.

        :return: a tuple storing h (horizontal) and z values:
        :type: tuple of two floats

        Examples:
          >>> apprFTuple(Plunge(0).toHZ())
          (1.0, 0.0)
          >>> apprFTuple(Plunge(90).toHZ())
          (0.0, -1.0)
          >>> apprFTuple(Plunge(-90).toHZ())
          (0.0, 1.0)
          >>> apprFTuple(Plunge(-45).toHZ(), ndec=6)
          (0.707107, 0.707107)
          >>> apprFTuple(Plunge(45).toHZ(), ndec=6)
          (0.707107, -0.707107)
        """

        return cos(self.p), -sin(self.p)

    @property
    def is_upward(self) -> bool:
        """
        Check whether the instance is pointing upward or horizontal.

        Examples:
          >>> Plunge(10).is_upward
          False
          >>> Plunge(0.0).is_upward
          False
          >>> Plunge(-45).is_upward
          True
        """

        return self.r < 0.0

    @property
    def is_downward(self) -> bool:
        """
        Check whether the instance is pointing downward or horizontal.

        Examples:
          >>> Plunge(15).is_downward
          True
          >>> Plunge(0.0).is_downward
          False
          >>> Plunge(-45).is_downward
          False
        """

        return self.r > 0.0


class Direct:
    """
    Class describing a direction, expressed as a polar direction in degrees.
    """

    def __init__(self,
                 az: numbers.Real,
                 pl: numbers.Real
                 ):
        """
        Creates a polar direction instance.

        :param az: the azimuth value in decimal degrees
        :param pl: the plunge value in decimal degrees
        """

        self._az = Azim(az)
        self._pl = Plunge(pl)

    @property
    def d(self) -> Tuple[numbers.Real, numbers.Real]:
        """
        Returns azimuth and plunge in decimal degrees as a tuple.

        :return: tuple of azimuth and plunge in decimal degrees

        Example:
          >>> Direct(100, 20).d
          (100.0, 20.0)
        """

        return self.az.d, self.pl.d

    @property
    def r(self) -> Tuple[numbers.Real, numbers.Real]:
        """
        Returns azimuth and plunge in radians as a tuple.

        :return: tuple of azimuth and plunge in radians

        Example:
          >>> Direct(90, 45).r
          (1.5707963267948966, 0.7853981633974483)
        """

        return self.az.r, self.pl.r

    @property
    def az(self) -> Azim:
        """
        Returns the Azim instance.

        :return: Azim
        """

        return self._az

    @property
    def pl(self) -> Plunge:
        """
        Returns the plunge instance.

        :return: Plunge
        """

        return self._pl

    '''
    @classmethod
    def fromAzPl(cls,
                 az: numbers.Real,
                 pl: numbers.Real,
                 unit='d'
                 ):
        """
        Class constructor from trend and plunge.

        :param az: trend value
        :param pl: plunge value
        :param unit: measurement unit, in degrees ('d') or radians ('r')
        :return: Orientation instance

        Examples:
          >>> Direct(30, 40)
          Direct(az: 30.00°, pl: 40.00°)
          >>> Direct(370, 80)
          Direct(az: 10.00°, pl: 80.00°)
          >>> Direct(pi/2, pi/4, unit='r')
          Direct(az: 90.00°, pl: 45.00°)
          >>> Direct(280, -100)
          Traceback (most recent call last):
          ...
          Exception: Input value in degrees must be between -90° and 90°
          >>> Direct("10", 0)
          Traceback (most recent call last):
          ...
          Exception: Input azimuth value must be int/float
          >>> Direct(100, np.nan)
          Traceback (most recent call last):
          ...
          Exception: Input plunge value must be finite
        """

        azim = Azim(az, unit=unit)
        plng = Plunge(pl, unit=unit)

        return cls(azim, plng)
    '''

    @classmethod
    def _from_xyz(cls,
                  x: numbers.Real,
                  y: numbers.Real,
                  z: numbers.Real
    ) -> 'Direct':
        """
        Private class constructor from three Cartesian values. Note: norm of components is unit.

        :param x: x component
        :param y: y component
        :param z: z component
        :return: Orientation instance
        """

        h = sqrt(x*x + y*y)

        az = Azim.fromXY(x, y)
        pl = Plunge.fromHZ(h, z)

        return cls(az.d, pl.d)

    @classmethod
    def fromXYZ(cls,
                x: numbers.Real,
                y: numbers.Real,
                z: numbers.Real
    ) -> 'Direct':
        """
        Class constructor from three generic Cartesian values.

        :param x: x component
        :param y: y component
        :param z: z component
        :return: Orientation instance

        Examples:
          >>> Direct.fromXYZ(1, 0, 0)
          Direct(az: 90.00°, pl: -0.00°)
          >>> Direct.fromXYZ(0, 1, 0)
          Direct(az: 0.00°, pl: -0.00°)
          >>> Direct.fromXYZ(0, 0, 1)
          Direct(az: 0.00°, pl: -90.00°)
          >>> Direct.fromXYZ(0, 0, -1)
          Direct(az: 0.00°, pl: 90.00°)
          >>> Direct.fromXYZ(1, 1, 0)
          Direct(az: 45.00°, pl: -0.00°)
          >>> Direct.fromXYZ(0.5, -0.5, -0.7071067811865476)
          Direct(az: 135.00°, pl: 45.00°)
          >>> Direct.fromXYZ(-0.5, 0.5, 0.7071067811865476)
          Direct(az: 315.00°, pl: -45.00°)
          >>> Direct.fromXYZ(0, 0, 0)
          Traceback (most recent call last):
          ...
          Exception: Input components have near-zero values
        """

        mag, norm_xyz = normXYZ(x, y, z)

        if norm_xyz is None:
            raise Exception("Input components have near-zero values")

        return cls._from_xyz(*norm_xyz)

    @classmethod
    def fromVect(cls,
                 vect: Vect3D
    ) -> [None, 'Direct', 'Axis']:
        """
        Calculate the polar direction parallel to the Vect instance.
        Trend range: [0°, 360°[
        Plunge range: [-90°, 90°], with negative values for upward-pointing
        geological axes and positive values for downward-pointing axes.

        Examples:
          >>> Direct.fromVect(Vect3D(1, 1, 1))
          Direct(az: 45.00°, pl: -35.26°)
          >>> Direct.fromVect(Vect3D(0, 1, 1))
          Direct(az: 0.00°, pl: -45.00°)
          >>> Direct.fromVect(Vect3D(1, 0, 1))
          Direct(az: 90.00°, pl: -45.00°)
          >>> Direct.fromVect(Vect3D(0, 0, 1))
          Direct(az: 0.00°, pl: -90.00°)
          >>> Direct.fromVect(Vect3D(0, 0, -1))
          Direct(az: 0.00°, pl: 90.00°)
          >>> Direct.fromVect(Vect3D(-1, 0, 0))
          Direct(az: 270.00°, pl: -0.00°)
          >>> Direct.fromVect(Vect3D(0, -1, 0))
          Direct(az: 180.00°, pl: -0.00°)
          >>> Direct.fromVect(Vect3D(-1, -1, 0))
          Direct(az: 225.00°, pl: -0.00°)
          >>> Direct.fromVect(Vect3D(0, 0, 0))
          Traceback (most recent call last):
          ...
          Exception: Input components have near-zero values
        """

        x, y, z = vect.toXYZ()
        return cls.fromXYZ(x, y, z)

    def __repr__(self) -> str:

        return "Direct(az: {:.2f}°, pl: {:.2f}°)".format(*self.d)

    def toXYZ(self) -> Tuple[numbers.Real, numbers.Real, numbers.Real]:
        """
        Converts a direction to a tuple of x, y and z cartesian components (with unit norm).

        :return: tuple of x, y and z components.

        Examples:
          >>> apprFTuple(Direct(90, 0).toXYZ())
          (1.0, 0.0, 0.0)
          >>> apprFTuple(Direct(135, 45).toXYZ(), ndec=6)
          (0.5, -0.5, -0.707107)
          >>> apprFTuple(Direct(135, 0).toXYZ(), ndec=6)
          (0.707107, -0.707107, 0.0)
          >>> apprFTuple(Direct(180, 45).toXYZ(), ndec=6)
          (0.0, -0.707107, -0.707107)
          >>> apprFTuple(Direct(225, -45).toXYZ(), ndec=6)
          (-0.5, -0.5, 0.707107)
          >>> apprFTuple(Direct(270, 90).toXYZ(), ndec=6)
          (0.0, 0.0, -1.0)
        """

        x, y = self.az.toXY()
        h, z = self.pl.toHZ()

        return x*h, y*h, z

    def copy(self):
        """
        Return a copy of the instance.

        Example:
          >>> Direct(10, 20).copy()
          Direct(az: 10.00°, pl: 20.00°)
        """

        return self.__class__(self.az.d, self.pl.d)

    def opposite(self):
        """
        Return the opposite direction.

        Example:
          >>> Direct(0, 30).opposite()
          Direct(az: 180.00°, pl: -30.00°)
          >>> Direct(315, 10).opposite()
          Direct(az: 135.00°, pl: -10.00°)
          >>> Direct(135, 0).opposite()
          Direct(az: 315.00°, pl: -0.00°)
        """

        az, pl = self.r

        az = (az + pi) % (2*pi)
        pl = -pl

        return self.__class__(degrees(az), degrees(pl))

    def mirrorHoriz(self):
        """
        Return the mirror Orientation using a horizontal plane.

        Example:
          >>> Direct(0, 30).mirrorHoriz()
          Direct(az: 0.00°, pl: -30.00°)
          >>> Direct(315, 10).mirrorHoriz()
          Direct(az: 315.00°, pl: -10.00°)
          >>> Direct(135, 0).mirrorHoriz()
          Direct(az: 135.00°, pl: -0.00°)
        """

        az = self.az.r
        pl = -self.pl.r

        return self.__class__(degrees(az), degrees(pl))

    @property
    def colatNorth(self) -> numbers.Real:
        """
        Calculates the colatitude from the North (top).

        :return: an angle between 0 and 180 (in degrees).
        :rtype: numbers.Real.

        Examples:
          >>> Direct(320, 90).colatNorth
          180.0
          >>> Direct(320, 45).colatNorth
          135.0
          >>> Direct(320, 0).colatNorth
          90.0
          >>> Direct(320, -45).colatNorth
          45.0
          >>> Direct(320, -90).colatNorth
          0.0
        """

        return plng2colatTop(self.pl.d)

    @property
    def colatSouth(self) -> numbers.Real:
        """
        Calculates the colatitude from the South (bottom).

        :return: an angle between 0 and 180 (in degrees).
        :rtype: numbers.Real.

        Examples:
          >>> Direct(320, 90).colatSouth
          0.0
          >>> Direct(320, 45).colatSouth
          45.0
          >>> Direct(320, 0).colatSouth
          90.0
          >>> Direct(320, -45).colatSouth
          135.0
          >>> Direct(320, -90).colatSouth
          180.0
        """

        return plng2colatBottom(self.pl.d)

    def as_versor(self) -> Vect3D:
        """
        Return the unit vector corresponding to the Direct instance.

        Examples:
          >>> Direct(0, 90).as_versor()
          Vect3D(0.0000, 0.0000, -1.0000)
          >>> Direct(0, -90).as_versor()
          Vect3D(0.0000, 0.0000, 1.0000)
          >>> Direct(90, 90).as_versor()
          Vect3D(0.0000, 0.0000, -1.0000)
        """

        az, pl = self.r
        cos_az, cos_pl = cos(az), cos(pl)
        sin_az, sin_pl = sin(az), sin(pl)
        north_coord = cos_pl * cos_az
        east_coord = cos_pl * sin_az
        down_coord = sin_pl

        return Vect3D(
            x=east_coord,
            y=north_coord,
            z=-down_coord
        )

    @property
    def is_upward(self) -> bool:
        """
        Check whether the instance is pointing upward or horizontal.

        Examples:
          >>> Direct(10, 15).is_upward
          False
          >>> Direct(257.4, 0.0).is_upward
          False
          >>> Direct(90, -45).is_upward
          True
        """

        return self.pl.is_upward

    @property
    def is_downward(self) -> bool:
        """
        Check whether the instance is pointing downward or horizontal.

        Examples:
          >>> Direct(10, 15).is_downward
          True
          >>> Direct(257.4, 0.0).is_downward
          False
          >>> Direct(90, -45).is_downward
          False
        """

        return self.pl.is_downward

    def upward(self) -> 'Direct':
        """
        Return upward-point geological vector.

        Examples:
          >>> Direct(90, -45).upward().is_sub_parallel(Direct(90.0, -45.0))
          True
          >>> Direct(180, 45).upward().is_sub_parallel(Direct(0.0, -45.0))
          True
          >>> Direct(0, 0).upward().is_sub_parallel(Direct(0.0, 0.0))
          True
          >>> Direct(0, 90).upward().is_sub_parallel(Direct(180.0, -90.0))
          True
          >>> Direct(90, -45).upward().is_sub_parallel(Direct(90.0, -35.0))
          False
          >>> Direct(180, 45).upward().is_sub_parallel(Direct(10.0, -45.0))
          False
          >>> Direct(0, 0).upward().is_sub_parallel(Direct(170.0, 0.0))
          False
          >>> Direct(0, 90).upward().is_sub_parallel(Direct(180.0, -80.0))
          False
        """

        if not self.is_downward:
            return self.copy()
        else:
            return self.opposite()

    def downward(self) -> 'Direct':
        """
        Return downward-pointing geological vector.

        Examples:
          >>> Direct(90, -45).downward().is_sub_parallel(Direct(270.0, 45.0))
          True
          >>> Direct(180, 45).downward().is_sub_parallel(Direct(180.0, 45.0))
          True
          >>> Direct(0, 0).downward().is_sub_parallel(Direct(180.0, 0.0))
          False
          >>> Direct(0, 90).downward().is_sub_parallel(Direct(0.0, 90.0))
          True
          >>> Direct(90, -45).downward().is_sub_parallel(Direct(270.0, 35.0))
          False
          >>> Direct(180, 45).downward().is_sub_parallel(Direct(170.0, 45.0))
          False
          >>> Direct(0, 0).downward().is_sub_parallel(Direct(180.0, 10.0))
          False
          >>> Direct(0, 90).downward().is_sub_parallel(Direct(0.0, 80.0))
          False
        """

        if not self.is_upward:
            return self.copy()
        else:
            return self.opposite()

    def is_abs_dip_within(self,
                          min_val: numbers.Real,
                          max_val: numbers.Real,
                          min_val_incl: bool = False,
                          max_value_incl: bool = True
                          ) -> bool:
        """
        Check whether the absolute value of the dip angle of an Direct instance is intersect a given range
        (default: minimum value is not included, maximum value is included).

        :param min_val: the minimum dip angle, positive, domain: 0-90°.
        :param max_val: the maximum dip angle, positive, domain: 0-90°.
        :param min_val_incl: is minimum value included, boolean.
        :param max_value_incl: is maximum value included, boolean.
        :return: Boolean

        Examples:
          >>> Direct(90, -45).is_abs_dip_within(30, 60)
          True
          >>> Direct(120, 0).is_abs_dip_within(0, 60)
          False
          >>> Direct(120, 0).is_abs_dip_within(0, 60, min_val_incl=True)
          True
          >>> Direct(120, 60).is_abs_dip_within(0, 60)
          True
        """

        abs_dip = abs(self.pl.d)

        if abs_dip < min_val or abs_dip > max_val:
            return False
        elif abs_dip == min_val:
            if min_val_incl:
                return True
            else:
                return False
        elif abs_dip == max_val:
            if max_value_incl:
                return True
            else:
                return False
        else:
            return True

    def is_sub_horizontal(self,
                          max_dip_angle=DIP_ANGLE_THRESHOLD
                          ) -> bool:
        """
        Check whether the instance is almost horizontal.

        Examples:
          >>> Direct(10, 15).is_sub_horizontal()
          False
          >>> Direct(257, 2).is_sub_horizontal()
          True
          >>> Direct(90, -5).is_sub_horizontal()
          False
        """

        return abs(self.pl.d) < max_dip_angle

    def is_sub_vertical(self,
                        min_dip_angle=90.0 - DIP_ANGLE_THRESHOLD
                        ) -> bool:
        """
        Check whether the instance is almost vertical.

        Examples:
          >>> Direct(10, 15).is_sub_vertical()
          False
          >>> Direct(257, 89).is_sub_vertical()
          True
        """

        return abs(self.pl.d) > min_dip_angle

    def angle_as_degrees(self,
                         another: 'Direct'
                         ) -> numbers.Real:
        """
        Calculate angle (in degrees) between the two Direct instances.
        Range is 0°-180°.

        Examples:
          >>> areClose(Direct(0, 90).angle_as_degrees(Direct(90, 0)), 90)
          True
          >>> areClose(Direct(0, 0).angle_as_degrees(Direct(270, 0)), 90)
          True
          >>> areClose(Direct(0, 0).angle_as_degrees(Direct(0, 0)), 0)
          True
          >>> areClose(Direct(0, 0).angle_as_degrees(Direct(180, 0)), 180)
          True
          >>> areClose(Direct(90, 0).angle_as_degrees(Direct(270, 0)), 180)
          True
        """

        angle_vers = self.as_versor().angle_as_degrees(another.as_versor())

        return angle_vers

    def is_sub_parallel(self,
                        another,
                        angle_tolerance=VECTOR_ANGLE_THRESHOLD
        ):
        """
        Check that two Direct instances are sub-parallel,

        :param another: an Direct instance
        :param angle_tolerance: the maximum allowed divergence angle (in degrees)
        :return: Boolean

        Examples:
          >>> Direct(0, 90).is_sub_parallel(Direct(90, 0))
          False
          >>> Direct(0, 0).is_sub_parallel(Direct(0, 1e-6))
          True
          >>> Direct(0, 90).is_sub_parallel(Direct(180, 0))
          False
          >>> Direct(0, 90).is_sub_parallel(Direct(0, -90))
          False
        """

        fst_gvect = self

        snd_geoelem = another

        angle = fst_gvect.angle_as_degrees(snd_geoelem)

        if isinstance(another, Plane):
            return angle > (90.0 - angle_tolerance)
        else:
            return angle <= angle_tolerance

    def is_sub_antiparallel(self,
                            another,
                            angle_tolerance=VECTOR_ANGLE_THRESHOLD
                            ) -> bool:
        """
        Check that two Vect instances are almost anti-parallel,

        :param another: a Vect instance
        :param angle_tolerance: the maximum allowed divergence angle (in degrees)
        :return: Boolean

        Examples:
          >>> Direct(0, 90).is_sub_antiparallel(Direct(90, -89.5))
          True
          >>> Direct(0, 0).is_sub_antiparallel(Direct(180, 1e-6))
          True
          >>> Direct(90, 45).is_sub_antiparallel(Direct(270, -45.5))
          True
          >>> Direct(45, 90).is_sub_antiparallel(Direct(0, -90))
          True
          >>> Direct(45, 72).is_sub_antiparallel(Direct(140, -38))
          False
        """

        return self.angle_as_degrees(another) > (180.0 - angle_tolerance)

    def is_sub_orthogonal(self,
                          another: 'Direct',
                          angle_tolerance: numbers.Real = VECTOR_ANGLE_THRESHOLD
                          ) -> bool:
        """
        Check that two Direct instance are sub-orthogonal

        :param another: a Direct instance
        :param angle_tolerance: the maximum allowed divergence angle (in degrees) from orthogonality
        :return: Boolean

         Examples:
          >>> Direct(0, 90).is_sub_orthogonal(Direct(90, 0))
          True
          >>> Direct(0, 0).is_sub_orthogonal(Direct(0, 1.e-6))
          False
          >>> Direct(0, 0).is_sub_orthogonal(Direct(180, 0))
          False
          >>> Direct(90, 0).is_sub_orthogonal(Direct(270, 89.5))
          True
          >>> Direct(0, 90).is_sub_orthogonal(Direct(0, 0.5))
          True
        """

        return 90.0 - angle_tolerance <= self.angle_as_degrees(another) <= 90.0 + angle_tolerance

    def normal_versor(self,
                      another: 'Direct'
                      ) -> Optional[Vect3D]:
        """
        Calculate the versor (Vect) defined by the vector product of two Direct instances.

        Examples:
          >>> Direct(0, 0).normal_versor(Direct(90, 0))
          Vect3D(0.0000, 0.0000, -1.0000)
          >>> Direct(45, 0).normal_versor(Direct(310, 0))
          Vect3D(0.0000, 0.0000, 1.0000)
          >>> Direct(0, 0).normal_versor(Direct(90, 90))
          Vect3D(-1.0000, 0.0000, -0.0000)
          >>> Direct(315, 45).normal_versor(Direct(315, 44.5)) is None
          True
        """

        if self.is_sub_parallel(another):
            return None
        else:
            return self.as_versor().cross_product(another.as_versor()).versor()

    def normal_plane(self) -> 'Plane':
        """
        Return the geological plane that is normal to the direction.

        Examples:
          >>> Direct(0, 45).normal_plane()
          Plane(180.00, +45.00)
          >>> Direct(0, -45).normal_plane()
          Plane(000.00, +45.00)
          >>> Direct(0, 90).normal_plane()
          Plane(180.00, +00.00)
        """

        down_orien = self.downward()
        dipdir = (down_orien.az.d + 180.0) % 360.0
        dipangle = 90.0 - down_orien.pl.d

        return Plane(dipdir, dipangle)

    def common_plane(self,
                     another
    ):
        """
        Calculate Plane instance defined by the two Vect instances.

        Examples:
          >>> Direct(0, 0).common_plane(Direct(90, 0)).is_sub_parallel(Plane(180.0, 0.0))
          True
          >>> Direct(0, 0).common_plane(Direct(90, 90)).is_sub_parallel(Plane(90.0, 90.0))
          True
          >>> Direct(45, 0).common_plane(Direct(135, 45)).is_sub_parallel(Plane(135.0, 45.0))
          True
          >>> Direct(315, 45).common_plane(Direct(135, 45)).is_sub_parallel(Plane(225.0, 90.0))
          True
          >>> Direct(0, 0).common_plane(Direct(90, 0)).is_sub_parallel(Plane(180.0, 10.0))
          False
          >>> Direct(0, 0).common_plane(Direct(90, 90)).is_sub_parallel(Plane(90.0, 80.0))
          False
          >>> Direct(45, 0).common_plane(Direct(135, 45)).is_sub_parallel(Plane(125.0, 45.0))
          False
          >>> Direct(315, 45).common_plane(Direct(135, 45)).is_sub_parallel(Plane(225.0, 80.0))
          False
          >>> Direct(315, 45).common_plane(Direct(315, 44.5)) is None
          True
        """

        normal_versor = self.normal_versor(another)
        if normal_versor is None:
            return None
        else:
            return Direct.fromVect(normal_versor).normal_plane()

    def normal_direction(self,
                         another: 'Direct'
    ) -> Optional['Direct']:
        """
        Calculate the instance that is normal to the two provided sources.
        Angle between sources must be larger than MIN_ANGLE_DEGR_DISORIENTATION,

        Example:
          >>> Direct(0, 0).normal_direction(Direct(0.5, 0)) is None
          True
          >>> Direct(0, 0).normal_direction(Direct(179.5, 0)) is None
          True
          >>> Direct(0, 0).normal_direction(Direct(5.1, 0))
          Direct(az: 0.00°, pl: 90.00°)
          >>> Direct(90, 45).normal_direction(Direct(90, 0))
          Direct(az: 180.00°, pl: -0.00°)
        """

        if self.is_sub_antiparallel(another):
            return None
        elif self.is_sub_parallel(another):
            return None
        else:
            return self.__class__.fromVect(self.normal_versor(another))


class Axis(Direct):
    """
    Polar Axis. Inherits from Orientation
    """

    def __init__(self,
                 az: numbers.Real,
                 pl: numbers.Real
                 ):

        super().__init__(az, pl)

    def __repr__(self):

        return "Axis(az: {:.2f}°, pl: {:.2f}°)".format(*self.d)

    @classmethod
    def from_direction(cls,
                       direction: Direct
    ) -> 'Axis':
        """
        Create Axis instance from a direction.

        Example:
          >>> Axis.from_direction(Direct(220, 32))
          Axis(az: 220.00°, pl: 32.00°)
        """

        return Axis(
            az=direction.az.d,
            pl=direction.pl.d
        )

    def as_direction(self) -> Direct:
        """
        Create Direct instance with the same attitude as the self instance.

        Example:
          >>> Axis(220, 32).as_direction()
          Direct(az: 220.00°, pl: 32.00°)
        """

        return Direct(
            az=self.az.d,
            pl=self.pl.d
        )

    def normal_axis(self,
                    another: 'Axis'
    ) -> Optional['Axis']:
        """
        Calculate the Axis instance that is perpendicular to the two provided.
        The two source Axis must not be subparallel (threshold is MIN_ANGLE_DEGR_DISORIENTATION),
        otherwise a SubparallelLineationException will be raised.

        Example:
          >>> Axis(0, 0).normal_axis(Axis(0.5, 0)) is None
          True
          >>> Axis(0, 0).normal_axis(Axis(180, 0)) is None
          True
          >>> Axis(90, 0).normal_axis(Axis(180, 0))
          Axis(az: 0.00°, pl: 90.00°)
          >>> Axis(90, 45).normal_axis(Axis(180, 0))
          Axis(az: 270.00°, pl: 45.00°)
          >>> Axis(270, 45).normal_axis(Axis(180, 90)).is_sub_parallel(Axis(180, 0))
          True
        """

        norm_orien = self.normal_direction(another)
        if norm_orien is None:
            return None
        else:
            return Axis.from_direction(norm_orien)

    def angle_as_degrees(self,
                         another
    ):
        """
        Calculate angle (in degrees) between the two Axis instances.
        Range is 0°-90°.

        Examples:
          >>> areClose(Axis(0, 90).angle_as_degrees(Axis(90, 0)), 90)
          True
          >>> areClose(Axis(0, 0).angle_as_degrees(Axis(270, 0)), 90)
          True
          >>> areClose(Axis(0, 0).angle_as_degrees(Axis(0, 0)), 0)
          True
          >>> areClose(Axis(0, 0).angle_as_degrees(Axis(180, 0)), 0)
          True
          >>> areClose(Axis(0, 0).angle_as_degrees(Axis(179, 0)), 1)
          True
          >>> areClose(Axis(0, -90).angle_as_degrees(Axis(0, 90)), 0)
          True
          >>> areClose(Axis(90, 0).angle_as_degrees(Axis(315, 0)), 45)
          True
        """

        angle_vers = self.as_versor().angle_as_degrees(another.as_versor())

        return min(angle_vers, 180.0 - angle_vers)


class RotationAxis(object):
    """
    Rotation axis, expressed by an Orientation and a rotation angle.
    """

    def __init__(self,
                 trend: numbers.Real,
                 plunge: numbers.Real,
                 rot_ang: numbers.Real
    ):
        """
        Constructor.

        :param trend: Float/Integer
        :param plunge: Float/Integer
        :param rot_ang: Float/Integer

        Example:
        >> RotationAxis(0, 90, 120)
        RotationAxis(0.0000, 90.0000, 120.0000)
        """

        self.dr = Direct(trend, plunge)
        self.a = float(rot_ang)

    @classmethod
    def fromQuater(cls,
                   quat: Quaternion
    ):
        """
        Calculates the Rotation Axis expressed by a quaternion.
        The resulting rotation asVect is set to point downward.
        Examples are taken from Kuipers, 2002, chp. 5.

        :return: RotationAxis instance.

        Examples:
          >>> RotationAxis.fromQuater(Quaternion(0.5, 0.5, 0.5, 0.5))
          RotationAxis(45.0000, -35.2644, 120.0000)
          >>> RotationAxis.fromQuater(Quaternion(sqrt(2)/2, 0.0, 0.0, sqrt(2)/2))
          RotationAxis(0.0000, -90.0000, 90.0000)
          >>> RotationAxis.fromQuater(Quaternion(sqrt(2)/2, sqrt(2)/2, 0.0, 0.0))
          RotationAxis(90.0000, -0.0000, 90.0000)
        """

        if abs(quat) < QUAT_MAGN_THRESH:

            rot_ang = 0.0
            rot_direct = Direct(0.0, 0.0)

        elif areClose(quat.scalar, 1):

            rot_ang = 0.0
            rot_direct = Direct(0.0, 0.0)

        else:

            unit_quat = quat.normalize()
            rot_ang = unit_quat.rotAngle()
            rot_direct = Direct.fromVect(unit_quat.vector())

        return RotationAxis(*rot_direct.d, rot_ang)

    @classmethod
    def fromDirect(cls,
                   direct: Direct,
                   angle: numbers.Real
    ):
        """
        Class constructor from a Direct instance and an angle value.

        :param direct: a Direct instance
        :param angle: numbers.Real.
        :return: RotationAxis instance

        Example:
          >>> RotationAxis.fromDirect(Direct(320, 12), 30)
          RotationAxis(320.0000, 12.0000, 30.0000)
          >>> RotationAxis.fromDirect(Direct(315.0, -0.0), 10)
          RotationAxis(315.0000, -0.0000, 10.0000)
        """

        return RotationAxis(*direct.d, angle)

    @classmethod
    def fromVect(cls,
                 vector: Vect3D,
                 angle: numbers.Real
                 ):
        """
        Class constructor from a Vect instance and an angle value.

        :param vector: a Vect instance
        :param angle: float value
        :return: RotationAxis instance

        Example:
          >>> RotationAxis.fromVect(Vect3D(0, 1, 0), 30)
          RotationAxis(0.0000, -0.0000, 30.0000)
          >>> RotationAxis.fromVect(Vect3D(1, 0, 0), 30)
          RotationAxis(90.0000, -0.0000, 30.0000)
          >>> RotationAxis.fromVect(Vect3D(0, 0, -1), 30)
          RotationAxis(0.0000, 90.0000, 30.0000)
        """

        direct = Direct.fromVect(vector)

        return RotationAxis.fromDirect(direct, angle)

    def __repr__(self):

        return "RotationAxis({:.4f}, {:.4f}, {:.4f})".format(*self.dr.d, self.a)

    @property
    def rotAngle(self) -> float:
        """
        Returns the rotation angle of the rotation axis.

        :return: rotation angle (Float)

        Example:
          >>> RotationAxis(10, 15, 230).rotAngle
          230.0
        """

        return self.a

    @property
    def rotDirect(self) -> Direct:
        """
        Returns the rotation axis, expressed as a Direct.

        :return: Direct instance

        Example:
          >>> RotationAxis(320, 40, 15).rotDirect
          Direct(az: 320.00°, pl: 40.00°)
          >>> RotationAxis(135, 0, -10).rotDirect
          Direct(az: 135.00°, pl: 0.00°)
          >>> RotationAxis(45, 10, 10).rotDirect
          Direct(az: 45.00°, pl: 10.00°)
        """

        return self.dr

    @property
    def versor(self) -> Vect3D:
        """
        Return the versor equivalent to the Rotation geological asVect.

        :return: Vect
        """

        return self.dr.as_versor()

    def specular(self):
        """
        Derives the rotation axis with opposite asVect direction
        and rotation angle that is the complement to 360°.
        The resultant rotation is equivalent to the original one.

        :return: RotationAxis instance.

        Example
          >>> RotationAxis(90, 45, 320).specular()
          RotationAxis(270.0000, -45.0000, 40.0000)
          >>> RotationAxis(135, 0, -10).specular()
          RotationAxis(315.0000, -0.0000, 10.0000)
          >>> RotationAxis(45, 10, 10).specular()
          RotationAxis(225.0000, -10.0000, 350.0000)
        """

        gvect_opp = self.rotDirect.opposite()
        opposite_angle = (360.0 - self.rotAngle) % 360.0

        return RotationAxis.fromDirect(gvect_opp, opposite_angle)

    def compl180(self):
        """
        Creates a new rotation axis that is the complement to 180 of the original one.

        :return: RotationAxis instance.

        Example:
          >>> RotationAxis(90, 45, 120).compl180()
          RotationAxis(90.0000, 45.0000, 300.0000)
          >>> RotationAxis(117, 34, 18).compl180()
          RotationAxis(117.0000, 34.0000, 198.0000)
          >>> RotationAxis(117, 34, -18).compl180()
          RotationAxis(117.0000, 34.0000, 162.0000)
        """

        rot_ang = - (180.0 - self.rotAngle) % 360.0
        return RotationAxis.fromDirect(self.dr, rot_ang)

    def strictlyEquival(self,
                        another,
                        angle_tolerance: numbers.Real=VECTOR_ANGLE_THRESHOLD
    ) -> bool:
        """
        Checks if two RotationAxis are almost equal, based on a strict checking
        of the Direct component and of the rotation angle.

        :param another: another RotationAxis instance, to be compared with
        :type another: RotationAxis
        :parameter angle_tolerance: the tolerance as the angle (in degrees)
        :type angle_tolerance: numbers.Real.
        :return: the equivalence (true/false) between the two compared RotationAxis
        :rtype: bool

        Examples:
          >>> ra_1 = RotationAxis(180, 10, 10)
          >>> ra_2 = RotationAxis(180, 10, 10.5)
          >>> ra_1.strictlyEquival(ra_2)
          True
          >>> ra_3 = RotationAxis(180.2, 10, 10.4)
          >>> ra_1.strictlyEquival(ra_3)
          True
          >>> ra_4 = RotationAxis(184.9, 10, 10.4)
          >>> ra_1.strictlyEquival(ra_4)
          False
        """

        if not self.dr.is_sub_parallel(another.dr, angle_tolerance):
            return False

        if not areClose(self.a, another.a, atol=1.0):
            return False

        return True

    def toRotQuater(self) -> Quaternion:
        """
        Converts the rotation axis to the equivalent rotation quaternion.

        :return: the rotation quaternion.
        :rtype: Quaternion
        """

        rotation_angle_rad = radians(self.a)
        rotation_vector = self.dr.as_versor()

        w = cos(rotation_angle_rad / 2.0)
        x, y, z = rotation_vector.scale(sin(rotation_angle_rad / 2.0)).toXYZ()

        return Quaternion(w, x, y, z).normalize()

    def toRotMatrix(self):
        """
        Derives the rotation matrix from the RotationAxis instance.

        :return: 3x3 numpy array
        """

        rotation_versor = self.versor
        phi = radians(self.a)

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

    def toMinRotAxis(self):
        """
        Calculates the minimum rotation axis from the given quaternion.

        :return: RotationAxis instance.
        """

        return self if abs(self.rotAngle) <= 180.0 else self.specular()

    @classmethod
    def randomNaive(cls):
        """
        Naive method for creating a random RotationAxis instance.
        :return: random rotation axis (not uniformly distributed in the space)
        :rtype: RotationAxis
        """

        random_trend = random.uniform(0, 360)
        random_dip = random.uniform(-90, 90)
        random_rotation = random.uniform(0, 360)

        return cls(
            trend=random_trend,
            plunge=random_dip,
            rot_ang=random_rotation
        )


def sortRotations(
        rotation_axes: List[RotationAxis]
) -> List[RotationAxis]:
    """
    Sorts a list or rotation axes, based on the rotation angle (absolute value),
    in an increasing order.

    :param rotation_axes: o list of RotationAxis objects.
    :return: the sorted list of RotationAxis

    Example:
      >>> rots = [RotationAxis(110, 14, -23), RotationAxis(42, 13, 17), RotationAxis(149, 87, 13)]
      >>> sortRotations(rots)
      [RotationAxis(149.0000, 87.0000, 13.0000), RotationAxis(42.0000, 13.0000, 17.0000), RotationAxis(110.0000, 14.0000, -23.0000)]
    """

    return sorted(rotation_axes, key=lambda rot_ax: abs(rot_ax.rotAngle))


def rotVectByAxis(
    v: Vect3D,
    rot_axis: RotationAxis
) -> Vect3D:
    """
    Rotates a vector.

    Implementation as in:
    https://en.wikipedia.org/wiki/Conversion_between_quaternions_and_Euler_angles
    Faster formula:
    t = 2 q x v
    v' = v + q0 t + q x t
    cited as:
    Janota, A; Šimák, V; Nemec, D; Hrbček, J (2015).
    "Improving the Precision and Speed of Euler Angles Computation from Low-Cost Rotation Sensor Data".
    Sensors. 15 (3): 7016–7039. doi:10.3390/s150307016. PMC 4435132. PMID 25806874.

    :param v: the vector to rotate
    :type v: Vect3D
    :param rot_axis: the rotation axis
    :type rot_axis: RotationAxis
    :return: the rotated vector
    :rtype: Vect3D

    Examples:
      >>> v = Vect3D(1,0,1)
      >>> rotation = RotationAxis(0, -90, 90)
      >>> rotVectByAxis(v, rotation)
      Vect3D(0.0000, 1.0000, 1.0000)
      >>> rotation = RotationAxis(0, 90, 90)
      >>> rotVectByAxis(v, rotation)
      Vect3D(0.0000, -1.0000, 1.0000)
      >>> rotation = RotationAxis(0, -90, 180)
      >>> rotVectByAxis(v, rotation)
      Vect3D(-1.0000, 0.0000, 1.0000)
      >>> rotation = RotationAxis(0, -90, 270)
      >>> rotVectByAxis(v, rotation)
      Vect3D(-0.0000, -1.0000, 1.0000)
      >>> rotation = RotationAxis(90, 0, 90)
      >>> rotVectByAxis(v, rotation)
      Vect3D(1.0000, -1.0000, 0.0000)
      >>> rotation = RotationAxis(90, 0, 180)
      >>> rotVectByAxis(v, rotation)
      Vect3D(1.0000, 0.0000, -1.0000)
      >>> rotation = RotationAxis(90, 0, 270)
      >>> rotVectByAxis(v, rotation)
      Vect3D(1.0000, 1.0000, -0.0000)
      >>> rotation = RotationAxis(90, 0, 360)
      >>> rotVectByAxis(v, rotation)
      Vect3D(1.0000, 0.0000, 1.0000)
      >>> rotation = RotationAxis(0, -90, 90)
      >>> v = Vect3D(0,0,3)
      >>> rotVectByAxis(v, rotation)
      Vect3D(0.0000, 0.0000, 3.0000)
      >>> rotation = RotationAxis(90, -45, 180)
      >>> rotVectByAxis(v, rotation)
      Vect3D(3.0000, -0.0000, -0.0000)
      >>> v = Vect3D(0,0,3)
      >>> rotVectByAxis(v, rotation)
      Vect3D(3.0000, -0.0000, -0.0000)
    """

    rot_quat = rot_axis.toRotQuater()
    q = rot_quat.vector()

    t = q.scale(2).cross_product(v)
    rot_v = v + t.scale(rot_quat.scalar) + q.cross_product(t)

    return rot_v


def rotVectByQuater(
        quat: Quaternion,
        vect: Vect3D
) -> Vect3D:
    """
    Calculates a rotated solution of a Vect instance given a normalized quaternion.
    Original formula in Ref. [1].
    Eq.6: R(qv) = q qv q(-1)

    :param quat: a Quaternion instance
    :param vect: a Vect instance
    :return: a rotated Vect instance

    Example:
      >>> q = Quaternion.i()  # rotation of 180° around the x axis
      >>> rotVectByQuater(q, Vect3D(0, 1, 0))
      Vect3D(0.0000, -1.0000, 0.0000)
      >>> rotVectByQuater(q, Vect3D(0, 1, 1))
      Vect3D(0.0000, -1.0000, -1.0000)
      >>> q = Quaternion.k()  # rotation of 180° around the z axis
      >>> rotVectByQuater(q, Vect3D(0, 1, 1))
      Vect3D(0.0000, -1.0000, 1.0000)
      >>> q = Quaternion.j()  # rotation of 180° around the y axis
      >>> rotVectByQuater(q, Vect3D(1, 0, 1))
      Vect3D(-1.0000, 0.0000, -1.0000)
    """

    q = quat.normalize()
    qv = Quaternion.fromVect(vect)

    rotated_v = q * (qv * q.inverse)

    return rotated_v.vector()


class Plane:
    """
    Geological plane.
    Defined by dip direction and dip angle (both in degrees):
     - dip direction: [0.0, 360.0[ clockwise, from 0 (North);
     - dip angle: [0, 90.0]: downward-pointing.
    """

    def __init__(self,
                 azim: numbers.Real,
                 dip_ang: numbers.Real,
                 is_rhr_strike: bool = False
    ):
        """
        Geological plane constructor.

        :param  azim:  azimuth of the plane (RHR strike or dip direction).
        :type  azim:  number or string convertible to float.
        :param  dip_ang:  Dip angle of the plane (0-90°).
        :type  dip_ang:  number or string convertible to float.
        :param is_rhr_strike: if the source azimuth is RHR strike (default is False, i.e. it is dip direction)
        :return: the instantiated geological plane.
        :rtype: Plane.

        Example:
          >>> Plane(0, 90)
          Plane(000.00, +90.00)
          >>> Plane(0, 90, is_rhr_strike=True)
          Plane(090.00, +90.00)
          >>> Plane(0, 90, True)
          Plane(090.00, +90.00)
          >>> Plane(0, 900)
          Traceback (most recent call last):
          ...
          Exception: Dip angle must be between 0° and 90°
        """

        def rhrstrk2dd(rhr_strk):
            """Converts RHR strike value to dip direction value.

            Example:
                >>> rhrstrk2dd(285.5)
                15.5
            """

            return (rhr_strk + 90.0) % 360.0

        if not isinstance(azim, numbers.Real):
            raise Exception("Source azimuth must be number")
        if not isinstance(dip_ang, numbers.Real):
            raise Exception("Source dip angle must be number")
        if not isinstance(is_rhr_strike, bool):
            raise Exception("Source azimuth type must be boolean")

        if not (0.0 <= dip_ang <= 90.0):
            raise Exception("Dip angle must be between 0° and 90°")

        if is_rhr_strike:
            self._dipdir = rhrstrk2dd(azim)
        else:
            self._dipdir = azim % 360.0
        self._dipangle = float(dip_ang)

    @property
    def dd(self):
        """
        Return the dip direction of the geological plane.

        Example:
          >>> Plane(34.2, 89.7).dd
          34.2
        """

        return self._dipdir

    @property
    def da(self):
        """
        Return the dip angle of the geological plane.

        Example:
          >>> Plane(183, 77).da
          77.0

        """

        return self._dipangle

    @property
    def dda(self):
        """
        Return a tuple storing the dip direction and dip angle values of a geological plane.

        Example:
          >>> gp = Plane(89.4, 17.2)
          >>> gp.dda
          (89.4, 17.2)
        """

        return self.dd, self.da

    @property
    def rhr_strike(self):
        """
        Return the strike according to the right-hand-rule.

        Examples:
          >>> Plane(90, 45).rhr_strike
          0.0
          >>> Plane(45, 89).rhr_strike
          315.0
          >>> Plane(275, 38).rhr_strike
          185.0
          >>> Plane(0, 38).rhr_strike
          270.0
        """

        return (self.dd - 90.0) % 360.0

    @property
    def srda(self):
        """
        Return a tuple storing the right-hand-rule strike and dip angle values of a geological plane.

        Example:
          >>> Plane(100, 17.2).srda
          (10.0, 17.2)
          >>> Plane(10, 87).srda
          (280.0, 87.0)
        """

        return self.rhr_strike, self.da

    @property
    def lhrStrike(self):
        """
        Return the strike according to the left-hand-rule.

        Examples:
          >>> Plane(90, 45).lhrStrike
          180.0
          >>> Plane(45, 89).lhrStrike
          135.0
          >>> Plane(275, 38).lhrStrike
          5.0
          >>> Plane(0, 38).lhrStrike
          90.0
        """

        return (self.dd + 90.0) % 360.0

    @property
    def slda(self):
        """
        Return a tuple storing the left-hand-rule strike and dip angle values of a geological plane.

        Example:
          >>> Plane(100, 17.2).slda
          (190.0, 17.2)
          >>> Plane(10, 87).slda
          (100.0, 87.0)
        """

        return self.lhrStrike, self.da

    def __repr__(self):

        return "Plane({:06.2f}, {:+06.2f})".format(*self.dda)

    def rhr_strike_direction(self):
        """
        Creates a direction instance that is parallel to the right-hand rule strike.

        Examples:
          >>> Plane(90, 45).rhr_strike_direction()
          Direct(az: 0.00°, pl: 0.00°)
          >>> Plane(45, 17).rhr_strike_direction()
          Direct(az: 315.00°, pl: 0.00°)
          >>> Plane(90, 0).rhr_strike_direction()
          Direct(az: 0.00°, pl: 0.00°)
        """

        return Direct(
            az=self.rhr_strike,
            pl=0.0)

    def lhrStrikeOrien(self):
        """
        Creates an Orientation instance that is parallel to the left-hand rule strike.

        :return: OrienM instance.

        Examples:
          >>> Plane(90, 45).lhrStrikeOrien()
          Direct(az: 180.00°, pl: 0.00°)
          >>> Plane(45, 17).lhrStrikeOrien()
          Direct(az: 135.00°, pl: 0.00°)
        """

        return Direct(
            az=self.lhrStrike,
            pl=0.0)

    def dipDirOrien(self):
        """
        Creates a OrienM instance that is parallel to the dip direction.

        :return: OrienM instance.

        Examples:
          >>> Plane(90, 45).dipDirOrien()
          Direct(az: 90.00°, pl: 45.00°)
          >>> Plane(45, 17).dipDirOrien()
          Direct(az: 45.00°, pl: 17.00°)
        """

        return Direct(
            az=self.dd,
            pl=self.da)

    def dipDirOppOrien(self):
        """
        Creates a OrienM instance that is anti-parallel to the dip direction.

        :return: OrienM instance.

        Examples:
          >>> Plane(90, 45).dipDirOppOrien()
          Direct(az: 270.00°, pl: -45.00°)
          >>> Plane(45, 17).dipDirOppOrien()
          Direct(az: 225.00°, pl: -17.00°)
        """

        return self.dipDirOrien().opposite()

    def mirrorVertPPlane(self):
        """
        Mirror a geological plane around a vertical plane
        creating a new one that has a dip direction opposite
        to the original one but with downward plunge.

        :return: geological plane
        :rtype: Plane

        Examples:
          >>> Plane(0, 45).mirrorVertPPlane()
          Plane(180.00, +45.00)
          >>> Plane(225, 80).mirrorVertPPlane()
          Plane(045.00, +80.00)
          >>> Plane(90, 90).mirrorVertPPlane()
          Plane(270.00, +90.00)
          >>> Plane(270, 0).mirrorVertPPlane()
          Plane(090.00, +00.00)
        """

        return Plane(
            azim=opposite_trend(self.dd),
            dip_ang=self.da)

    def normDirectFrwrd(self):
        """
        Return the direction normal to the geological plane,
        pointing in the same direction as the geological plane.

        Example:
            >>> Plane(90, 55).normDirectFrwrd()
            Direct(az: 90.00°, pl: -35.00°)
            >>> Plane(90, 90).normDirectFrwrd()
            Direct(az: 90.00°, pl: 0.00°)
            >>> Plane(90, 0).normDirectFrwrd()
            Direct(az: 90.00°, pl: -90.00°)
        """

        tr = self.dd % 360.0
        pl = self.da - 90.0

        return Direct(
            az=tr,
            pl=pl)

    def normDirectBckwrd(self):
        """
        Return the direction normal to the geological plane,
        pointing in the opposite direction to the geological plane.

        Example:
            >>> Plane(90, 55).normDirectBckwrd()
            Direct(az: 270.00°, pl: 35.00°)
            >>> Plane(90, 90).normDirectBckwrd()
            Direct(az: 270.00°, pl: -0.00°)
            >>> Plane(90, 0).normDirectBckwrd()
            Direct(az: 270.00°, pl: 90.00°)
        """

        return self.normDirectFrwrd().opposite()

    def normDirectDown(self):
        """
        Return the direction normal to the geological plane and
        pointing downward.

        Example:
            >>> Plane(90, 55).normDirectDown()
            Direct(az: 270.00°, pl: 35.00°)
            >>> Plane(90, 90).normDirectDown()
            Direct(az: 90.00°, pl: 0.00°)
            >>> Plane(90, 0).normDirectDown()
            Direct(az: 270.00°, pl: 90.00°)
        """

        return self.normDirectFrwrd().downward()

    def normDirectUp(self):
        """
        Return the direction normal to the polar plane,
        pointing upward.

        Example:
            >>> Plane(90, 55).normDirectUp()
            Direct(az: 90.00°, pl: -35.00°)
            >>> Plane(90, 90).normDirectUp()
            Direct(az: 90.00°, pl: 0.00°)
            >>> Plane(90, 0).normDirectUp()
            Direct(az: 90.00°, pl: -90.00°)
        """

        return self.normDirectFrwrd().upward()

    def normal_direction(self) -> 'Direct':
        """
        Wrapper to down_normal_gv.

        :return: downward-pointing Direct instance normal to the Plane self instance
        """

        return self.normDirectDown()

    def normal_axis(self):
        """
        Normal Axis.

        :return: Axis normal to the Plane self instance
        """

        return Axis.from_direction(self.normDirectDown())

    def angle_degr(self, another: 'Plane'):
        """
        Calculate angle (in degrees) between two geoplanes.
        Range is 0°-90°.

        Examples:
          >>> Plane(100.0, 50.0).angle_degr(Plane(100.0, 50.0))
          0.0
          >>> Plane(300.0, 10.0).angle_degr(Plane(300.0, 90.0))
          80.0
          >>> Plane(90.0, 90.0).angle_degr(Plane(270.0, 90.0))
          0.0
          >>> areClose(Plane(90.0, 90.0).angle_degr(Plane(130.0, 90.0)), 40)
          True
          >>> areClose(Plane(90, 70).angle_degr(Plane(270, 70)), 40)
          True
          >>> areClose(Plane(90.0, 10.0).angle_degr(Plane(270.0, 10.0)), 20.0)
          True
          >>> areClose(Plane(90.0, 10.0).angle_degr(Plane(270.0, 30.0)), 40.0)
          True
        """

        if not isinstance(another, Plane):
            raise Exception("Second instance for angle is of {} type".format(type(another)))

        gpl_axis = Axis.from_direction(self.normDirectFrwrd())
        an_axis = Axis.from_direction(another.normDirectFrwrd())

        return gpl_axis.angle_as_degrees(an_axis)

    def is_sub_parallel(self,
                        another,
                        angle_tolerance: numbers.Real = PLANE_ANGLE_THRESHOLD
    ):
        """
        Check that two GPlanes are sub-parallel

        :param another: a Plane instance
        :param angle_tolerance: the maximum allowed divergence angle (in degrees)
        :return: Boolean

         Examples:
          >>> Plane(0, 90).is_sub_parallel(Plane(270, 90))
          False
          >>> Plane(0, 90).is_sub_parallel(Plane(180, 90))
          True
          >>> Plane(0, 90).is_sub_parallel(Plane(0, 0))
          False
          >>> Plane(0, 0).is_sub_parallel(Plane(0, 1e-6))
          True
          >>> Plane(0, 0).is_sub_parallel(Plane(0, 1.1))
          False
        """

        return self.angle_degr(another) < angle_tolerance

    def contains(self,
        direct: 'Direct',
        angle_tolerance: numbers.Real = PLANE_ANGLE_THRESHOLD
    ) -> bool:
        """
        Check that a plane contains a direction instance.

        :param direct: a Direct instance
        :param angle_tolerance: the tolerance angle
        :return: True or False

        Examples:
          >>> Plane(90, 0).contains(Direct(60, 0))
          True
          >>> Plane(90, 0).contains(Axis(60, 0))
          True
          >>> Plane(90, 0).contains(Direct(60, 10))
          False
        """

        plane_norm = self.normal_axis()

        return direct.is_sub_orthogonal(plane_norm, angle_tolerance)

    def is_sub_orthogonal(self,
                          another,
                          angle_tolerance: numbers.Real = PLANE_ANGLE_THRESHOLD
    ):
        """
        Check that two GPlanes are sub-orthogonal.

        :param another: a Plane instance
        :param angle_tolerance: the maximum allowed divergence angle (in degrees)
        :return: Boolean

         Examples:
          >>> Plane(0, 90).is_sub_orthogonal(Plane(270, 90))
          True
          >>> Plane(0, 90).is_sub_orthogonal(Plane(180, 90))
          False
          >>> Plane(0, 90).is_sub_orthogonal(Plane(0, 0))
          True
          >>> Plane(0, 0).is_sub_orthogonal(Plane(0, 88))
          False
          >>> Plane(0, 0).is_sub_orthogonal(Plane(0, 45))
          False
        """

        fst_axis = Axis.from_direction(self.normal_direction())

        if isinstance(another, Plane):
            snd_gaxis = Axis.from_direction(another.normal_direction())
        else:
            raise Exception("Not accepted argument type for isSubOrthog method")

        angle = fst_axis.angle_as_degrees(snd_gaxis)

        if isinstance(another, Plane):
            return angle > 90.0 - angle_tolerance
        else:
            return angle < angle_tolerance

    def rakeToDirect(self,
                     rake
    ):
        """
        Calculate the Direct instance given a Plane instance and a rake value.
        The rake is defined according to the Aki and Richards, 1980 conventions:
        rake = 0° -> left-lateral
        rake = 90° -> reverse
        rake = +/- 180° -> right-lateral
        rake = -90° -> normal

        Examples:
          >>> Plane(180, 45).rakeToDirect(0.0)
          Direct(az: 90.00°, pl: -0.00°)
          >>> Plane(180, 45).rakeToDirect(90.0)
          Direct(az: 0.00°, pl: -45.00°)
          >>> Plane(180, 45).rakeToDirect(-90.0)
          Direct(az: 180.00°, pl: 45.00°)
          >>> Plane(180, 45).rakeToDirect(180.0).is_sub_parallel(Direct(270.00, 0.00))
          True
          >>> Plane(180, 45).rakeToDirect(-180.0)
          Direct(az: 270.00°, pl: 0.00°)
        """

        rk = radians(rake)
        strk = radians(self.rhr_strike)
        dip = radians(self.da)

        x = cos(rk) * sin(strk) - sin(rk) * cos(dip) * cos(strk)
        y = cos(rk) * cos(strk) + sin(rk) * cos(dip) * sin(strk)
        z = sin(rk) * sin(dip)

        return Direct.fromXYZ(x, y, z)

    def isVLowAngle(self,
                    dip_angle_threshold: numbers.Real = angle_gplane_thrshld
    ):
        """
        Checks if a geological plane is very low angle.

        :param dip_angle_threshold: the limit for the plane angle, in degrees
        :type dip_angle_threshold: numbers.Real.
        :return: bool flag indicating if it is very low angle

        Examples:
          >>> Plane(38.9, 1.2).isVLowAngle()
          True
          >>> Plane(38.9, 7.4).isVLowAngle()
          False
        """

        return self.da < dip_angle_threshold

    def isVHighAngle(self,
                     dip_angle_threshold: numbers.Real = angle_gplane_thrshld
    ):
        """
        Checks if a geological plane is very high angle.

        :param dip_angle_threshold: the limit for the plane angle, in degrees
        :type dip_angle_threshold: numbers.Real.
        :return: bool flag indicating if it is very high angle

        Examples:
          >>> Plane(38.9, 11.2).isVHighAngle()
          False
          >>> Plane(38.9, 88.4).isVHighAngle()
          True
        """

        return self.da > (90.0 - dip_angle_threshold)

    '''
    def to_cartesian_plane(self,
                           pt: Point
                           ) -> CPlane:
        """
        Given a Plane instance and a provided Point instance,
        calculate the corresponding Plane instance.

        Example:
          >>> Plane(0, 0).to_cartesian_plane(Point(0, 0, 0))
          CPlane(0.0000, 0.0000, 1.0000, -0.0000, -1)
          >>> Plane(90, 45).to_cartesian_plane(Point(0, 0, 0))
          CPlane(0.7071, 0.0000, 0.7071, -0.0000, -1)
          >>> Plane(0, 90).to_cartesian_plane(Point(0, 0, 0))
          CPlane(0.0000, 1.0000, -0.0000, -0.0000, -1)
        """

        normal_versor = self.normDirectFrwrd().as_versor()
        a, b, c = normal_versor.x, normal_versor.y, normal_versor.z
        d = - (a * pt.x + b * pt.y + c * pt.z)
        return CPlane(a, b, c, d)
    '''

    def slope_x_dir(self) -> numbers.Real:
        """
        Calculate the slope of a given plane along the x direction.
        The plane orientation  is expressed following the geological convention.

        :return: the slope along the x direction
        :rtype: numbers.Real.

        Example:
        """
        return - sin(radians(self.dd)) * tan(radians(self.da))

    def slope_y_dir(self) -> numbers.Real:
        """
        Calculate the slope of a given plane along the y direction.
        The plane orientation  is expressed following the geological convention.

        :return: the slope along the y direction
        :rtype: numbers.Real.

        Example:
        """
        return - cos(radians(self.dd)) * tan(radians(self.da))


if __name__ == "__main__":

    import doctest
    doctest.testmod()

