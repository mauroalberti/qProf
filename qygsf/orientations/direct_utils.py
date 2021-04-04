
from typing import Optional

import numbers

from math import atan2, degrees


def mod360(
        val: numbers.Real
) -> numbers.Real:
    """
    Return the module 360 of the originila value.

    :param val: value.
    :type val: numbers.Real.
    :return: value divided by mod 360.
    :rtype: numbers.Real.

    Examples:
    >>> mod360(10)
    10.0
    >>> mod360(360)
    0.0
    >>> mod360(-50)
    310.0
    >>> mod360(400)
    40.0

    """

    return val % 360.0


def opposite_trend(
        tr: numbers.Real
) -> numbers.Real:
    """
    Calculate the trend opposite to the original one.

    :return: the opposite trend.

    Examples:
    >>> opposite_trend(0)
    180.0
    >>> opposite_trend(45)
    225.0
    >>> opposite_trend(90)
    270.0
    >>> opposite_trend(180)
    0.0
    >>> opposite_trend(270)
    90.0
    """

    return mod360(tr + 180.0)


def angle_east_anticlock(
        x: numbers.Real,
        y: numbers.Real
) -> Optional[numbers.Real]:
    """
    Calculates the angle (in uDegrees, positive anti-clockwise) in the 2D plane
    given x (East) and y (North) Cartesian xyz.

    :param x: x component
    :param y: y component
    :return: angle (in decimal degrees) or None

    Examples:
      >>> angle_east_anticlock(0, 0) is None
      True
      >>> angle_east_anticlock(1, 0)
      0.0
      >>> angle_east_anticlock(1, 1)
      45.0
      >>> angle_east_anticlock(0, 5)
      90.0
      >>> angle_east_anticlock(-1, 0)
      180.0
      >>> angle_east_anticlock(0, -4)
      270.0
      >>> angle_east_anticlock(1, -1)
      315.0
    """
    if x == 0.0 and y == 0.0:
        return None
    else:
        return mod360(degrees(atan2(y, x)))


def angle_east_clock(
        x: numbers.Real,
        y: numbers.Real
) -> Optional[numbers.Real]:
    """
    Calculate the angle from East (in uDegrees, positive clockwise) in the 2D plane
    given x (East) and y (North) Cartesian xyz.

    :param x: x component
    :param y: y component
    :return: angle (in decimal degrees) or None

    Examples:
      >>> angle_east_clock(0, 0) is None
      True
      >>> angle_east_clock(1, 0)
      0.0
      >>> angle_east_clock(1, -1)
      45.0
      >>> angle_east_clock(0, -4)
      90.0
      >>> angle_east_clock(-1, 0)
      180.0
      >>> angle_east_clock(0, 5)
      270.0
      >>> angle_east_clock(1, 1)
      315.0
    """

    ang = angle_east_anticlock(x, y)
    if ang is None:
        return None
    else:
        return mod360(360.0 - ang)


def angle_north_clock(
        x: numbers.Real,
        y: numbers.Real
) -> Optional[numbers.Real]:
    """
    Calculate the angle from North (in degrees, positive clockwise) in the 2D plane
    given x (East) and y (North) Cartesian xyz.

    Examples:
      >>> angle_north_clock(0, 0) is None
      True
      >>> angle_north_clock(0, 5)
      0.0
      >>> angle_north_clock(1, 1)
      45.0
      >>> angle_north_clock(1, 0)
      90.0
      >>> angle_north_clock(1, -1)
      135.0
      >>> angle_north_clock(0, -1)
      180.0
      >>> angle_north_clock(-1, 0)
      270.0
      >>> angle_north_clock(0, -4)
      180.0
    """

    ang = angle_east_clock(x, y)
    if ang is None:
        return None
    else:
        return mod360(90.0 + ang)


def plng2colatTop(
        plunge: numbers.Real
) -> numbers.Real:
    """
    Calculates the colatitude angle from the top.

    :param plunge: an angle from -90° (upward-pointing) to 90° (downward-pointing)
    :type plunge: numbers.Real.
    :return: the colatitude angle
    :rtype: numbers.Real.

    Examples:
      >>> plng2colatTop(90)
      180.0
      >>> plng2colatTop(45)
      135.0
      >>> plng2colatTop(0)
      90.0
      >>> plng2colatTop(-45)
      45.0
      >>> plng2colatTop(-90)
      0.0
    """

    return 90.0 + plunge


def plng2colatBottom(
        plunge: numbers.Real
) -> numbers.Real:
    """
    Calculates the colatitude angle from the bottom.

    :param plunge: an angle from -90° (upward-pointing) to 90° (downward-pointing)
    :type plunge: numbers.Real.
    :return: the colatitude angle
    :rtype: numbers.Real.

    Examples:
      >>> plng2colatBottom(90)
      0.0
      >>> plng2colatBottom(45)
      45.0
      >>> plng2colatBottom(0)
      90.0
      >>> plng2colatBottom(-45)
      135.0
      >>> plng2colatBottom(-90)
      180.0
    """

    return 90.0 - plunge


def slope(
        h: numbers.Real,
        v: numbers.Real
) -> Optional[numbers.Real]:
    """
    Slope (in decimal degrees) given horizontal and vertical lengths
    both input are assumed positive.

    :param h: the horizontal distance
    :param v: the vertical offset
    :return: the slope in decimal degrees or None

    Examples:
      >>> slope(0, 0) is None
      True
      >>> slope(1, 1)
      45.0
      >>> slope(1, 0)
      0.0
    """

    if h == 0.0 and v == 0.0:
        return None
    elif h == 0.0:
        return 90.0
    else:
        return degrees(atan2(v, h))

