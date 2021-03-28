
from typing import Tuple

import math

import numbers


def areClose(
    a: numbers.Real,
    b: numbers.Real,
    rtol: numbers.Real = 1e-012,
    atol: numbers.Real = 1e-12,
    equal_nan: bool = False,
    equal_inf: bool = False
) -> bool:
    """
    Mimics math.isclose from Python 3.5 (see: https://src_docs.python.org/3.5/library/math.html)

    Example:
      >>> areClose(1.0, 1.0)
      True
      >>> areClose(1.0, 1.000000000000001)
      True
      >>> areClose(1.0, 1.0000000001)
      False
      >>> areClose(0.0, 0.0)
      True
      >>> areClose(0.0, 0.000000000000001)
      True
      >>> areClose(0.0, 0.0000000001)
      False
      >>> areClose(100000.0, 100000.0)
      True
      >>> areClose(100000.0, 100000.0000000001)
      True
      >>> areClose(float('nan'), float('nan'))
      False
      >>> areClose(float('nan'), 1000000)
      False
      >>> areClose(1.000000000001e300, 1.0e300)
      False
      >>> areClose(1.0000000000001e300, 1.0e300)
      True
      >>> areClose(float('nan'), float('nan'), equal_nan=True)
      True
      >>> areClose(float('inf'), float('inf'))
      False
      >>> areClose(float('inf'), 1.0e300)
      False
      >>> areClose(float('inf'), float('inf'), equal_inf=True)
      True
    """

    # nan cases
    if equal_nan and math.isnan(a) and math.isnan(b):
        return True
    elif math.isnan(a) or math.isnan(b):
        return False

    # inf cases
    if equal_inf and math.isinf(a) and a > 0 and math.isinf(b) and b > 0:
        return True
    elif equal_inf and math.isinf(a) and a < 0 and math.isinf(b) and b < 0:
        return True
    elif math.isinf(a) or math.isinf(b):
        return False

    # regular case
    return abs(a - b) <= max(rtol * max(abs(a), abs(b)), atol)


def apprFloat(
    val: numbers.Real,
    ndec: numbers.Integral = 1
) -> numbers.Real:
    """
    Rounds a numeric value to ndec.

    :param val: value to round
    :param ndec: number of decimals used
    :return: rounded float value

    Examples:
      >>> apprFloat(0.00001)
      0.0
      >>> apprFloat(1.425324e-7)
      0.0
    """

    rval = round(val, ndec)
    if rval == 0.0:
        rval = round(0.0, ndec)

    return rval


def apprFTuple(
    tup: Tuple[numbers.Real, ...],
    ndec=1
) -> Tuple[numbers.Real, ...]:
    """
    Rounds numeric values inside a tuple to ndec decimals

    :param tup: tuple of numbers.Real values
    :param ndec: number of decimals used
    :return: tuple with rounded numbers

    Examples:
      >>> apprFTuple((-2.4492935982947064e-16, 1.0))
      (0.0, 1.0)
      >>> apprFTuple((-1.0, -1.8369701987210297e-16))
      (-1.0, 0.0)
    """

    return tuple(map(lambda val: apprFloat(val, ndec), tup))


if __name__ == '__main__':

    import doctest
    doctest.testmod()
