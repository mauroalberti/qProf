
import math
import numbers
from typing import Tuple

import numpy as np

from .defaults import MIN_VECTOR_MAGNITUDE


def normXYZ(
        x: numbers.Real,
        y: numbers.Real,
        z: numbers.Real
) -> Tuple:
    """
    Normalize numeric values.

    :param x: x numeric value
    :param y: y numeric value
    :param z: z numeric value
    :return: the magnitude and a tuple of three float values
    """

    # input vals checks

    vals = [x, y, z]
    if not all(map(lambda val: isinstance(val, numbers.Real), vals)):
        raise Exception("Input values must be integer or float")
    elif not all(map(math.isfinite, vals)):
        raise Exception("Input values must be finite (#01)")

    mag = math.sqrt(x*x + y*y + z*z)

    if mag <= MIN_VECTOR_MAGNITUDE:
        norm_xyz = None
    else:
        norm_xyz = x/mag, y/mag, z/mag

    return mag, norm_xyz


def isclose(
        a,
        b,
        rtol=1e-012,
        atol=1e-12,
        equal_nan=False,
        equal_inf=False
):
    """
    Mimics math.isclose from Python 3.5 (see: https://docs.python.org/3.5/library/math.html)

    Example:
      >>> isclose(1.0, 1.0)
      True
      >>> isclose(1.0, 1.000000000000001)
      True
      >>> isclose(1.0, 1.0000000001)
      False
      >>> isclose(0.0, 0.0)
      True
      >>> isclose(0.0, 0.000000000000001)
      True
      >>> isclose(0.0, 0.0000000001)
      False
      >>> isclose(100000.0, 100000.0)
      True
      >>> isclose(100000.0, 100000.0000000001)
      True
      >>> isclose(np.nan, np.nan)
      False
      >>> isclose(np.nan, 1000000)
      False
      >>> isclose(1.000000000001e300, 1.0e300)
      False
      >>> isclose(1.0000000000001e300, 1.0e300)
      True
      >>> isclose(np.nan, np.nan, equal_nan=True)
      True
      >>> isclose(np.inf, np.inf)
      False
      >>> isclose(np.inf, 1.0e300)
      False
      >>> isclose(np.inf, np.inf, equal_inf=True)
      True
    """

    if equal_nan and a is np.nan and b is np.nan:
        return True
    elif equal_inf and a is np.inf and b is np.inf:
        return True
    elif a is np.inf or b is np.inf:
        return False
    else:
        return abs(a - b) <= max(rtol * max(abs(a), abs(b)), atol)