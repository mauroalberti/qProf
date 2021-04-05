
import numbers

import math
from math import floor

import numpy as np


def interp_linear(
        frac_s: numbers.Real,
        v0: numbers.Real,
        v1: numbers.Real
) -> numbers.Real:
    """
    Interpolate a number in a linear way.

    :param frac_s: the fractional distance between the start and end point. Range 0-1.
    :type frac_s: numbers.Real.
    :param v0: the value at the start point.
    :type v0: numbers.Real.
    :param v1: the value at the end point.
    ;:type v1: numbers.Real.
    :return: the interpolated value.

    Examples:
      >>> interp_linear(0, 10, 20)
      10
      >>> interp_linear(1, 10, 20)
      20
      >>> interp_linear(-1, 10, 20)
      0
      >>> interp_linear(2, 10, 20)
      30
      >>> interp_linear(0.3, 0, 10)
      3.0
      >>> interp_linear(0.75, 0, 10)
      7.5
    """

    delta_z = v1 - v0
    return v0 + frac_s * delta_z


def scalars_bilin_interp(
        i: numbers.Real,
        j: numbers.Real,
        v00: numbers.Real,
        v01: numbers.Real,
        v10: numbers.Real,
        v11: numbers.Real
) -> numbers.Real:
    """
    Return an interpolated number based on a bilinear interpolation.

    :param i: the delta i relative to the preceding cell center.
    :param j: the delta j relative to the preceding cell center.
    :param v00: the z value of the (i=0, j=0) cell center.
    :param v01: the z value of the (i=0, j=1) cell center.
    :param v10: the z value of the (i=1, j=0) cell center.
    :param v11: the z value of the (i=1, j=1) cell center.
    :return: the interpolated z value.
    """

    grid_val_y0 = v00 + (v10 - v00) * i
    grid_val_y1 = v01 + (v11 - v01) * i

    return grid_val_y0 + (grid_val_y1 - grid_val_y0) * j


if __name__ == "__main__":

    import doctest

    doctest.testmod()


def array_bilin_interp(
        arr: np.ndarray,
        i: numbers.Real,
        j: numbers.Real
) -> numbers.Real:
    """
    Interpolate the z value at a given i,j values couple.
    Interpolation method: bilinear.

    0, 0   0, 1

    1, 0,  1, 1

    :param arr: array with values for which the interpolation will be made.
    :type arr: Numpy array.
    :param i: i array index of the point (may be fractional).
    :type i: numbers.Real.
    :param j: j array index of the point (may be fractional).
    :type j: numbers.Real.
    :return: interpolated z value (may be math.nan).
    :rtype: numbers.Real.
    """

    i_max, j_max = arr.shape
    di = i - floor(i)
    dj = j - floor(j)

    if i < 0.0 or j < 0.0:
        return math.nan
    elif i > i_max - 1 or j > j_max - 1:
        return math.nan
    elif i == i_max - 1 and j == j_max - 1:
        return arr[i, j]
    elif i == i_max - 1:
        v0 = arr[i, int(floor(j))]
        v1 = arr[i, int(floor(j + 1))]
        return interp_linear(
            frac_s=dj,
            v0=v0,
            v1=v1)
    elif j == j_max - 1:
        v0 = arr[int(floor(i)), j]
        v1 = arr[int(floor(i + 1)), j]
        return interp_linear(
            frac_s=di,
            v0=v0,
            v1=v1)
    else:
        v00 = arr[int(floor(i)), int(floor(j))]
        v01 = arr[int(floor(i)), int(floor(j + 1))]
        v10 = arr[int(floor(i + 1)), int(floor(j))]
        v11 = arr[int(floor(i + 1)), int(floor(j + 1))]
        return scalars_bilin_interp(di, dj, v00, v01, v10, v11)