
import numbers


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
    :type i: numbers.Real.
    :param j: the delta j relative to the preceding cell center.
    :type j: numbers.Real.
    :param v00: the z value of the (i=0, j=0) cell center.
    :type v00: numbers.Real.
    :param v01: the z value of the (i=0, j=1) cell center.
    :type v01: numbers.Real.
    :param v10: the z value of the (i=1, j=0) cell center.
    :type v10: numbers.Real.
    :param v11: the z value of the (i=1, j=1) cell center.
    :type v11: numbers.Real.
    :return: the interpolated z value.
    :rtype: numbers.Real.
    """

    grid_val_y0 = v00 + (v10 - v00) * i
    grid_val_y1 = v01 + (v11 - v01) * i

    return grid_val_y0 + (grid_val_y1 - grid_val_y0) * j


if __name__ == "__main__":

    import doctest

    doctest.testmod()
