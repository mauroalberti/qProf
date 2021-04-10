
from typing import Callable

from math import sin, cos

import numbers


from .geotransform import *


def ij_transfer_func(
        i: numbers.Real,
        j: numbers.Real,
        geotransform: GeoTransform,
        z_transfer_func: Callable,
        i_shift=0.5,
        j_shift=0.5
) -> numbers.Real:
    """
    Return a z value as the result of a function (transfer_func_z) applied to a
    (i+i_shift,j+j_shift) point (i.e., with defaultvalues, the cell center, not the cell top-left corner)
    given a geotransform.

    :param  i:  array i (-y) coordinate of a single point.
    :type  i:  numbers.Real.
    :param  j:  array j (x) coordinate of a single point.
    :type  j:  numbers.Real.
    :param  geotransform:  geotransform
    :type  geotransform:  GeoTransform.
    :param  z_transfer_func:  function that calculates the z value given x and y input
    :type  z_transfer_func:  function.
    :param i_shift: cell unit shift in the i direction with respect to the cell top-left corner. Default is 0.5, i.e. half the cell size
    :type i_shift: numbers.Real
    :param j_shift: cell unit shift in the j direction with respect to the cell top-left corner. Default is 0.5, i.e. half the cell size
    :type j_shift: numbers.Real
    :return: z value
    :rtype: numbers.Real.

    Examples:
    """

    return z_transfer_func(*ijPixToxyGeogr(geotransform, i + i_shift, j + j_shift))


def array_from_geotransform_function(
        row_num: numbers.Integral,
        col_num: numbers.Integral,
        geotransform: GeoTransform,
        z_transfer_func: Callable
) -> np.ndarray:
    """
    Creates an array of z values based on functions that map (i,j) indices (to be created)
    into (x, y) values and then z values.

    :param  row_num:  row number of the array to be created.
    :type  row_num:  numbers.Integral.
    :param  col_num:  column number of the array to be created.
    :type  col_num:  numbers.Integral.
    :param  geotransform:  the used geotransform.
    :type  geotransform:  GeoTransform.
    :param  z_transfer_func:  function that derives z given a (x, y) point.
    :type  z_transfer_func:  Callable.

    :return:  array of z values
    :rtype: np.ndarray of numbers.Real numbers.

    Examples:
    """

    array = np.fromfunction(
        function=ij_transfer_func,
        shape=(row_num, col_num),
        dtype=np.float64,
        geotransform=geotransform,
        z_transfer_func=z_transfer_func)

    return np.asarray(array)


def grad_j(
        fld: np.ndarray,
        cell_size_j: numbers.Real,
        edge_order: numbers.Integral = 2
) -> np.ndarray:
    """
    Calculates the array gradient along the j axis.

    :param fld: array.
    :type fld: np.ndarray.
    :param cell_size_j: the cell spacing in the x direction.
    :type cell_size_j: numbers.Real.
    :param edge_order: the type of edge order used in the Numpy gradient method.
    :type edge_order: numbers.Integral.
    :return: gradient field.
    :rtype: np.ndarray.

    Examples:
    """

    return np.gradient(
        fld,
        edge_order=edge_order,
        axis=1
    ) / cell_size_j


def grad_i(
        fld: np.ndarray,
        cell_size_i: numbers.Real,
        edge_order: numbers.Integral = 2
) -> np.ndarray:
    """
    Calculates the array gradient along the i axis.

    :param fld: array.
    :type fld: np.ndarray.
    :param cell_size_i: the cell spacing in the y direction.
    :type cell_size_i: numbers.Real.
    :param edge_order: the type of edge order used in the Numpy gradient method.
    :type edge_order: numbers.Integral.
    :return: gradient field.
    :rtype: np.ndarray.

    Examples:
    """

    return np.gradient(fld, edge_order=edge_order, axis=0) / cell_size_i


def grad_iminus(
        fld: np.ndarray,
        cell_size_i: numbers.Real,
        edge_order: numbers.Integral = 2
) -> np.ndarray:
    """
    Calculates the array gradient along the -i axis.

    :param fld: array.
    :type fld: np.ndarray.
    :param cell_size_i: the cell spacing in the y direction.
    :type cell_size_i: numbers.Real.
    :param edge_order: the type of edge order used in the Numpy gradient method.
    :type edge_order: numbers.Integral.
    :return: gradient field.
    :rtype: np.ndarray.

    Examples:
    """

    return - np.gradient(fld, edge_order=edge_order, axis=0) / cell_size_i


def dir_deriv(
        fld: np.ndarray,
        cell_size_x: numbers.Real,
        cell_size_y: numbers.Real,
        direct_rad: numbers.Real,
        dx_edge_order: numbers.Integral = 2,
        dy_edge_order: numbers.Integral = 2
) -> np.ndarray:
    """
    Calculates the directional derivative in the provided direction.

    :param fld: the field.
    :type fld: Numpy array.
    :param cell_size_x: the cell size along the x axis.
    :type cell_size_x: numbers.Real.
    :param cell_size_y: the cell size along the y
    :param direct_rad: the direction, expressed as radians.
    :type direct_rad: numbers.Real.
    :param dx_edge_order: the edge order of the gradient along x.
    :type dx_edge_order: numbers.Integral.
    :param dy_edge_order: the edge order of the gradient along y.
    :type dy_edge_order: numbers.Integral.
    :return: the directional derivative array.
    :rtype: Numpy array.
    """

    df_dx = grad_j(
        fld=fld,
        cell_size_j=cell_size_x,
        edge_order=dx_edge_order)

    df_dy = grad_iminus(
        fld=fld,
        cell_size_i=cell_size_y,
        edge_order=dy_edge_order)

    return df_dx * sin(direct_rad) + df_dy * cos(direct_rad)


def magnitude(
        fld_x: np.ndarray,
        fld_y: np.ndarray
) -> np.ndarray:
    """
    Calculates the magnitude given two 2D arrays:
    the first represents the vector field x component, the second the vector field y component.

    :param fld_x: vector field x component.
    :type fld_x: np.ndarray.
    :param fld_y: vector field y component.
    :type fld_y: np.ndarray.
    :return: magnitude field.
    :rtype: np.ndarray.

    Examples:
    """

    return np.sqrt(fld_x ** 2 + fld_y ** 2)


def orients_r(
        fld_x: np.ndarray,
        fld_y: np.ndarray
) -> np.ndarray:
    """
    Calculates the orientations (as radians) given two 2D arrays:
    the first represents the vector field x component, the second the vector field y component.

    :param fld_x: vector field x component.
    :type fld_x: np.ndarray.
    :param fld_y: vector field y component.
    :type fld_y: np.ndarray.
    :return: orientation field, in radians.
    :rtype: np.ndarray.

    Examples:
    """

    azimuth_rad = np.arctan2(fld_x, fld_y)
    azimuth_rad = np.where(azimuth_rad < 0.0, azimuth_rad + 2*np.pi, azimuth_rad)

    return azimuth_rad


def orients_d(
        fld_x: np.ndarray,
        fld_y: np.ndarray
) -> np.ndarray:
    """
    Calculates the orientations (as decimal degrees) given two 2D arrays:
    the first represents the vector field x component, the second the vector field y component.

    :param fld_x: vector field x component.
    :type fld_x: np.ndarray.
    :param fld_y: vector field y component.
    :type fld_y: np.ndarray.
    :return: orientation field, in decimal degrees.
    :rtype: np.ndarray.

    Examples:
    """

    return np.degrees(orients_r(fld_x, fld_y))


def divergence(
        fld_x: np.ndarray,
        fld_y: np.ndarray,
        cell_size_x: numbers.Real,
        cell_size_y: numbers.Real
) -> np.ndarray:
    """
    Calculates the divergence from two 2D arrays:
    the first represents the vector field x component, the second the vector field y component.

    :param fld_x: vector field x component.
    :type fld_x: np.ndarray.
    :param fld_y: vector field y component.
    :type fld_y: np.ndarray.
    :param cell_size_x: the cell spacing in the x direction.
    :type cell_size_x: numbers.Real.
    :param cell_size_y: the cell spacing in the y direction.
    :type cell_size_y: numbers.Real.
    :return: divergence field.
    :rtype: np.ndarray.

    Examples:
    """

    dfx_dx = grad_j(fld_x, cell_size_x)
    dfy_dy = grad_iminus(fld_y, cell_size_y)

    return dfx_dx + dfy_dy


def curl_module(
        fld_x: np.ndarray,
        fld_y: np.ndarray,
        cell_size_x: numbers.Real,
        cell_size_y: numbers.Real
) -> np.ndarray:
    """
    Calculates the curl module from two 2D arrays:
    the first represents the vector field x component, the second the vector field y component.

    :param fld_x: vector field x component.
    :type fld_x: np.ndarray.
    :param fld_y: vector field y component.
    :type fld_y: np.ndarray.
    :param cell_size_x: the cell spacing in the x direction.
    :type cell_size_x: numbers.Real.
    :param cell_size_y: the cell spacing in the y direction.
    :type cell_size_y: numbers.Real.
    :return: curl field.
    :rtype: np.ndarray.

    Examples:
    """

    dfx_dy = grad_iminus(fld_x, cell_size_y, edge_order=2)
    dfy_dx = grad_j(fld_y, cell_size_x, edge_order=1)

    return dfy_dx - dfx_dy


def magn_grads(
        fld_x: np.ndarray,
        fld_y: np.ndarray,
        dir_cell_sizes: List[numbers.Real],
        axis: str = ''
) -> List[np.ndarray]:
    """
    Calculates the magnitude gradient along the given direction, based on the field-defining two 2D arrays:
    the first representing the x component, the second the y component.

    :param fld_x: vector field x component.
    :type fld_x: np.ndarray.
    :param fld_y: vector field y component.
    :type fld_y: np.ndarray.
    :param dir_cell_sizes: list of cell spacing(s) in the considered direction(s).
    :type dir_cell_sizes: list of numbers.Real(s).
    :param axis: declares the axis ('x' or 'y') or the axes('', i.e., empty string) for both x and y directions.
    :type axis: str.
    :return: magnitude gradient field(s) along the considered direction.
    :rtype: list of np.ndarray.
    :raises: Exception.

    Examples:
    """

    magn = magnitude(fld_x, fld_y)
    if axis == 'x':
        return [grad_j(magn, dir_cell_sizes[0])]
    elif axis == 'y':
        return [grad_iminus(magn, dir_cell_sizes[0])]
    elif axis == '':
        return [grad_j(magn, dir_cell_sizes[0]), grad_iminus(magn, dir_cell_sizes[1])]
    else:
        raise Exception("Axis must be 'x' or 'y' or '' (for both x and y). '{}' given".format(axis))


def magn_grad_along_flowlines(
        fld_x: np.ndarray,
        fld_y: np.ndarray,
        cell_size_x: numbers.Real,
        cell_size_y: numbers.Real
) -> np.ndarray:
    """
    Calculates gradient along flow lines.

    :param fld_x: vector field x component.
    :type fld_x: np.ndarray.
    :param fld_y: vector field y component.
    :type fld_y: np.ndarray.
    :param cell_size_x: the cell spacing in the x direction.
    :type cell_size_x: numbers.Real.
    :param cell_size_y: the cell spacing in the y direction.
    :type cell_size_y: numbers.Real.
    :return: the flowline gradient field
    :rtype: np.ndarray.
    """

    orien_rad = orients_r(fld_x, fld_y)

    dm_dx, dm_dy = magn_grads(
        fld_x=fld_x,
        fld_y=fld_y,
        dir_cell_sizes=[cell_size_x, cell_size_y])

    velocity_gradient = dm_dx * np.sin(orien_rad) + dm_dy * np.cos(orien_rad)

    return velocity_gradient
