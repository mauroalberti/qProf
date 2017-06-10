# -*- coding: utf-8 -*-

from __future__ import division

from numpy import *  # general import for compatibility with formula input
from numpy.linalg import svd

from .errors import AnaliticSurfaceCalcException


def point_solution(a_array, b_array):
    """
    finds a non-unique solution
    for a set of linear equations
    """

    try:
        return linalg.lstsq(a_array, b_array)[0]
    except:
        return None, None, None


def xyz_svd(xyz_array):
    # modified after: 
    # http://stackoverflow.com/questions/15959411/best-fit-plane-algorithms-why-different-results-solved

    try:
        return dict(result=svd(xyz_array))
    except:
        return dict(result=None)


def formula_to_grid(array_range, array_size, formula):

    a_min, a_max, b_max, b_min = array_range  # note: b range reversed for conventional j order in arrays
    array_rows, array_cols = array_size

    a_array = linspace(a_min, a_max, num=array_cols)
    b_array = linspace(b_max, b_min, num=array_rows)  # note: reversed for conventional j order in arrays

    try:
        a_list, b_list = [a for a in a_array for b in b_array], [b for a in a_array for b in b_array]
    except:
        raise AnaliticSurfaceCalcException("Error in a-b values")

    try:
        z_list = [eval(formula) for a in a_array for b in b_array]
    except:
        raise AnaliticSurfaceCalcException("Error in applying formula to a and b array values")

    return a_list, b_list, z_list


def is_number(s):
    """
    Check if string can be converted to number.

    @param  s:  parameter to check.
    @type  s:  string

    @return:  boolean, whether string can be converted to a number (float).

    """

    try:
        float(s)
    except:
        return False
    else:
        return True


def to_float(iterable):

    return [float(item) for item in iterable]


def almost_zero(val):

    tolerance = 1e-10
    if abs(val) > tolerance:
        return False
    else:
        return True


def ij_transfer_func(i, j, transfer_funcs):
    """
    Return a p_z value as the result of a function (transfer_func_z) applied to a (p_x,p_y) point.
    This point is derived from a (i,j) point given two "transfer" functions (transfer_func_y, transfer_func_x).
    All three functions are stored into a tuple (transfer_funcs).

    @param  i:  array i (-p_y) coordinate of a single point.
    @type  i:  float.
    @param  j:  array j (p_x) coordinate of a single point.
    @type  j:  float.
    @param  transfer_funcs:  tuple storing three functions (transfer_func_x, transfer_func_y, transfer_func_z)
                            that derives p_y from i (transfer_func_y), p_x from j (transfer_func_x)
                            and p_z from (p_x,p_y) (transfer_func_z).
    @type  transfer_funcs:  Tuple of Functions.

    @return:  p_z value - float.

    """

    transfer_func_x, transfer_func_y, transfer_func_z = transfer_funcs

    return transfer_func_z(transfer_func_x(j), transfer_func_y(i))


def array_from_function(row_num, col_num, x_transfer_func, y_transfer_func, z_transfer_func):
    """
    Creates an array of p_z values based on functions that map (i,j) indices (to be created)
    into (p_x, p_y) values and then p_z values.

    @param  row_num:  row number of the array to be created.
    @type  row_num:  int.
    @param  col_num:  column number of the array to be created.
    @type  col_num:  int.
    @param  x_transfer_func:  function that derives p_x given a j array index.
    @type  x_transfer_func:  Function.
    @param  y_transfer_func:  function that derives p_y given an i array index.
    @type  y_transfer_func:  Function.
    @param  z_transfer_func:  function that derives p_z given a (p_x,p_y) point.
    @type  z_transfer_func:  Function.

    @return:  array of p_z value - array of float numbers.

    """

    transfer_funcs = (x_transfer_func, y_transfer_func, z_transfer_func)

    return fromfunction(ij_transfer_func, (row_num, col_num), transfer_funcs=transfer_funcs)

