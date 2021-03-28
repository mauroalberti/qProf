# coding: utf-8

# # Divergence result check

# Theoretical example based on:
# 
# https://www.khanacademy.org/math/multivariable-calculus/multivariable-derivatives/divergence-and-curl-articles/a/divergence
# 


import unittest

from pygsf.geometries.rasters.fields import *


# Z transfer functions

# These analytical functions define the value of the cells, from the provided x and y geographic coordinates.

def z_func_fx(x, y):

    return 2 * x - y


def z_func_fy(x, y):

    return y * y


def z_func_div(x, y):

    return 2 + 2 * y


class TestDivergence(unittest.TestCase):

    def setUp(self):

        pass

    def test_divergence_1(self):
        """
        Test the divergence calculation.

        :return:
        """

        # The geotransform defines the raster-to-geographic coordinates mapping.

        gt1 = GeoTransform(1500, 3000, 10, 10)

        # Z fields (x and y components)

        # The gridded field values are calculated for the vector field x- and y- components,
        # as well as the expected teorethical divergence field.

        rows = 100
        cols = 50

        # Vector field x-component

        fx = array_from_function(row_num=rows, col_num=cols, geotransform=gt1, z_transfer_func=z_func_fx)

        # Vector field y-component

        fy = array_from_function(row_num=rows, col_num=cols, geotransform=gt1, z_transfer_func=z_func_fy)

        # Theoretical divergence

        div_theor = array_from_function(row_num=rows, col_num=cols, geotransform=gt1, z_transfer_func=z_func_div)

        # Divergence as resulting from pygsf calculation:

        div_pygsf = divergence(fx, fy, 10, 10)

        assert np.allclose(div_theor, div_pygsf)

    def test_divergence_2(self):
        """
        Test the divergence calculation.

        :return:
        """

        def z_func_fx(x, y):
            return 0.0001 * x * y ** 3

        def z_func_fy(x, y):
            return - 0.0002 * x ** 2 * y

        def z_func_div(x, y):
            return 0.0001 * y ** 3 - 0.0002 * x ** 2

        rows = 100
        cols = 200

        size_x = 10
        size_y = 10

        tlx = 500.0
        tly = 250.0

        gt1 = GeoTransform(
            inTopLeftX=tlx,
            inTopLeftY=tly,
            inPixWidth=size_x,
            inPixHeight=size_y)

        fx = array_from_function(
            row_num=rows,
            col_num=cols,
            geotransform=gt1,
            z_transfer_func=z_func_fx)

        fy = array_from_function(
            row_num=rows,
            col_num=cols,
            geotransform=gt1,
            z_transfer_func=z_func_fy)

        # theoretical divergence

        theor_div = array_from_function(
            row_num=rows,
            col_num=cols,
            geotransform=gt1,
            z_transfer_func=z_func_div)

        # estimated divergence

        div = divergence(
            fld_x=fx,
            fld_y=fy,
            cell_size_x=size_x,
            cell_size_y=size_y)

        assert np.allclose(theor_div, div)

    def test_curlmodule_1(self):
        """
        Test the curl module calculation.

        :return:
        """

        def z_func_fx(x, y):

            return y

        def z_func_fy(x, y):

            return - x

        # ### geotransform and grid definitions

        rows = 200
        cols = 200

        size_x = 10
        size_y = 10

        tlx = -1000.0
        tly = 1000.0

        gt = GeoTransform(
            inTopLeftX=tlx,
            inTopLeftY=tly,
            inPixWidth=size_x,
            inPixHeight=size_y)

        fx = array_from_function(
            row_num=rows,
            col_num=cols,
            geotransform=gt,
            z_transfer_func=z_func_fx)

        fy = array_from_function(
            row_num=rows,
            col_num=cols,
            geotransform=gt,
            z_transfer_func=z_func_fy)

        curl_mod = curl_module(
            fld_x=fx,
            fld_y=fy,
            cell_size_x=size_x,
            cell_size_y=size_y)

        assert np.allclose(-2.0, curl_mod)


