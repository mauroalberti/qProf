# coding: utf-8

import unittest

from pygsf.georeferenced.rasters import *


class TestRKFInterpolation(unittest.TestCase):

    def test_rkf_interpolation(self):

        # ## Velocity field with circular motion

        # The example circular motion vector field has components:
        #
        # **v** = y **i** - x **j**

        # as deriving from the equation:
        #
        # **v** = - **z** x **r**

        # where **z** is the vertical vector, **r** the position vector and *x* the vector product.

        k = 2 * math.pi

        def z_func_fx(x, y):

            return k*y

        def z_func_fy(x, y):

            return -k*x

        # The velocity field parameters for testing the results are:
        #
        # v = w * r
        #
        # w = v / r
        #
        # |v| = sqrt(k^2*y^2 + k^2*x^2) = k * r

        # 1 cycle -> 2 pi r
        #
        # v = ds / dt -> ds = v * dt
        #
        # 2 pi r = v dt
        #
        # 2 pi r = v T -> T = 2 pi r / v = 2 pi / k

        # geotransform and grid definitions

        rows = 100
        cols = 100

        size_x = 1
        size_y = 1

        tlx = -50.0
        tly = 50.0

        gt1 = GeoTransform(
            inTopLeftX=tlx,
            inTopLeftY=tly,
            inPixWidth=size_x,
            inPixHeight=size_y)

        # vector field x-component

        fx1 = array_from_function(
            row_num=rows,
            col_num=cols,
            geotransform=gt1,
            z_transfer_func=z_func_fx)

        # vector field y-component

        fy1 = array_from_function(
            row_num=rows,
            col_num=cols,
            geotransform=gt1,
            z_transfer_func=z_func_fy)

        # geoarray definition

        ga = GeoArray(
            inGeotransform=gt1,
            epsg_code=32633,
            inLevels=[fx1, fy1])

        time_increm = 1.0e-4

        period = 2 * math.pi / k

        number_of_cycles = 10

        steps = number_of_cycles * (period / time_increm)

        first_pt = Point2D(
            x=0,
            y=20)

        str_pt = first_pt
        pts_x, pts_y = [first_pt.x], [first_pt.y]

        for n in range(int(steps)):

            end_pt, error = interpolate_rkf(
                geoarray=ga,
                delta_time=time_increm,
                start_pt=str_pt)

            if end_pt is None:
                break

            pts_x.append(end_pt.x)
            pts_y.append(end_pt.y)
            str_pt = end_pt

        assert end_pt.is_coincident(first_pt)
        # After 10 cycles the calculated point position is in the expected (initial) position: x=0, y=200.
