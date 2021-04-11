# coding: utf-8

# # Check bilinear interpolation


import unittest
import math
import numpy as np

from ..mathematics.scalars import areClose
from ..geometries.grids.interpolations import scalars_bilin_interp as s_bilinear_interp
from ..geometries.grids.interpolations import array_bilin_interp as a_bilinear_interp


class TestBilinearInterpolations(unittest.TestCase):

    def test_scala_interpolation(self):

        assert areClose(s_bilinear_interp(
            i=0,
            j=0,
            v00=42.0,
            v01=17.2,
            v10=-44.2,
            v11=19.4),
            42.0)

        assert areClose(s_bilinear_interp(
            i=0,
            j=1,
            v00=42.0,
            v01=17.2,
            v10=-44.2,
            v11=19.4),
            17.2)

        assert areClose(s_bilinear_interp(
            i=1,
            j=0,
            v00=42.0,
            v01=17.2,
            v10=-44.2,
            v11=19.4),
            -44.2)

        assert areClose(s_bilinear_interp(
            i=1,
            j=1,
            v00=42.0,
            v01=17.2,
            v10=-44.2,
            v11=19.4),
            19.4)

        assert areClose(s_bilinear_interp(
            i=0,
            j=0.5,
            v00=10.0,
            v01=30.0,
            v10=-30,
            v11=-10),
            20.0)

        assert areClose(s_bilinear_interp(
            i=0.5,
            j=0.0,
            v00=10.0,
            v01=20.2,
            v10=-34.2,
            v11=19.4),
            -12.1)

        assert areClose(s_bilinear_interp(
            i=0.5,
            j=0.5,
            v00=0.0,
            v01=40.0,
            v10=-30.0,
            v11=-10.0),
            0.0)

    def test_array_interpolation(self):

        arr = np.array([[0, 1],
                        [1, 2]])

        assert a_bilinear_interp(
            arr=arr,
            i=-1,
            j=-1) is math.nan

        assert a_bilinear_interp(
            arr=arr,
            i=2,
            j=2) is math.nan

        assert a_bilinear_interp(
            arr=arr,
            i=0.5,
            j=1.5) is math.nan

        assert areClose(a_bilinear_interp(
            arr=arr,
            i=0,
            j=0),
            0.0)

        assert areClose(a_bilinear_interp(
            arr=arr,
            i=0,
            j=1),
            1.0)

        assert areClose(a_bilinear_interp(
            arr=arr,
            i=1,
            j=0),
            1.0)

        assert areClose(a_bilinear_interp(
            arr=arr,
            i=1,
            j=1),
            2.0)

        assert areClose(a_bilinear_interp(
            arr=arr,
            i=0.5,
            j=0),
            0.5)

        assert areClose(a_bilinear_interp(
            arr=arr,
            i=0.5,
            j=1),
            1.5)

        assert areClose(a_bilinear_interp(
            arr=arr,
            i=0.75,
            j=1),
            1.75)

        assert areClose(a_bilinear_interp(
            arr=arr,
            i=0.5,
            j=0.5),
            1.0)

        assert areClose(a_bilinear_interp(
            arr=arr,
            i=1,
            j=0.5),
            1.5)

        assert areClose(a_bilinear_interp(
            arr=arr,
            i=1,
            j=0.75),
            1.75)


if __name__ == '__main__':

    unittest.main()

