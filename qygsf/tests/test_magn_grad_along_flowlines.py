# Testing magnitude gradient along flowlines

from math import sqrt, degrees, atan

import unittest

from pygsf.geometries.rasters.fields import *
from pygsf.geometries.rasters.geotransform import *


class TestMagnGradFlwlns(unittest.TestCase):

    def setUp(self):

        pass

    def test_magn_grad_flns_1(self):
        """
        Test the gradients calculations.

        :return:
        """

        fx = np.array([
            [1, 1, 1],
            [1, 1, 1],
            [1, 1, 1]
        ])

        fy = np.array([
            [2, 2, 2],
            [2, 2, 2],
            [2, 2, 2]
        ])

        magn = magnitude(fx, fy)

        assert np.allclose(
            magn,
            sqrt(1 + 2**2))

        oriens_d = orients_d(fx, fy)

        assert np.allclose(
            oriens_d,
            degrees(atan(0.5)))

        mgflwln = magn_grad_along_flowlines(fx, fy, 10, 10)

        assert np.allclose(
            mgflwln,
            0.0)

    def test_magn_grad_flns_2(self):
        """
        Test the gradients calculations.

        :return:
        """

        fx = np.array([
            [1, 2, 3],
            [1, 2, 3],
            [1, 2, 3]
        ])

        fy = np.array([
            [0, 0, 0],
            [0, 0, 0],
            [0, 0, 0]
        ])

        oriens_d = orients_d(fx, fy)

        assert np.allclose(
            oriens_d,
            90.0)

        mag = magnitude(fx, fy)

        assert np.allclose(
            mag,
            np.sqrt(fx**2 + fy**2))

        mgflwlns = magn_grad_along_flowlines(
                fld_x=fx,
                fld_y=fy,
                cell_size_x=10,
                cell_size_y=10)

        assert np.allclose(
            mgflwlns,
            0.1)

    def test_magn_grad_flns_3(self):
        """
        Test the gradients calculations.

        :return:
        """

        fx = np.array([
            [1, 2, 3],
            [2, 3, 4],
            [3, 4, 5]
        ])

        fy = -fx

        mgflwlns = magn_grad_along_flowlines(
                fld_x=fx,
                fld_y=fy,
                cell_size_x=1,
                cell_size_y=1)

        assert np.allclose(
            mgflwlns,
            2.0)

        mgflwlns = magn_grad_along_flowlines(
                fld_x=fx,
                fld_y=fy,
                cell_size_x=10,
                cell_size_y=10)

        assert np.allclose(
            mgflwlns,
            0.2)
