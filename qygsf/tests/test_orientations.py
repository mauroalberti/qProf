
import unittest

from pygsf.orientations.orientations import *
from pygsf.geometries.shapes.space3d import *
from pygsf.orientations.direct_utils import *
from pygsf.geometries.shapes.space3d import *
from pygsf.geometries.shapes.space3d import CPlane3D


class TestOrientations(unittest.TestCase):

    def setUp(self):

        pass

    def test_direct_general(self):
        """
        Check expected OrienM results for downward dip.
        """

        assert Direct(90, 90).is_downward
        assert Direct(90, -45).is_upward
        assert areClose(Direct(90, 90).as_versor().z, -1.0)
        assert areClose(Direct(90, -90).as_versor().z, 1.0)
        assert areClose(Direct(0, 90).upward().as_versor().z, 1.0)
        assert areClose(Direct(0, -90).downward().as_versor().z, -1.0)

    def test_direct_angle(self):

        assert areClose(Direct(90, 45).angle_as_degrees(
            Direct(90, 55)), 10.)
        assert areClose(Direct(90, 45).angle_as_degrees(
            Direct(270, 10)), 125.)
        assert areClose(Direct(90, 90).angle_as_degrees(
            Direct(135, 90)), 0.)
        assert areClose(Direct(0, 0).angle_as_degrees(
            Direct(135, 0)), 135.)
        assert areClose(Direct(0, 80).angle_as_degrees(
            Direct(180, 80)), 20.)

    def test_axis_angle(self):

        assert areClose(Axis.fromAzPl(90, 0).angle_as_degrees(
            Axis.fromAzPl(270, 0)), 0.)

    def test_plane_normal(self):

        assert areClose(
            Plane(90, 45).normDirectFrwrd().angle_as_degrees(
                Direct(90, -45)), 0.)

    def test_plane2cplane(self):

        pl = CPlane3D.from_geological_plane(Plane(90, 45), Point2D(0, 0, 0))
        assert areClose(pl.angle_as_degrees(CPlane3D(1, 0, 1, 0)), 0.0)

    def test_plane_angle(self):

        assert areClose(
            Plane(90, 45).angle_degr(
                Plane(90, 45)), 0.)
        assert areClose(
            Plane(90, 45).angle_degr(
                Plane(90, 55)), 10.)
        assert areClose(
            Plane(90, 5).angle_degr(
                Plane(270, 5)), 10.)
        assert areClose(
            Plane(90, 85).angle_degr(
                Plane(270, 85)), 10.)
        assert areClose(
            Plane(0, 0).angle_degr(
                Plane(0, 10)), 10.)
        assert areClose(
            Plane(0, 0).angle_degr(
                Plane(180, 0)), 0.)

    def tearDown(self):

        pass


if __name__ == '__main__':

    unittest.main()
