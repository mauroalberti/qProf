
import unittest

from pygsf.geometries.shapes.space3d import *
from pygsf.orientations.orientations import *


class TestSegmentIntersections(unittest.TestCase):

    def test_segment_intersections(self, n=100):

        for _ in range(n):

            random_segment = Segment2D.random()

            random_interval_value = random.uniform(0, 1)

            rotation_axis = RotationAxis.randomNaive()

            center_point = random_segment.pointAt(random_interval_value)

            rotated_segment = random_segment.rotate(
                rotation_axis=rotation_axis,
                center_point=center_point
            )

            intersection_pt = intersect_segments(random_segment, rotated_segment)

            if not intersection_pt or not intersection_pt.is_coincident(center_point, tolerance=1e-2):
                print(n)
                print(random_segment)
                print(random_interval_value)
                print(rotation_axis)
                print(center_point)
                print(rotated_segment)
                print(intersection_pt)
                if intersection_pt:
                    print(center_point.distance(intersection_pt))


if __name__ == '__main__':
    unittest.main()

