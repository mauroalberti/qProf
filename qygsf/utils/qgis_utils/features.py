
from typing import List, Tuple

from qgis.core import *


from ...geometries.shapes.space3d import *


def extract_from_linestring_3d(
        linestring: QgsLineString
) -> Line3D:

    pts = []

    for pt in linestring:

        pts.append(
            Point3D(
                x=pt.x,
                y=pt.y,
                z=pt.z
            )
        )

    return Line3D(pts)


def extract_from_multilinestring_3d(
        multilinestring: QgsMultiLineString
) -> MultiLine3D:

    lines = []

    for line in multilinestring:

        lines.append(
            extract_from_linestring_3d(line)
        )

    return MultiLine3D(lines)
