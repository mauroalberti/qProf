
from typing import Optional

from qgis.core import *


def vector_type(
        layer: QgsMapLayer.VectorLayer
) -> Optional[str]:

    if layer.geometryType() == QgsWkbTypes.PointGeometry:
        return "point"
    elif layer.geometryType() == QgsWkbTypes.LineGeometry:
        return "line"
    elif layer.geometryType() == QgsWkbTypes.PolygonGeometry:
        return "polygon"
    else:
        return None

