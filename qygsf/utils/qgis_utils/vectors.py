import numbers
from typing import Union, Optional, List, Tuple

from qgis.core import *

from ...geometries.shapes.abstract import *
from ...geometries.shapes.space2d import *
from ...geometries.shapes.space3d import *
from ...geometries.shapes.space4d import *

from .lines import *


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


def try_geoms_attrs(
        layer: QgsVectorLayer,
        field_indices: Optional[List[numbers.Integral]] = None
) -> Tuple[bool, Union[str, List[Tuple[QgsGeometry, Tuple]]]]:
    """
    Returns geometry (unchanged) and attributes of (selected) layer records.
    """

    try:

        assert isinstance(layer, QgsVectorLayer)

        if field_indices is None:
            field_indices = []

        records = []

        if layer.selectedFeatureCount() > 0:
            features = layer.selectedFeatures()
        else:
            features = layer.getFeatures()

        for feature in features:

            geom = feature.geometry()

            attrs = feature.fields().toList()
            rec_data = tuple([feature.attribute(attrs[field_ndx].name()) for field_ndx in field_indices])

            records.append((geom, rec_data))

        return True, records

    except Exception as e:

        return False, str(e)


def try_extract_lines_infos(
        layer: QgsVectorLayer,
        field_indices: Optional[List[numbers.Integral]] = None
) -> Tuple[bool, Union[str, List[Tuple[Line, Tuple]]]]:
    """
    Extract geometry as line and attributes from line layer.
    """

    try:

        success, result = try_geoms_attrs(
            layer=layer,
            field_indices=field_indices
        )

        if not success:
            msg = result
            return False, msg

        records = result

        lines_infos = []

        for line, *res in records:

            if not line.isMultipart():

                success, result = try_extract_line(line)

                if not success:
                    msg = result
                    return False, msg

                ln = result

                lines_infos.append((ln, res))

            else:

                for subline in line:

                    success, result = try_extract_line(subline)

                    if not success:
                        msg = result
                        return False, msg

                    ln = result

                    lines_infos.append((ln, res))

        return True, lines_infos

    except Exception as e:

        return False, str(e)

