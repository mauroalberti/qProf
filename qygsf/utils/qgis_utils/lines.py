
from qgis.core import QgsVectorLayer, QgsCoordinateReferenceSystem,QgsPointXY

from ..qgis_utils.points import project_point, project_qgs_point
from ...geometries.shapes.space2d import *


def polyline_to_xytuple_list(
        qgsline
):

    return [(qgspoint.x(), qgspoint.y()) for qgspoint in qgsline]


def multipolyline_to_xytuple_list2(
        qgspolyline
):

    return [polyline_to_xytuple_list(qgsline) for qgsline in qgspolyline]


def try_get_line_traces(
    line_shape,
    order_field_ndx: Optional[numbers.Integral] = None
) -> Tuple[bool, Union[str, Tuple]]:

    try:

        success, result = try_line_geoms_with_order_infos(
            line_shape,
            order_field_ndx
        )

        if not success:
            msg = result
            return False, msg

        profile_orig_lines, order_values = result

        return True, (profile_orig_lines, order_values)

    except Exception as e:

        return False, str(e)


def project_line2d(
    src_line2d: Line2D,
    src_crs: QgsCoordinateReferenceSystem,
    dest_crs: QgsCoordinateReferenceSystem
) -> Line2D:

    projected_pts = []

    for pt in src_line2d.pts():
        projected_pt = project_point(
            pt=pt,
            srcCrs=src_crs,
            destCrs=dest_crs
        )
        projected_pts.append(projected_pt)

    return Line2D(pts=projected_pts)


def try_load_line_layer(
    line_layer: QgsVectorLayer,
    project_crs,
    line_order_fld_ndx: Optional[numbers.Integral],
    invert_direction: bool
) -> Tuple[bool, Union[str, List[Line2D]]]:

    try:

        areLinesToReorder = False if line_order_fld_ndx is None else True

        # get profile path from input line layer

        success, result = try_get_line_traces(
            line_layer,
            line_order_fld_ndx
        )

        if not success:
            raise Exception(result)

        profile_orig_lines, order_values = result

        processed_lines = []

        if areLinesToReorder:

            processed_lines.append(merge_lines2d(profile_orig_lines, order_values))

        else:

            for orig_line in profile_orig_lines:
                processed_lines.append(merge_line2d(orig_line))

        # process input line layer

        projected_lines = []
        for ndx, processed_line in enumerate(processed_lines):

            projected_lines.append(
                project_line2d(
                    processed_line,
                    line_layer.crs(),
                    project_crs
                )
            )

        profiles = [line.remove_coincident_points() for line in projected_lines]

        if invert_direction:
            profiles = [line.invert_direction() for line in profiles]

        return True, profiles

    except Exception as e:

        return False, str(e)


def try_line_geoms_with_order_infos(
    line_layer,
    order_field_ndx: Optional[numbers.Integral] = None
) -> Tuple[bool, Union[str, Tuple[List, List]]]:

    try:

        lines = []
        order_values = []

        if line_layer.selectedFeatureCount() > 0:
            features = line_layer.selectedFeatures()
        else:
            features = line_layer.getFeatures()

        dummy_progressive = 0

        for feature in features:

            dummy_progressive += 1

            order_val = feature[order_field_ndx] if order_field_ndx is not None else dummy_progressive

            order_values.append(order_val)

            geom = feature.geometry()

            if geom.isMultipart():

                lines.append(
                    (
                        'multiline',
                        multipolyline_to_xytuple_list2(geom.asMultiPolyline())  # geom is QVector<QgsPolyline>
                     )
                )
                # now it's a list of list of (x,y) tuples

            else:

                lines.append(
                    (
                        'line',
                        polyline_to_xytuple_list(geom.asPolyline())  # geom is QVector<QgsPointXY>
                    )
                )

        return True, (lines, order_values)

    except Exception as e:

        return False, str(e)


def line_geoms_attrs(
    line_layer,
    field_list=None
):

    if field_list is None:
        field_list = []

    lines = []

    if line_layer.selectedFeatureCount() > 0:
        features = line_layer.selectedFeatures()
    else:
        features = line_layer.getFeatures()

    provider = line_layer.dataProvider()
    field_indices = [provider.fieldNameIndex(field_name) for field_name in field_list]

    for feature in features:
        geom = feature.geometry()
        if geom.isMultipart():
            rec_geom = multipolyline_to_xytuple_list2(geom.asMultiPolyline())
        else:
            rec_geom = [polyline_to_xytuple_list(geom.asPolyline())]

        attrs = feature.fields().toList()
        rec_data = [str(feature.attribute(attrs[field_ndx].name())) for field_ndx in field_indices]

        lines.append([rec_geom, rec_data])

    return lines


def project_xy_list(
    src_crs_xy_list,
    srcCrs,
    destCrs
):

    pt_list_dest_crs = []
    for x, y in src_crs_xy_list.pts:
        srcPt = QgsPointXY(x, y)
        destPt = project_qgs_point(srcPt, srcCrs, destCrs)
        pt_list_dest_crs = pt_list_dest_crs.append([destPt.x(), destPt.y()])

    return pt_list_dest_crs

