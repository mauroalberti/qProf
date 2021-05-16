
from operator import itemgetter

from qgis.core import *

from .points import *
from .vectors import *
from ...geometries.shapes.collections import *
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

        success, result = try_line_geoms_with_field_infos(
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
    line_name_fld_ndx: Optional[numbers.Integral],
    invert_direction: bool
) -> Tuple[bool, Union[str, NamedLines]]:

    try:

        areLinesToReorder = False if line_order_fld_ndx is None else True

        # get profile path from input line layer

        if line_order_fld_ndx is None:
            if line_name_fld_ndx is None:
                field_indices = None
            else:
                field_indices = [line_name_fld_ndx]
        else:
            if line_name_fld_ndx is None:
                field_indices = [line_order_fld_ndx]
            else:
                field_indices = [line_order_fld_ndx, line_name_fld_ndx]

        success, result = try_extract_lines_infos(
            layer=line_layer,
            field_indices=field_indices
        )

        if not success:
            msg = result
            return False, msg

        lines_infos = result

        print(f"DEBUG: lines_infos: {lines_infos}")

        if areLinesToReorder:

            lines_infos.sort(key=itemgetter(1))

        # process input line layer

        projected_lines = []
        for ndx, (line, *other) in enumerate(lines_infos):

            projected_lines.append(
                project_line2d(
                    line,
                    line_layer.crs(),
                    project_crs
                )
            )

        profiles = [line.remove_coincident_points() for line in projected_lines]

        if invert_direction:
            profiles = [line.reversed() for line in profiles]

        """
        if line_order_fld_ndx is None:
            if line_name_fld_ndx is None:
                field_indices = None
            else:
                field_indices = [line_name_fld_ndx]
        else:
            if line_name_fld_ndx is None:
                field_indices = [line_order_fld_ndx]
            else:
                field_indices = [line_order_fld_ndx, line_name_fld_ndx]
                """

        print(f"DEBUG: calculating names")

        if line_name_fld_ndx is None:
            names = [""]*len(lines_infos)
        elif line_order_fld_ndx is None:
            names = [name for _, name in lines_infos]
        else:
            names = [name for _, _, name in lines_infos]

        return True, NamedLines(list(zip(names, profiles)))

    except Exception as e:

        return False, str(e)


def try_line_geoms_with_field_infos(
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
    """
    Deprecated: use 'try_line_geoms_attrs'
    """

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


def try_extract_line(
        line: QgsLineString
) -> Tuple[bool, Union[str, Line]]:

    try:

        pts = []
        dim = set()

        print(f"DEBUG: pre-for point in line:")
        print(f"DEBUG: type(line): {type(line)}")
        for point in line.points():
            print(f"DEBUG: type(point): {type(point)}")
            x = point.x()
            y = point.y()
            z = point.z()

            print(f"DEBUG: x, y, z: {x} ,{y}, {z}")

            if np.isfinite(z):
                pt = Point3D(x, y, z)
                dim.add(3)
            else:
                pt = Point2D(x, y)
                dim.add(2)

            pts.append(pt)

        if len(dim) != 1:
            return False, f"Dim is {len(dim)}"

        dim = list(dim)[0]

        print(f"DEBUG: dim: {dim}")
        if dim == 2:
            return True, Line2D(pts)
        elif dim == 3:
            return True, Line3D(pts)
        else:
            return False, f"Got a dim of {dim}"

    except Exception as e:

        return False, str(e)


def try_extract_lines_infos(
        layer: QgsVectorLayer,
        field_indices: Optional[List[numbers.Integral]] = None
) -> Tuple[bool, Union[str, List[Tuple]]]:
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

            print(f"DEBUG: *res: {res}")

            if not line.isMultipart():

                success, result = try_extract_line(line)

                if not success:
                    msg = result
                    return False, msg

                ln = result

                new_rec = (ln, *res)
                print(f"DEBUG: new_rec: {new_rec}")

                lines_infos.append((ln, *res))

            else:

                print(f"DEBUG: pre-for subline in line:")
                for subline in line.parts():
                    print(f"DEBUG: inside line loop")
                    success, result = try_extract_line(subline)

                    if not success:
                        msg = result
                        return False, msg

                    ln = result

                    lines_infos.append((ln, *res))

        return True, lines_infos

    except Exception as e:

        return False, str(e)

