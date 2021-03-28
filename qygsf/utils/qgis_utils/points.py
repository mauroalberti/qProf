
from qgis.core import *


def distance_projected_pts(
        x,
        y,
        delta_x,
        delta_y,
        src_crs,
        dest_crs
):
    qgspt_start_src_crs = qgs_pt(x, y)
    qgspt_end_src_crs = qgs_pt(x + delta_x, y + delta_y)

    qgspt_start_dest_crs = project_qgs_point(qgspt_start_src_crs, src_crs, dest_crs)
    qgspt_end_dest_crs = project_qgs_point(qgspt_end_src_crs, src_crs, dest_crs)

    pt2_start_dest_crs = Point4D(qgspt_start_dest_crs.x(), qgspt_start_dest_crs.y())
    pt2d_end_dest_crs = Point4D(qgspt_end_dest_crs.x(), qgspt_end_dest_crs.y())

    return pt2_start_dest_crs.dist_2d(pt2d_end_dest_crs)


def project_point(
        pt: Point4D,
        srcCrs: QgsCoordinateReferenceSystem,
        destCrs: QgsCoordinateReferenceSystem
) -> Point4D:

    qgs_pt = QgsPointXY(pt.x, pt.y)

    proj_qgs_pt = project_qgs_point(qgs_pt, srcCrs, destCrs)
    proj_x, proj_y = proj_qgs_pt.x(), proj_qgs_pt.y()

    return Point4D(
        x=proj_x,
        y=proj_y,
        z=pt.z,
        t=pt.t
    )


def project_qgs_point(
        qgsPt: QgsPointXY,
        srcCrs: QgsCoordinateReferenceSystem,
        destCrs: QgsCoordinateReferenceSystem
) -> QgsPointXY:

    return QgsCoordinateTransform(
        srcCrs,
        destCrs,
        QgsProject.instance()
    ).transform(qgsPt)


def pt_geoms_attrs(
    pt_layer,
    field_list=None
):

    if field_list is None:
        field_list = []

    if pt_layer.selectedFeatureCount() > 0:
        features = pt_layer.selectedFeatures()
    else:
        features = pt_layer.getFeatures()

    provider = pt_layer.dataProvider()
    field_indices = [provider.fieldNameIndex(field_name) for field_name in field_list if field_name]

    # retrieve selected features with their geometry and relevant attributes
    rec_list = []
    for feature in features:

        # fetch point geometry
        pt = feature.geometry().asPoint()

        attrs = feature.fields().toList()

        # creates feature attribute list
        feat_list = [pt.x(), pt.y()]
        for field_ndx in field_indices:
            feat_list.append(str(feature.attribute(attrs[field_ndx].name())))

        # add to result list
        rec_list.append(feat_list)

    return rec_list


def qgs_pt(x, y):

    return QgsPointXY(x, y)

