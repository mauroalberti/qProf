
import datetime

from qgis.core import QgsCoordinateReferenceSystem, QgsPointXY, QgsCoordinateTransform, QgsProject

from ...geometries.shapes.space3d import *
from ...geometries.shapes.space4d import Point4D
from ...georeferenced.geodetic import geodetic2ecef
from ..time import standard_gpstime_to_seconds


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

    pt2_start_dest_crs = Point2D(qgspt_start_dest_crs.x(), qgspt_start_dest_crs.y())
    pt2d_end_dest_crs = Point2D(qgspt_end_dest_crs.x(), qgspt_end_dest_crs.y())

    return pt2_start_dest_crs.distance(pt2d_end_dest_crs)


def project_point(
        pt: Point2D,
        srcCrs: QgsCoordinateReferenceSystem,
        destCrs: QgsCoordinateReferenceSystem
) -> Point2D:

    qgs_pt = QgsPointXY(pt.x, pt.y)

    proj_qgs_pt = project_qgs_point(qgs_pt, srcCrs, destCrs)
    proj_x, proj_y = proj_qgs_pt.x(), proj_qgs_pt.y()

    return Point2D(
        x=proj_x,
        y=proj_y
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


def calculate_pts_in_projection(pts_in_orig_crs, srcCrs, destCrs):

    pts_in_prj_crs = []
    for pt in pts_in_orig_crs:
        qgs_pt = QgsPointXY(pt.x, pt.y)
        qgs_pt_prj_crs = project_qgs_point(qgs_pt, srcCrs, destCrs)
        pts_in_prj_crs.append(Point2D(qgs_pt_prj_crs.x(), qgs_pt_prj_crs.y()))
    return pts_in_prj_crs


class TrackPointGPX(object):

    def __init__(self,
                 lat: numbers.Real,
                 lon: numbers.Real,
                 elev: numbers.Real,
                 time: datetime.datetime):

        self.lat = float(lat)
        self.lon = float(lon)
        self.elev = float(elev)
        self.time = time

    def as_pt4d(self):

        x, y, _ = geodetic2ecef(self.lat, self.lon, self.elev)
        t = standard_gpstime_to_seconds(self.time)

        return Point4D(x, y, self.elev, t)

    def project(self,
                dest_crs: QgsCoordinateReferenceSystem
                ) -> Point4D:

        pt = Point2D(
            x=self.lon,
            y=self.lat
        )

        crs = QgsCoordinateReferenceSystem("EPSG:4326")

        projected_pt = project_point(
                pt=pt,
                srcCrs=crs,
                destCrs=dest_crs)

        return Point4D(
            x=projected_pt.x,
            y=projected_pt.y,
            z=self.elev,
            t=self.time
        )