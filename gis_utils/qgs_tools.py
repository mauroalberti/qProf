from __future__ import division

from math import isnan, sin, cos, asin, radians, degrees, floor, ceil, sqrt

import numpy as np

from osgeo import ogr, osr

from qgis.core import QgsMapLayerRegistry, QgsMapLayer, QGis, QgsCoordinateTransform, QgsPoint, QgsRaster
from qgis.gui import *

from PyQt4.QtCore import *
from PyQt4.QtGui import *

from .errors import VectorIOException
from ..gsf.geometry import Point


def get_on_the_fly_projection(canvas):

    on_the_fly_projection = True if canvas.hasCrsTransformEnabled() else False

    if on_the_fly_projection:
        project_crs = canvas.mapRenderer().destinationCrs()
    else:
        project_crs = None

    return on_the_fly_projection, project_crs


def vector_type(layer):

    if not layer.type() == QgsMapLayer.VectorLayer:
        raise VectorIOException("Layer is not vector")

    if layer.geometryType() == QGis.Point:
        return "point"
    elif layer.geometryType() == QGis.Line:
        return "line"
    elif layer.geometryType() == QGis.Polygon:
        return "polygon"
    else:
        raise VectorIOException("Unknown vector type")


def loaded_layers():

    return QgsMapLayerRegistry.instance().mapLayers().values()


def loaded_vector_layers():

    return filter(lambda layer: layer.type() == QgsMapLayer.VectorLayer,
                  loaded_layers())


def loaded_polygon_layers():

    return filter(lambda layer: layer.geometryType() == QGis.Polygon,
                  loaded_vector_layers())


def loaded_line_layers():

    return filter(lambda layer: layer.geometryType() == QGis.Line,
                  loaded_vector_layers())


def loaded_point_layers():

    return filter(lambda layer: layer.geometryType() == QGis.Point,
                  loaded_vector_layers())


def loaded_raster_layers():

    return filter(lambda layer: layer.type() == QgsMapLayer.RasterLayer,
                  loaded_layers())


def loaded_monoband_raster_layers():

    return filter(lambda layer: layer.bandCount() == 1,
                  loaded_raster_layers())


def pt_geoms_attrs(pt_layer, field_list=None):

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


def line_geoms_attrs(line_layer, field_list=None):

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


def line_geoms_with_id(line_layer, curr_field_ndx):

    lines = []
    progress_ids = []

    if line_layer.selectedFeatureCount() > 0:
        features = line_layer.selectedFeatures()
    else:
        features = line_layer.getFeatures()

    dummy_progressive = 0
    for feature in features:
        try:
            progress_ids.append(int(feature[curr_field_ndx]))
        except:
            dummy_progressive += 1
            progress_ids.append(dummy_progressive)

        geom = feature.geometry()
        if geom.isMultipart():
            lines.append(
                ('multiline', multipolyline_to_xytuple_list2(geom.asMultiPolyline())))  # typedef QVector<QgsPolyline>
            # now is a list of list of (x,y) tuples
        else:
            lines.append(('line', polyline_to_xytuple_list(geom.asPolyline())))  # typedef QVector<QgsPoint>

    return lines, progress_ids


def polyline_to_xytuple_list(qgsline):

    assert len(qgsline) > 0
    return [(qgspoint.x(), qgspoint.y()) for qgspoint in qgsline]


def multipolyline_to_xytuple_list2(qgspolyline):

    return [polyline_to_xytuple_list(qgsline) for qgsline in qgspolyline]


def field_values(layer, curr_field_ndx):

    values = []

    if layer.selectedFeatureCount() > 0:
        features = layer.selectedFeatures()
    else:
        features = layer.getFeatures()

    for feature in features:
        values.append(feature.attributes()[curr_field_ndx])

    return values


def vect_attrs(layer, field_list):

    if layer.selectedFeatureCount() > 0:
        features = layer.selectedFeatures()
    else:
        features = layer.getFeatures()

    provider = layer.dataProvider()
    field_indices = [provider.fieldNameIndex(field_name) for field_name in field_list]

    # retrieve (selected) attributes features
    data_list = []
    for feature in features:
        attrs = feature.fields().toList()
        data_list.append([feature.attribute(attrs[field_ndx].name()) for field_ndx in field_indices])

    return data_list


def raster_qgis_params(raster_layer):

    name = raster_layer.name()

    rows = raster_layer.height()
    cols = raster_layer.width()

    extent = raster_layer.extent()

    xMin = extent.xMinimum()
    xMax = extent.xMaximum()
    yMin = extent.yMinimum()
    yMax = extent.yMaximum()

    cellsizeEW = (xMax - xMin) / float(cols)
    cellsizeNS = (yMax - yMin) / float(rows)

    # TODO: get real no data value from QGIS
    if raster_layer.dataProvider().srcHasNoDataValue(1):
        nodatavalue = raster_layer.dataProvider().srcNoDataValue(1)
    else:
        nodatavalue = np.nan

    try:
        crs = raster_layer.crs()
    except:
        crs = None

    return name, cellsizeEW, cellsizeNS, rows, cols, xMin, xMax, yMin, yMax, nodatavalue, crs


def qgs_pt(x, y):

    return QgsPoint(x, y)


def project_qgs_point(qgsPt, srcCrs, destCrs):

    return QgsCoordinateTransform(srcCrs, destCrs).transform(qgsPt)


def project_point(pt, srcCrs, destCrs):

    qgs_pt = QgsPoint(pt.x, pt.y)
    proj_qgs_pt = project_qgs_point(qgs_pt, srcCrs, destCrs)
    proj_x, proj_y = proj_qgs_pt.x(), proj_qgs_pt.y()

    return Point(proj_x, proj_y)


def project_xy_list(src_crs_xy_list, srcCrs, destCrs):

    pt_list_dest_crs = []
    for x, y in src_crs_xy_list.pts:
        srcPt = QgsPoint(x, y)
        destPt = project_qgs_point(srcPt, srcCrs, destCrs)
        pt_list_dest_crs = pt_list_dest_crs.append([destPt.x(), destPt.y()])

    return pt_list_dest_crs


def qcolor2rgbmpl(qcolor):

    red = qcolor.red() / 255.0
    green = qcolor.green() / 255.0
    blue = qcolor.blue() / 255.0
    return red, green, blue


"""
Modified from: profiletool, script: tools/ptmaptool.py

#-----------------------------------------------------------
# 
# Profile
# Copyright (C) 2008  Borys Jurgiel
# Copyright (C) 2012  Patrice Verchere
#-----------------------------------------------------------
# 
# licensed under the terms of GNU GPL 2
# 
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License along
# with this program; if not, print to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#
#---------------------------------------------------------------------
"""


class PointMapToolEmitPoint(QgsMapToolEmitPoint):

    def __init__(self, canvas, button):

        super(PointMapToolEmitPoint, self).__init__(canvas)
        self.canvas = canvas
        self.cursor = QCursor(Qt.CrossCursor)
        self.button = button

    def setCursor(self, cursor):

        self.cursor = QCursor(cursor)


class MapDigitizeTool(QgsMapTool):

    def __init__(self, canvas):

        QgsMapTool.__init__(self, canvas)
        self.canvas = canvas
        self.cursor = QCursor(Qt.CrossCursor)

    def canvasMoveEvent(self, event):

        self.emit(SIGNAL("moved"), {'x': event.pos().x(), 'y': event.pos().y()})

    def canvasReleaseEvent(self, event):

        if event.button() == Qt.RightButton:
            button_type = "rightClicked"
        elif event.button() == Qt.LeftButton:
            button_type = "leftClicked"
        else:
            return

        self.emit(SIGNAL(button_type), {'x': event.pos().x(), 'y': event.pos().y()})

    def canvasDoubleClickEvent(self, event):

        self.emit(SIGNAL("doubleClicked"), {'x': event.pos().x(), 'y': event.pos().y()})

    def activate(self):

        QgsMapTool.activate(self)
        self.canvas.setCursor(self.cursor)

    def deactivate(self):

        QgsMapTool.deactivate(self)

    def isZoomTool(self):

        return False

    def setCursor(self, cursor):

        self.cursor = QCursor(cursor)


class QGisRasterParameters(object):

    def __init__(self, name, cellsizeEW, cellsizeNS, rows, cols, xMin, xMax, yMin, yMax, nodatavalue, crs):

        self.name = name
        self.cellsizeEW = cellsizeEW
        self.cellsizeNS = cellsizeNS
        self.rows = rows
        self.cols = cols
        self.xMin = xMin
        self.xMax = xMax
        self.yMin = yMin
        self.yMax = yMax
        self.nodatavalue = nodatavalue
        self.crs = crs

    def point_in_dem_area(self, point):
        """
        Check that a point is within or on the boundary of the grid area.
        Assume grid has no rotation.
        
        :param point: qProf.gsf.geometry.Point
        :return: bool
        """

        if self.xMin <= point.x <= self.xMax and \
                                self.yMin <= point.y <= self.yMax:
            return True
        else:
            return False

    def point_in_interpolation_area(self, point):
        """
        Check that a point is within or on the boundary of the area defined by
        the extreme cell center values.
        Assume grid has no rotation.
        
        :param point: qProf.gsf.geometry.Point
        :return: bool
        """

        if self.xMin + self.cellsizeEW / 2.0 <= point.x <= self.xMax - self.cellsizeEW / 2.0 and \
           self.yMin + self.cellsizeNS / 2.0 <= point.y <= self.yMax - self.cellsizeNS / 2.0:
            return True
        else:
            return False

    def geogr2raster(self, point):
        """
        Convert from geographic to raster-based coordinates.
        Assume grid has no rotation.

        :param point: qProf.gsf.geometry.Point
        :return: dict
        """

        x = (point.x - (self.xMin + self.cellsizeEW / 2.0)) / self.cellsizeEW
        y = (point.y - (self.yMin + self.cellsizeNS / 2.0)) / self.cellsizeNS

        return dict(x=x, y=y)

    def raster2geogr(self, array_dict):
        """
        Convert from raster-based to geographic coordinates.
        Assume grid has no rotation.

        :param array_dict: dict
        :return: qProf.gsf.geometry.Point instance
        """

        assert 'x' in array_dict
        assert 'y' in array_dict

        x = self.xMin + (array_dict['x'] + 0.5) * self.cellsizeEW
        y = self.yMin + (array_dict['y'] + 0.5) * self.cellsizeNS

        return Point(x, y)


def get_z(dem_layer, point):

    identification = dem_layer.dataProvider().identify(QgsPoint(point.x, point.y), QgsRaster.IdentifyFormatValue)
    if not identification.isValid():
        return np.nan
    else:
        try:
            result_map = identification.results()
            return float(result_map[1])
        except:
            return np.nan


def interpolate_bilinear(dem, qrpDemParams, point):
    """        
    :param dem: qgis._core.QgsRasterLayer
    :param qrpDemParams: qProf.gis_utils.qgs_tools.QGisRasterParameters
    :param point: qProf.gis_utils.features.Point
    :return: float
    """

    dArrayCoords = qrpDemParams.geogr2raster(point)

    floor_x_raster = floor(dArrayCoords["x"])
    ceil_x_raster = ceil(dArrayCoords["x"])
    floor_y_raster = floor(dArrayCoords["y"])
    ceil_y_raster = ceil(dArrayCoords["y"])

    # bottom-left center
    p1 = qrpDemParams.raster2geogr(dict(x=floor_x_raster,
                                        y=floor_y_raster))
    # bottom-right center
    p2 = qrpDemParams.raster2geogr(dict(x=ceil_x_raster,
                                        y=floor_y_raster))
    # top-left center
    p3 = qrpDemParams.raster2geogr(dict(x=floor_x_raster,
                                        y=ceil_y_raster))
    # top-right center
    p4 = qrpDemParams.raster2geogr(dict(x=ceil_x_raster,
                                        y=ceil_y_raster))

    z1 = get_z(dem, p1)
    z2 = get_z(dem, p2)
    z3 = get_z(dem, p3)
    z4 = get_z(dem, p4)

    delta_x = point.x - p1.x
    delta_y = point.y - p1.y

    z_x_a = z1 + (z2 - z1) * delta_x / qrpDemParams.cellsizeEW
    z_x_b = z3 + (z4 - z3) * delta_x / qrpDemParams.cellsizeEW

    return z_x_a + (z_x_b - z_x_a) * delta_y / qrpDemParams.cellsizeNS


def interpolate_z(dem, dem_params, point):
    """
        dem_params: type qProf.gis_utils.qgs_tools.QGisRasterParameters
        point: type qProf.gis_utils.features.Point
    """

    if dem_params.point_in_interpolation_area(point):
        return interpolate_bilinear(dem, dem_params, point)
    elif dem_params.point_in_dem_area(point):
        return get_z(dem, point)
    else:
        return np.nan


def get_zs_from_dem(struct_pts_2d, demObj):

    z_list = []
    for point_2d in struct_pts_2d:
        interp_z = interpolate_z(demObj.layer, demObj.params, point_2d)
        z_list.append(interp_z)

    return z_list


def xy_from_canvas(canvas, position):

    mapPos = canvas.getCoordinateTransform().toMapCoordinates(position["x"], position["y"])

    return mapPos.x(), mapPos.y()


def get_prjcrs_as_proj4str(canvas):

    # get project CRS information
    hasOTFP, project_crs = get_on_the_fly_projection(canvas)
    if hasOTFP:
        proj4_str = str(project_crs.toProj4())
        project_crs_osr = osr.SpatialReference()
        project_crs_osr.ImportFromProj4(proj4_str)
        return project_crs_osr
    else:
        return None