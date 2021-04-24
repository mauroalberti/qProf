
from typing import Tuple, List, Union

import numbers
from copy import deepcopy

from math import ceil, floor

from collections import namedtuple

import numpy as np

from qgis.core import QgsRaster, QgsPointXY, QgsRasterLayer


from .points import *  # distance_projected_pts, calculate_pts_in_projection
from .project import *  # projectCrs

from ...geometries.shapes.space2d import Point2D
from ...geometries.shapes.space3d import Point3D

raster_parameters_fields = [
    'name',
    'cellsizeEW',
    'cellsizeNS',
    'rows',
    'cols',
    'xMin',
    'xMax',
    'yMin',
    'yMax',
    'nodatavalue',
    'crs'
]

RasterParameters = namedtuple(
    'RasterParameters',
    raster_parameters_fields
)


class QGisRasterParameters(object):

    def __init__(self,
                 name,
                 cellsizeEW,
                 cellsizeNS,
                 rows,
                 cols,
                 xMin,
                 xMax,
                 yMin,
                 yMax,
                 nodatavalue,
                 crs
                 ):

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

    def point_in_dem_area(self,
                          point
                          ):
        """
        Check that a point is within or on the boundary of the grid area.
        Assume grid has no rotation.

        :return: bool
        """

        if self.xMin <= point.x <= self.xMax and \
                self.yMin <= point.y <= self.yMax:
            return True
        else:
            return False

    def point_in_interpolation_area(self,
                                    point
                                    ):
        """
        Check that a point is within or on the boundary of the area defined by
        the extreme cell center values.
        Assume grid has no rotation.

        :return: bool
        """

        if self.xMin + self.cellsizeEW / 2.0 <= point.x <= self.xMax - self.cellsizeEW / 2.0 and \
                self.yMin + self.cellsizeNS / 2.0 <= point.y <= self.yMax - self.cellsizeNS / 2.0:
            return True
        else:
            return False

    def geogr2raster(self,
                     point
                     ):
        """
        Convert from geographic to raster-based coordinates.
        Assume grid has no rotation.

        :return: dict
        """

        x = (point.x - (self.xMin + self.cellsizeEW / 2.0)) / self.cellsizeEW
        y = (point.y - (self.yMin + self.cellsizeNS / 2.0)) / self.cellsizeNS

        return dict(x=x, y=y)

    def raster2geogr(self,
         array_dict
         ) -> Point2D:
        """
        Convert from raster-based to geographic coordinates.
        Assume grid has no rotation.

        :param array_dict: dict
        :return: the point in geographic planar coordinates
        """

        assert 'x' in array_dict
        assert 'y' in array_dict

        x = self.xMin + (array_dict['x'] + 0.5) * self.cellsizeEW
        y = self.yMin + (array_dict['y'] + 0.5) * self.cellsizeNS

        return Point2D(x, y)


def get_z(
    dem_layer,
    point
):

    identification = dem_layer.dataProvider().identify(QgsPointXY(point.x, point.y), QgsRaster.IdentifyFormatValue)
    if not identification.isValid():
        return np.nan
    else:
        try:
            result_map = identification.results()
            return float(result_map[1])
        except:
            return np.nan


def get_zs_from_dem(
    struct_pts_2d,
    demObj
):

    z_list = []
    for point_2d in struct_pts_2d:
        interp_z = interpolate_z(demObj.layer, demObj.params, point_2d)
        z_list.append(interp_z)

    return z_list


def interpolate_bilinear(
    dem: QgsRasterLayer,
    qrpDemParams: QGisRasterParameters,
    point: Point2D
) -> numbers.Real:
    """
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


def interpolate_z(
    qgs_raster_layer: QgsRasterLayer,
    qgis_raster_parameters: QGisRasterParameters,
    point2d: Point2D
):
    """
    """

    if qgis_raster_parameters.point_in_interpolation_area(point2d):
        return interpolate_bilinear(qgs_raster_layer, qgis_raster_parameters, point2d)
    elif qgis_raster_parameters.point_in_dem_area(point2d):
        return get_z(qgs_raster_layer, point2d)
    else:
        return np.nan


def try_raster_qgis_params(
        raster_layer
) -> Tuple[bool, Union[str, Tuple]]:

    try:

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
        if raster_layer.dataProvider().sourceHasNoDataValue(1):
            nodatavalue = raster_layer.dataProvider().sourceNoDataValue(1)
        else:
            nodatavalue = np.nan

        try:
            crs = raster_layer.crs()
        except:
            crs = None

        return True, (name, cellsizeEW, cellsizeNS, rows, cols, xMin, xMax, yMin, yMax, nodatavalue, crs)

    except Exception as e:

        return False, str(e)


def raster_qgis_params(
        raster_layer
):
    """
    Deprecated: use 'try_extract_raster_qgis_params'
    """

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
    if raster_layer.dataProvider().sourceHasNoDataValue(1):
        nodatavalue = raster_layer.dataProvider().sourceNoDataValue(1)
    else:
        nodatavalue = np.nan

    try:
        crs = raster_layer.crs()
    except:
        crs = None

    return name, cellsizeEW, cellsizeNS, rows, cols, xMin, xMax, yMin, yMax, nodatavalue, crs


def get_dem_resolution_in_prj_crs(
        dem,
        dem_params,
        prj_crs
):

    cellsizeEW, cellsizeNS = dem_params.cellsizeEW, dem_params.cellsizeNS
    xMin, yMin = dem_params.xMin, dem_params.yMin

    if dem.crs() != prj_crs:
        cellsizeEW_prj_crs = distance_projected_pts(
            xMin,
            yMin,
            cellsizeEW,
            0,
            dem.crs(),
            prj_crs
        )
        cellsizeNS_prj_crs = distance_projected_pts(
            xMin,
            yMin,
            0,
            cellsizeNS,
            dem.crs(),
            prj_crs
        )
    else:
        cellsizeEW_prj_crs = cellsizeEW
        cellsizeNS_prj_crs = cellsizeNS

    return 0.5 * (cellsizeEW_prj_crs + cellsizeNS_prj_crs)


def get_min_dem_resolution(
    selected_dems: List,
    selected_dem_parameters: List
) -> numbers.Real:

    dem_resolutions_prj_crs_list = []

    for dem, dem_params in zip(selected_dems, selected_dem_parameters):
        dem_resolutions_prj_crs_list.append(
            get_dem_resolution_in_prj_crs(
                dem,
                dem_params,
                projectCrs())
        )

    min_dem_resolution = np.nanmin(dem_resolutions_prj_crs_list)

    if min_dem_resolution > 1:
        sample_distance = round(min_dem_resolution)
    else:
        sample_distance = min_dem_resolution

    return sample_distance


def calculate_projected_3d_pts(
    struct_pts,
    structural_pts_crs,
    demObj
):

    demCrs = demObj.params.crs

    # check if on-the-fly-projection is set on
    project_crs = projectCrs()

    # set points in the project crs
    if structural_pts_crs != project_crs:
        struct_pts_in_prj_crs = calculate_pts_in_projection(struct_pts, structural_pts_crs, project_crs)
    else:
        struct_pts_in_prj_crs = deepcopy(struct_pts)

        # project the source points from point layer crs to DEM crs
    # if the two crs are different
    if structural_pts_crs != demCrs:
        struct_pts_in_dem_crs = calculate_pts_in_projection(struct_pts, structural_pts_crs, demCrs)
    else:
        struct_pts_in_dem_crs = deepcopy(struct_pts)

        # - 3D structural points, with x, y, and z extracted from the current DEM
    struct_pts_z = get_zs_from_dem(struct_pts_in_dem_crs, demObj)

    assert len(struct_pts_in_prj_crs) == len(struct_pts_z)

    return [Point3D(pt.x, pt.y, z) for (pt, z) in zip(struct_pts_in_prj_crs, struct_pts_z)]

