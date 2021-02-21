
from typing import Tuple

from osgeo import osr

from qgis.core import *


def projectCrs():

    return QgsProject.instance().crs()


def projectHasValidCrs() -> bool:

    return QgsProject.instance().crs().isValid()


def projectHasGeographicCrs() -> bool:

    return QgsProject.instance().crs().isGeographic()


def check_project_planar_crs() -> Tuple[bool, str]:

    project = QgsProject.instance()

    if project.count() == 0:
        msg = "Is a project open or has layers loaded?"
        return False, msg

    curr_proj_crs = project.crs()

    if not curr_proj_crs.isValid():
        msg = "Current project crs is not valid.\nPlease apply a valid crs to the current project."
        return False, msg

    if curr_proj_crs.isGeographic():
        msg = "Current project crs is geographic.\nPlease apply a planar crs to the current project."
        return False, msg

    msg = "Project open and with planar crs defined."
    return True, msg


def loaded_layers():

    return list(QgsProject.instance().mapLayers().values())


def loaded_vector_layers():

    return [layer for layer in loaded_layers() if layer.type() == QgsMapLayer.VectorLayer]


def loaded_polygon_layers():

    return [layer for layer in loaded_vector_layers() if layer.geometryType() == QgsWkbTypes.PolygonGeometry]


def loaded_line_layers():

    return [layer for layer in loaded_vector_layers() if layer.geometryType() == QgsWkbTypes.LineGeometry]


def loaded_point_layers():

    return [layer for layer in loaded_vector_layers() if layer.geometryType() == QgsWkbTypes.PointGeometry]


def loaded_raster_layers():

    return [layer for layer in loaded_layers() if layer.type() == QgsMapLayer.RasterLayer]


def loaded_monoband_raster_layers():

    return [layer for layer in loaded_raster_layers() if layer.bandCount() == 1]


def proj4str():

    project_crs = QgsProject.instance().crs()
    proj4_str = str(project_crs.toProj4())
    project_crs_osr = osr.SpatialReference()
    project_crs_osr.ImportFromProj4(proj4_str)

    return project_crs_osr
