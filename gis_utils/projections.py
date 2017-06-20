
from qgis.core import QgsPoint

from .features import CartesianLine2DT
from ..gsf.geometry import Point
from .qgs_tools import project_qgs_point


def line2d_change_crs(source_line, source_crs, destination_crs):
    
    destination_line = CartesianLine2DT()
    for pt in source_line.pts:
        source_qgspt = QgsPoint(pt.x, pt.y)
        destination_qgspt = project_qgs_point(source_qgspt, source_crs, destination_crs)
        destination_line = destination_line.add_pt(Point(destination_qgspt.x(), destination_qgspt.y()))
        
    return destination_line

