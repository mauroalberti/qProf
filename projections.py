
from geosurf.spatial import CartesianPoint2DT, Line2DT

from geosurf.qgs_tools import project_qgs_point

from qgis.core import QgsPoint

   
def line2d_change_crs(source_line, source_crs, destination_crs):
    
    destination_line = Line2DT()
    for pt in source_line._pts:
        source_pt = QgsPoint(pt.x(), pt.y())
        destination_pt = project_qgs_point(source_pt, source_crs, destination_crs)
        destination_line = destination_line.add_pt(CartesianPoint2DT(destination_pt.x(), destination_pt.y()))
        
    return destination_line

