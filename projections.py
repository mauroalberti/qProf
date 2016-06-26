
from geosurf.spatial import CartesianPoint2DT, CartesianLine2DT

from geosurf.qgs_tools import project_qgs_point

from qgis.core import QgsPoint

   
def line2d_change_crs(source_line, source_crs, destination_crs):
    
    destination_line = CartesianLine2DT()
    for pt in source_line.pts:
        source_qgspt = QgsPoint(pt.p_x, pt.p_y)
        destination_qgspt = project_qgs_point(source_qgspt, source_crs, destination_crs)
        destination_line = destination_line.add_pt(CartesianPoint2DT(destination_qgspt.x(), destination_qgspt.y()))
        
    return destination_line

