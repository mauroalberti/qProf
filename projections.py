
from geosurf.spatial import Point_2D, Line_2D

from geosurf.qgs_tools import project_qgs_point

from qgis.core import QgsPoint

   
def project_line_2d( srcLine, srcCrs, destCrs ):
    
    destLine = Line_2D()    
    for pt in srcLine._pts:        
        srcPt = QgsPoint (pt._x, pt._y)
        destPt = project_qgs_point( srcPt, srcCrs, destCrs )
        destLine = destLine.add_pt( Point_2D( destPt.x(), destPt.y() ) )
        
    return destLine
        
    
        