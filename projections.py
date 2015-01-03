
from geosurf.spatial import Point2D, Line2D

from geosurf.qgs_tools import project_qgs_point

from qgis.core import QgsPoint

   
def project_line_2d( srcLine, srcCrs, destCrs ):
    
    destLine = Line2D()    
    for pt in srcLine._pts:        
        srcPt = QgsPoint (pt._x, pt._y)
        destPt = project_qgs_point( srcPt, srcCrs, destCrs )
        destLine = destLine.add_pt( Point2D( destPt.x(), destPt.y() ) )
        
    return destLine
        
    
        