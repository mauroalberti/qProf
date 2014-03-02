
from geosurf.spatial import Point, Line

from qgs_tools.tools import project_point

from qgis.core import QgsPoint

   
def project_line( srcLine, srcCrs, destCrs ):
    
    destLine = Line()
    
    for pt in srcLine.points:        
        srcPt = QgsPoint (pt.x, pt.y)
        destPt = project_point( srcPt, srcCrs, destCrs )
        destLine = destLine.add_point( Point( destPt.x(), destPt.y() ) )
        
    return destLine
        
    
        