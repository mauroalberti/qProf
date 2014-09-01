
import numpy as np

from PyQt4.QtCore import *
from PyQt4.QtGui import *

from qgis.core import *
from qgis.core import QgsMapLayerRegistry, QgsMapLayer, QGis, QgsFeature, \
                      QgsCoordinateTransform, QgsPoint
from qgis.gui import *

from .spatial import Point_2D, Line_2D
from .errors import VectorIOException


def vector_type( layer ):
    
    if not layer.type() == QgsMapLayer.VectorLayer:
        raise VectorIOException, "Layer is not vector"
    
    if layer.geometryType() == QGis.Point:
        return "point"
    elif layer.geometryType() == QGis.Line:
        return "line"        
    elif layer.geometryType() == QGis.Polygon:
        return "polygon"
    else: 
        raise VectorIOException, "Unknown vector type" 
       
    
def loaded_layers():
    
    return QgsMapLayerRegistry.instance().mapLayers().values()

    
def loaded_vector_layers():
 
    return filter( lambda layer: layer.type() == QgsMapLayer.VectorLayer, 
                   loaded_layers() )
    
            
def loaded_line_layers():        
    
    return filter( lambda layer: layer.geometryType() == QGis.Line, 
                          loaded_vector_layers() )


def loaded_point_layers():

    return filter( lambda layer: layer.geometryType() == QGis.Point, 
                          loaded_vector_layers() )
    
 
def loaded_raster_layers( ):        
          
    return filter( lambda layer: layer.type() == QgsMapLayer.RasterLayer, 
                          loaded_layers() )


def loaded_monoband_raster_layers( ):        
          
    return filter( lambda layer: layer.bandCount() == 1, 
                   loaded_raster_layers() )
       
           
def pt_geoms_attrs( pt_layer, field_list = [] ):
    
    if pt_layer.selectedFeatureCount() > 0:
        features = pt_layer.selectedFeatures()
    else:
        features = pt_layer.getFeatures() 
    
    provider = pt_layer.dataProvider()    
    field_indices = [ provider.fieldNameIndex( field_name ) for field_name in field_list ]

    # retrieve selected features with their geometry and relevant attributes
    rec_list = [] 
    for feature in features:
             
        # fetch point geometry
        pt = feature.geometry().asPoint()

        attrs = feature.fields().toList() 

        # creates feature attribute list
        feat_list = [ pt.x(), pt.y() ]
        for field_ndx in field_indices:
            feat_list.append( str( feature.attribute( attrs[ field_ndx ].name() ) ) )

        # add to result list
        rec_list.append( feat_list )
        
    return rec_list


def line_geoms_attrs( line_layer, field_list = [] ):
    """
    return: lines, list
    """
    
    lines = []
    
    if line_layer.selectedFeatureCount() > 0:
        features = line_layer.selectedFeatures()
    else:
        features = line_layer.getFeatures()

    provider = line_layer.dataProvider() 
    field_indices = [ provider.fieldNameIndex( field_name ) for field_name in field_list ]
                
    for feature in features:
        geom = feature.geometry()
        if geom.isMultipart():
            rec_geom = multipolyline_to_xytuple_list2( geom.asMultiPolyline() )
        else:
            rec_geom = [ polyline_to_xytuple_list( geom.asPolyline() ) ]
            
        attrs = feature.fields().toList()
        rec_data = [ str( feature.attribute( attrs[ field_ndx ].name() ) ) for field_ndx in field_indices ]
            
        lines.append( [ rec_geom, rec_data ] )
            
    return lines
           
       
def line_geoms_with_id( line_layer, curr_field_ndx ):    
        
    lines = []
    progress_ids = [] 
    dummy_progressive = 0 
      
    line_iterator = line_layer.getFeatures()
   
    for feature in line_iterator:
        try:
            progress_ids.append( int( feature[ curr_field_ndx ] ) )
        except:
            dummy_progressive += 1
            progress_ids.append( dummy_progressive )
             
        geom = feature.geometry()         
        if geom.isMultipart():
            lines.append( 'multiline', multipolyline_to_xytuple_list2( geom.asMultiPolyline() ) ) # typedef QVector<QgsPolyline>
            # now is a list of list of (x,y) tuples
        else:           
            lines.append( ( 'line', polyline_to_xytuple_list( geom.asPolyline() ) ) ) # typedef QVector<QgsPoint>
                         
    return lines, progress_ids
              
                   
def polyline_to_xytuple_list( qgsline):
    
    assert len( qgsline ) > 0    
    return [ ( qgspoint.x(), qgspoint.y() ) for qgspoint in qgsline ]


def multipolyline_to_xytuple_list2( qgspolyline ):
    
    return [ polyline_to_xytuple_list( qgsline) for qgsline in qgspolyline ] 


def field_values( layer, curr_field_ndx ): 
    
    values = []
    iterator = layer.getFeatures()
    
    for feature in iterator:
        values.append( feature.attributes()[ curr_field_ndx ] ) 
            
    return values
    
    
def vect_attrs( layer, field_list):
    
    if layer.selectedFeatureCount() > 0:
        features = layer.selectedFeatures()
    else:
        features = layer.getFeatures()
        
    provider = layer.dataProvider()   
    field_indices = [ provider.fieldNameIndex( field_name ) for field_name in field_list ]

    # retrieve (selected) attributes features
    data_list = [] 
    for feature in features:        
        attrs = feature.fields().toList()     
        data_list.append( [ feature.attribute( attrs[ field_ndx ].name() ) for field_ndx in field_indices ] )
        
    return data_list    
    
    
def raster_qgis_params( raster_layer ):
    
    name = raster_layer.name()
                  
    rows = raster_layer.height()
    cols = raster_layer.width()
    
    extent = raster_layer.extent()
    
    xMin = extent.xMinimum()
    xMax = extent.xMaximum()        
    yMin = extent.yMinimum()
    yMax = extent.yMaximum()
        
    cellsizeEW = (xMax-xMin) / float(cols)
    cellsizeNS = (yMax-yMin) / float(rows)
    
    #TODO: get real no data value from QGIS
    if raster_layer.dataProvider().srcHasNoDataValue(1):
        nodatavalue = raster_layer.dataProvider().srcNoDataValue ( 1 )
    else:
        nodatavalue = np.nan
    
    try:
        crs = raster_layer.crs()
    except:
        crs = None
    
    return name, cellsizeEW, cellsizeNS, rows, cols, xMin, xMax, yMin, yMax, nodatavalue, crs    


def qgs_point( x, y ):
    
    return QgsPoint(x, y)
        

def project_qgs_point( qgsPt, srcCrs, destCrs ):
    
    return QgsCoordinateTransform( srcCrs, destCrs ).transform( qgsPt )

    
def project_line_2d( srcLine, srcCrs, destCrs ):
    
    destLine = Line_2D()    
    for pt in srcLine._pts:        
        srcPt = QgsPoint (pt._x, pt._y)
        destPt = project_qgs_point( srcPt, srcCrs, destCrs )
        destLine = destLine.add_pt( Point_2D( destPt.x(), destPt.y() ) )
        
    return destLine
        
        

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

class PointMapToolEmitPoint( QgsMapToolEmitPoint ):

    def __init__( self, canvas, button ):
        
        super( PointMapToolEmitPoint, self ).__init__( canvas )
        self.canvas = canvas
        self.cursor = QCursor( Qt.CrossCursor )
        self.button = button


    def setCursor( self, cursor ):
        
        self.cursor = QCursor( cursor )
        
                

class MapDigitizeTool( QgsMapTool ):
    

    def __init__(self, canvas ):
        
        QgsMapTool.__init__(self,canvas)
        self.canvas = canvas
        self.cursor = QCursor( Qt.CrossCursor )
                

    def canvasMoveEvent( self, event ):
        
        self.emit( SIGNAL("moved"), {'x': event.pos().x(), 'y': event.pos().y()} )


    def canvasReleaseEvent( self, event ):
        
        if event.button() == Qt.RightButton:
            button_type =  "rightClicked"
        elif event.button() == Qt.LeftButton:
            button_type =  "leftClicked"
        else:
            return      
        
        self.emit( SIGNAL( button_type ), {'x': event.pos().x(), 'y': event.pos().y()} )


    def canvasDoubleClickEvent( self, event ):
        
        self.emit( SIGNAL("doubleClicked"), {'x': event.pos().x(), 'y': event.pos().y()} )


    def activate( self ):
        
        QgsMapTool.activate( self )
        self.canvas.setCursor(self.cursor)


    def deactivate( self ):
        
        #self.emit( SIGNAL("deactivate") )
        QgsMapTool.deactivate( self )


    def isZoomTool( self ):
        
        return False


    def setCursor( self, cursor ):
        
        self.cursor = QCursor(cursor)



    
    
    
    
