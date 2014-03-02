
import numpy as np

from qgis.core import QgsMapLayerRegistry, QgsMapLayer, QGis, QgsFeature, \
                      QgsCoordinateTransform, QgsPoint


def get_map_layers_from_qgis():
    
    return QgsMapLayerRegistry.instance().mapLayers().values()
    

def get_current_raster_layers( ):        
          
    return filter( lambda layer: layer.type() == QgsMapLayer.RasterLayer, 
                          get_map_layers_from_qgis() )


def get_current_singleband_raster_layers( ):        
          
    return filter( lambda layer: layer.type() == QgsMapLayer.RasterLayer and layer.bandCount() == 1, 
                          get_map_layers_from_qgis() )
    
    
def get_current_vector_layers():
 
    return filter( lambda layer: layer.type() == QgsMapLayer.VectorLayer, 
                          get_map_layers_from_qgis() )
    
            
def get_current_line_layers():        
    
    return filter( lambda layer: layer.geometryType() == QGis.Line, 
                          get_current_vector_layers() )


def get_current_point_layers():

    return filter( lambda layer: layer.geometryType() == QGis.Point, 
                          get_current_vector_layers() )
    
    
    
def get_selected_features_attr( pt_layer, field_list ):
    
    num_selected_features = pt_layer.selectedFeatureCount() 
    if num_selected_features == 0:
        return False, "No geological layer feature selected" 
    
    provider = pt_layer.dataProvider()
    
    field_indices = [ provider.fieldNameIndex( field_name ) for field_name in field_list ]

    # retrieve selected features with their geometry and relevant attributes
    data_list = [] 
    for feat in pt_layer.selectedFeatures():
             
        # fetch point geometry
        pt = feat.geometry().asPoint()
        
        # fetch map of attributes as a dictionary: key = field index, value = QgsFeatureAttribute
        attrs = feat.fields().toList() 

        # creates feature attribute list
        feat_list = [ pt.x(), pt.y() ]
        for field_ndx in field_indices:
            feat_list.append( str( feat.attribute( attrs[ field_ndx ].name() ) ) )

        # add to result list
        data_list.append( feat_list )
        
    return True, data_list


def read_vector_line_qgs( line_layer, curr_field_ndx ):
    
    
    line_crs = line_layer.crs()
    
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
            lines.append( 'multiline', convert_qgspolyline_to_lines_list( geom.asMultiPolyline() ) ) # typedef QVector<QgsPolyline>
            # now is a list of list of (x,y,z) tuples
        else:           
            lines.append( ( 'line', convert_qgsline( geom.asPolyline() ) ) ) # typedef QVector<QgsPoint>
                         
    return lines, progress_ids, line_crs
              
                   
def convert_qgsline( qgsline):
    
    return [ ( qgspoint.x(), qgspoint.y(), 0.0 ) for qgspoint in qgsline ]


def convert_qgspolyline_to_lines_list( qgspolyline ):
    
    return [ convert_qgsline( qgsline) for qgsline in qgspolyline ]                   
                
    
def get_raster_params_via_qgis( raster_layer ):
                  
    rows = raster_layer.height()
    cols = raster_layer.width()
    
    extent = raster_layer.extent()
    
    xMin = extent.xMinimum()
    xMax = extent.xMaximum()        
    yMin = extent.yMinimum()
    yMax = extent.yMaximum()
        
    cellsizeEW = (xMax-xMin) / float(cols)
    cellsizeNS = (yMax-yMin) / float(rows)
    
    nodatavalue = None
    
    try:
        crs = raster_layer.crs()
    except:
        crs = None
    
    return cellsizeEW, cellsizeNS, rows, cols, xMin, xMax, yMin, yMax, nodatavalue, crs    


def make_qgs_point( x, y):
    
    return QgsPoint(x, y)
        

def project_point( qgsPt, srcCrs, destCrs ):
    
    return QgsCoordinateTransform( srcCrs, destCrs ).transform( qgsPt )

   
    
    
    
    
    
