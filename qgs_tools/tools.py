
import numpy as np

from qgis.core import QgsMapLayerRegistry, QgsMapLayer, QGis, QgsFeature, \
                      QgsCoordinateTransform, QgsPoint


def get_map_layers_from_qgis():
    
    return QgsMapLayerRegistry.instance().mapLayers().values()

    
def get_current_vector_layers():
 
    return filter( lambda layer: layer.type() == QgsMapLayer.VectorLayer, 
                   get_map_layers_from_qgis() )
    
            
def get_current_line_layers():        
    
    return filter( lambda layer: layer.geometryType() == QGis.Line, 
                          get_current_vector_layers() )


def get_current_point_layers():

    return filter( lambda layer: layer.geometryType() == QGis.Point, 
                          get_current_vector_layers() )
    
 
def get_current_raster_layers( ):        
          
    return filter( lambda layer: layer.type() == QgsMapLayer.RasterLayer, 
                          get_map_layers_from_qgis() )


def get_current_singleband_raster_layers( ):        
          
    return filter( lambda layer: layer.bandCount() == 1, 
                   get_current_raster_layers() )
       

def get_layer_attributes( layer, field_list):
    
    if layer.selectedFeatureCount() > 0:
        features = layer.selectedFeatures()
    else:
        features = layer.getFeatures()
        
    provider = layer.dataProvider()   
    field_indices = [ provider.fieldNameIndex( field_name ) for field_name in field_list ]

    # retrieve (selected) attributes features
    data_list = [] 
    for feat in features:        
        attrs = feat.fields().toList()     
        data_list.append( [ feat.attribute( attrs[ field_ndx ].name() ) for field_ndx in field_indices ] )
        
    return data_list    
    
            
def get_pt_layer_attrs_from_field_list( pt_layer, field_list ):
    
    if pt_layer.selectedFeatureCount() > 0:
        features = pt_layer.selectedFeatures()
    else:
        features = pt_layer.getFeatures() 
    
    provider = pt_layer.dataProvider()
    
    field_indices = [ provider.fieldNameIndex( field_name ) for field_name in field_list ]

    # retrieve selected features with their geometry and relevant attributes
    data_list = [] 
    for feat in features:
             
        # fetch point geometry
        pt = feat.geometry().asPoint()

        attrs = feat.fields().toList() 

        # creates feature attribute list
        feat_list = [ pt.x(), pt.y() ]
        for field_ndx in field_indices:
            feat_list.append( str( feat.attribute( attrs[ field_ndx ].name() ) ) )

        # add to result list
        data_list.append( feat_list )
        
    return data_list


def read_line_layer_geometries( line_layer):
    
    lines = []
    line_iterator = line_layer.getFeatures()     
    for feature in line_iterator:
        geom = feature.geometry()
        if geom.isMultipart():
            lines.append( qgis_multipolyline_to_list_of_xy_list( geom.asMultiPolyline() ) )
        else:
            lines.append( [ qgis_polyline_to_xy_list( geom.asPolyline() ) ] )
            
    return lines
    

def read_layer_field( layer, curr_field_ndx ): 
    
    values = []
    iterator = layer.getFeatures()
    
    for feature in iterator:
        values.append( feature.attributes()[ curr_field_ndx ] ) 
            
    return values
            
       
def read_vector_line_qgs( line_layer, curr_field_ndx ):    
        
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
            lines.append( 'multiline', qgis_multipolyline_to_list_of_xy_list( geom.asMultiPolyline() ) ) # typedef QVector<QgsPolyline>
            # now is a list of list of (x,y) tuples
        else:           
            lines.append( ( 'line', qgis_polyline_to_xy_list( geom.asPolyline() ) ) ) # typedef QVector<QgsPoint>
                         
    return lines, progress_ids
              
                   
def qgis_polyline_to_xy_list( qgsline):
    
    assert len( qgsline ) > 0    
    return [ ( qgspoint.x(), qgspoint.y() ) for qgspoint in qgsline ]


def qgis_multipolyline_to_list_of_xy_list( qgspolyline ):
    
    return [ qgis_polyline_to_xy_list( qgsline) for qgsline in qgspolyline ]                   
                
    
def get_qgis_raster_params( raster_layer ):
    
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


def make_qgs_point( x, y):
    
    return QgsPoint(x, y)
        

def project_point( qgsPt, srcCrs, destCrs ):
    
    return QgsCoordinateTransform( srcCrs, destCrs ).transform( qgsPt )

   
    
    
    
    
    
