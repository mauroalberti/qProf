
from qgis.core import QgsMapLayerRegistry, QgsMapLayer, QGis, QgsFeature


def get_map_layers_from_qgis():
    
    return QgsMapLayerRegistry.instance().mapLayers().values()
    

def get_current_raster_layers( ):        
          
    return filter( lambda layer: layer.type() == QgsMapLayer.RasterLayer, 
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
        return False, "No feature selected" 
    
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





