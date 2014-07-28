
import os
from osgeo import ogr, osr


from .errors import OGR_IO_Errors


def create_def_field( field_def ):
    
    fieldDef = ogr.FieldDefn( field_def['name'], field_def['ogr_type'] )
    if field_def['ogr_type'] == ogr.OFTString:
        fieldDef.SetWidth(field_def['width'] )
        
    return fieldDef
 

def create_shapefile( path, geom_type, fields_dict_list, crs = None, layer_name="layer" ):
    
    """
    crs_prj4: projection in Proj4 text format
    geom_type = OGRwkbGeometryType: ogr.wkbPoint, ....  
    list of:  
        field dict: 'name', 
                    'type': ogr.OFTString,
                            ogr.wkbLineString, 
                            ogr.wkbLinearRing,
                            ogr.wkbPolygon,
                        
                    'width',    
    """
        
    driver = ogr.GetDriverByName("ESRI Shapefile")
    
    outShapefile = driver.CreateDataSource( str( path ) )
    if outShapefile is None:
        raise OGR_IO_Errors, 'Unable to save shapefile in provided path'

    """
    if crs is not None:
        spatialReference = osr.SpatialReference()
        spatialReference.ImportFromProj4( crs )    
        outShapelayer = outShapefile.CreateLayer("layer", geom_type, spatialReference)
    else:
    """ 
    outShapelayer = outShapefile.CreateLayer("layer", geom_type=geom_type )
        
    map( lambda field_def_params : outShapelayer.CreateField( create_def_field( field_def_params ) ), fields_dict_list ) 

    return outShapefile, outShapelayer


def open_shapefile( path, fields_dict_list ):
    
    driver = ogr.GetDriverByName("ESRI Shapefile")    
    
    dataSource = driver.Open( str( path ), 0 )
    
    if dataSource is None:
        raise OGR_IO_Errors, 'Unable to open shapefile in provided path'  
       
    point_shapelayer = dataSource.GetLayer()

    prev_solution_list = []
    in_point = point_shapelayer.GetNextFeature()
    while in_point:
        rec_id = int( in_point.GetField('id') )
        x = in_point.GetField('x')
        y = in_point.GetField('y')
        z = in_point.GetField('z')
        dip_dir = in_point.GetField('dip_dir') 
        dip_ang = in_point.GetField('dip_ang')
        descript = in_point.GetField('descript') 
        prev_solution_list.append( [ rec_id,x,y,z,dip_dir,dip_ang,descript ] )         
        in_point.Destroy()
        in_point = point_shapelayer.GetNextFeature()    
    
    #point_shapelayer.Destroy()
    dataSource.Destroy()
    
    if os.path.exists( path ):
        driver.DeleteDataSource( str( path ) )
    
    outShapefile, outShapelayer = create_shapefile( path, ogr.wkbPoint, fields_dict_list, crs = None, layer_name="layer" )
    return outShapefile, outShapelayer, prev_solution_list
        

def write_point_result( point_shapefile, point_shapelayer, recs_list2 ):

    outshape_featdef = point_shapelayer.GetLayerDefn()

    for curr_Pt in recs_list2: 

        # pre-processing for new feature in output layer
        curr_Pt_geom = ogr.Geometry(ogr.wkbPoint)
        curr_Pt_geom.AddPoint( curr_Pt[1], curr_Pt[2], curr_Pt[3] )
            
        # create a new feature
        curr_Pt_shape = ogr.Feature( outshape_featdef )
        curr_Pt_shape.SetGeometry( curr_Pt_geom )
        curr_Pt_shape.SetField( 'id', curr_Pt[0] )
                                
        curr_Pt_shape.SetField( 'x', curr_Pt[1] )
        curr_Pt_shape.SetField( 'y', curr_Pt[2] ) 
        curr_Pt_shape.SetField( 'z', curr_Pt[3] ) 

        curr_Pt_shape.SetField('dip_dir', curr_Pt[4] )
        curr_Pt_shape.SetField('dip_ang', curr_Pt[5] )             

        curr_Pt_shape.SetField('descript', curr_Pt[6] ) 
        
        # add the feature to the output layer
        point_shapelayer.CreateFeature(curr_Pt_shape)            
        
        # destroy no longer used objects
        curr_Pt_geom.Destroy(); curr_Pt_shape.Destroy()






    
    
    
    
    
        