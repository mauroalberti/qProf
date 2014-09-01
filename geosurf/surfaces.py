
from numpy import * # general import for compatibility with formula input

from .deformations import calculate_geographic_scale_matrix, calculate_geographic_rotation_matrix, calculate_geographic_offset
from .deformations import define_deformation_matrices
from .errors import AnaliticSurfaceCalcException
 

def calculate_abz_lists( array_range, array_size, formula ):
    
    a_min, a_max, b_max, b_min = array_range
    array_rows, array_cols = array_size
    
    a_array = linspace( a_min, a_max, num=array_cols )
    b_array = linspace( b_max, b_min, num=array_rows ) # note: reversed for conventional j order in arrays

    try:
        a_list, b_list = [ a for a in a_array for b in b_array ], [ b for a in a_array for b in b_array ]
    except: 
        raise AnaliticSurfaceCalcException, "Error in a-b values"
    
    try:
        z_list = [ eval( formula ) for a in a_array for b in b_array ]
    except:
        raise AnaliticSurfaceCalcException, "Error in formula application to array"
         
    return a_list, b_list, z_list
            
                                   
def calculate_geosurface( analytical_params, geogr_params, deform_params ):

    array_range, array_size, formula = analytical_params
    (geog_x_min, geog_y_min), (area_length, area_width), area_rot_ang_deg = geogr_params    

    a_min, a_max, b_min, b_max = array_range
    a_range, b_range = a_max-a_min, b_max-b_min
    
    # calculate array from formula    
    try:
        X, Y, Z = calculate_abz_lists( array_range, array_size, formula ) 
    except AnaliticSurfaceCalcException, msg:
        raise AnaliticSurfaceCalcException, msg

    # calculate geographic transformations to surface
    geographic_scale_matrix = calculate_geographic_scale_matrix( a_range, b_range, area_length, area_width )
    geographic_rotation_matrix = calculate_geographic_rotation_matrix( area_rot_ang_deg )

    geographic_transformation_matrix = dot( geographic_rotation_matrix, geographic_scale_matrix )

    geographic_offset_matrix = calculate_geographic_offset( geographic_transformation_matrix, 
                                                            array( [a_min, b_min, 0.0] ),
                                                            array( [geog_x_min, geog_y_min, 0.0] ) )

    # apply total transformations to grid points 
    deformations = define_deformation_matrices( deform_params )   
    geosurface_X = []; geosurface_Y = []; geosurface_Z = []
    for x, y, z in zip( X, Y, Z ):
        pt = dot( geographic_transformation_matrix, array([x,y,z]) ) + geographic_offset_matrix
        for deformation in deformations:
            if deformation['increment'] == 'additive':
                pt = pt + deformation['matrix']
            elif deformation['increment'] == 'multiplicative':
                pt = dot( deformation['matrix'], pt )
        geosurface_X.append( pt[0] )
        geosurface_Y.append( pt[1] )
        geosurface_Z.append( pt[2] )
        
    return geosurface_X, geosurface_Y, geosurface_Z 
      
            
