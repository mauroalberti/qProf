
from math import radians, sin, cos, tan
import numpy as np

from .spatial import GeolAxis


def calculate_geographic_scale_matrix( a_range, b_range, grid_length, grid_width ):
    
    assert a_range > 0.0
    assert b_range > 0.0
    assert grid_length > 0.0
    assert grid_width > 0.0
            
    sx = grid_length / a_range
    sy = grid_width / b_range
    sz = 1
    
    return np.array( [ ( sx , 0.0, 0.0 ), ( 0.0, sy, 0.0 ), ( 0.0, 0.0, sz ) ] )


def calculate_geographic_rotation_matrix( grid_rot_angle_degr ):
    
    grid_rot_angle_rad = radians( grid_rot_angle_degr )
    sin_rot_angle = sin( grid_rot_angle_rad )
    cos_rot_angle = cos( grid_rot_angle_rad )
    
    return np.array( [ ( cos_rot_angle, -sin_rot_angle, 0.0 ), 
                       ( sin_rot_angle,  cos_rot_angle, 0.0 ), 
                       ( 0.0,            0.0,           1.0 ) ] )


def calculate_geographic_offset( transformation_matrix, llc_point_matr, llc_point_geog):
                                           
    return llc_point_geog - np.dot( transformation_matrix, llc_point_matr )


def calculate_rotation_matrix( rot_axis_trend, rot_axis_plunge, rot_angle ):
    
    phi = radians( rot_angle )

    rotation_versor = GeolAxis(rot_axis_trend, rot_axis_plunge).to_versor()

    l = rotation_versor._x
    m = rotation_versor._y
    n = rotation_versor._z

    cos_phi = cos( phi )
    sin_phi = sin( phi )
    
    a11 = cos_phi + ((l*l)*(1-cos_phi))
    a12 = ((l*m)*(1-cos_phi))-(n*sin_phi)
    a13 = ((l*n)*(1-cos_phi))+(m*sin_phi)
    
    a21 = ((l*m)*(1-cos_phi))+(n*sin_phi)
    a22 = cos_phi + ((m*m)*(1-cos_phi))
    a23 = ((m*n)*(1-cos_phi))-(l*sin_phi)
    
    a31 = ((l*n)*(1-cos_phi))-(m*sin_phi)
    a32 = ((m*n)*(1-cos_phi))+(l*sin_phi)
    a33 = cos_phi + ((n*n)*(1-cos_phi))
    
    return np.array( [ ( a11, a12, a13 ),
                       ( a21, a22, a23 ),
                       ( a31, a32, a33 ) ] )
        

def calculate_scaling_matrix( scale_factor_x, scale_factor_y, scale_factor_z ):
    
    return np.array( [ ( scale_factor_x, 0.0, 0.0 ),
                       ( 0.0, scale_factor_y, 0.0 ),
                       ( 0.0, 0.0, scale_factor_z ) ] )
    

def calculate_simple_shear_horiz_matrix( phi_angle_degr, alpha_angle_degr ):
    
    phi_angle_rad = radians(phi_angle_degr)    
    alpha_angle_rad = radians(alpha_angle_degr)    
    
    gamma = tan( phi_angle_rad )
    sin_a = sin( alpha_angle_rad )
    cos_a = cos( alpha_angle_rad )
            
    return np.array( [ (1.0-gamma*sin_a*cos_a, gamma*cos_a*cos_a, 0.0 ),
                       ( -gamma*sin_a*sin_a, 1.0+gamma*sin_a*cos_a, 0.0),
                       (0.0, 0.0, 1.0) ])        
            

def calculate_simple_shear_vert_matrix( phi_angle_degr, alpha_angle_degr ):
    
    phi_angle_rad = radians(phi_angle_degr)    
    alpha_angle_rad = radians(alpha_angle_degr)    
    
    gamma = tan( phi_angle_rad )
    sin_a = sin( alpha_angle_rad )
    cos_a = cos( alpha_angle_rad )
            
    return np.array( [ ( 1.0, 0.0, gamma*cos_a ),
                       ( 0.0, 1.0, gamma*sin_a ),
                       ( 0.0, 0.0, 1.0 ) ])        

            
def define_deformation_matrices( deform_params ):
      
    deformation_matrices = []
    
    for deform_param in deform_params :
        if deform_param['type'] == 'displacement':
            displ_x = deform_param['parameters']['delta_x']
            displ_y = deform_param['parameters']['delta_y']
            displ_z = deform_param['parameters']['delta_z']            
            deformation = { 'increment': 'additive', 'matrix': np.array( [ displ_x, displ_y, displ_z ] ) }
        elif deform_param['type'] == 'rotation':
            rotation_matrix = calculate_rotation_matrix( deform_param['parameters']['rotation axis trend'],
                                                         deform_param['parameters']['rotation axis plunge'],
                                                         deform_param['parameters']['rotation angle'] )
            deformation = { 'increment': 'multiplicative', 'matrix': rotation_matrix }
        elif deform_param['type'] == 'scaling':
            scaling_matrix = calculate_scaling_matrix( deform_param['parameters']['x factor'],
                                                       deform_param['parameters']['y factor'],
                                                       deform_param['parameters']['z factor'] )            
            deformation = { 'increment': 'multiplicative', 'matrix': scaling_matrix }
        elif deform_param['type'] == 'simple shear - horizontal':
            simple_shear_horiz_matrix = calculate_simple_shear_horiz_matrix( deform_param['parameters']['phi angle (degr.)'],
                                                                             deform_param['parameters']['alpha angle (degr.)'] )            
            deformation = { 'increment': 'multiplicative', 'matrix': simple_shear_horiz_matrix }
        elif deform_param['type'] == 'simple shear - vertical':
            simple_shear_vert_matrix = calculate_simple_shear_vert_matrix( deform_param['parameters']['phi angle (degr.)'],
                                                                             deform_param['parameters']['alpha angle (degr.)'] )            
            deformation = { 'increment': 'multiplicative', 'matrix': simple_shear_vert_matrix }
       
                               
        deformation_matrices.append( deformation )       


    return deformation_matrices
    
    
