
import numpy as np

from .utils import almost_zero


class RefFrame( object ):
    
    def __init__( self, versor_x, versor_y, versor_z ):        
        
        assert almost_zero( versor_x.scalar_product( versor_y ) )
        assert almost_zero( versor_x.scalar_product( versor_z ) )
        assert almost_zero( versor_y.scalar_product( versor_z ) )
   
        self.axes = [ versor_x, versor_y, versor_z ]

        
    def rotation_matrix( self, rotated_frame ):
       
        for frame_axis in self.axes:        
            assert almost_zero( frame_axis.lenght_3d() - 1.0 )
        for frame_axis in rotated_frame.axes:
            assert almost_zero( frame_axis.lenght_3d() - 1.0 )
                
        rot_matrix = np.zeros( (3,3) ) * np.nan
        
        for i, rot_frame_versor in enumerate( rotated_frame.axes ):
            for j, init_frame_versor in enumerate( self.axes ):
                rot_matrix[ i, j ] = init_frame_versor.scalar_product( rot_frame_versor )
        
        return rot_matrix
    
    