
import numpy as np
import json

   
def save_surface_vtk( output_filepath, geodata ):

    geosurface, grid_dims = geodata
    X, Y, Z = geosurface
    
    X_arr = np.array(X, dtype=float)
    Y_arr = np.array(Y, dtype=float)
    Z_arr = np.array(Z, dtype=float)

    n_points = np.size( X_arr )
    
    n_rows, n_cols = grid_dims        
                      
    with open( output_filepath, 'w' ) as outfile:
    
        outfile.write( '# vtk DataFile Version 2.0\n' )
        outfile.write( 'Geosurface - qgSurf vers. 0.3.0\n' ) 
        outfile.write( 'ASCII\n' )   
        outfile.write( '\nDATASET POLYDATA\n' )
        
        outfile.write( 'POINTS %d float\n' % n_points )
        for n in xrange( n_points ):
            outfile.write( '%.4f %.4f %.4f\n' % ( X_arr[n], Y_arr[n], Z_arr[n] ) )
        
        outfile.write( '\n' )      
        
        outfile.write( 'TRIANGLE_STRIPS %d %d\n' % ( n_cols-1, (n_cols-1)*(1+n_rows*2) ) )
    
        num_points_strip = n_rows * 2   
        for l in xrange( n_cols - 1 ):
            triangle_strip_string = "%d " % num_points_strip
            for p in xrange( n_rows ):            
                triangle_strip_string += "%d %d " % ( (l+1)*n_rows+p, l*n_rows+p )
            triangle_strip_string += "\n" 
            outfile.write( triangle_strip_string )      
    
            
def save_surface_grass( output_filepath, geodata ):
    # Save in Grass format

    geosurface, grid_dims = geodata
    X, Y, Z = geosurface

    X_arr = np.array(X, dtype=float)
    Y_arr = np.array(Y, dtype=float)
    Z_arr = np.array(Z, dtype=float)
    
    n_rows, n_cols = grid_dims   
                            
    with open( output_filepath, 'w' ) as outfile:
        outfile.write( 'VERTI:\n' ) 
        for l in xrange( n_cols - 1 ):
            for p in xrange( n_rows - 1 ):
                start_point_ndx = l*n_rows + p
                forward_line_point_ndx = start_point_ndx + n_rows
                outfile.write('F 4\n')
                outfile.write(' %.4f %.4f %.4f\n' % ( X_arr[start_point_ndx], Y_arr[start_point_ndx], Z_arr[start_point_ndx] ) )
                outfile.write(' %.4f %.4f %.4f\n' % ( X_arr[start_point_ndx+1], Y_arr[start_point_ndx+1], Z_arr[start_point_ndx+1] ) )
                outfile.write(' %.4f %.4f %.4f\n' % ( X_arr[forward_line_point_ndx], Y_arr[forward_line_point_ndx], Z_arr[forward_line_point_ndx] ) )        
                outfile.write(' %.4f %.4f %.4f\n' % ( X_arr[start_point_ndx], Y_arr[start_point_ndx], Z_arr[start_point_ndx] ) )            
                outfile.write('F 4\n')
                outfile.write(' %.4f %.4f %.4f\n' % ( X_arr[forward_line_point_ndx], Y_arr[forward_line_point_ndx], Z_arr[forward_line_point_ndx] ) ) 
                outfile.write(' %.4f %.4f %.4f\n' % ( X_arr[start_point_ndx+1], Y_arr[start_point_ndx+1], Z_arr[start_point_ndx+1] ) )
                outfile.write(' %.4f %.4f %.4f\n' % ( X_arr[forward_line_point_ndx+1], Y_arr[forward_line_point_ndx+1], Z_arr[forward_line_point_ndx+1] ) )
                outfile.write(' %.4f %.4f %.4f\n' % ( X_arr[forward_line_point_ndx], Y_arr[forward_line_point_ndx], Z_arr[forward_line_point_ndx] ) ) 


def save_surface_gas( output_filepath, geodata ):
    
    with open( output_filepath, 'w' ) as outfile:
        json.dump( geodata, outfile )
    
        
        
                         
            
