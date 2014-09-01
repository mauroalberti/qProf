# -*- coding: utf-8 -*-

import numpy as np
 
 
def is_number(s):
    """
    Check if string can be converted to number.

    @param  s:  parameter to check.
    @type  s:  string
    
    @return:  boolean, whether string can be converted to a number (float).
    
    """
    try:
        float(s)
    except:
        return False
    else:
        return True
 

def almost_zero( val ):
    
    tolerance = 1e-10    
    if abs( val ) > tolerance: return False
    else: return True
    
       
def ij_transfer_func( i, j, transfer_funcs ):
    """
    Return a z value as the result of a function (transfer_func_z) applied to a (x,y) point.
    This point is derived from a (i,j) point given two "transfer" functions (transfer_func_y, transfer_func_x).
    All three functions are stored into a tuple (transfer_funcs).
    
    @param  i:  array i (-y) coordinate of a single point.
    @type  i:  float.
    @param  j:  array j (x) coordinate of a single point.
    @type  j:  float. 
    @param  transfer_funcs:  tuple storing three functions (transfer_func_x, transfer_func_y, transfer_func_z)
                            that derives y from i (transfer_func_y), x from j (transfer_func_x)
                            and z from (x,y) (transfer_func_z).
    @type  transfer_funcs:  Tuple of Functions.
    
    @return:  z value - float.   
    
    """    
    transfer_func_x, transfer_func_y, transfer_func_z = transfer_funcs

    return transfer_func_z( transfer_func_x(j), transfer_func_y(i) )
    

def array_from_function( row_num, col_num, x_transfer_func, y_transfer_func, z_transfer_func ):
    """
    Creates an array of z values based on functions that map (i,j) indices (to be created) 
    into (x, y) values and then z values.
    
    @param  row_num:  row number of the array to be created.
    @type  row_num:  int.
    @param  col_num:  column number of the array to be created.
    @type  col_num:  int.
    @param  x_transfer_func:  function that derives x given a j array index.
    @type  x_transfer_func:  Function.    
    @param  y_transfer_func:  function that derives y given an i array index.
    @type  y_transfer_func:  Function.    
    @param  z_transfer_func:  function that derives z given a (x,y) point.
    @type  z_transfer_func:  Function.
    
    @return:  array of z value - array of float numbers.    
                
    """
    
    transfer_funcs = ( x_transfer_func, y_transfer_func, z_transfer_func )    
    
    return np.fromfunction( ij_transfer_func, ( row_num, col_num ), transfer_funcs=transfer_funcs )
        
 