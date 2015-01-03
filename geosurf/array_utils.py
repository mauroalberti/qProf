
from __future__  import division

from numpy import *  # general import for compatibility with formula input
from numpy.linalg import svd

from .errors import AnaliticSurfaceCalcException


def point_solution( a_array, b_array ):
    """
    finds a non-unique solution
    for a set of linear equations
    """
    
    try:
        return linalg.lstsq( a_array, b_array )[ 0 ]
    except:
        return None, None, None


def xyz_svd( xyz_array ):
    # modified after: 
    # http://stackoverflow.com/questions/15959411/best-fit-plane-algorithms-why-different-results-solved

    try:
        return dict( result = svd( xyz_array ) )
    except:
        return dict( result = None )


def formula_to_grid( array_range, array_size, formula ):
    
    a_min, a_max, b_max, b_min = array_range # note: b range reversed for conventional j order in arrays
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
        raise AnaliticSurfaceCalcException, "Error in applying formula to a and b array values"
         
    return a_list, b_list, z_list
            

 

