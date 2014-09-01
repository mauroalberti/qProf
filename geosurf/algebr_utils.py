
import numpy as np


def point_solution( a_array, b_array ):
    """
    finds a non-unique solution
    for a set of linear equations
    """
    
    try:
        return np.linalg.lstsq( a_array, b_array )[ 0 ]
    except:
        return None, None, None
