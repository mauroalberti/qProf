
from numpy.linalg import svd


def xyz_svd( xyz_array ):
    # modified after: 
    # http://stackoverflow.com/questions/15959411/best-fit-plane-algorithms-why-different-results-solved

    try:
        return dict( result = svd( xyz_array ) )
    except:
        return dict( result = None )

