
from typing import Tuple, Union, List

from functools import singledispatch
import numpy as np

from pygsf.mathematics.arrays import *
from pygsf.mathematics.vectors import *

from pygsf.orientations.orientations import *

#from pygsf.geometries.shapes.space2d import *
from pygsf.geometries.shapes.space3d import *

"""
@singledispatch
def mean(
        shapes: list
) -> Shape2D:

    return None


@mean.register(list)
def mean(
        shapes: list
) -> Point2D:
    '''Mean points center'''

    return Point2D(
        x=np.mean(shapes.xs()),
        y=np.mean(shapes.ys())
    )
"""

def try_derive_bestfitplane(
    points: Points3D
) -> Tuple[bool, Union[str, CPlane3D]]:

    npaXyz = points.asXyzArray()

    #print(points.asXyzArray())

    xyz_mean = np.mean(npaXyz, axis=0)

    svd = xyzSvd(npaXyz - xyz_mean)

    if svd['result'] is None:
        return False, "Unable to calculate result"

    _, _, eigenvectors = svd['result']

    lowest_eigenvector = eigenvectors[-1, : ]  # Solution is last row

    normal = lowest_eigenvector[: 3 ] / np.linalg.norm(lowest_eigenvector[: 3 ])
    normal_vector = Vect(normal[0], normal[1], normal[2])
    normal_direct = Direct.fromVect(normal_vector)

    return True, normal_direct.normal_plane()