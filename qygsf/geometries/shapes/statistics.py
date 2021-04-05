
from .space3d import *

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

    svd_result = svd(npaXyz - xyz_mean)

    if svd_result is None:
        return False, "Unable to calculate result"

    _, _, eigenvectors = svd_result

    lowest_eigenvector = eigenvectors[-1, :]  # Solution is last row

    normal = lowest_eigenvector[: 3 ] / np.linalg.norm(lowest_eigenvector[: 3 ])
    normal_vector = Vect3D(normal[0], normal[1], normal[2])
    normal_direct = Direct.fromVect(normal_vector)

    return True, normal_direct.normal_plane()