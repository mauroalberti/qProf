from qygsf import Vect, Point4D
from qygsf.mathematics.defaults import MIN_SCALAR_VALUE, MIN_ANGLE_DEGR_VALUE
from qygsf.utils.arrays import point_solution

MIN_2D_SEPARATION_THRESHOLD = 1e-10
MINIMUM_SEPARATION_THRESHOLD = 1e-10
MINIMUM_VECTOR_MAGNITUDE = 1e-10


class ParamLine3D(object):
    """
    parametric line
    srcPt: source Point
    l, m, n: .....
    """

    def __init__(self, srcPt, l, m, n):

        assert -1.0 <= l <= 1.0
        assert -1.0 <= m <= 1.0
        assert -1.0 <= n <= 1.0

        self._srcPt = srcPt
        self._l = l
        self._m = m
        self._n = n

    def intersect_cartes_plane(self, cartes_plane):
        """
        Return intersection point between parametric line and Cartesian plane
        """

        # line parameters
        x1, y1, z1 = self._srcPt.x, self._srcPt.y, self._srcPt.z
        l, m, n = self._l, self._m, self._n

        # Cartesian plane parameters
        a, b, c, d = cartes_plane.a, cartes_plane.b, cartes_plane.c, cartes_plane.d

        try:
            k = (a * x1 + b * y1 + c * z1 + d) / (a * l + b * m + c * n)
        except ZeroDivisionError:
            return None

        return Point4D(x1 - l * k,
                       y1 - m * k,
                       z1 - n * k)


def eq_xy_pair(xy_pair_1, xy_pair_2):

    if xy_pair_1[0] == xy_pair_2[0] and xy_pair_1[1] == xy_pair_2[1]:
        return True

    return False

'''
def remove_equal_consecutive_xypairs(xy_list):

    out_xy_list = [xy_list[0]]

    for n in range(1, len(xy_list)):
        if not eq_xy_pair(xy_list[n], out_xy_list[-1]):
            out_xy_list.append(xy_list[n])

    return out_xy_list
'''


class CPlane3D:
    """
    Cartesian plane.
    Expressed by equation:
    ax + by + cz + d = 0

    Note: Plane is locational - its position in space is defined.
    This contrast with GPlane, defined just by its attitude, but with undefined position

    """

    def __init__(self, a, b, c, d):

        self._a = float(a)
        self._b = float(b)
        self._c = float(c)
        self._d = float(d)

    @property
    def a(self):
        """
        Return a coefficient of a Plane instance.

        Example:
          >>> CPlane3D(1, 0, 0, 2).a
          1.0
        """

        return self._a

    @property
    def b(self):
        """
        Return b coefficient of a Plane instance.

        Example:
          >>> CPlane3D(1, 4, 0, 2).b
          4.0
        """

        return self._b

    @property
    def c(self):
        """
        Return a coefficient of a Plane instance.

        Example:
          >>> CPlane3D(1, 0, 5.4, 2).c
          5.4
        """

        return self._c

    @property
    def d(self):
        """
        Return a coefficient of a Plane instance.

        Example:
          >>> CPlane3D(1, 0, 0, 2).d
          2.0
        """

        return self._d

    @property
    def v(self):
        """
        Return coefficients of a Plane instance.

        Example:
          >>> CPlane3D(1, 1, 7, -4).v
          (1.0, 1.0, 7.0, -4.0)
        """
        return self.a, self.b, self.c, self.d

    @classmethod
    def from_points(cls, pt1, pt2, pt3):
        """
        Create a Plane from three given Point instances.

        Example:
          >>> CPlane3D.from_points(Point4D(0, 0, 0), Point4D(1, 0, 0), Point4D(0, 1, 0))
          Plane(0.0000, 0.0000, 1.0000, 0.0000)
          >>> CPlane3D.from_points(Point4D(0, 0, 0), Point4D(0, 1, 0), Point4D(0, 0, 1))
          Plane(1.0000, 0.0000, 0.0000, 0.0000)
        """

        matr_a = np.array([[pt1.y, pt1.z, 1],
                           [pt2.y, pt2.z, 1],
                           [pt3.y, pt3.z, 1]])

        matr_b = - np.array([[pt1.x, pt1.z, 1],
                             [pt2.x, pt2.z, 1],
                             [pt3.x, pt3.z, 1]])

        matr_c = np.array([[pt1.x, pt1.y, 1],
                           [pt2.x, pt2.y, 1],
                           [pt3.x, pt3.y, 1]])

        matr_d = - np.array([[pt1.x, pt1.y, pt1.z],
                             [pt2.x, pt2.y, pt2.z],
                             [pt3.x, pt3.y, pt3.z]])

        return cls(np.linalg.det(matr_a),
                   np.linalg.det(matr_b),
                   np.linalg.det(matr_c),
                   np.linalg.det(matr_d))

    def __repr__(self):

        return "Plane({:.4f}, {:.4f}, {:.4f}, {:.4f})".format(*self.v)

    @property
    def nversor(self):
        """
        Return the versor normal to the cartesian plane.

        Examples:
          >>> CPlane3D(0, 0, 5, -2).nversor
          Vect(0.0000, 0.0000, 1.0000)
          >>> CPlane3D(0, 7, 0, 5).nversor
          Vect(0.0000, 1.0000, 0.0000)
        """

        return Vect(self.a, self.b, self.c).versor()

    def gplane_point(self):
        """
        Converts a cartesian plane into a geological plane
        and a point lying in the plane (non-unique solution).

        Examples:
          >>> gpl, pt = CPlane3D(0, 0, 1, -1).gplane_point()
          >>> gpl
          GPlane(000.00, +00.00)
          >>> pt
          Point(0.0000, 0.0000, 1.0000, nan)
        """

        geol_plane = self.nversor.gvect.normal_gplane
        point = Point4D(*point_solution(np.array([[self.a, self.b, self.c]]),
                                        np.array([-self.d])))
        return geol_plane, point

    def inters_versor(self, another):
        """
        Return intersection versor for two intersecting planes.

        >>> a = CPlane3D(1, 0, 0, 0)
        >>> b = CPlane3D(0, 0, 1, 0)
        >>> a.inters_versor(b)
        Vect(0.0000, -1.0000, 0.0000)
        """

        return self.nversor.vp(another.nversor).versor()

    def inters_point(self,
                     another
                     ):
        """
        Return point on intersection line (obviously non-unique solution)
        for two planes.

        >>> pa = CPlane3D(1, 0, 0, 0)
        >>> pb = CPlane3D(0, 0, 1, 0)
        >>> pa.inters_point(pb)
        Point(0.0000, 0.0000, 0.0000, nan)
        """

        # find a point lying on the intersection line (this is a non-unique solution)
        a = np.array([[self.a, self.b, self.c], [another.a, another.b, another.c]])
        b = np.array([-self.d, -another.d])
        x, y, z = point_solution(a, b)

        return Point4D(x, y, z)

    def is_point_inplane(self, pt):
        """
          Check whether a point lie in a plane.

          >>> pl = CPlane3D(0, 0, 1, 0)
          >>> pt = Point4D(0, 1, 0)
          >>> pl.is_point_inplane(pt)
          True
        """

        if abs(self.a * pt.x + self.b * pt.y + self.c * pt.z + self.d) < MIN_SCALAR_VALUE:
            return True
        else:
            return False

    def angle(self, another):
        """
        Calculate angle (in degrees) between two planes.

        Examples:
          >>> CPlane3D(1,0,0,0).angle(CPlane3D(0,1,0,0))
          90.0
          >>> CPlane3D(1,0,0,0).angle(CPlane3D(1,0,1,0))
          45.0
          >>> CPlane3D(1,0,0,0).angle(CPlane3D(1,0,0,0))
          0.0
        """

        angle_degr = self.nversor.angle(another.nversor)
        if abs(angle_degr) < MIN_ANGLE_DEGR_VALUE:
            angle_degr = 0.0
        elif angle_degr > 90.0:
            angle_degr = 180.0 - angle_degr

        return angle_degr