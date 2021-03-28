
from pygsf.mathematics.scalars import *
from pygsf.mathematics.vectors import *


class RefFrame(object):

    def __init__(self,
                 versor_x: Vect,
                 versor_y: Vect):
        """
        Reference frame constructor.

        :param versor_x: Vect instance representing the x axis orientation
        :type versor_x: Vect
        :param versor_y: Vect instance representing the y axis orientation
        :type versor_y: Vect

        Examples:
        """

        if not (versor_x.is_close_to_1 and versor_y.is_close_to_1):
            raise Exception("Input vectors must be near unit")

        if not areClose(versor_x.angle_as_degrees(versor_y), 90.0):
            raise Exception("Input vectors must be sub-orthogonal")

        self._x = versor_x
        self._y = versor_y

    @property
    def x(self) -> Vect:
        """
        Return the x axis as a vector.

        :return: x axis
        :rtype: Vect

        Examples:
          >>> RefFrame(Vect(1,0,0), Vect(0,1,0)).x
          Vect(1.0000, 0.0000, 0.0000, EPSG: -1)
        """

        return self._x

    @property
    def y(self) -> Vect:
        """
        Return the y as a vector.

        :return: y axis
        :rtype: Vect

        Examples:
          >>> RefFrame(Vect(1,0,0), Vect(0,1,0)).y
          Vect(0.0000, 1.0000, 0.0000, EPSG: -1)
        """

        return self._y

    @property
    def z(self) -> Vect:
        """
        Return the z as a vector.

        :return: z axis
        :rtype: Vect

        Examples:
          >>> RefFrame(Vect(1,0,0), Vect(0,1,0)).z
          Vect(0.0000, 0.0000, 1.0000, EPSG: -1)
        """

        return self.x.cross_product(self.y)


if __name__ == "__main__":

    import doctest
    doctest.testmod()
