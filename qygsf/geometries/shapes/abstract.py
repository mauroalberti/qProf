import abc


class Shape(object, metaclass=abc.ABCMeta):

    @abc.abstractmethod
    def area(self):
        """Calculate shape area"""

    @abc.abstractmethod
    def length(self):
        """Calculate shape area"""


class Point(Shape, metaclass=abc.ABCMeta):

    def area(self):
        """Calculate shape area"""

        return 0.0

    def length(self):
        """Calculate shape area"""

        return 0.0


class Segment(Shape, metaclass=abc.ABCMeta):

    def area(self):
        """Calculate shape area"""

        return 0.0


class Line(Shape, metaclass=abc.ABCMeta):

    def area(self):
        """Calculate shape area"""

        return 0.0


class Polygon(Shape, metaclass=abc.ABCMeta):

    @abc.abstractmethod
    def num_side(self):
        """Return number of sides"""


