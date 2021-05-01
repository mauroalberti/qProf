from .elements import *
from ..utils.arrays import *
from ..geometries.shapes.space4d import *


class AttitudesProfile(list):

    def __init__(self, atts: List[ProfileAttitude]):

        check_type(atts, "Attitude georeferenced", List)
        for el in atts:
            check_type(el, "Attitude projection", ProfileAttitude)

        super(AttitudesProfile, self).__init__(atts)


class IntersectionsProfile(list):
    """
    A list of intersections (as arrays with ids).

    """

    def __init__(self, intersections: List[ArrayList]):

        check_type(intersections, "Profiles intersections", List)
        for el in intersections:
            check_type(el, "Profile intersections", ArrayList)

        super(IntersectionsProfile, self).__init__(intersections)


