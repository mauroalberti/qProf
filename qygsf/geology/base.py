
from collections import namedtuple

from pygsf.geology.faults import *
from pygsf.geometries.shapes.space3d import *

georef_att_flds = [
    'id',
    'posit',
    'attitude'
]

GeorefAttitude = namedtuple(
    'GeorefAttitude',
    georef_att_flds
)


class StructuralSet:

    def __init__(self,
                 location: Point3D,
                 stratifications: Optional[List[CPlane3D]] = None,
                 foliations: Optional[List[CPlane3D]] = None,
                 faults: Optional[List[Fault]] = None
                 ):
        """
        Creates a structural set.

        :param location:
        :type location: Point.
        :param stratifications:
        :type stratifications: Optional[List[CPlane]].
        :param foliations:
        :type foliations: Optional[List[CPlane]].
        :param faults:
        :type faults: Optional[List[CPlane]].
        """

        if not isinstance(location, Point3D):
            raise Exception("Location should be Point3D but is {}".format(type(location)))

        checks = [
            (stratifications, "Stratification"),
            (foliations, "Foliations"),
            (faults, "Faults")
        ]

        for var, name in checks:
            if var:
                if not isinstance(var, List):
                    raise Exception("{} should be a List but is {}".format(name, type(var)))
                for el in var:
                    if not isinstance(el, CPlane3D):
                        raise Exception("{} should be CPlane but is {}".format(name, type(el)))

        self._location = location

        if not stratifications:
            self._strats = []
        else:
            self._strats = stratifications

        if not foliations:
            self._foliats = []
        else:
            self._foliats = foliations

        if not faults:
            self._faults = []
        else:
            self._faults = faults








