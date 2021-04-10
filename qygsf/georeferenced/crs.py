
from typing import Any
import numbers


min_epsg_crs_code = 2000  # checked 2019-06-14 in EPSG database


class Crs(object):
    """
    CRS class.
    Currently it is in a basic form,
    just managing simple comparisons and validity checks.

    """

    def __init__(self, epsg_cd: numbers.Integral = -1):

        self._epsg = int(epsg_cd)

    @property
    def epsg_code(self) -> int:

        return self._epsg

    def valid(self):

        return self.epsg_code >= min_epsg_crs_code

    def __repr__(self):

        return "EPSG:{}".format(self.epsg_code)

    def __eq__(self, another) -> bool:
        """
        Checks for equality between Crs instances.
        Currently it considers equal two Crs instances when they have the
        same EPSG code, even an invalid one (i.e., -1).

        :param another: the Crs instance to compare with.
        :type another: Crs.
        :return: whether the input Crs instance is equal to the current one.
        :rtype: bool.
        :raise: Exception.
        """

        if not (isinstance(another, Crs)):
            raise Exception("Input instance should be Crs but is {}".format(type(another)))

        return self.epsg_code == another.epsg_code


def check_crs(
    template_element: Any,
    checked_element: Any
) -> None:
    """
    Check whether two spatial elements have the same georeferenced, raising an exception when not equal.
    The two elements should both implement the georeferenced property.

    :param template_element: first spatial element.
    :type template_element: Any
    :param checked_element: second spatial element.
    :type checked_element: Any
    :return: nothing
    :rtype: None
    :raise: Exception
    """

    if checked_element.crs != template_element.crs:
        raise Exception("checked {} instance has {} EPSG code but {} expected".format(
            type(checked_element).__name__,
            checked_element.epsg_code,
            template_element.epsg_code
        )
    )


def check_epsg(
    spatial_element: Any,
    epsg_code: numbers.Integral
) -> None:
    """
    Check whether a spatial element has a given EPSG code, raising an exception when not true.
    The spatial element should implement the epsg_code method.

    :param spatial_element: spatial element
    :type spatial_element: Any
    :param epsg_code: the EPSG code
    :type epsg_code: numbers.Integral
    :return: nothing
    :rtype: None
    :raise: Exception
    """

    if spatial_element.epsg_code != epsg_code:
        raise Exception("checked {} instance has {} EPSG code but {} expected".format(
            type(spatial_element).__name__,
            spatial_element.epsg_code,
            epsg_code
        )
    )
