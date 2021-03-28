
from typing import Any


def check_type(
        var: Any,
        name: str,
        expected_types
):
    """
    Checks the type of the variable, raising an exception when not equal.

    :param var:
    :param name:
    :param expected_types:
    :return: None
    :raise: Exception
    """

    if not (isinstance(var, expected_types)):
        raise Exception("{} should be {} but {} got".format(name, expected_types, type(var)))


def check_optional_type(var, name, expected_type):
    """
    Checks the type of the optional variable, raising an exception when not equal.

    :param var:
    :param name:
    :param expected_type: Any
    :return: None
    :rtype: None
    :raise: Exception
    """

    if var:
        if not (isinstance(var, expected_type)):
            raise Exception("{} should be {} but got {}".format(name, expected_type, type(var)))


