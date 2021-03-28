# -*- coding: utf-8 -*-


from .ptbaxes import *
from ..orientations.orientations import Axis, Direct


def isDirect(obj) -> bool:

    return isinstance(obj, Direct) and not isinstance(obj, Axis)


def isAxis(obj) -> bool:

    return isinstance(obj, Axis)


def isPlane(obj) -> bool:

    return isinstance(obj, Plane)


def isSlick(obj) -> bool:

    return isinstance(obj, Slick)


def isFault(obj) -> bool:

    return isinstance(obj, Fault)


def isPTBAxes(obj) -> bool:

    return isinstance(obj, PTBAxes)


def isUpward(obj: Direct) -> bool:

    return obj.is_upward


def isNotUpward(obj: Direct) -> bool:

    return not obj.is_upward



