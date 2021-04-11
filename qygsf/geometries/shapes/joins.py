from enum import Enum
from typing import Union, List, Optional

from .space2d import *
from .space3d import *


class JoinTypes(Enum):
    """
    Enumeration for Line and Segment type.
    """

    START_START = 1  # start point coincident with start point
    START_END   = 2  # start point coincident with end point
    END_START   = 3  # end point coincident with start point
    END_END     = 4  # end point coincident with end point


def analizeJoins2D(
        first: Union[Line2D, Segment2D],
        second: Union[Line2D, Segment2D]
) -> List[Optional[JoinTypes]]:
    """
    Analyze join types between two lines/segments.

    :param first: a line or segment.
    :param second: a line or segment.
    :return: a list of join types.

    Examples:
      >>> first = Segment2D(Point2D(x=0,y=0), Point2D(x=1,y=0))
      >>> second = Segment2D(Point2D(x=1,y=0), Point2D(x=0,y=0))
      >>> analizeJoins2D(first, second)
      [<JoinTypes.START_END: 2>, <JoinTypes.END_START: 3>]
      >>> first = Segment2D(Point2D(x=0,y=0), Point2D(x=1,y=0))
      >>> second = Segment2D(Point2D(x=2,y=0), Point2D(x=3,y=0))
      >>> analizeJoins2D(first, second)
      []
    """

    join_types = []

    if first.start_pt.is_coincident(second.start_pt):
        join_types.append(JoinTypes.START_START)

    if first.start_pt.is_coincident(second.end_pt):
        join_types.append(JoinTypes.START_END)

    if first.end_pt.is_coincident(second.start_pt):
        join_types.append(JoinTypes.END_START)

    if first.end_pt.is_coincident(second.end_pt):
        join_types.append(JoinTypes.END_END)

    return join_types


def analizeJoins3D(
        first: Union[Line3D, Segment3D],
        second: Union[Line3D, Segment3D]
) -> List[Optional[JoinTypes]]:
    """
    Analyze join types between two lines/segments.

    :param first: a line or segment.
    :type first: Line or Segment.
    :param second: a line or segment.
    :param second: Line or Segment.
    :return: a list of join types.
    :rtype: List[Optional[JoinTypes]].

    Examples:
    """

    join_types = []

    if first.start_pt.is_coincident(second.start_pt):
        join_types.append(JoinTypes.START_START)

    if first.start_pt.is_coincident(second.end_pt):
        join_types.append(JoinTypes.START_END)

    if first.end_pt.is_coincident(second.start_pt):
        join_types.append(JoinTypes.END_START)

    if first.end_pt.is_coincident(second.end_pt):
        join_types.append(JoinTypes.END_END)

    return join_types