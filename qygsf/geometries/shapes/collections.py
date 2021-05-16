import numbers
from typing import Optional, List, Tuple

import numpy as np

from ...utils.types import *
from .abstract import *
from .space2d import XYArrayPair


class NamedLines(list):
    """
    Class representing multiple grid profiles
    derived from a single trace,
    that may be a straight or a non-straight line.
    """

    def __init__(self,
                 named_lines: Optional[List[Tuple[str, Line]]] = None):

        if named_lines is None:
            named_lines = []  # list of name and Line
        else:
            print(f"DEBUG: type(named_lines) -> {type(named_lines)}")
            check_type(named_lines, "Named lines", List)
            for name, line in named_lines:
                check_type(name, "Name", (str, numbers.Integral))
                check_type(line, "Line", Line)

        print(f"DEBUG: intializing NamedLines with {named_lines}")
        super(NamedLines, self).__init__(named_lines)

    def __repr__(self):

        return f"""NamedLines with {len(self)} Line(s):
        {self}
        """

    def num_lines(self) -> numbers.Integral:
        """
        Returns the number of lines present in the collection.
        """

        return len(self)

    """
    def clear_lines(self):

        self = []
    """

    """
    def set_named_lines(self,
                        named_grid_profiles
                        ):

        self = named_grid_profiles
    """

    """
    def add_named_lines(self, named_line):

        self += named_line
    """

    def s_max(self):

        return np.nanmax([line.length_2d() for _, line in self])

    def z_min(self):

        return np.nanmin([line.z_min() for _, line in self])

    def z_max(self):

        return np.nanmax([line.z_max() for _, line in self])

    def to_sz_arrays(self,
                     ndx: numbers.Integral) -> XYArrayPair:
        """
        Transform the indexed line in an XYArrayPair instance
        (distance and z along the line profile).
        """

        name, line = self[ndx]
        s_array = line.incremental_length_2d()
        z_array = line.z_array()

        return XYArrayPair(
            x_array=s_array,
            y_array=z_array,
            id=name
        )

