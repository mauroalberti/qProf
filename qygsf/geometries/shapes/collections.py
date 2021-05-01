import numbers
from typing import Optional, List, Tuple, Union

import numpy as np

from .space2d import XYArrayPair
from .space3d import Line3D
from .space4d import Line4D


class NamedLines:
    """
    Class representing multiple grid profiles
    derived from a single trace,
    that may be a straight or a non-straight line.
    """

    def __init__(self,
                 named_lines: Optional[List[Tuple[str, Union[Line3D, Line4D]]]] = None):

        if named_lines is None:
            self._named_lines = []  # list of name and Line3/4D
        else:
            self._named_lines = named_lines

    def __repr__(self):

        return f"""NamedLines with {len(self._named_lines)} Line(s):
        {self._named_lines}
        """

    def num_lines(self) -> numbers.Integral:
        """
        Returns the number of lines present in the collection.
        """

        return len(self._named_lines)

    def clear_lines(self):

        self._named_lines = None

    def set_named_lines(self,
                        named_grid_profiles
                        ):

        self._named_lines = named_grid_profiles

    @property
    def named_lines(self):
        """

        :return:
        """

        return self._named_lines

    def add_named_lines(self, named_line):

        self._named_lines += named_line

    def s_max(self):

        return np.nanmax([line.length_2d() for _, line in self._named_lines])

    def z_min(self):

        return np.nanmin([line.z_min() for _, line in self._named_lines])

    def z_max(self):

        return np.nanmax([line.z_max() for _, line in self._named_lines])

    def to_sz_arrays(self,
                     ndx: numbers.Integral) -> XYArrayPair:
        """
        Transform the indexed line in an XYArrayPair instance
        (distance and z along the line profile).
        """

        name, line = self._named_lines[ndx]
        s_array = line.incremental_length_2d()
        z_array = line.z_array()

        return XYArrayPair(
            x_array=s_array,
            y_array=z_array,
            id=name
        )

