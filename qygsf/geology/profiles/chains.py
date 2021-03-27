
from typing import Optional

from array import array

from pygsf.geology.profiles.elements import *
from pygsf.utils.arrays import *


'''
class ProfileSubpartsSet(list):
    """
    Describes subelements (point or segment-type) of a profile
    """

    def __init__(self,
        subelements: Optional[List] = None):
        """
        """

        if subelements is None:
            subelements = []

        super(ProfileSubpartsSet, self).__init__(subelements)
'''

class TopographicProfile:
    """

    """

    def __init__(self, s_array: array, z_array: array):

        check_type(s_array, "Distances array", array)
        if s_array.typecode != 'd':
            raise Exception("s array must be of type double")

        num_steps = len(s_array)

        check_type(z_array, "Scalar values array", array)
        if z_array.typecode != 'd':
            raise Exception("z array must be of type double")
        if len(z_array) != num_steps:
            raise Exception("z array must have the same length of s array")

        self._num_steps = num_steps
        self._s = s_array
        self._z = z_array

    def s_arr(self) -> array:
        """
        Return the s array.

        :return: the s array.
        :rtype: array.
        """

        return self._s

    def z_arr(self) -> array:
        """
        Return the z arrays.

        :return: the z array.
        :rtype: array.
        """

        return self._z

    def z(self,
        ndx: numbers.Integral
    ) -> numbers.Real:
        """
        Returns the z value with the index ndx.

        :param ndx: the index in the z array
        :type ndx: numbers.Integral
        :return: the z value corresponding to the ndx index
        :rtype: numbers.Real
        """

        return self.z_arr()[ndx]

    def s(self,
        ndx: numbers.Integral
    ) -> numbers.Real:
        """
        Returns the s value with the index ndx.

        :param ndx: the index in the s array
        :type ndx: numbers.Integral
        :return: the s value corresponding to the ndx index
        :rtype: numbers.Real
        """

        return self.s_arr()[ndx]

    def s_min(self) -> numbers.Real:
        """
        Returns the minimum s value.

        :return: the minimum s value.
        :rtype: numbers.Real.
        """

        return min(self._s)

    def s_max(self) -> numbers.Real:
        """
        Returns the maximum s value.

        :return: the maximum s value.
        :rtype: numbers.Real.
        """

        return max(self._s)

    def z_min(self) -> numbers.Real:
        """
        Returns the minimum z value.

        :return: the minimum z value.
        :rtype: numbers.Real.
        """

        return min(self._z)

    def z_max(self) -> numbers.Real:
        """
        Returns the maximum z value.

        :return: the maximum z value.
        :rtype: numbers.Real.
        """

        return max(self._z)

    def num_steps(self) -> int:
        """
        Return the number of steps of the profiles.

        :return: number of steps of the profiles.
        :rtype: int.
        """

        return self._num_steps

    def profile_length(self) -> numbers.Real:
        """
        Returns the length of the profile.

        :return: length of profile.
        :rtype: numbers.Real.
        """

        return self.s_arr()[-1]

    def s_upper_index(self,
        s_val: numbers.Real
    ) -> Optional[numbers.Integral]:
        """
        Returns the optional index in the s array of the provided value.

        :param s_val: the value to search the index for in the s array
        :type s_val: numbers.Real
        :return: the optional index in the s array of the provided value
        :rtype: Optional[numbers.Integral]

        Examples:
          >>> p = TopographicProfile(array('d', [ 0.0,  1.0,  2.0,  3.0, 3.14]), array('d', [10.0, 20.0, 0.0, 14.5, 17.9]))
          >>> p.s_upper_index(-1) is None
          True
          >>> p.s_upper_index(5) is None
          True
          >>> p.s_upper_index(0.5)
          1
          >>> p.s_upper_index(0.75)
          1
          >>> p.s_upper_index(1.0)
          1
          >>> p.s_upper_index(2.0)
          2
          >>> p.s_upper_index(2.5)
          3
          >>> p.s_upper_index(3.11)
          4
          >>> p.s_upper_index(3.14)
          4
          >>> p.s_upper_index(0.0)
          0
        """

        check_type(s_val, "Input value", numbers.Real)

        if s_val < self.s_min() or s_val > self.s_max():
            return None

        return np.searchsorted(self.s_arr(), s_val)

    def z_for_s(self,
                s_val: numbers.Real
                ) -> Optional[numbers.Integral]:
        """
        Returns the optional interpolated z value in the z array of the provided s value.

        :param s_val: the value to search the index for in the s array
        :type s_val: numbers.Real
        :return: the optional index in the s array of the provided value
        :rtype: Optional[numbers.Integral]

        Examples:
          >>> p = TopographicProfile(array('d', [ 0.0,  1.0,  2.0,  3.0, 3.14]), array('d', [10.0, 20.0, 0.0, 14.5, 17.9]))
          >>> p.z_for_s(-1) is None
          True
          >>> p.z_for_s(5) is None
          True
          >>> p.z_for_s(0.5)
          15.0
          >>> p.z_for_s(0.75)
          17.5
          >>> p.z_for_s(2.5)
          7.25
          >>> p.z_for_s(3.14)
          17.9
          >>> p.z_for_s(0.0)
          10.0
          >>> p.z_for_s(1.0)
          20.0
        """

        check_type(s_val, "Input value", numbers.Real)

        ndx = self.s_upper_index(s_val)

        if ndx is not None:

            if ndx == 0:
                return self.z(0)

            val_z_i = self.z(ndx-1)
            val_z_i_next = self.z(ndx)
            delta_val_z = val_z_i_next - val_z_i

            if delta_val_z == 0.0:
                return val_z_i

            val_s_i = self.s(ndx-1)
            val_s_i_next = self.s(ndx)
            delta_val_s = val_s_i_next - val_s_i

            if delta_val_s == 0.0:
                return val_z_i

            d_s = s_val - val_s_i

            return val_z_i + d_s*delta_val_z/delta_val_s

        else:

            return None

    def s_subset(self,
        s_start: numbers.Real,
        s_end: Optional[numbers.Real] = None
    ) -> array:
        """
        Return the s array values defined by the provided s range (extremes included).
        When the end value is not provided, a single-valued array is returned.

        :param s_start: the start s value (distance along the profile)
        :type s_start: numbers.Real
        :param s_end: the optional end s value (distance along the profile)
        :type s_end: numbers.Real
        :return: the s array subset, possibly with just a value
        :rtype: array

        Examples:
          >>> p = TopographicProfile(array('d', [ 0.0,  1.0,  2.0,  3.0, 3.14]), array('d', [10.0, 20.0, 0.0, 14.5, 17.9]))
          >>> p.s_subset(1.0)
          array('d', [1.0])
          >>> p.s_subset(0.0)
          array('d', [0.0])
          >>> p.s_subset(0.75)
          array('d', [0.75])
          >>> p.s_subset(3.14)
          array('d', [3.14])
          >>> p.s_subset(1.0, 2.0)
          array('d', [1.0, 2.0])
          >>> p.s_subset(0.75, 2.0)
          array('d', [0.75, 1.0, 2.0])
          >>> p.s_subset(0.75, 2.5)
          array('d', [0.75, 1.0, 2.0, 2.5])
          >>> p.s_subset(0.75, 3.0)
          array('d', [0.75, 1.0, 2.0, 3.0])
          >>> p.s_subset(0.75, 0.5)
          NotImplemented
          >>> p.s_subset(-1, 1)
          NotImplemented
          >>> p.s_subset(-1)
          NotImplemented
          >>> p.s_subset(0.0, 10)
          NotImplemented
          >>> p.s_subset(0.0, 10)
          NotImplemented
          >>> p.s_subset(0.0, 3.14)
          array('d', [0.0,  1.0,  2.0,  3.0, 3.14])
        """

        if not s_end and not self.s_min() <= s_start <= self.s_max():

            return NotImplemented

        if s_end and not self.s_min() <= s_start <= s_end <= self.s_max():

            return NotImplemented

        if s_end is None or s_end == s_start:

            return array('d', [s_start])

        result = array('d', [])

        s_start_upper_index_value = self.s_upper_index(s_start)

        if s_start < self.s(s_start_upper_index_value):
            result.append(s_start)

        s_end_upper_index_value = self.s_upper_index(s_end)

        for ndx in range(s_start_upper_index_value, s_end_upper_index_value):
            result.append(self.s(ndx))

        if s_end > self.s(s_end_upper_index_value - 1):
            result.append(s_end)

        return result

    def zs_from_s_range(self,
                        s_start: numbers.Real,
                        s_end: Optional[numbers.Real] = None
                        ) -> array:
        """
        Return the z array values defined by the provided s range (extremes included).
        When the end value is not provided, a single-valued array is returned.

        :param s_start: the start s value (distance along the profile)
        :type s_start: numbers.Real
        :param s_end: the optional end s value (distance along the profile)
        :type s_end: numbers.Real
        :return: the z array, possibly with just a value
        :rtype: array

        Examples:
          >>> p = TopographicProfile(array('d', [ 0.0,  1.0,  2.0,  3.0, 3.14]), array('d', [10.0, 20.0, 0.0, 14.5, 17.9]))
          >>> p.zs_from_s_range(1.0)
          array('d', [20.0])
          >>> p.zs_from_s_range(0.0)
          array('d', [10.0])
          >>> p.zs_from_s_range(0.75)
          array('d', [17.5])
          >>> p.zs_from_s_range(3.14)
          array('d', [17.9])
          >>> p.zs_from_s_range(1.0, 2.0)
          array('d', [20.0, 0.0])
          >>> p.zs_from_s_range(0.75, 2.0)
          array('d', [17.5, 20.0, 0.0])
          >>> p.zs_from_s_range(0.75, 2.5)
          array('d', [17.5, 20.0, 0.0, 7.25])
          >>> p.zs_from_s_range(0.75, 3.0)
          array('d', [17.5, 20.0, 0.0, 14.5])
          >>> p.zs_from_s_range(0.75, 0.5)
          NotImplemented
          >>> p.zs_from_s_range(-1, 1)
          NotImplemented
          >>> p.zs_from_s_range(-1)
          NotImplemented
          >>> p.zs_from_s_range(0.0, 10)
          NotImplemented
          >>> p.zs_from_s_range(0.0, 10)
          NotImplemented
          >>> p.zs_from_s_range(0.0, 3.14)
          array('d', [10.0, 20.0, 0.0, 14.5, 17.9])
        """

        s_subset = self.s_subset(s_start, s_end)

        if s_subset is NotImplemented:
            return NotImplemented

        return array('d', map(self.z_for_s, s_subset))


class ProfileAttitudes(list):

    def __init__(self, atts: List[ProfileAttitude]):

        check_type(atts, "Attitude georeferenced", List)
        for el in atts:
            check_type(el, "Attitude projection", ProfileAttitude)

        super(ProfileAttitudes, self).__init__(atts)


class ProfilesIntersections(list):

    def __init__(self, intersections: List[ArrayList]):

        check_type(intersections, "Profiles intersections", List)
        for el in intersections:
            check_type(el, "Profile intersections", ArrayList)

        super(ProfilesIntersections, self).__init__(intersections)


