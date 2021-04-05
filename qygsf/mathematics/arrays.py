
from typing import Optional, Sequence, List

import numpy as np

from .scalars import *


def arrToTuple(
        arr1D: np.ndarray
) -> Tuple[numbers.Real, ...]:
    """
    Converts a 1D arrays into a tuple of floats.
    It assumes a 1D np.ndarray as input.
    Modified from: https://stackoverflow.com/questions/10016352/convert-numpy-array-to-tuple

    :param arr1D: the 1D-arrays whose components have to be extracted.
    :type arr1D: numpy array.
    :return: a tuple derived from the array values extraction.
    :rtype: tuple of numbers.Real.

    Examples:
      >>> levels = np.array([1,2,3,4,5])
      >>> arrToTuple(levels)
      (1.0, 2.0, 3.0, 4.0, 5.0)
    """

    return tuple(map(float, arr1D))


def toFloats(
        iterable_obj: Sequence[numbers.Real]
) -> List[numbers.Real]:
    """
    Converts an iterable object storing float-compatible values to a list of floats.

    :param iterable_obj: iterable object storing float-compatible values
    :type iterable_obj: iterable storing float-compatible values
    :return:
    :rtype: list of Floats.

    Examples:
      >>> toFloats([1, 2, 3])
      [1.0, 2.0, 3.0]
    """

    return [float(item) for item in iterable_obj]


def arraysAreClose(
        a_array: np.ndarray,
        b_array: np.ndarray,
        rtol: numbers.Real = 1e-012,
        atol: numbers.Real = 1e-12,
        equal_nan: bool = False,
        equal_inf: bool = False
) -> bool:
    """
    Check for equivalence between two numpy arrays.

    :param a_array: first array to be compared.
    :type a_array: numpy array.
    :param b_array: second array to be compared with the first one.
    :type b_array: numpy array.
    :param rtol: relative tolerance.
    :type rtol:
    :param atol: absolute tolerance.
    :type atol:
    :param equal_nan: consider nan values equivalent or not.
    :type equal_nan:
    :param equal_inf: consider inf values equivalent or not.
    :type equal_inf:
    :return: whether the two arrays are close as component values.
    :rtype: bool.

    Examples:
      >>> arraysAreClose(np.array([1,2,3]), np.array([1,2,3]))
      True
      >>> arraysAreClose(np.array([[1,2,3], [4, 5, 6]]), np.array([1,2,3]))
      False
      >>> arraysAreClose(np.array([[1,2,3], [4,5,6]]), np.array([[1,2,3], [4,5,6]]))
      True
      >>> arraysAreClose(np.array([[1,2,np.nan], [4,5,6]]), np.array([[1,2,np.nan], [4,5,6]]))
      False
      >>> arraysAreClose(np.array([[1,2,np.nan], [4,5,6]]), np.array([[1,2,np.nan], [4,5,6]]), equal_nan=True)
      True
    """
    if a_array.shape != b_array.shape:
        return False

    are_close = []
    for a, b in np.nditer([a_array, b_array]):
        are_close.append(areClose(a.item(0), b.item(0), rtol=rtol, atol=atol, equal_nan=equal_nan, equal_inf=equal_inf))

    return all(are_close)


def arraysSameShape(
        a_array: np.ndarray,
        b_array: np.ndarray
) -> bool:
    """
    Checks that two arrays have the same shape.

    :param a_array: first array
    :type a_array: Numpy array.
    :param b_array: second array
    :type b_array: Numpy array.
    :return: whether the two arrays have the same shape.
    :rtype: bool.

    Examples:
      >>> arraysSameShape(np.ones((2,2)), np.ones(4))
      False
      >>> arraysSameShape(np.ones((2,2,3)), np.zeros((2,2,3)))
      True
    """

    return a_array.shape == b_array.shape


def pointSolution(
        a_array: np.ndarray,
        b_array: np.ndarray
) -> Tuple[Optional[numbers.Real], Optional[numbers.Real], Optional[numbers.Real]]:
    """
    Finds a non-unique solution for a set of linear equations.

    :param a_array:
    :type a_array: numpy array.
    :param b_array:
    :type b_array: numpy array.
    :return: an optional tuple of solutions
    :rtype: Tuple[Optional[numbers.Real], Optional[numbers.Real], Optional[numbers.Real]]

    Examples:
    """

    try:

        return np.linalg.lstsq(a_array, b_array, rcond=None)[0]

    except Exception:

        return None, None, None


def svd(xyz_array
        ) -> Optional[Tuple[np.ndarray, np.ndarray, np.ndarray]]:
    """
    Calculates the SVD solution given a Numpy array.

    # modified after:
    # http://stackoverflow.com/questions/15959411/best-fit-plane-algorithms-why-different-results-solved

    :param xyz_array:
    :type xyz_array: numpy array.
    :return:
    :rtype:

    Examples:
    """

    try:
        return np.linalg.svd(xyz_array)
    except Exception:
        return None


if __name__ == "__main__":
    import doctest

    doctest.testmod()
