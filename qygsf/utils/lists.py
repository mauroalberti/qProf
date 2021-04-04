
from typing import List, Callable, Any


def find_val(
        func: Callable[[List], Any],
        lst: List
) -> Any:
    """
    Applies a function to a list when not empty,
    otherwise return None.

    :param func: function to apply to the list.
    :type func: function.
    :param lst: list to be processed.
    :type lst: List. May be empty.
    :return: result of function application or None.
    :rtype: Any.
    """

    if lst:
        return func(lst)
    else:
        return None


def list2_to_list(list2):
    """
    input: a list of list
    output: a list
    """

    out_list = []
    for list1 in list2:
        for el in list1:
            out_list.append(el)

    return out_list


def list3_to_list(list3):
    """
    input: a list of list of list
    output: a list
    """

    out_list = []
    for list2 in list3:
        for list1 in list2:
            out_list += list1

    return out_list
