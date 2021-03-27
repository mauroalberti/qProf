
from typing import List, Any


def sort_by_external_key(
        unsorted_list: List[Any],
        sorting_values: List
):

    zipped = zip(sorting_values, unsorted_list)
    return map(lambda v: v[1], sorted(zipped, key=lambda v: v[0]))