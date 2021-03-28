
from typing import List, Dict, Union

import numpy as np


def get_statistics(
        vals: Union[List, np.ndarray]
) -> Dict:
    """

    :param vals: the values, as a list or a numpy array
    :type vals: list or numpy array.
    :return: the statistics values.
    :rtype: a dictionary.
    """

    array = np.asarray(vals)

    amin = np.nanmin(array)
    amax = np.nanmax(array)
    mean = np.nanmean(array)
    var = np.nanvar(array)
    std = np.nanstd(array)

    stats = dict(min=amin,
                 max=amax,
                 mean=mean,
                 var=var,
                 std=std)

    return stats

