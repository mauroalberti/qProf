
import numpy as np


def get_statistics(topo_array):

    min = np.nanmin(topo_array)
    max = np.nanmax(topo_array)
    mean = np.nanmean(topo_array)
    var = np.nanvar(topo_array)
    std = np.nanstd(topo_array)

    stats = dict(min=min,
                 max=max,
                 mean=mean,
                 var=var,
                 std=std)

    return stats