
import numpy as np


def get_statistics(topo_array):

    stat_min = np.nanmin(topo_array)
    stat_max = np.nanmax(topo_array)
    stat_mean = np.nanmean(topo_array)
    stat_var = np.nanvar(topo_array)
    stat_std = np.nanstd(topo_array)

    stats = dict(min=stat_min,
                 max=stat_max,
                 mean=stat_mean,
                 var=stat_var,
                 std=stat_std)

    return stats