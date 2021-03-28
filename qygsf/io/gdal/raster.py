
from typing import Any, Tuple, Dict, Optional, Union

import os
import gdal
import numpy as np

from pygsf.defaults import *
from pygsf.utils.types import *
from pygsf.defaults import *
from pygsf.mathematics.scalars import *
from pygsf.georeferenced.rasters import *
#from pygsf.georeferenced.crs import *
from pygsf.geometries.rasters.geotransform import *


GRID_NULL_VALUE = -99999  # should already be imported but it isn't


def read_raster(
        file_ref: Any
) -> Tuple[gdal.Dataset, Optional[GeoTransform], int, str]:
    """
    Read a raster layer.

    :param file_ref: the reference to the raster
    :type file_ref: Any
    :return: the dataset, its geotransform, the number of bands, the projection.
    :rtype: tuple made up by a gdal.Dataset instance, an optional Geotransform object, and int and a string.
    :raises: RasterIOException

    Examples:
    """

    # open raster file and check operation success

    dataset = gdal.Open(file_ref, gdal.GA_ReadOnly)
    if not dataset:
        raise RasterIOException("No input data open")

    # get raster descriptive infos

    gt = dataset.GetGeoTransform()
    if gt:
        geotransform = GeoTransform.fromGdalGt(gt)
    else:
        geotransform = None

    num_bands = dataset.RasterCount

    projection = dataset.GetProjection()

    return dataset, geotransform, num_bands, projection


def read_band(
        dataset: gdal.Dataset,
        bnd_ndx: int = 1
) -> Tuple[dict, 'np.array']:
    """
    Read data and metadata of a rasters band based on GDAL.

    :param dataset: the source raster dataset
    :type dataset: gdal.Dataset
    :param bnd_ndx: the index of the band (starts from 1)
    :type bnd_ndx: int
    :return: the band parameters and the data values
    :rtype: dict of data parameters and values as a numpy.array
    :raises: RasterIOException

    Examples:

    """

    band = dataset.GetRasterBand(bnd_ndx)
    data_type = gdal.GetDataTypeName(band.DataType)

    unit_type = band.GetUnitType()

    stats = band.GetStatistics(False, False)
    if stats is None:
        dStats = dict(
            min=None,
            max=None,
            mean=None,
            std_dev=None)
    else:
        dStats = dict(
            min=stats[0],
            max=stats[1],
            mean=stats[2],
            std_dev=stats[3])

    noDataVal = band.GetNoDataValue()

    nOverviews = band.GetOverviewCount()

    colorTable = band.GetRasterColorTable()

    if colorTable:
        nColTableEntries = colorTable.GetCount()
    else:
        nColTableEntries = 0

    # read data from band

    grid_values = band.ReadAsArray()
    if grid_values is None:
        raise RasterIOException("Unable to read data from rasters")

    # transform data into numpy array

    data = np.asarray(grid_values)

    # if nodatavalue exists, set null values to NaN in numpy array
    if noDataVal is not None and np.isfinite(noDataVal):
        data = np.where(abs(data - noDataVal) > 1e-10, data, np.NaN)

    band_params = dict(
        dataType=data_type,
        unitType=unit_type,
        stats=dStats,
        noData=noDataVal,
        numOverviews=nOverviews,
        numColorTableEntries=nColTableEntries)

    return band_params, data


def try_read_raster_band(
    raster_source: str,
    bnd_ndx: int = 1
) -> Tuple[bool, Union[str, Tuple[GeoTransform, str, Dict, 'np.array']]]:

    # get raster parameters and data
    try:
        dataset, geotransform, num_bands, projection = read_raster(raster_source)
    except (IOError, TypeError, RasterIOException) as err:
        return False, "Exception with reading {}: {}".format(raster_source, err)

    band_params, data = read_band(dataset, bnd_ndx)

    return True, (geotransform, projection, band_params, data)


def try_write_esrigrid(
    geoarray: GeoArray,
    outgrid_fn: str,
    esri_nullvalue: numbers.Integral = GRID_NULL_VALUE,
    level_ndx: int = 0
) -> Tuple[bool, str]:
    """
    Writes ESRI ascii grid.
    
    :param geoarray: 
    :param outgrid_fn: 
    :param esri_nullvalue: 
    :param level_ndx: index of the level array to write.
    :type level_ndx: int.
    :return: success and descriptive message
    :rtype: tuple made up by a boolean and a string
    """
    
    outgrid_fn = str(outgrid_fn)

    # checking existence of output slope grid

    if os.path.exists(outgrid_fn):
        return False, "Output grid '{}' already exists".format(outgrid_fn)

    try:
        outputgrid = open(outgrid_fn, 'w')  # create the output ascii file
    except Exception:
        return False, "Unable to create output grid '{}'".format(outgrid_fn)

    if outputgrid is None:
        return False, "Unable to create output grid '{}'".format(outgrid_fn)

    if geoarray.has_rotation:
        return False, "Grid has axes rotations defined"

    cell_size_x = geoarray.src_cellsize_j
    cell_size_y = geoarray.src_cellsize_i

    if not areClose(cell_size_x, cell_size_y):
        return False, "Cell sizes in the x- and y- directions are not similar"

    arr = geoarray.level(level_ndx)
    if arr is None:
        return False, "Array with index {} does not exist".format((level_ndx))

    num_rows, num_cols = arr.shape
    llc_x, llc_y = geoarray.level_llc(level_ndx)

    # writes header of grid ascii file

    outputgrid.write("NCOLS %d\n" % num_cols)
    outputgrid.write("NROWS %d\n" % num_rows)
    outputgrid.write("XLLCORNER %.8f\n" % llc_x)
    outputgrid.write("YLLCORNER %.8f\n" % llc_y)
    outputgrid.write("CELLSIZE %.8f\n" % cell_size_x)
    outputgrid.write("NODATA_VALUE %f\n" % esri_nullvalue)

    esrigrid_outvalues = np.where(np.isnan(arr), esri_nullvalue, arr)

    # output of results

    for i in range(0, num_rows):
        for j in range(0, num_cols):
            outputgrid.write("%.8f " % (esrigrid_outvalues[i, j]))
        outputgrid.write("\n")

    outputgrid.close()

    return True, "Data saved in {}".format(outgrid_fn)


if __name__ == "__main__":

    import doctest
    doctest.testmod()



