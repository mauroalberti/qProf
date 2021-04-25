
import numbers

import affine
import numpy as np

from ..geometries.grids.fields import *
from ..geometries.grids.geotransform import *

from ..geometries.shapes.space2d import *
from ..geometries.shapes.space3d import *
from ..geometries.grids.interpolations import *

from .crs import *


class GeoArray:
    """
    GeoArray class.
    Stores and process georeferenced raster data.
    """

    def __init__(self,
                 inGeotransform: GeoTransform,
                 epsg_code: numbers.Integral = -1,
                 inLevels: Optional[List[np.ndarray]] = None
                 ):
        """
        GeoArray class constructor.

        :param  inGeotransform:  the geotransform
        :type  inGeotransform:  GeoTransform.
        :param epsg_code: the projection EPSG code.
        :type epsg_code: numbers.Integral
        :param  inLevels:  the nd-array storing the data.
        :type  inLevels:  np.ndarray.

        :return:  None.

        Examples:
        """

        self._gt = inGeotransform
        self._crs = Crs(epsg_code)
        if inLevels is None:
            self._levels = []
        else:
            self._levels = inLevels

    def geotransform(self):
        """
        Returns geotransform.

        :return: the geotransform.
        :rtype: GeoTransform.
        """

        return self._gt

    @classmethod
    def fromRasterio(cls,
                     array: np.ndarray,
                     affine_transform: affine.Affine,
                     epsg_code: numbers.Integral
                     ):
        """
        Create a GeoArray instance from RasterIO-derived input data.

        :param array: the numpy array with the Geoarray values
        :type array: np.ndarray
        :param affine_transform: the affine transformation from image to geographic coordinates
        :type affine_transform: affine.Affine
        :param epsg_code: the EPSG code
        :type epsg_code: numbers.Integral
        """

        gt = GeoTransform.fromAffine(affine_transform)

        return GeoArray(
            inGeotransform=gt,
            epsg_code=epsg_code,
            inLevels=[array]
        )

    @property
    def crs(self) -> Crs:
        """
        Return the geoarray georeferenced.

        :return: the georeferenced.
        :rtype: Crs.
        """

        return self._crs

    @property
    def epsg_code(self) -> numbers.Integral:
        """
        Return the geoarray georeferenced EPSG code.

        :return: the georeferenced EPSG  code.
        :rtype: numbers.Integral.
        """

        return self.crs.epsg_code

    def define_epsg(self, epsg_cd: numbers.Integral):
        """
        Overwrite the geoarray EPSG code.

        :return:
        """

        if not isinstance(epsg_cd, numbers.Integral):
            raise Exception("Provided EPSG code must be integer")

        self._crs = Crs(epsg_cd)

    def __repr__(self) -> str:
        """
        Represents a GeoArray instance as a shortened text.

        :return: a textual shortened representation of a GeoArray instance.
        :rtype: basestring.
        """

        num_bands = self.levels_num
        epsg_code = self.epsg_code
        bands_txt = ""
        for band_ndx in range(num_bands):
            band = self.level(level_ndx=band_ndx)
            rows, cols = band.shape
            bmin, bmax = band.min(), band.max()
            bands_txt += "\nBand {}: {} rows x {} cols; min: {},  max: {}".format(band_ndx+1, rows, cols, bmin, bmax)

        txt = "GeoArray with {} band(s) - CRS: EPSG: {}\n{}".format(num_bands, epsg_code, bands_txt)

        return txt

    @property
    def src_cellsize_j(self) -> numbers.Real:
        """
        Get the cell size of the geoarray in the x direction.

        :return: cell size in the x (j) direction.
        :rtype: numbers.Real.

        Examples:
        """

        return abs(self._gt.pixWidth)

    @property
    def src_cellsize_i(self) -> numbers.Real:
        """
        Get the cell size of the geoarray in the y direction.

        :return: cell size in the y (-i) direction.
        :rtype: numbers.Real.

        Examples:
        """

        return abs(self._gt.pixHeight)

    @property
    def mean_cellsize(self) -> numbers.Real:
        """
        Get the mean cell size of the geoarray.

        :return: mean cell size.

        Examples:
        """

        return (self.src_cellsize_i + self.src_cellsize_j) / 2.0

    @property
    def levels_num(self) -> numbers.Integral:
        """
        Returns the number of levels (dimensions) of the geoarray.

        :return: number of levels.
        :rtype: numbers.Integral.

        Examples:
          >>> gt = GeoTransform(0, 0, 10, 10)
          >>> GeoArray(gt, -1, [np.array([[1, 2], [3, 4]])]).levels_num
          1
          >>> GeoArray(gt, -1, [np.array([[1, 2], [3, 4]]), np.ones((4, 3, 2))]).levels_num
          2
        """

        return len(self._levels)

    def level(self, level_ndx: numbers.Integral=0):
        """
        Return the array corresponding to the requested level
        if existing else None.

        :param level_ndx: the index of the requested level.
        :type level_ndx: numbers.Integral.
        :return: the array or None.
        :rtype: optional array.

        Examples:
        """

        if 0 <= level_ndx < self.levels_num:
            return self._levels[level_ndx]
        else:
            return None

    def level_shape(self, level_ndx: numbers.Integral=0) -> Optional[Tuple[numbers.Integral, numbers.Integral]]:
        """
        Returns the shape (num. rows and num. columns) of the considered level grid.

        :param level_ndx: index of the level (grid) to consider.
        :type level_ndx: numbers.Integral.
        :return: number of rows and columns of the specific grid.
        :rtype: optional tuple of two int values.

        Examples:
          >>> gt = GeoTransform(0, 0, 10, 10)
          >>> GeoArray(gt, -1, [np.array([[1, 2], [3, 4]])]).level_shape()
          (2, 2)
          >>> GeoArray(gt, -1, [np.array([[1, 2], [3, 4]]), np.ones((4, 3, 2))]).level_shape(1)
          (4, 3, 2)
        """

        if 0 <= level_ndx < self.levels_num:
            return self._levels[level_ndx].shape
        else:
            return None

    def level_llc(self, level_ndx: numbers.Integral = 0) -> Optional[Tuple[numbers.Integral, numbers.Integral]]:
        """
        Deprecated. Use "band_corners_pixcoords" instead.

        Returns the coordinates of the lower-left corner.

        :param level_ndx: index of the level (grid) to consider.
        :type level_ndx: numbers.Integral.
        :return: x and y values of the lower-left corner of the specific grid.
        :rtype: optional tuple of two numbers.Integral values.

        Examples:
        """

        shape = self.level_shape(level_ndx)
        if not shape:
            return None

        llc_i_pix, llc_j_pix = shape[0], 0

        return self.ijPixToxy(llc_i_pix, llc_j_pix)

    def band_corners_pixcoords(self, level_ndx: numbers.Integral = 0) -> \
            Tuple[Tuple[numbers.Real, numbers.Real], Tuple[numbers.Real, numbers.Real], Tuple[numbers.Real, numbers.Real], Tuple[numbers.Real, numbers.Real]]:
        """
        Returns the pixel coordinates of the top-left, top-right, bottom-right and bottom-left band corners.

        :param level_ndx: index of the level (grid) to consider.
        :type level_ndx: numbers.Integral.
        :return: pixel coordinates of the top-left, top-right, bottom-right and bottom-left band corners.
        :rtype: four tuples of numbers.Real pairs.

        Examples:
          >>> gt = GeoTransform(0, 0, 10, 10)
          >>> ga = GeoArray(gt, -1, [np.array([[1, 2, 3], [4, 5, 6]])])
          >>> ga.band_corners_pixcoords()
          ((0.0, 0.0), (0.0, 3.0), (2.0, 3.0), (2.0, 0.0))
        """

        shape = self.level_shape(level_ndx)
        num_rows, num_cols = shape

        top_left_ijpix = (0.0, 0.0)
        top_right_ijpix = (0.0, float(num_cols))
        btm_right_ijpix = (float(num_rows), float(num_cols))
        btm_left_ijpix = (float(num_rows), 0.0)

        return top_left_ijpix, top_right_ijpix, btm_right_ijpix, btm_left_ijpix

    def band_corners_geogcoords(self, level_ndx: numbers.Integral = 0) -> \
            Tuple[Tuple[numbers.Real, numbers.Real], Tuple[numbers.Real, numbers.Real], Tuple[numbers.Real, numbers.Real], Tuple[numbers.Real, numbers.Real]]:
        """
        Returns the geographic coordinates of the top-left, top-right, bottom-right and bottom-left band corners.

        :param level_ndx: index of the level (grid) to consider.
        :type level_ndx: numbers.Integral.
        :return: geographic coordinates of the top-left, top-right, bottom-right and bottom-left band corners.
        :rtype: four tuples of numbers.Real pairs.

        Examples:
          >>> gt = GeoTransform(1500, 3000, 10, 10)
          >>> ga = GeoArray(gt, -1, [np.array([[1, 2, 3], [4, 5, 6]])])
          >>> ga.band_corners_geogcoords()
          ((1500.0, 3000.0), (1530.0, 3000.0), (1530.0, 2980.0), (1500.0, 2980.0))
        """

        top_left_ijpix, top_right_ijpix, btm_right_ijpix, btm_left_ijpix = self.band_corners_pixcoords(level_ndx=level_ndx)

        top_left_geogcoord = self.ijPixToxy(*top_left_ijpix)
        top_right_geogcoord = self.ijPixToxy(*top_right_ijpix)
        btm_right_geogcoord = self.ijPixToxy(*btm_right_ijpix)
        btm_left_geogcoord = self.ijPixToxy(*btm_left_ijpix)

        return top_left_geogcoord, top_right_geogcoord, btm_right_geogcoord, btm_left_geogcoord

    def xyToijArr(self, x: numbers.Real, y: numbers.Real) -> Tuple[numbers.Real, numbers.Real]:
        """
        Converts from geographic to array coordinates.

        :param x: x geographic component.
        :type x: numbers.Real.
        :param y: y geographic component.
        :type y: numbers.Real.
        :return: i and j values referred to array.
        :type: tuple of two numbers.Real values.

        Examples:
        """

        return ijPixToijArray(*xyGeogrToijPix(self._gt, x, y))

    def xyToijPix(self, x: numbers.Real, y: numbers.Real) -> Tuple[numbers.Real, numbers.Real]:
        """
        Converts from geographic to pixel coordinates.

        :param x: x geographic component
        :type x: numbers.Real
        :param y: y geographic component
        :type y: numbers.Real
        :return: i and j values referred to grid.
        :type: tuple of two numbers.Real values

        Examples:
        """

        return xyGeogrToijPix(self._gt, x, y)

    def ijArrToxy(self, i: numbers.Real, j: numbers.Real) -> Tuple[numbers.Real, numbers.Real]:
        """
        Converts from array indices to geographic coordinates.

        :param i: i array component.
        :type i: numbers.Real.
        :param j: j array component.
        :type j: numbers.Real.
        :return: x and y geographic coordinates.
        :type: tuple of two numbers.Real values.

        Examples:
        """

        i_pix, j_pix = ijArrToijPix(i, j)

        return ijPixToxyGeogr(self._gt, i_pix, j_pix)

    def ijPixToxy(self, i: numbers.Real, j: numbers.Real) -> Tuple[numbers.Real, numbers.Real]:
        """
        Converts from grid indices to geographic coordinates.

        :param i: i pixel component.
        :type i: numbers.Real.
        :param j: j pixel component.
        :type j: numbers.Real.
        :return: x and y geographic coordinates.
        :type: tuple of two numbers.Real values.

        Examples:
        """

        return ijPixToxyGeogr(self._gt, i, j)

    @property
    def has_rotation(self) -> bool:
        """
        Determines if a geoarray has axis rotations defined.

        :return: true if there are rotations, false otherwise.
        :rtype: bool.

        Examples:
        """

        return self._gt.has_rotation

    def geotransf_cell_sizes(self) -> Tuple[numbers.Real, numbers.Real]:
        """
        Calculates the geotransformed cell sizes.

        :return: a pair of numbers.Real values, representing the cell sizes in the j and i directions.
        """

        factor = 100

        start_pt = Point2D(*self.ijArrToxy(0, 0))
        end_pt_j = Point2D(*self.ijArrToxy(0, factor))
        end_pt_i = Point2D(*self.ijArrToxy(factor, 0))

        return end_pt_j.distance(start_pt) / factor, end_pt_i.distance(start_pt) / factor

    def xy(self, level_ndx: numbers.Integral=0) -> Optional[Tuple[np.ndarray, np.ndarray]]:
        """
        Returns the two arrays storing respectively the x and the y coordinates
        of the grid cell centers for the chosen level (default is first level).

        :param level_ndx: the index of the
        :return: two arrays storing the geographical coordinates of the grid centers.
        :rtype: tuple made up by two numbers.Real arrays.

        Examples:
        """

        res = self.level_shape(level_ndx)

        if not res:
            return None
        else:
            num_rows, num_cols = res
            return gtToxyCellCenters(self._gt, num_rows, num_cols)

    def interpolate_bilinear(self,
         x: numbers.Real,
         y: numbers.Real,
         level_ndx=0
    ) -> Optional[numbers.Real]:
        """
        Interpolate the z value at a point, given its geographic coordinates.
        Interpolation method: bilinear.

        :param x: x geographic coordinate.
        :type x: numbers.Real.
        :param y: y geographic coordinate.
        :type y: numbers.Real.
        :param level_ndx: the index of the used array.
        :type level_ndx: numbers.Integral.
        :return: the interpolated z value.
        :rtype: optional numbers.Real.

        Examples:
        """

        i, j = self.xyToijArr(x, y)

        return array_bilin_interp(self._levels[level_ndx], i, j)

    def  interpolate_bilinear_point(self,
                                   pt: Point3D,
                                   level_ndx=0
    ) -> Optional[Point3D]:
        """
        Interpolate the z value at a point, returning a Point with elevation extracted from the DEM.
        Interpolation method: bilinear.

        :param pt: the positional point.
        :type pt: Point.
        :param level_ndx: the index of the used array.
        :type level_ndx: numbers.Integral.
        :return: a point with the same x-y position of the input point and with z equal to the interpolated z value.
        :rtype: optional Point.

        Examples:
        """

        check_type(pt, "Input point", Point2D)

        #check_crs(self, pt)

        x, y = pt.x, pt.y

        z = self.interpolate_bilinear(x=x, y=y, level_ndx=level_ndx)

        if z and isfinite(z):
            return Point3D(x, y, z)
        else:
            return None

    def magnitude_field(self, ndx_fx=0, ndx_fy=1) -> 'GeoArray':
        """
        Calculates magnitude field as a geoarray.

        :param ndx_fx: index of x field.
        :type ndx_fx: integer.
        :param ndx_fy: index of y field.
        :type ndx_fy: integer.
        :return: a geoarray storing the magnitude field.
        :rtype: GeoArray.

        Examples:
        """

        magn = magnitude(
            fld_x=self._levels[ndx_fx],
            fld_y=self._levels[ndx_fy])

        return GeoArray(
            inGeotransform=self._gt,
            epsg_code=self.epsg_code,
            inLevels=[magn]
        )

    def orientations(self, ndx_fx=0, ndx_fy=1) -> 'GeoArray':
        """
        Calculates orientations field as a geoarray.

        :param ndx_fx: index of x field.
        :type ndx_fx: integer.
        :param ndx_fy: index of y field.
        :type ndx_fy: integer.
        :return: a geoarray storing the orientation field.
        :rtype: GeoArray.

        Examples:
        """

        orient = orients_d(
            fld_x=self._levels[ndx_fx],
            fld_y=self._levels[ndx_fy])

        return GeoArray(
            inGeotransform=self._gt,
            epsg_code=self.epsg_code,
            inLevels=[orient]
        )

    def divergence_2D(self, ndx_fx=0, ndx_fy=1) -> 'GeoArray':
        """
        Calculates divergence of a 2D field as a geoarray.

        :param ndx_fx: index of x field.
        :type ndx_fx: integer.
        :param ndx_fy: index of y field.
        :type ndx_fy: integer.
        :return: a geoarray storing the divergence field.
        :rtype: GeoArray.

        Examples:
        """

        div = divergence(
            fld_x=self._levels[ndx_fx],
            fld_y=self._levels[ndx_fy],
            cell_size_x=self.src_cellsize_j,
            cell_size_y=self.src_cellsize_i)

        return GeoArray(
            inGeotransform=self._gt,
            epsg_code=self.epsg_code,
            inLevels=[div]
        )

    def curl_module(self, ndx_fx=0, ndx_fy=1) -> 'GeoArray':
        """
        Calculates curl module of a 2D field as a geoarray.

        :param ndx_fx: index of x field.
        :type ndx_fx: integer.
        :param ndx_fy: index of y field.
        :type ndx_fy: integer.
        :return: a geoarray storing the curl module field.
        :rtype: GeoArray.

        Examples:
        """

        curl_m = curl_module(
            fld_x=self._levels[ndx_fx],
            fld_y=self._levels[ndx_fy],
            cell_size_x=self.src_cellsize_j,
            cell_size_y=self.src_cellsize_i)

        return GeoArray(
            inGeotransform=self._gt,
            epsg_code=self.epsg_code,
            inLevels=[curl_m])

    def magnitude_grads(self, axis: str= '', ndx_fx: numbers.Integral=0, ndx_fy: numbers.Integral=1) -> 'GeoArray':
        """
        Calculates the magnitude gradient along the x, y axis or both, of a 2D field as a geoarray.

        :param axis: axis along wich to calculate the gradient, 'x' or 'y', or '' (predefined) for both x and y.
        :type axis: str.
        :param ndx_fx: index of x field.
        :type ndx_fx: integer.
        :param ndx_fy: index of y field.
        :type ndx_fy: integer.
        :return: a geoarray storing the magnitude gradient along the x, y axis (or both) field.
        :rtype: GeoArray.
        :raises: Exception.

        Examples:
        """

        if axis == 'x':
            cell_sizes = [self.src_cellsize_j]
        elif axis == 'y':
            cell_sizes = [self.src_cellsize_i]
        elif axis == '':
            cell_sizes = [self.src_cellsize_j, self.src_cellsize_i]
        else:
            raise Exception("Axis must be 'x' or 'y. '{}' given".format(axis))

        magnitude_gradients = magn_grads(
            fld_x=self._levels[ndx_fx],
            fld_y=self._levels[ndx_fy],
            dir_cell_sizes=cell_sizes,
            axis=axis)

        return GeoArray(
            inGeotransform=self._gt,
            epsg_code=self.epsg_code,
            inLevels=magnitude_gradients)

    def grad_flowlines(self, ndx_fx: numbers.Integral=0, ndx_fy: numbers.Integral=1) -> 'GeoArray':
        """
        Calculates gradient along flow lines.

        :param ndx_fx: index of x field.
        :type ndx_fx: integer.
        :param ndx_fy: index of y field.
        :type ndx_fy: integer.
        :return: a geoarray storing the flowline gradient field
        :rtype: GeoArray
        """

        flowln_grad = magn_grad_along_flowlines(
            fld_x=self._levels[ndx_fx],
            fld_y=self._levels[ndx_fy],
            cell_size_x=self.src_cellsize_j,
            cell_size_y=self.src_cellsize_i)

        return GeoArray(
            inGeotransform=self._gt,
            epsg_code=self.epsg_code,
            inLevels=[flowln_grad])


def point_velocity(
        geoarray: GeoArray,
        pt: Point2D
) -> Tuple[Optional[numbers.Real], Optional[numbers.Real]]:
    """
    Return the velocity components of a 2D-flow field at a point location, based on bilinear interpolation.

    :param geoarray: the flow field expressed as a GeoArray.
    :type geoarray: GeoArray.
    :param pt: the point for which the velocity comnponents are extracted.
    :type pt: Point.
    :return: the x and y velocity components of the flow field at the point location.
    :rtype: tuple of two numbers.Real values.

    Examples:
    """

    x, y = pt.toXY()
    vx = geoarray.interpolate_bilinear(
        x=x,
        y=y,
        level_ndx=0)
    vy = geoarray.interpolate_bilinear(
        x=x,
        y=y,
        level_ndx=1)

    return vx, vy


def interpolate_rkf(
        geoarray: GeoArray,
        delta_time: numbers.Real,
        start_pt: Point2D
    ) -> Tuple[Optional[Point2D], Optional[numbers.Real]]:
    """
    Interpolate point-like object position according to the Runge-Kutta-Fehlberg method.

    :param geoarray: the flow field expressed as a GeoArray.
    :type geoarray: GeoArray.
    :param delta_time: the flow field expressed as a GeoArray.
    :type delta_time: GeoArray.
    :param start_pt: the initial point.
    :type start_pt: Point.
    :return: the estimated point-like object position at the incremented time, with the estimation error.
    :rtype: tuple of optional point and optional numbers.Real.

    Examples:
    """

    check_type(geoarray, "Geoarray", GeoArray)

    check_type(delta_time, "Delta time", numbers.Real)

    check_type(start_pt, "Start point", Point2D)

    k1_vx, k1_vy = point_velocity(geoarray, start_pt)

    if k1_vx is None or k1_vy is None:
        return None, None

    k2_pt = Point2D(
        x=start_pt.x + 0.25 * delta_time * k1_vx,
        y=start_pt.y + 0.25 * delta_time * k1_vy
    )

    k2_vx, k2_vy = point_velocity(geoarray, k2_pt)

    if k2_vx is None or k2_vy is None:
        return None, None

    k3_pt = Point2D(
        x=start_pt.x + (3.0 / 32.0) * delta_time * k1_vx + (9.0 / 32.0) * delta_time * k2_vx,
        y=start_pt.y + (3.0 / 32.0) * delta_time * k1_vy + (9.0 / 32.0) * delta_time * k2_vy
    )

    k3_vx, k3_vy = point_velocity(geoarray, k3_pt)

    if k3_vx is None or k3_vy is None:
        return None, None

    k4_pt = Point2D(
        x=start_pt.x + (1932.0 / 2197.0) * delta_time * k1_vx - (7200.0 / 2197.0) * delta_time * k2_vx + (7296.0 / 2197.0) * delta_time * k3_vx,
        y=start_pt.y + (1932.0 / 2197.0) * delta_time * k1_vy - (7200.0 / 2197.0) * delta_time * k2_vy + (7296.0 / 2197.0) * delta_time * k3_vy
    )

    k4_vx, k4_vy = point_velocity(geoarray, k4_pt)

    if k4_vx is None or k4_vy is None:
        return None, None

    k5_pt = Point2D(
        x=start_pt.x + (439.0 / 216.0) * delta_time * k1_vx - 8.0 * delta_time * k2_vx + (3680.0 / 513.0) * delta_time * k3_vx - (845.0 / 4104.0) * delta_time * k4_vx,
        y=start_pt.y + (439.0 / 216.0) * delta_time * k1_vy - 8.0 * delta_time * k2_vy + (3680.0 / 513.0) * delta_time * k3_vy - (845.0 / 4104.0) * delta_time * k4_vy
    )

    k5_vx, k5_vy = point_velocity(geoarray, k5_pt)

    if k5_vx is None or k5_vy is None:
        return None, None

    k6_pt = Point2D(
        x=start_pt.x - (8.0 / 27.0) * delta_time * k1_vx + 2.0 * delta_time * k2_vx - (3544.0 / 2565.0) * delta_time * k3_vx + (1859.0 / 4104.0) * delta_time * k4_vx - (
                          11.0 / 40.0) * delta_time * k5_vx,
        y=start_pt.y - (8.0 / 27.0) * delta_time * k1_vy + 2.0 * delta_time * k2_vy - (3544.0 / 2565.0) * delta_time * k3_vy + (1859.0 / 4104.0) * delta_time * k4_vy - (
                          11.0 / 40.0) * delta_time * k5_vy
    )

    k6_vx, k6_vy = point_velocity(geoarray, k6_pt)

    if k6_vx is None or k6_vy is None:
        return None, None

    rkf_4o_x = start_pt.x + delta_time * (
            (25.0 / 216.0) * k1_vx + (1408.0 / 2565.0) * k3_vx + (2197.0 / 4104.0) * k4_vx - (
            1.0 / 5.0) * k5_vx)
    rkf_4o_y = start_pt.y + delta_time * (
            (25.0 / 216.0) * k1_vy + (1408.0 / 2565.0) * k3_vy + (2197.0 / 4104.0) * k4_vy - (
            1.0 / 5.0) * k5_vy)
    temp_pt = Point2D(
        x=rkf_4o_x,
        y=rkf_4o_y
    )

    interp_x = start_pt.x + delta_time * (
            (16.0 / 135.0) * k1_vx + (6656.0 / 12825.0) * k3_vx + (28561.0 / 56430.0) * k4_vx - (
            9.0 / 50.0) * k5_vx + (2.0 / 55.0) * k6_vx)
    interp_y = start_pt.y + delta_time * (
            (16.0 / 135.0) * k1_vy + (6656.0 / 12825.0) * k3_vy + (28561.0 / 56430.0) * k4_vy - (
            9.0 / 50.0) * k5_vy + (2.0 / 55.0) * k6_vy)
    interp_pt = Point2D(
        x=interp_x,
        y=interp_y
    )

    interp_pt_error_estim = interp_pt.distance(temp_pt)

    return interp_pt, interp_pt_error_estim


def line_on_grid(
        ga: GeoArray,
        profile_line: Line2D
) -> Optional[Line2D]:
    """
    Calculates a line draped on a grid.

    :param ga: geoarray
    :type ga: GeoArray.
    :param profile_line: the profile line.
    :type profile_line: Line2D
    :return: the profile.
    :rtype: Optional[Line].
    """

    lnProfile = Line2D()

    for point in profile_line.pts():

        z = ga.interpolate_bilinear(point.x, point.y)
        if z:
            lnProfile.add_pt(
                Point2D(
                    x=point.x,
                    y=point.y
                )
            )

    return lnProfile


def ijarr2xyz(
        ijarr2xy_func: Callable,
        xy2z_func: Callable,
        i: numbers.Real,
        j: numbers.Real
) -> Tuple[numbers.Real, numbers.Real, numbers.Real]:
    """
    Return a tuple of (x, y, z) values, starting by array indices.

    :param ijarr2xy_func: a function converting from array to geographic coordinates.
    :param xy2z_func: a callable converting from x, y geographic coordinates to a z value.
    :param i: i index.
    :param j: j index.
    :return: Point
    """

    x, y = ijarr2xy_func(i, j)
    z = xy2z_func(x, y)
    return x, y, z


def xyarr2segmentslope(
        xy2z_func: Callable,
        arrij2xy_func: Callable,
        i: numbers.Real,
        j: numbers.Real,
        i_start=0.0,
        j_start=0.0
) -> numbers.Real:
    """
    Calculates the segment slope along a gridded direction defined by its end point i, j array coordinates.
    Assumed start point is array coordinates 0, 0.

    :param xy2z_func: a callable deriving a z value from geographic x-y coordinates..
    :param arrij2xy_func: a function converting from array coordinates to geographic coordinates.
    :param i: i index of end point.
    :param j: j index of end point.
    :param i_start: i index of start point. Default is 0.0.
    :param j_start:j index of start point. Default is 0.0.
    :return: segment slope.
    :rtype: numbers.Real.
    """

    start_point = Point3D(*ijarr2xyz(
        ijarr2xy_func=arrij2xy_func,
        xy2z_func=xy2z_func,
        i=i_start,
        j=j_start))

    end_point = Point3D(*ijarr2xyz(
        ijarr2xy_func=arrij2xy_func,
        xy2z_func=xy2z_func,
        i=i,
        j=j))

    return Segment3D(start_point, end_point).ratio_delta_zs()


def segment_intersections_array(
        m_arr1: np.ndarray,
        m_arr2: np.ndarray,
        q_arr1: np.ndarray,
        q_arr2: np.ndarray,
        cell_size: numbers.Real,
        m_delta_tol: Optional[numbers.Real] = 1e-6,
        q_delta_tol: Optional[numbers.Real] = 1e-6
) -> np.ndarray:
    """
    Creates array that gives the residual index [0-1[ of the intersection between segments along the considered
    array axis (i or j) whose m (slope) and q (y-axis intersection) values along the considered array axis (i or j)
    are defined in the two pairs of input arrays.

    :param m_arr1: array storing values of grid 1 segment slopes.
    :param m_arr2: array storing values of grid 2 segment slopes.
    :param q_arr1: array storing values of grid 1 segment y-axis intersections.
    :param q_arr2: array storing values of grid 2 segment y-axis intersections.
    :param cell_size: cell size of the two io along the considered direction. Required the same in the two io.
    :param m_delta_tol: optional tolerance for delta between grid 1 and grid 2 segment slopes.
    :param q_delta_tol: optional tolerance for delta between grid 1 and grid 2 segment y-axis intersections.
    :return: array with values of intersection residual indices [0 - 1[
    """

    # if segments slope are not sub-equal, we calculate the intersection residual slope using the required formula

    inters_residual_indices = np.where(abs(m_arr1 - m_arr2) < m_delta_tol, np.NaN,
                                       (q_arr2 - q_arr1) / (cell_size * (m_arr1 - m_arr2)))

    # if the elevations at the left cell center are sub-equal,
    # the residual index is set to 0.0, i.e. there is an intersection at the left cell

    inters_with_coincident_starts = np.where(abs(q_arr1 - q_arr2) < q_delta_tol, 0.0, inters_residual_indices)

    # we filter out residual indices that are not intersect cell span, i.e., not between 0.0 (included) and 1.0 (excluded)

    inters_intracells_residuals = np.where(
        np.logical_and(inters_with_coincident_starts >= 0.0, inters_with_coincident_starts < 1.0),
        inters_with_coincident_starts, np.NaN)

    return inters_intracells_residuals


def arrayTo3DPts(
        direction: str,
        arr: np.ndarray,
        ij2xy_func: Callable,
        xy2z_func: Callable
) -> List[Point3D]:
    """
    Converts an array of along-direction (i- or j-) intra-cell segments [0 -> 1[ into
    a list of 3D points.

    :param direction: considered intersection direction: 'i' (for i axis) or 'j' (for j axis).
    :param arr: array of along-direction (i- or j-) intra-cell segments [0 -> 1[.
    :param ij2xy_func: function to convert from array indices to x-y geographic coordinates.
    :param xy2z_func: function that calculates z value given x and y coordinates.
    :return: list of 3D points.
    :raise: Exception when direction is not 'i' or 'j'
    """

    pts = []
    for i in range(arr.shape[0]):
        for j in range(arr.shape[1]):
            val = arr[i, j]
            if np.isfinite(val):
                if direction == 'i':
                    i_int, j_int = i + val, j
                elif direction == 'j':
                    i_int, j_int = i, j + val
                else:
                    raise Exception('Unexpected array direction value: {}'.format(direction))
                x, y = ij2xy_func(i_int, j_int)
                z = xy2z_func(x, y)
                pts.append(Point3D(x, y, z))

    return pts


def plane_dem_intersection(
        srcPlaneAttitude: Plane,
        srcPt: Point3D,
        geo_array: GeoArray,
        level_ndx: numbers.Integral = 0) -> List[Point3D]:
    """
    Calculates the intersections (as points) between the grid and a planar analytical surface.

    :param srcPlaneAttitude: orientation of the surface (currently only planes).
    :param srcPt: point that the plane must contain.
    :param geo_array: the input GeoArray storing the used grid.
    :param level_ndx: the grid level to use from the provided geoarray. Default is first (index equal to zero).
    :return: list of unique intersecting points.

    Examples:
    """

    # dem values as a Numpy array

    q_d = geo_array.level(
        level_ndx=level_ndx)

    # row and column numbers of the dem

    row_num, col_num = q_d.shape

    # plane closure that, given (x, y), derive z

    plane_z_closure = closure_plane_from_geo(
        srcPlaneAttitude,
        srcPt
    )

    # plane elevations at grid cell centers

    q_p = array_from_geotransform_function(
        row_num=row_num,
        col_num=col_num,
        geotransform=geo_array.geotransform(),
        z_transfer_func=plane_z_closure)

    index_multiplier = 100  # sufficiently large value to ensure a precise slope values

    mi_p = xyarr2segmentslope(
        xy2z_func=plane_z_closure,
        arrij2xy_func=geo_array.ijArrToxy,
        i=index_multiplier,
        j=0) * np.ones((row_num, col_num))

    mj_p = xyarr2segmentslope(
        xy2z_func=plane_z_closure,
        arrij2xy_func=geo_array.ijArrToxy,
        i=0,
        j=index_multiplier) * np.ones((row_num, col_num))

    # 2D array of DEM segment parameters

    cell_size_j, cell_size_i = geo_array.geotransf_cell_sizes()

    mj_d = grad_j(
        fld=q_d,
        cell_size_j=cell_size_j)

    mi_d = grad_iminus(
        fld=q_d,
        cell_size_i=cell_size_i)

    # intersection points

    intersection_pts_j = segment_intersections_array(
        m_arr1=mj_d,
        m_arr2=mj_p,
        q_arr1=q_d,
        q_arr2=q_p,
        cell_size=cell_size_j)

    intersection_pts_j = arrayTo3DPts(
        direction='j',
        arr=intersection_pts_j,
        ij2xy_func=geo_array.ijArrToxy,
        xy2z_func=plane_z_closure)

    intersection_pts_i = segment_intersections_array(
        m_arr1=mi_d,
        m_arr2=mi_p,
        q_arr1=q_d,
        q_arr2=q_p,
        cell_size=cell_size_i)

    intersection_pts_i = arrayTo3DPts(
        direction='i',
        arr=intersection_pts_i,
        ij2xy_func=geo_array.ijArrToxy,
        xy2z_func=plane_z_closure)

    unique_pts = intersection_pts_j + intersection_pts_i

    return unique_pts