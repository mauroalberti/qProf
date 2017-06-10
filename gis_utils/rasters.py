
from math import ceil, floor
import copy

import numpy as np

from ..gsf.geometry import MIN_SEPARATION_THRESHOLD, Point


class ArrCoord(object):
    """
    2D Array coordinates.
    Manages coordinates in the raster (array) space.

    """

    def __init__(self, ival=0.0, jval=0.0):
        """
        @param  ival:  the i (-y) array coordinate of the point.
        @type  ival:  number or string convertible to float.
        @param  jval:  the j (x) array coordinate of the point.
        @type  jval:  number or string convertible to float.

        @return:  self.
        """
        self._i = float(ival)
        self._j = float(jval)

    def g_i(self):
        """
        Get i (row) coordinate value.

        @return:  the i (-y) array coordinate of the point - float.
        """
        return self._i

    def s_i(self, ival):
        """
        Set i (row) coordinate value.

        @param  ival:  the i (-y) array coordinate of the point.
        @type  ival:  number or string convertible to float.

        @return:  self.
        """
        self._i = float(ival)

    # set property for i
    i = property(g_i, s_i)

    def g_j(self):
        """
        Get j (column) coordinate value.

        @return:  the j (x) array coordinate of the point - float.
        """
        return self._j

    def s_j(self, jval):
        """
        Set j (column) coordinate value.

        @param  jval:  the j (x) array coordinate of the point.
        @type  jval:  number or string convertible to float.

        @return:  self.
        """
        self._j = jval

    j = property(g_j, s_j)

    def grid2geogcoord(self, currGeoGrid):
        currPt_geogr_y = currGeoGrid.domain.trcorner.p_y - self.i * currGeoGrid.cellsize_y
        currPt_geogr_x = currGeoGrid.domain.llcorner.p_x + self.j * currGeoGrid.cellsize_x

        return Point(currPt_geogr_x, currPt_geogr_y)


class RectangularDomain(object):
    """
    Rectangular spatial domain class.

    """

    def __init__(self, pt_llc=None, pt_trc=None):
        """
        Class constructor.

        @param  pt_llc:  lower-left corner of the domain.
        @type  pt_llc:  Point.
        @param  pt_trc:  top-right corner of the domain.
        @type  pt_trc:  Point.

        @return:  RectangularDomain instance.
        """
        self._llcorner = pt_llc
        self._trcorner = pt_trc

    @property
    def llcorner(self):
        """
        Get lower-left corner of the spatial domain.

        @return:  lower-left corner of the spatial domain - Point.
        """
        return self._llcorner

    @property
    def trcorner(self):
        """
        Get top-right corner of the spatial domain.

        @return:  top-right corner of the spatial domain - Point.
        """
        return self._trcorner

    @property
    def xrange(self):
        """
        Get x range of spatial domain.

        @return:  x range - float.
        """
        return self.trcorner.p_x - self.llcorner.p_x

    @property
    def yrange(self):
        """
        Get y range of spatial domain.

        @return:  y range - float.
        """
        return self.trcorner.p_y - self.llcorner.p_y

    @property
    def zrange(self):
        """
        Get z range of spatial domain.

        @return:  z range - float.
        """
        return self.trcorner.p_z - self.llcorner.p_z

    @property
    def horiz_area(self):
        """
        Get horizontal area of spatial domain.

        @return:  area - float.
        """
        return self.xrange * self.yrange


class Grid(object):
    """
    Grid class.
    Stores and manages the most of data and processing.

    """

    def __init__(self, source_filename=None, grid_params=None, grid_data=None):
        """
        Grid class constructor.

        @param  source_filename:  name of file from which data and geo-parameters derive.
        @type  source_filename:  string.
        @param  grid_params:  the geo-parameters of the grid.
        @type  grid_params:  class GDALParameters.
        @param  grid_data:  the array storing the data.
        @type  grid_data:  2D np.array.

        @return:  self.
        """
        self._sourcename = source_filename

        if grid_params is not None:
            pt_llc = grid_params.llcorner
            pt_trc = grid_params.trcorner
        else:
            pt_llc = None
            pt_trc = None

        self._grid_domain = RectangularDomain(pt_llc, pt_trc)

        if grid_data is not None:
            self._grid_data = grid_data.copy()
        else:
            self._grid_data = None

    def s_domain(self, domain):
        """
        Set spatial domain.

        @param  domain:  Spatial domain to be attributed to the current Grid instance.
        @type  domain:  class RectangularDomain.

        @return: self
        """

        del self._grid_domain
        self._grid_domain = copy.deepcopy(domain)

    def g_domain(self):
        """
        Get spatial domain.

        @return: the spatial domain of the current Grid instance - class RectangularDomain.
        """

        return self._grid_domain

    def d_domain(self):
        """
        Delete current spatial domain of the Grid instance.

        @return: self
        """

        del self._grid_domain

    # set property for spatial domain
    domain = property(g_domain, s_domain, d_domain)

    def s_grid_data(self, data_array):
        """
        Set grid data array.

        @param data_array: numpy.array of data values.
        @param type: 2D numpy.array.

        @return: self.
        """

        if self._grid_data is not None:
            del self._grid_data

        self._grid_data = data_array.copy()

    def g_grid_data(self):
        """
        Get grid data array.

        @return: 2D numpy.array.
        """

        return self._grid_data

    def d_grid_data(self):
        """
        Delete grid data array.

        @return: self.
        """

        del self._grid_data

    data = property(g_grid_data, s_grid_data, d_grid_data)

    def grid_extent(self):
        """
        Return the xmin, xmax and ymin, ymax values as a dictionary
        """

        return dict(xmin=self.domain.llcorner.p_x,
                    xmax=self.domain.trcorner.p_x,
                    ymin=self.domain.llcorner.p_y,
                    ymax=self.domain.trcorner.p_y)

    @property
    def xmin(self):

        return self.grid_extent()['xmin']

    @property
    def xmax(self):

        return self.grid_extent()['xmax']

    @property
    def ymin(self):

        return self.grid_extent()['ymin']

    @property
    def ymax(self):

        return self.grid_extent()['ymax']

    @property
    def row_num(self):
        """
        Get row number of the grid domain.

        @return: number of rows of data array - int.
        """

        return np.shape(self.data)[0]

    @property
    def col_num(self):
        """
        Get column number of the grid domain.

        @return: number of columns of data array - int.
        """

        return np.shape(self.data)[1]

    @property
    def cellsize_x(self):
        """
        Get the cell size of the grid in the x direction.

        @return: cell size in the x (j) direction - float.
        """

        return self.domain.xrange / float(self.col_num)

    @property
    def cellsize_y(self):
        """
        Get the cell size of the grid in the y direction.

        @return: cell size in the y (-i) direction - float.
        """

        return self.domain.yrange / float(self.row_num)

    @property
    def cellsize_h(self):
        """
        Get the mean horizontal cell size.

        @return: mean horizontal cell size - float.
        """

        return (self.cellsize_x + self.cellsize_y) / 2.0

    def geog2array_coord(self, curr_Pt):
        """
        Converts from geographic to raster (array) coordinates.

        @param curr_Pt: point whose geographical coordinates will be converted to raster (array) ones.
        @type curr_Pt: Point.

        @return: point coordinates in raster (array) frame - class ArrCoord.
        """
        currArrCoord_grid_i = (self.domain.trcorner.p_y - curr_Pt.p_y) / self.cellsize_y
        currArrCoord_grid_j = (curr_Pt.p_x - self.domain.llcorner.p_x) / self.cellsize_x

        return ArrCoord(currArrCoord_grid_i, currArrCoord_grid_j)

    def x(self):
        """
        Creates an array storing the geographical coordinates of the cell centers along the x axis.
        Direction is from left to right.

        @return: numpy.array, shape: 1 x col_num.
        """

        x_values = self.domain.llcorner.p_x + self.cellsize_x * (0.5 + np.arange(self.col_num))

        return x_values[np.newaxis, :]

    def y(self):
        """
        Creates an array storing the geographical coordinates of the cell centers along the y axis.
        Direction is from top to bottom.

        @return: numpy.array, shape: row_num x 1.
        """

        y_values = self.domain.trcorner.p_y - self.cellsize_y * (0.5 + np.arange(self.row_num))

        return y_values[:, np.newaxis]

    def grad_forward_y(self):
        """
        Return an array representing the forward gradient in the y direction (top-wards), with values scaled by cell size.

        @return: numpy.array, same shape as current Grid instance
        """

        gf = np.zeros(np.shape(self.data)) * np.NaN
        gf[1:, :] = self.data[:-1, :] - self.data[1:, :]

        return gf / float(self.cellsize_y)

    def grad_forward_x(self):
        """
        Return an array representing the forward gradient in the x direction (right-wards), with values scaled by cell size.

        @return: numpy.array, same shape as current Grid instance
        """

        gf = np.zeros(np.shape(self.data), ) * np.NaN
        gf[:, :-1] = self.data[:, 1:] - self.data[:, :-1]

        return gf / float(self.cellsize_x)

    def interpolate_bilinear(self, curr_Pt_array_coord):
        """
        Interpolate the z value at a point, given its array coordinates.
        Interpolation method: bilinear.

        @param curr_Pt_array_coord: array coordinates of the point for which the interpolation will be made.
        @type curr_Pt_array_coord: class ArrCoord.

        @return: interpolated z value - float.
        """

        currPt_cellcenter_i = curr_Pt_array_coord.i - 0.5
        currPt_cellcenter_j = curr_Pt_array_coord.j - 0.5

        assert currPt_cellcenter_i > 0, currPt_cellcenter_j > 0

        grid_val_00 = self.data[int(floor(currPt_cellcenter_i)), int(floor(currPt_cellcenter_j))]
        grid_val_01 = self.data[int(floor(currPt_cellcenter_i)), int(ceil(currPt_cellcenter_j))]
        grid_val_10 = self.data[int(ceil(currPt_cellcenter_i)), int(floor(currPt_cellcenter_j))]
        grid_val_11 = self.data[int(ceil(currPt_cellcenter_i)), int(ceil(currPt_cellcenter_j))]

        delta_i = currPt_cellcenter_i - floor(currPt_cellcenter_i)
        delta_j = currPt_cellcenter_j - floor(currPt_cellcenter_j)

        grid_val_y0 = grid_val_00 + (grid_val_10 - grid_val_00) * delta_i
        grid_val_y1 = grid_val_01 + (grid_val_11 - grid_val_01) * delta_i

        grid_val_interp = grid_val_y0 + (grid_val_y1 - grid_val_y0) * delta_j

        return grid_val_interp

    def intersection_with_surface(self, surf_type, srcPt, srcPlaneAttitude):
        """
        Calculates the intersections (as points) between DEM (the self object) and an analytical surface.
        Currently it works only with planes.

        @param surf_type: type of considered surface (i.e., plane, the only case implemented at present).
        @type surf_type: String.
        @param srcPt: point, expressed in geographical coordinates, that the plane must contain.
        @type srcPt: Point.
        @param srcPlaneAttitude: orientation of the surface (currently only planes).
        @type srcPlaneAttitude: class GPlane.

        @return: tuple of four arrays
        """

        if surf_type == 'plane':

            # closures to compute the geographic coordinates (in x- and y-) of a cell center
            # the grid coordinates of the cell center are expressed by i and j
            grid_coord_to_geogr_coord_x_closure = lambda j: self.domain.llcorner.p_x + self.cellsize_x * (0.5 + j)
            grid_coord_to_geogr_coord_y_closure = lambda i: self.domain.trcorner.p_y - self.cellsize_y * (0.5 + i)

            # arrays storing the geographical coordinates of the cell centers along the x- and y- axes
            cell_center_x_array = self.x()
            cell_center_y_array = self.y()

            ycoords_x, xcoords_y = np.broadcast_arrays(cell_center_x_array, cell_center_y_array)

            #### x-axis direction intersections

            # 2D array of DEM segment parameters
            x_dem_m = self.grad_forward_x()
            x_dem_q = self.data - cell_center_x_array * x_dem_m

            # closure for the planar surface that, given (x,y), will be used to derive z
            plane_z_closure = srcPlaneAttitude.plane_from_geo(srcPt)

            # 2D array of plane segment parameters
            x_plane_m = srcPlaneAttitude.plane_x_coeff()
            x_plane_q = np.array_from_function(self.row_num(), 1, lambda j: 0, grid_coord_to_geogr_coord_y_closure,
                                               plane_z_closure)

            # 2D array that defines denominator for intersections between local segments
            x_inters_denomin = np.where(x_dem_m != x_plane_m, x_dem_m - x_plane_m, np.NaN)

            coincident_x = np.where(x_dem_q != x_plane_q, np.NaN, ycoords_x)

            xcoords_x = np.where(x_dem_m != x_plane_m, (x_plane_q - x_dem_q) / x_inters_denomin, coincident_x)
            xcoords_x = np.where(xcoords_x < ycoords_x, np.NaN, xcoords_x)
            xcoords_x = np.where(xcoords_x >= ycoords_x + self.cellsize_x, np.NaN, xcoords_x)

            #### y-axis direction intersections

            # 2D array of DEM segment parameters
            y_dem_m = self.grad_forward_y()
            y_dem_q = self.data - cell_center_y_array * y_dem_m

            # 2D array of plane segment parameters
            y_plane_m = srcPlaneAttitude.plane_y_coeff()
            y_plane_q = np.array_from_function(1, self.col_num, grid_coord_to_geogr_coord_x_closure, lambda i: 0,
                                               plane_z_closure)

            # 2D array that defines denominator for intersections between local segments
            y_inters_denomin = np.where(y_dem_m != y_plane_m, y_dem_m - y_plane_m, np.NaN)
            coincident_y = np.where(y_dem_q != y_plane_q, np.NaN, xcoords_y)

            ycoords_y = np.where(y_dem_m != y_plane_m, (y_plane_q - y_dem_q) / y_inters_denomin, coincident_y)

            # filter out cases where intersection is outside cell range
            ycoords_y = np.where(ycoords_y < xcoords_y, np.NaN, ycoords_y)
            ycoords_y = np.where(ycoords_y >= xcoords_y + self.cellsize_y, np.NaN, ycoords_y)

            for i in xrange(xcoords_x.shape[0]):
                for j in xrange(xcoords_x.shape[1]):
                    if abs(xcoords_x[i, j] - ycoords_x[i, j]) < MIN_SEPARATION_THRESHOLD and abs(
                                    ycoords_y[i, j] - xcoords_y[i, j]) < MIN_SEPARATION_THRESHOLD:
                        ycoords_y[i, j] = np.NaN

            return xcoords_x, xcoords_y, ycoords_x, ycoords_y

