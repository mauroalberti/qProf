# -*- coding: utf-8 -*-


from __future__ import division

import json
import numpy as np
import os
from osgeo import ogr, osr

from .errors import RasterParametersException, OGRIOException, AnaliticSurfaceIOException
from .spatial import CartesianPoint3DT


class GDALParameters(object):
    """
    Manage GDAL parameters from rasters.
    
    """

    # class constructor
    def __init__(self):
        """
        Class constructor.
        
        @return:  generic-case GDAL parameters.
        """
        self._nodatavalue = None
        self._topleftX = None
        self._topleftY = None
        self._pixsizeEW = None
        self._pixsizeNS = None
        self._rows = None
        self._cols = None
        self._rotation_GT_2 = 0.0
        self._rotation_GT_4 = 0.0

    def s_noDataValue(self, nodataval):
        """
        Set raster no data value.
        
        @param  nodataval:  the raster no-data value.
        @type  nodataval:  None, otherwise number or string convertible to float.
                
        @return:  self.         
        """

        try:
            self._nodatavalue = float(nodataval)
        except:
            self._nodatavalue = None

    def g_noDataValue(self):
        """
        Get raster no-data value.
        
        @return:  no-data value - float.
        """
        return self._nodatavalue

    # set property for no-data value
    noDataValue = property(g_noDataValue, s_noDataValue)

    def s_topLeftX(self, topleftX):
        """
        Set top-left corner x value of the raster.
        
        @param  topleftX:  the top-left corner x value, according to GDAL convention.
        @type  topleftX:  number or string convertible to float.
                
        @return:  self.
        """
        self._topleftX = float(topleftX)

    def g_topLeftX(self):
        """
        Get top-left corner x value of the raster.
        
        @return:  the top-left corner x value, according to GDAL convention - float.
        """
        return self._topleftX

    # set property for topleftX
    topLeftX = property(g_topLeftX, s_topLeftX)

    def s_topLeftY(self, topleftY):
        """
        Set top-left corner y value of the raster.
        
        @param  topleftY:  the top-left corner y value, according to GDAL convention.
        @type  topleftY:  number or string convertible to float.
                
        @return:  self.
        """
        self._topleftY = float(topleftY)

    def g_topLeftY(self):
        """
        Get top-left corner y value of the raster.
        
        @return:  the top-left corner y value, according to GDAL convention - float.
        """
        return self._topleftY

    # set property for topleftY
    topLeftY = property(g_topLeftY, s_topLeftY)

    def s_pixSizeEW(self, pixsizeEW):
        """
        Set East-West size of the raster cell.
        
        @param  pixsizeEW:  the top-left y value, according to GDAL convention.
        @type  pixsizeEW:  number or string convertible to float.
                
        @return:  self.
        """
        self._pixsizeEW = float(pixsizeEW)

    def g_pixSizeEW(self):
        """
        Get East-West size of the raster cell.
        
        @return:  the East-West size of the raster cell - float.        
        """
        return self._pixsizeEW

    # set property for topleftY
    pixSizeEW = property(g_pixSizeEW, s_pixSizeEW)

    # pixsizeNS

    def s_pixSizeNS(self, pixsizeNS):
        """
        Set North-South size of the raster cell.
        
        @param  pixsizeNS:  the North-South size of the raster cell.
        @type  pixsizeNS:  number or string convertible to float.
                
        @return:  self.
        """
        self._pixsizeNS = float(pixsizeNS)

    def g_pixSizeNS(self):
        """
        Get North-South size of the raster cell.
        
        @return:  the North-South size of the raster cell - float.        
        """
        return self._pixsizeNS

    # set property for topleftY
    pixSizeNS = property(g_pixSizeNS, s_pixSizeNS)

    def s_rows(self, rows):
        """
        Set row number.
        
        @param  rows:  the raster row number.
        @type  rows:  number or string convertible to int.
                
        @return:  self.
        """
        self._rows = int(rows)

    def g_rows(self):
        """
        Get row number.
        
        @return:  the raster row number - int.        
        """
        return self._rows

    # set property for rows
    rows = property(g_rows, s_rows)

    def s_cols(self, cols):
        """
        Set column number.
        
        @param  cols:  the raster column number.
        @type  cols:  number or string convertible to int.
                
        @return:  self.
        """
        self._cols = int(cols)

    def g_cols(self):
        """
        Get column number.
        
        @return:  the raster column number - int. 
        """
        return self._cols

    # set property for cols
    cols = property(g_cols, s_cols)

    def s_rotation_GT_2(self, rotation_GT_2):
        """
        Set rotation GT(2) (see GDAL documentation).
        
        @param  rotation_GT_2:  the raster rotation value GT(2).
        @type  rotation_GT_2:  number or string convertible to float.
                
        @return:  self.
        """
        self._rotation_GT_2 = float(rotation_GT_2)

    def g_rotation_GT_2(self):
        """
        Get rotation GT(2) (see GDAL documentation).
        
        @return:  the raster rotation value GT(2). - float. 
        """
        return self._rotation_GT_2

    # set property for rotation_GT_2
    rotGT2 = property(g_rotation_GT_2, s_rotation_GT_2)

    def s_rotation_GT_4(self, rotation_GT_4):
        """
        Set rotation GT(4) (see GDAL documentation)
        
        @param  rotation_GT_4:  the raster rotation value GT(4).
        @type  rotation_GT_4:  number or string convertible to float.
                
        @return:  self.
        """
        self._rotation_GT_4 = float(rotation_GT_4)

    def g_rotation_GT_4(self):
        """
        Get rotation GT(4) (see GDAL documentation).

        @return:  the raster rotation value GT(4) - float. 
        """
        return self._rotation_GT_4

    # set property for rotation_GT_4
    rotGT4 = property(g_rotation_GT_4, s_rotation_GT_4)

    def check_params(self, tolerance=1e-06):
        """
        Check absence of axis rotations or pixel size differences in the raster band.
        
        @param  tolerance:  the maximum threshold for both pixel N-S and E-W difference, or axis rotations.
        @type  tolerance:  float.
                
        @return:  None when successful, RasterParametersException when pixel differences or axis rotations.
        
        @raise: RasterParametersException - raster geometry incompatible with this module (i.e. different cell sizes or axis rotations).          
        """
        # check if pixel size can be considered the same in the two axis directions
        if abs(abs(self._pixsizeEW) - abs(self._pixsizeNS)) / abs(self._pixsizeNS) > tolerance:
            raise RasterParametersException, 'Pixel sizes in x and y directions are different in raster'

            # check for the absence of axis rotations
        if abs(self._rotation_GT_2) > tolerance or abs(self._rotation_GT_4) > tolerance:
            raise RasterParametersException, 'There should be no axis rotation in raster'

        return

    def llcorner(self):
        """
        Creates a point at the lower-left corner of the raster.
        
        @return:  new CartesianPoint3DT instance.
        """
        return CartesianPoint3DT(self.topLeftX, self.topLeftY - abs(self.pixSizeNS) * self.rows)

    def trcorner(self):
        """
        Create a point at the top-right corner of the raster.

        @return:  new CartesianPoint3DT instance.
        """
        return CartesianPoint3DT(self.topLeftX + abs(self.pixSizeEW) * self.cols, self.topLeftY)

    def geo_equiv(self, other, tolerance=1.0e-6):
        """
        Checks if two rasters are geographically equivalent.

        @param  other:  a grid to be compared with self.
        @type  other:  Grid instance.
        @param  tolerance:  the maximum threshold for pixel sizes, topLeftX or topLeftY differences.
        @type  tolerance:  float.
               
        @return:  Boolean.        
        """
        if 2 * (self.topLeftX - other.topLeftX) / (self.topLeftX + other.topLeftX) > tolerance or \
                                        2 * (self.topLeftY - other.topLeftY) / (
                                    self.topLeftY + other.topLeftY) > tolerance or \
                                        2 * (abs(self.pixSizeEW) - abs(other.pixSizeEW)) / (
                                    abs(self.pixSizeEW) + abs(other.pixSizeEW)) > tolerance or \
                                        2 * (abs(self.pixSizeNS) - abs(other.pixSizeNS)) / (
                                    abs(self.pixSizeNS) + abs(other.pixSizeNS)) > tolerance or \
                        self.rows != other.rows or self.cols != other.cols or self.projection != other.projection:
            return False
        else:
            return True


class QGisRasterParameters(object):
    # class constructor
    def __init__(self, name, cellsizeEW, cellsizeNS, rows, cols, xMin, xMax, yMin, yMax, nodatavalue, crs):

        self.name = name
        self.cellsizeEW = cellsizeEW
        self.cellsizeNS = cellsizeNS
        self.rows = rows
        self.cols = cols
        self.xMin = xMin
        self.xMax = xMax
        self.yMin = yMin
        self.yMax = yMax
        self.nodatavalue = nodatavalue
        self.crs = crs

    def point_in_dem_area(self, point):

        if self.xMin <= point.p_x <= self.xMax and \
                                self.yMin <= point.p_y <= self.yMax:
            return True
        else:
            return False

    def point_in_interpolation_area(self, point):

        if self.xMin + self.cellsizeEW / 2.0 <= point.p_x <= self.xMax - self.cellsizeEW / 2.0 and \
                                        self.yMin + self.cellsizeNS / 2.0 <= point.p_y <= self.yMax - self.cellsizeNS / 2.0:
            return True
        else:
            return False

    def geogr2raster(self, point):

        x = (point.p_x - (self.xMin + self.cellsizeEW / 2.0)) / self.cellsizeEW
        y = (point.p_y - (self.yMin + self.cellsizeNS / 2.0)) / self.cellsizeNS

        return dict(x=x, y=y)

    def raster2geogr(self, array_dict):

        point = CartesianPoint3DT()
        point._x = self.xMin + (array_dict['x'] + 0.5) * self.cellsizeEW
        point._y = self.yMin + (array_dict['y'] + 0.5) * self.cellsizeNS

        return point


def read_line_shapefile_via_ogr(line_shp_path):
    """
    Read line shapefile using OGR.

    @param  line_shp_path:  parameter to check.
    @type  line_shp_path:  QString or string
    
    """
    # reset layer parameters 

    if line_shp_path is None or line_shp_path == '':
        return dict(success=False, error_message='No input path')

        # open input vector layer
    shape_driver = ogr.GetDriverByName("ESRI Shapefile")

    line_shape = shape_driver.Open(str(line_shp_path), 0)

    # layer not read
    if line_shape is None:
        return dict(success=False, error_message='Unable to open input shapefile')

        # get internal layer
    lnLayer = line_shape.GetLayer(0)

    # set vector layer extent   
    layer_extent = lnLayer.GetExtent()
    lines_extent = {'xmin': layer_extent[0], 'xmax': layer_extent[1], 'ymin': layer_extent[2], 'ymax': layer_extent[3]}

    # initialize lists storing vertex coordinates of line
    lines_points = []

    # start reading layer features        
    curr_line = lnLayer.GetNextFeature()

    # loop in layer features                 
    while curr_line:

        line_points = []

        line_geom = curr_line.GetGeometryRef()

        if line_geom is None:
            line_shape.Destroy()
            return dict(success=False, error_message='No geometry ref')

        if line_geom.GetGeometryType() != ogr.wkbLineString and \
                        line_geom.GetGeometryType() != ogr.wkbMultiLineString:
            line_shape.Destroy()
            return dict(success=False, error_message='Not a linestring/multilinestring')

        for i in range(line_geom.GetPointCount()):
            x, y, z = line_geom.GetX(i), line_geom.GetY(i), line_geom.GetZ(i)

            line_points.append(CartesianPoint3DT(x, y, z))

        lines_points.append(line_points)

        curr_line = lnLayer.GetNextFeature()

    line_shape.Destroy()

    return dict(success=True, extent=lines_extent, vertices=lines_points)


def geosurface_export_vtk(output_filepath, geodata):
    geosurface_XYZ, grid_dims = geodata
    X, Y, Z = geosurface_XYZ

    X_arr = np.array(X, dtype=float)
    Y_arr = np.array(Y, dtype=float)
    Z_arr = np.array(Z, dtype=float)

    n_points = np.size(X_arr)

    n_rows, n_cols = grid_dims

    with open(output_filepath, 'w') as outfile:

        outfile.write('# vtk DataFile Version 2.0\n')
        outfile.write('Geosurface - qgSurf vers. 0.3.0\n')
        outfile.write('ASCII\n')
        outfile.write('\nDATASET POLYDATA\n')

        outfile.write('POINTS %d float\n' % n_points)
        for n in xrange(n_points):
            outfile.write('%.4f %.4f %.4f\n' % (X_arr[n], Y_arr[n], Z_arr[n]))

        outfile.write('\n')

        outfile.write('TRIANGLE_STRIPS %d %d\n' % (n_cols - 1, (n_cols - 1) * (1 + n_rows * 2)))

        num_points_strip = n_rows * 2
        for l in xrange(n_cols - 1):
            triangle_strip_string = "%d " % num_points_strip
            for p in xrange(n_rows):
                triangle_strip_string += "%d %d " % ((l + 1) * n_rows + p, l * n_rows + p)
            triangle_strip_string += "\n"
            outfile.write(triangle_strip_string)


def geosurface_export_grass(output_filepath, geodata):
    # Save in Grass format

    geosurface_XYZ, grid_dims = geodata
    X, Y, Z = geosurface_XYZ

    X_arr = np.array(X, dtype=float)
    Y_arr = np.array(Y, dtype=float)
    Z_arr = np.array(Z, dtype=float)

    n_rows, n_cols = grid_dims

    with open(output_filepath, 'w') as outfile:
        outfile.write('VERTI:\n')
        for l in xrange(n_cols - 1):
            for p in xrange(n_rows - 1):
                start_point_ndx = l * n_rows + p
                forward_line_point_ndx = start_point_ndx + n_rows
                outfile.write('F 4\n')
                outfile.write(
                    ' %.4f %.4f %.4f\n' % (X_arr[start_point_ndx], Y_arr[start_point_ndx], Z_arr[start_point_ndx]))
                outfile.write(' %.4f %.4f %.4f\n' % (
                    X_arr[start_point_ndx + 1], Y_arr[start_point_ndx + 1], Z_arr[start_point_ndx + 1]))
                outfile.write(' %.4f %.4f %.4f\n' % (
                    X_arr[forward_line_point_ndx], Y_arr[forward_line_point_ndx], Z_arr[forward_line_point_ndx]))
                outfile.write(
                    ' %.4f %.4f %.4f\n' % (X_arr[start_point_ndx], Y_arr[start_point_ndx], Z_arr[start_point_ndx]))
                outfile.write('F 4\n')
                outfile.write(' %.4f %.4f %.4f\n' % (
                    X_arr[forward_line_point_ndx], Y_arr[forward_line_point_ndx], Z_arr[forward_line_point_ndx]))
                outfile.write(' %.4f %.4f %.4f\n' % (
                    X_arr[start_point_ndx + 1], Y_arr[start_point_ndx + 1], Z_arr[start_point_ndx + 1]))
                outfile.write(' %.4f %.4f %.4f\n' % (
                    X_arr[forward_line_point_ndx + 1], Y_arr[forward_line_point_ndx + 1],
                    Z_arr[forward_line_point_ndx + 1]))
                outfile.write(' %.4f %.4f %.4f\n' % (
                    X_arr[forward_line_point_ndx], Y_arr[forward_line_point_ndx], Z_arr[forward_line_point_ndx]))


def geosurface_export_esri_generate(output_filepath, geodata):
    # Save geosurface (GAS) in Esri generate  format

    geosurface_XYZ, grid_dims = geodata
    X, Y, Z = geosurface_XYZ

    X_arr = np.array(X, dtype=float)
    Y_arr = np.array(Y, dtype=float)
    Z_arr = np.array(Z, dtype=float)

    n_rows, n_cols = grid_dims

    progr_id = 0
    with open(output_filepath, 'w') as outfile:
        outfile.write('VERTI:\n')
        for l in xrange(n_cols - 1):
            for p in xrange(n_rows - 1):
                start_point_ndx = l * n_rows + p
                forward_line_point_ndx = start_point_ndx + n_rows
                progr_id += 1
                outfile.write('%d\n' % progr_id)
                outfile.write(
                    ' %.4f %.4f %.4f\n' % (X_arr[start_point_ndx], Y_arr[start_point_ndx], Z_arr[start_point_ndx]))
                outfile.write(' %.4f %.4f %.4f\n' % (
                    X_arr[start_point_ndx + 1], Y_arr[start_point_ndx + 1], Z_arr[start_point_ndx + 1]))
                outfile.write(' %.4f %.4f %.4f\n' % (
                    X_arr[forward_line_point_ndx], Y_arr[forward_line_point_ndx], Z_arr[forward_line_point_ndx]))
                outfile.write(
                    ' %.4f %.4f %.4f\n' % (X_arr[start_point_ndx], Y_arr[start_point_ndx], Z_arr[start_point_ndx]))
                outfile.write('END\n')
                progr_id += 1
                outfile.write('%d\n' % progr_id)
                outfile.write(' %.4f %.4f %.4f\n' % (
                    X_arr[forward_line_point_ndx], Y_arr[forward_line_point_ndx], Z_arr[forward_line_point_ndx]))
                outfile.write(' %.4f %.4f %.4f\n' % (
                    X_arr[start_point_ndx + 1], Y_arr[start_point_ndx + 1], Z_arr[start_point_ndx + 1]))
                outfile.write(' %.4f %.4f %.4f\n' % (
                    X_arr[forward_line_point_ndx + 1], Y_arr[forward_line_point_ndx + 1],
                    Z_arr[forward_line_point_ndx + 1]))
                outfile.write(' %.4f %.4f %.4f\n' % (
                    X_arr[forward_line_point_ndx], Y_arr[forward_line_point_ndx], Z_arr[forward_line_point_ndx]))
                outfile.write('END\n')
        outfile.write('END\n')


def geosurface_save_gas(output_filepath, geodata):
    with open(output_filepath, 'w') as outfile:
        json.dump(geodata, outfile)


def shapefile_create_def_field(field_def):
    fieldDef = ogr.FieldDefn(field_def['name'], field_def['ogr_type'])
    if field_def['ogr_type'] == ogr.OFTString:
        fieldDef.SetWidth(field_def['width'])

    return fieldDef


def shapefile_create(path, geom_type, fields_dict_list, crs=None):
    """
    crs_prj4: projection in Proj4 text format
    geom_type = OGRwkbGeometryType: ogr.wkbPoint, ....  
    list of:  
        field dict: 'name', 
                    'type': ogr.OFTString,
                            ogr.wkbLineString, 
                            ogr.wkbLinearRing,
                            ogr.wkbPolygon,
                        
                    'width',    
    """

    driver = ogr.GetDriverByName("ESRI Shapefile")

    outShapefile = driver.CreateDataSource(str(path))
    if outShapefile is None:
        raise OGRIOException, 'Unable to save shapefile in provided path'

    if crs is not None:
        spatialReference = osr.SpatialReference()
        spatialReference.ImportFromProj4(crs)
        outShapelayer = outShapefile.CreateLayer("layer", geom_type, spatialReference)
    else:
        outShapelayer = outShapefile.CreateLayer("layer", geom_type=geom_type)

    map(lambda field_def_params: outShapelayer.CreateField(shapefile_create_def_field(field_def_params)),
        fields_dict_list)

    return outShapefile, outShapelayer


def ogr_get_solution_shapefile(path, fields_dict_list):
    driver = ogr.GetDriverByName("ESRI Shapefile")

    dataSource = driver.Open(str(path), 0)

    if dataSource is None:
        raise OGRIOException, 'Unable to open shapefile in provided path'

    point_shapelayer = dataSource.GetLayer()

    prev_solution_list = []
    in_point = point_shapelayer.GetNextFeature()
    while in_point:
        rec_id = int(in_point.GetField('id'))
        x = in_point.GetField('x')
        y = in_point.GetField('y')
        z = in_point.GetField('z')
        dip_dir = in_point.GetField('dip_dir')
        dip_ang = in_point.GetField('dip_ang')
        descript = in_point.GetField('descript')
        prev_solution_list.append([rec_id, x, y, z, dip_dir, dip_ang, descript])
        in_point.Destroy()
        in_point = point_shapelayer.GetNextFeature()

    dataSource.Destroy()

    if os.path.exists(path):
        driver.DeleteDataSource(str(path))

    outShapefile, outShapelayer = shapefile_create(path, ogr.wkbPoint25D, fields_dict_list, crs=None)
    return outShapefile, outShapelayer, prev_solution_list


def ogr_write_point_result(point_shapelayer, field_list, rec_values_list2):
    outshape_featdef = point_shapelayer.GetLayerDefn()

    for rec_value_list in rec_values_list2:

        # pre-processing for new feature in output layer
        curr_Pt_geom = ogr.Geometry(ogr.wkbPoint25D)
        curr_Pt_geom.AddPoint(rec_value_list[1], rec_value_list[2], rec_value_list[3])

        # create a new feature
        curr_Pt_shape = ogr.Feature(outshape_featdef)
        curr_Pt_shape.SetGeometry(curr_Pt_geom)

        for fld_name, fld_value in zip(field_list, rec_value_list):
            curr_Pt_shape.SetField(fld_name, fld_value)

        # add the feature to the output layer
        point_shapelayer.CreateFeature(curr_Pt_shape)

        # destroy no longer used objects
        curr_Pt_geom.Destroy()
        curr_Pt_shape.Destroy()


def geosurface_read_gas_input(infile_path):
    try:
        with open(infile_path, 'r') as infile:
            input_geosurface = json.load(infile)
    except:
        raise AnaliticSurfaceIOException, "Check input file name"

    src_analytical_params = input_geosurface['analytical surface']
    src_geographical_params = input_geosurface['geographical params']
    try:
        src_deformational_params = input_geosurface['deformational params']
    except:
        src_deformational_params = []

    return src_analytical_params, src_geographical_params, src_deformational_params


def geosurface_export_shapefile_pt3d(shapefile_path, geodata, fields_dict_list, crs=None):
    point_shapefile, point_shapelayer = shapefile_create(shapefile_path,
                                                         ogr.wkbPoint25D,
                                                         fields_dict_list,
                                                         crs)

    field_list = [field_dict["name"] for field_dict in fields_dict_list]

    geosurface_XYZ, _ = geodata
    X, Y, Z = geosurface_XYZ
    assert len(X) == len(Y)
    assert len(X) == len(Z)
    ids = range(len(X))

    rec_values_list2 = zip(ids, X, Y, Z)
    ogr_write_point_result(point_shapelayer, field_list, rec_values_list2)
