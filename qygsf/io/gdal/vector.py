
from typing import Tuple, Union

import os

try:
    from osgeo import ogr
except ImportError:
    import ogr

try:
    from osgeo import gdal
except ImportError:
    import gdal

try:
    from osgeo import osr
except ImportError:
    import osr

import geopandas as gpd

from ...georeferenced.geoshapes3d import *
from ...utils.types import *


ogr_simpleline_types = [
    ogr.wkbLineString,
    ogr.wkbLineString25D,
    ogr.wkbLineStringM,
    ogr.wkbLineStringZM
]

ogr_multiline_types = [
    ogr.wkbMultiLineString,
    ogr.wkbMultiLineString25D,
    ogr.wkbMultiLineStringM,
    ogr.wkbMultiLineStringZM
]


def try_read_as_geodataframe(
    path: str
) -> Tuple[bool, Union[str, gpd.GeoDataFrame]]:
    """
    Try reading a GIS vectorial dataset as a geopandas GeoDataFrame.

    :param path: the vectorial geodataset full path
    :type path: str
    :return: success status and (error message or results).
    :rtype: Tuple[bool, Union[str, gpd.GeoDataFrame]]
    """

    try:
        v = gpd.read_file(
            filename=path
        )
        return True, v
    except Exception as e:
        return False, str(e)


def try_open_shapefile(
        path: str
) -> Tuple[bool, Union[ogr.Layer, str]]:

    dataSource = ogr.Open(path)

    if dataSource is None:
        return False, "Unable to open shapefile in provided path"

    shapelayer = dataSource.GetLayer()

    return True, shapelayer


def reading_line_shapefile(
        shp_path: str,
        flds: Optional[List[str]] = None,
        read_z: bool = False
    ) -> Tuple[bool, Union[str, List[Tuple[list, tuple]]]]:
    """
    Read results geometries from a line shapefile using ogr.
    TODO: it could read also other formats, but it has to be checked.

    :param shp_path: line shapefile path.
    :param flds: the fields to extract values from.
    :return: success status and (error message or results).
    """

    try:

        if read_z:
            Point = Point3D
            Line = Line3D
        else:
            Point = Point2D
            Line = Line2D

        # check input path

        check_type(shp_path, "Shapefile path", str)
        if shp_path == '':
            return False, "Input shapefile path should not be empty"
        if not os.path.exists(shp_path):
            return False, "Input shapefile path does not exist"

        # open input vector layer

        ds = ogr.Open(shp_path, 0)

        if ds is None:
            return False, "Input shapefile path not read"

        # get internal layer

        layer = ds.GetLayer()

        # get projection

        srs = layer.GetSpatialRef()
        srs.AutoIdentifyEPSG()
        authority = srs.GetAuthorityName(None)
        if authority.upper() == "EPSG":
            epsg_cd = int(srs.GetAuthorityCode(None))
        else:
            epsg_cd = -1

        # initialize list storing results

        results = []

        # loop in layer features

        for feat in layer:

            # get attributes

            if flds:
                feat_attributes = tuple(map(lambda fld_nm: feat.GetField(fld_nm), flds))
            else:
                feat_attributes = ()

            # get geometries

            feat_geometries = []

            curr_geom = feat.GetGeometryRef()

            if curr_geom is None:
                del ds
                return False, "Input shapefile path not read"

            geometry_type = curr_geom.GetGeometryType()
            if geometry_type in ogr_simpleline_types:
                geom_type = "simpleline"
            elif geometry_type in ogr_multiline_types:
                geom_type = "multiline"
            else:
                del ds
                return False, "Geometry type is {}, line expected".format(geom_type)

            if geom_type == "simpleline":

                line = Line()

                for i in range(curr_geom.GetPointCount()):
                    x, y = curr_geom.GetX(i), curr_geom.GetY(i)
                    coords = [x, y]
                    if read_z:
                        z = curr_geom.GetZ(i)
                        coords.append(z)

                    line.add_pt(Point(*coords))

                feat_geometries = line

            else:  # multiline case

                for line_geom in curr_geom:

                    line = Line()

                    for i in range(line_geom.GetPointCount()):
                        x, y = line_geom.GetX(i), line_geom.GetY(i)
                        coords = [x, y]
                        if read_z:
                            z = curr_geom.GetZ(i)
                            coords.append(z)
                        line.add_pt(Point(*coords))

                    feat_geometries.append(line)

            results.append((feat_geometries, feat_attributes))

        del ds

        return True, results

    except Exception as e:

        return False, str(e)


def try_read_line_shapefile(
        shp_path: str,
        flds: Optional[List[str]] = None
    ) -> Tuple[bool, Union[str, List[Tuple[GeoMultiLine3D, tuple]]]]:
    """
    Deprecated. Use 'reading_line_shapefile'.

    Read results geometries from a line shapefile using ogr.
    TODO: it could read also other formats, but it has to be checked.

    :param shp_path: line shapefile path.
    :param flds: the fields to extract values from.
    :return: success status and (error message or results).
    """

    # check input path

    check_type(shp_path, "Shapefile path", str)
    if shp_path == '':
        return False, "Input shapefile path should not be empty"
    if not os.path.exists(shp_path):
        return False, "Input shapefile path does not exist"

    # open input vector layer

    ds = ogr.Open(shp_path, 0)

    if ds is None:
        return False, "Input shapefile path not read"

    # get internal layer

    lyr = ds.GetLayer()

    # get projection

    srs = lyr.GetSpatialRef()
    srs.AutoIdentifyEPSG()
    authority = srs.GetAuthorityName(None)
    if authority.upper() == "EPSG":
        epsg_cd = int(srs.GetAuthorityCode(None))
    else:
        epsg_cd = -1

    # initialize list storing results

    results = []

    # loop in layer features

    for feat in lyr:

        # get attributes

        if flds:
            feat_attributes = tuple(map(lambda fld_nm: feat.GetField(fld_nm), flds))
        else:
            feat_attributes = ()

        # get geometries

        # feat_geometries = []

        curr_geom = feat.GetGeometryRef()

        if curr_geom is None:
            del ds
            return False, "Input shapefile path not read"

        geometry_type = curr_geom.GetGeometryType()
        if geometry_type in ogr_simpleline_types:
            geom_type = "simpleline"
        elif geometry_type in ogr_multiline_types:
            geom_type = "multiline"
        else:
            del ds
            return False, "Geometry type is {}, line expected".format(geometry_type)

        multiline = GeoMultiLine3D(epsg_cd=epsg_cd)

        if geom_type == "simpleline":

            line = Line3D()

            for i in range(curr_geom.GetPointCount()):
                x, y, z = curr_geom.GetX(i), curr_geom.GetY(i), curr_geom.GetZ(i)

                line.add_pt(Point3D(x, y, z))

            multiline.add_line(line)

        else:  # multiline case

            for line_geom in curr_geom:

                line = Line3D()

                for i in range(line_geom.GetPointCount()):
                    x, y, z = line_geom.GetX(i), line_geom.GetY(i), line_geom.GetZ(i)

                    line.add_pt(Point3D(x, y, z))

                multiline.add_line(line)

        results.append((multiline, feat_attributes))

    del ds

    return True, results


def read_linestring_geometries(
        line_shp_path: str
) -> Optional[GeoMultiLine3D]:
    """
    Read linestring geometries from a shapefile using ogr.
    The geometry type of the input shapefile must be LineString (MultiLineString is not currently managed).

    It returns a MultiLine instance.

    :param line_shp_path:  parameter to check.
    :return: the result of data reading
    """

    # check input path

    if line_shp_path is None or line_shp_path == '':
        return None

    # open input vector layer

    shape_driver = ogr.GetDriverByName("ESRI Shapefile")

    datasource = shape_driver.Open(str(line_shp_path), 0)

    # layer not read

    if datasource is None:
        return None

    # get internal layer

    layer = datasource.GetLayer()

    # get projection

    try:
        srs = layer.GetSpatialRef()
        srs.AutoIdentifyEPSG()
        authority = srs.GetAuthorityName(None)
        if authority == "EPSG":
            epsg_cd = int(srs.GetAuthorityCode(None))
        else:
            epsg_cd = -1
    except:
        epsg_cd = -1

    # initialize list storing vertex coordinates of lines

    lines = []

    # start reading layer features

    feature = layer.GetNextFeature()

    # loop in layer features

    while feature:

        geometry = feature.GetGeometryRef()

        if geometry is None:
            datasource.Destroy()
            return None

        geometry_type = geometry.GetGeometryType()
        if geometry_type not in ogr_simpleline_types:
            datasource.Destroy()
            return None

        line = Line3D()

        for i in range(geometry.GetPointCount()):

            x, y, z = geometry.GetX(i), geometry.GetY(i), geometry.GetZ(i)

            line.add_pt(Point3D(x, y, z))

        feature.Destroy()

        lines.append(line)

        feature = layer.GetNextFeature()

    datasource.Destroy()

    multiline = GeoMultiLine3D(
        lines=lines,
        epsg_cd=epsg_cd
    )

    return multiline


def parse_ogr_type(ogr_type_str: str) -> 'ogr.OGRFieldType':
    """
    Parse the provided textual field type to return an actual OGRFieldType.

    :param ogr_type_str: the string referring to the ogr field type.
    :type ogr_type_str: str.
    :return: the actural ogr type.
    :rtype: OGRFieldType.
    :raise: Exception.
    """

    if ogr_type_str.endswith("OFTInteger"):
        return ogr.OFTInteger
    elif ogr_type_str.endswith("OFTIntegerList"):
        return ogr.OFTIntegerList
    elif ogr_type_str.endswith("OFTReal"):
        return ogr.OFTReal
    elif ogr_type_str.endswith("OFTRealList"):
        return ogr.OFTRealList
    elif ogr_type_str.endswith("OFTString"):
        return ogr.OFTString
    elif ogr_type_str.endswith("OFTStringList"):
        return ogr.OFTStringList
    elif ogr_type_str.endswith("OFTBinary"):
        return ogr.OFTBinary
    elif ogr_type_str.endswith("OFTDate"):
        return ogr.OFTDate
    elif ogr_type_str.endswith("OFTTime"):
        return ogr.OFTTime
    elif ogr_type_str.endswith("OFTDateTime"):
        return ogr.OFTDateTime
    elif ogr_type_str.endswith("OFTInteger64"):
        return ogr.OFTInteger64
    elif ogr_type_str.endswith("OFTInteger64List"):
        return ogr.OFTInteger64List
    else:
        raise Exception("Debug: not recognized ogr type")


class GDALParameters(object):
    """
    Manage GDAL parameters from grids.

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
            raise Exception('Pixel sizes in x and y directions are different in raster')

            # check for the absence of axis rotations
        if abs(self._rotation_GT_2) > tolerance or abs(self._rotation_GT_4) > tolerance:
            raise Exception('There should be no axis rotation in raster')

        return

    def llcorner(self):
        """
        Creates a point at the lower-left corner of the raster.

        @return:  new Point instance.
        """
        return Point2D(self.topLeftX, self.topLeftY - abs(self.pixSizeNS) * self.rows)

    def trcorner(self):
        """
        Create a point at the top-right corner of the raster.

        @return:  new Point instance.
        """
        return Point2D(self.topLeftX + abs(self.pixSizeEW) * self.cols, self.topLeftY)

    def geo_equiv(self, other, tolerance=1.0e-6):
        """
        Checks if two grids are geographically equivalent.

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
                        self.rows != other.rows or self.cols != other.cols:
            return False
        else:
            return True


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

            line_points.append(Point3D(x, y, z))

        lines_points.append(line_points)

        curr_line = lnLayer.GetNextFeature()

    line_shape.Destroy()

    return dict(success=True, extent=lines_extent, vertices=lines_points)


def shapefile_create_def_field(field_def):
    """

    :param field_def:
    :return:
    """

    name = field_def['name']
    ogr_type = parse_ogr_type(field_def['ogr_type'])

    fieldDef = ogr.FieldDefn(name, ogr_type)
    if ogr_type == ogr.OFTString:
        fieldDef.SetWidth(int(field_def['width']))

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
        raise Exception('Unable to save shapefile in provided path')

    if crs is not None:
        spatial_reference = osr.SpatialReference()
        spatial_reference.ImportFromProj4(crs)
        outShapelayer = outShapefile.CreateLayer("layer", spatial_reference, geom_type)
    else:
        outShapelayer = outShapefile.CreateLayer("layer", None, geom_type)

    if not outShapelayer:
        return None, None

    for field_def_params in fields_dict_list:
        field_def = shapefile_create_def_field(field_def_params)
        outShapelayer.CreateField(field_def)

    return outShapefile, outShapelayer


def try_write_pt_shapefile(point_layer, geoms: List[Tuple[numbers.Real, numbers.Real, numbers.Real]], field_names: List[str], attrs: List[Tuple]) -> Tuple[bool, str]:
    """
    Add point records in an existing shapefile, filling attribute values.

    :param point_layer: the existing shapefile layer in which to write.
    :param geoms: the geometric coordinates of the points.
    :type geoms: List of x, y, and z coordinates.
    :param field_names: the field names of the attribute table.
    :type field_names: list of strings.
    :param attrs: the values for each record.
    :type attrs: list of tuple.
    :return: success status and related messages.
    :rtype: tuple of a boolean and a string.
    """

    len_geoms = len(geoms)
    len_attrs = len(attrs)

    if len_geoms != len_attrs:
        return False, "Function error: geometries are {} while attributes are {}".format(len_geoms, len_attrs)

    if len_geoms == 0:
        return True, "No values to be added in shapefile"

    try:

        outshape_featdef = point_layer.GetLayerDefn()

        for ndx_rec in range(len_geoms):

            # pre-processing for new feature in output layer

            curr_Pt_geom = ogr.Geometry(ogr.wkbPoint25D)
            curr_Pt_geom.AddPoint(*geoms[ndx_rec])

            # create a new feature

            curr_pt_shape = ogr.Feature(outshape_featdef)
            curr_pt_shape.SetGeometry(curr_Pt_geom)

            rec_attrs = attrs[ndx_rec]

            for ndx_fld, fld_nm in enumerate(field_names):

                curr_pt_shape.SetField(fld_nm, rec_attrs[ndx_fld])

            # add the feature to the output layer
            point_layer.CreateFeature(curr_pt_shape)

            # destroy no longer used objects
            curr_Pt_geom.Destroy()
            curr_pt_shape.Destroy()

        del outshape_featdef

        return True, ""

    except Exception as e:

        return False, "Exception: {}".format(e)


def try_write_point_shapefile(path: str, field_names: List[str], values: List[Tuple], ndx_x_val: int) -> Tuple[bool, str]:
    """
    Note: candidate for future deprecation.

    Add point records in an existing shapefile, filling attribute values.
    The point coordinates, i.e. x, y, z start at ndx_x_val index (index is zero-based) and are
    assumed to be sequential in order (i.e., 0, 1, 2 or 3, 4, 5).

    :param path: the path of the existing shapefile in which to write.
    :type path: string.
    :param field_names: the field names of the attribute table.
    :type field_names: list of strings.
    :param values: the values for each record.
    :type values: list of tuple.
    :param ndx_x_val: the index of the x coordinate. Y and z should follow.
    :type ndx_x_val: int.
    :return: success status and related messages.
    :rtype: tuple of a boolean and a string.
    """

    success = False
    msg = ""

    try:

        dataSource = ogr.Open(path, 1)

        if dataSource is None:
            return False, "Unable to open shapefile in provided path"

        point_layer = dataSource.GetLayer()

        outshape_featdef = point_layer.GetLayerDefn()

        for pt_vals in values:

            # pre-processing for new feature in output layer
            curr_Pt_geom = ogr.Geometry(ogr.wkbPoint)
            curr_Pt_geom.AddPoint(pt_vals[ndx_x_val], pt_vals[ndx_x_val+1], pt_vals[ndx_x_val+2])

            # create a new feature
            curr_pt_shape = ogr.Feature(outshape_featdef)
            curr_pt_shape.SetGeometry(curr_Pt_geom)

            for ndx, fld_nm in enumerate(field_names):

                curr_pt_shape.SetField(fld_nm, pt_vals[ndx])

            # add the feature to the output layer
            point_layer.CreateFeature(curr_pt_shape)

            # destroy no longer used objects
            curr_Pt_geom.Destroy()
            curr_pt_shape.Destroy()

        del outshape_featdef
        del point_layer
        del dataSource

        success = True

    except Exception as e:

        msg = e

    finally:

        return success, msg


def try_write_line_shapefile(path: str, field_names: List[str], values: Dict) -> Tuple[bool, str]:
    """
    Add point records in an existing shapefile, filling attribute values.


    :param path: the path of the existing shapefile in which to write.
    :type path: string.
    :param field_names: the field names of the attribute table.
    :type field_names: list of strings.
    :param values: the values for each record.
    :type values: dict with values made up by two dictionaries.
    :return: success status and related messages.
    :rtype: tuple of a boolean and a string.
    """

    success = False
    msg = ""

    try:

        dataSource = ogr.Open(path, 1)

        if dataSource is None:
            return False, "Unable to open shapefile in provided path"

        line_layer = dataSource.GetLayer()

        outshape_featdef = line_layer.GetLayerDefn()

        for curr_id in sorted(values.keys()):

            # pre-processing for new feature in output layer
            line_geom = ogr.Geometry(ogr.wkbLineString)

            for id_xyz in values[curr_id]["pts"]:
                x, y, z = id_xyz
                line_geom.AddPoint(x, y, z)

            # create a new feature
            line_shape = ogr.Feature(outshape_featdef)
            line_shape.SetGeometry(line_geom)

            for ndx, fld_nm in enumerate(field_names):

                line_shape.SetField(fld_nm, values[curr_id]["vals"][ndx])

            # add the feature to the output layer
            line_layer.CreateFeature(line_shape)

            # destroy no longer used objects
            line_geom.Destroy()
            line_shape.Destroy()

        del outshape_featdef
        del line_layer
        del dataSource

        success = True

    except Exception as e:

        msg = str(e)

    finally:

        return success, msg


def ogr_get_solution_shapefile(path, fields_dict_list):
    """

    :param path:
    :param fields_dict_list:
    :return:
    """

    driver = ogr.GetDriverByName("ESRI Shapefile")

    dataSource = driver.Open(str(path), 0)

    if dataSource is None:
        raise Exception('Unable to open shapefile in provided path')

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


def ogr_write_point_result(
        point_shapelayer,
        field_list,
        rec_values_list2,
        geom_type=ogr.wkbPoint25D
):
    outshape_featdef = point_shapelayer.GetLayerDefn()

    for rec_value_list in rec_values_list2:

        # pre-processing for new feature in output layer
        curr_Pt_geom = ogr.Geometry(geom_type)
        if geom_type == ogr.wkbPoint25D:
            curr_Pt_geom.AddPoint(rec_value_list[1], rec_value_list[2], rec_value_list[3])
        else:
            curr_Pt_geom.AddPoint(rec_value_list[1], rec_value_list[2])

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

