
from typing import Set, Any, List
import numbers

import geopandas as gpd

from ...geometries.shapes.space2d import Point2D


def geodataframe_geom_types(
    geodataframe: gpd.GeoDataFrame
) -> Set[str]:
    # Side effects: none
    """
    Return a set storing the geometric types in a GeoDataFrame instance.

    :param geodataframe: the input geodataframe
    :type geodataframe: gpd.GeoDataFrame
    :return: the set of geometric types values
    :rtype: Set[str]
    """

    return set(geodataframe.geom_type)


def containsPoints(
    geodataframe: gpd.GeoDataFrame
) -> bool:
    # Side effects: none
    """
    Check if a GeoDataFrame instance contains points.

    :param geodataframe: the input geodataframe
    :type geodataframe: gpd.GeoDataFrame
    :return: if a GeoDataFrame instance contains points
    :rtype: bool
    """

    return 'Point' in geodataframe_geom_types(geodataframe)


def containsLines(
    geodataframe: gpd.GeoDataFrame
) -> bool:
    # Side effects: none
    """
    Check if a GeoDataFrame instance contains lines.

    :param geodataframe: the input geodataframe
    :type geodataframe: gpd.GeoDataFrame
    :return: if a GeoDataFrame instance contains lines
    :rtype: bool
    """

    return 'LineString' in geodataframe_geom_types(geodataframe)


def containsPolygons(
    geodataframe: gpd.GeoDataFrame
) -> bool:
    # Side effects: none
    """
    Check if a GeoDataFrame instance contains polygons.

    :param geodataframe: the input geodataframe
    :type geodataframe: gpd.GeoDataFrame
    :return: if a GeoDataFrame instance contains polygons
    :rtype: bool
    """

    return 'Polygon' in geodataframe_geom_types(geodataframe)


def extract_geometries(
    geodataframe: gpd.GeoDataFrame
) -> gpd.geoseries.GeoSeries:
    # Side effects: none
    """
    Extract geometries from a GeoDataFrame instance.

    :param geodataframe: the input geodataframe
    :type geodataframe: gpd.GeoDataFrame
    :return: the geometries stored in the GeoDataFrame instance
    :rtype: geopandas.geoseries.GeoSeries
    """

    return geodataframe.geometry


def extract_geometry(
    geodataframe: gpd.GeoDataFrame,
    ndx: numbers.Integral
) -> Any:
    # Side effects: none
    """
    Extract a geometry from a GeoDataFrame instance,
    given the geometry index.

    :param geodataframe: the input geodataframe
    :type geodataframe: gpd.GeoDataFrame
    :param ndx: the geometry index
    :type ndx: numbers.Integral
    :return: the geometry stored in the GeoDataFrame instance
    :rtype: Shapely geometry
    """

    return extract_geometries(geodataframe)[ndx]


def get_epsg(
    geodataframe: gpd.GeoDataFrame
) -> numbers.Integral:
    # Side effects: None
    """
    Extract the EPSG code of the data

    :param geodataframe: the input geodataframe
    :type geodataframe: gpd.GeoDataFrame
    :return: the EPSG code or -1
    :rtype: numbers.Integral
    """

    crs_dict = geodataframe.crs
    #print("Source georeferenced: {}".format(crs_dict))
    epsg = -1
    try:
        val = crs_dict["init"]
        if val.lower().startswith("epsg"):
            epsg = int(val.split(":", 1)[1])
    except Exception as e:
        pass

    #print("Extracted EPSG code: {}".format(epsg))
    return epsg


def extract_line_points(
    geodataframe: gpd.GeoDataFrame,
    ndx: numbers.Integral,
    epsg_code: numbers.Integral
) -> List[Point2D]:
    """
    Extract a geometry from a GeoDataFrame instance,
    given the geometry index.

    :param geodataframe: the input geodataframe
    :type geodataframe: gpd.GeoDataFrame
    :param ndx: the geometry index
    :type ndx: numbers.Integral
    :param epsg_code: the EPSG code
    :type epsg_code: numbers.Integral
    :return: the geometry stored in the GeoDataFrame instance
    :rtype: Shapely geometry
    """

    geometry = extract_geometry(
        geodataframe=geodataframe,
        ndx=ndx
    )

    xs, ys = geometry.xy

    pts = []

    for x, y in zip(xs, ys):
        pts.append(
            Point2D(
                x=x,
                y=y
            )
        )

    return pts


