
from geopandas import GeoDataFrame

from pygsf.orientations.orientations import *
from pygsf.geometries.shapes.space3d import *
from pygsf.io.geodataframes.geodataframes import *

from .base import GeorefAttitude


def try_extract_flat_georeferenced_attitudes(
        geodataframe: GeoDataFrame,
        azim_fldnm: str,
        dip_ang_fldnm: str,
        id_fldnm: Optional[str] = None,
        is_rhrstrike: bool = False
) -> Tuple[bool, Union[str, List[GeorefAttitude]]]:
    """
    Try extracting the georeferenced _attitudes from a geopandas GeoDataFrame instance representing point records.

    :param geodataframe: the source geodataframe
    :param azim_fldnm: the name of the azimuth field in the geodataframe
    :param dip_ang_fldnm: the name of the dip angle field in the geodataframe
    :param id_fldnm: the name of the id field in the geodataframe
    :param is_rhrstrike: whether the dip azimuth is strike RHR
    :return: the success status and an error message or a collection of georeferenced attitudes, one for each source record
    """

    try:

        #epsg = get_epsg(geodataframe)

        attitudes = []

        for ndx, row in geodataframe.iterrows():

            pt = row['geometry']
            x, y = pt.x, pt.y

            if id_fldnm:
                azimuth, dip_ang, rec_id = row[azim_fldnm], row[dip_ang_fldnm], row[id_fldnm]
            else:
                azimuth, dip_ang, rec_id = row[azim_fldnm], row[dip_ang_fldnm], ndx + 1

            if is_rhrstrike:
                azimuth = (azimuth + 90.0) % 360.0

            attitudes.append(
                GeorefAttitude(
                    rec_id,
                    Point2D(x, y),
                    Plane(azimuth, dip_ang)
                )
            )

        return True, attitudes

    except Exception as e:

        return False, str(e)


