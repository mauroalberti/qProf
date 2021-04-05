
import numbers

from ..geometries.shapes.space3d import *
from ..geometries.shapes.space4d import *

# Earth WGS84 parameters

WGS84 = {'semi-major axis': 6378137.0,
         'first eccentricity squared': 6.69437999014e-3}


def n_phi(phi_rad: numbers.Real) -> numbers.Real:
    """
    It return the N(phi) parameter.
    See: https://en.wikipedia.org/wiki/Geographic_coordinate_conversion#From_geodetic_to_ECEF_coordinates
    N(phi) =  a / sqrt( 1 - e^2 * sin(phi)^2
    where phi is the latitude (in radians) and e^2 is the first eccentricity squared.

    :param phi_rad: the latitude expressed in radians.
    :type phi_rad: numbers.Real.
    :return: the N(phi) value.
    :rtype: numbers.Real.
    """

    a = WGS84['semi-major axis']
    e_squared = WGS84['first eccentricity squared']
    return a / sqrt(1.0 - e_squared * sin(phi_rad) ** 2)


def geodetic2ecef(lat: numbers.Real, lon: numbers.Real, height: numbers.Real) -> Tuple[numbers.Real, numbers.Real, numbers.Real]:
    """
    Converts from geodetic (lat-long-height) to Cartesian ECEF reference system.
    See: https://en.wikipedia.org/wiki/Geographic_coordinate_conversion#From_geodetic_to_ECEF_coordinates

    :param lat: latitude.
    :type lat: numbers.Real.
    :param lon: longitude.
    ;:type lon: numbers.Real.
    :param height: height.
    :type height: numbers.Real.
    :return: x, y and z coordinates.
    :rtype: tuple of three numbers.Real values.
    """

    e_squared = WGS84['first eccentricity squared']

    lat_rad, lon_rad = radians(lat), radians(lon)

    nphi = n_phi(lat_rad)

    x = (nphi + height) * cos(lat_rad) * cos(lon_rad)
    y = (nphi + height) * cos(lat_rad) * sin(lon_rad)
    z = (nphi * (1 - e_squared) + height) * sin(lat_rad)

    return x, y, z


latitude_one_degree_45degr_meters = 111131.745  #  source: http://www.csgnetwork.com/degreelenllavcalc.html, consulted on 2018-12-23


def projectionType(id_code):
    """
    NOTE: currently it is a stub code.
    TODO: make more general.

    Determines if the provided projection type is polar, planar or unknown.

    :param id_code: string.
    :return: True if polar, False if planar
    :rtype: bool.
    """

    if id_code == "EPSG:4326":
        return "polar"

    return "unknown"


def latLengthOneMinutePrime() -> numbers.Real:
    """
    Approximate length (in meters) of one minute prime at latitude 45°.
    :return: length in meters.
    :rtype: numbers.Real
    """

    return latitude_one_degree_45degr_meters / 60.0


def latLengthOneMinuteSecond() -> numbers.Real:
    """
    Approximate length (in meters) of one minute second at latitude 45°.
    :return: length in meters.
    :rtype: numbers.Real
    """

    return latitude_one_degree_45degr_meters / 3600.0


def pt_4326_ecef(pt: Point3D) -> Optional[Point3D]:
    """
    Project a point from EPSG:4326 to ECEF

    :param pt: Point instance.
    :type pt: Point.
    :return: the projected Point instance.
    :rtype: Point.
    """

    lon, lat, height = pt.x, pt.y, pt.z
    x, y, z = geodetic2ecef(lat, lon, height)

    return Point3D(
        x=x,
        y=y,
        z=z
    )


def line_4326_ecef(line: Line3D) -> Optional[Line3D]:
    """
    Converts from WGS84 to ECEF reference system, provided its CRS is EPSG:4326.

    :return: a line with ECEF coordinates (EPSG:4978).
    :rtype: optional Line.
    """

    pts = [pt_4326_ecef(pt) for pt in line.pts()]

    return Line3D(
        pts=pts
    )


if __name__ == "__main__":

    print("Approximate length of one minute prime at 45° latitude: {} meters".format(latLengthOneMinutePrime()))
    print("Approximate length of one minute second at 45° latitude: {} meters".format(latLengthOneMinuteSecond()))

