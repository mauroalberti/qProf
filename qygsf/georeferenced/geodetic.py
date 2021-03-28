from qygsf.utils.time import *
from qygsf.utils.qgis_utils.points import *


WGS84 = {'semi-major axis': 6378137.0,
         'first eccentricity squared': 6.69437999014e-3}


def n_phi(phi_rad):
    a = WGS84['semi-major axis']
    e_squared = WGS84['first eccentricity squared']
    return a / sqrt(1.0 - e_squared * sin(phi_rad) ** 2)


def geodetic2ecef(lat, lon, height):

    e_squared = WGS84['first eccentricity squared']

    lat_rad, lon_rad = radians(lat), radians(lon)

    nphi = n_phi(lat_rad)

    x = (nphi + height) * cos(lat_rad) * cos(lon_rad)
    y = (nphi + height) * cos(lat_rad) * sin(lon_rad)
    z = (nphi * (1 - e_squared) + height) * sin(lat_rad)

    return x, y, z


class TrackPointGPX(object):

    def __init__(self,
                 lat: numbers.Real,
                 lon: numbers.Real,
                 elev: numbers.Real,
                 time: datetime.datetime):

        self.lat = float(lat)
        self.lon = float(lon)
        self.elev = float(elev)
        self.time = time

    def as_pt3dt(self):

        x, y, _ = geodetic2ecef(self.lat, self.lon, self.elev)
        t = standard_gpstime_to_seconds(self.time)

        return Point4D(x, y, self.elev, t)

    def project(self,
                dest_crs: QgsCoordinateReferenceSystem):

        pt = Point4D(
            x=self.lon,
            y=self.lat,
            z=self.elev,
            t=self.time
        )
        crs = QgsCoordinateReferenceSystem("EPSG:4326")

        projected_pt = project_point(
                pt=pt,
                srcCrs=crs,
                destCrs=dest_crs)

        return projected_pt
