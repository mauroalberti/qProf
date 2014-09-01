
from math import sqrt, radians, sin, cos

from .spatial import Point_4D

WGS84 = { 'semi-major axis': 6378137.0, 
          'first eccentricity squared': 6.69437999014e-3}
   
        
def N_Phi( phi_rad ):

    a = WGS84['semi-major axis']
    e_squared = WGS84['first eccentricity squared']
    return a / sqrt( 1.0 - e_squared*sin(phi_rad)**2)


def geodetic2ecef ( lat, lon, height ):
    
    e_squared = WGS84['first eccentricity squared']
    
    lat_rad, lon_rad = radians( lat ), radians( lon )
    
    n_phi = N_Phi( lat_rad )
    
    X = ( n_phi + height ) * cos( lat_rad ) * cos( lon_rad )
    Y = ( n_phi + height ) * cos( lat_rad ) * sin( lon_rad ) 
    Z = ( n_phi*(1-e_squared) + height ) * sin( lat_rad )
    
    return X, Y, Z 
    

class TrackPointGPX( object ):
    
    def __init__(self, lat, lon, elev, time):
        
        self.lat = float( lat )
        self.lon = float( lon )
        self.time = time      
        self.elev = float( elev )     
        
        
    def toPoint4D ( self ):
        
        x, y, z = geodetic2ecef ( self.lat, self.lon, self.elev )
        t = self.time
        
        return Point_4D( x, y, z, t )
        
        
    