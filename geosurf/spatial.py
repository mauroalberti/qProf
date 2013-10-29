# -*- coding: utf-8 -*-


from __future__  import division

from math import sqrt
from math import floor, ceil
from math import sin, cos, tan, radians, asin, atan, atan2, degrees

import numpy as np

import copy

from .utils import array_from_function
from .algebr_utils import point_solution


MINIMUM_SEPARATION_THRESHOLD = 1e-10
MINIMUM_VECTOR_MAGNITUDE = 1e-10 
       
 

class Point(object):
    """
    3D points.
    
    """
        
    # class constructor
    def __init__( self, x = np.nan, y = np.nan, z = np.nan ):
        """ 
        Point class constructor.
        
        @param  x:  x coordinate of the point.
        @type  x:  number or string convertible to float.
        @param  y:  y coordinate of the point.
        @type  y:  number or string convertible to float.
        @param  z:  z coordinate of the point.
        @type  z:  number or string convertible to float.
                   
        @return:  Point instance located at (x,y,z).        
        """        
        self._x = float( x )
        self._y = float( y )
        self._z = float( z )

    
    def g_x( self ):
        """
        Get point x value.
        
        @return:  x coordinate of Point instance - float.
        """        
        return self._x
    

    def s_x( self, xval ):
        """
        Set point x value.
        
        @param  xval:  x coordinate to be set to the point.
        @type  xval:  number or string convertible to float.
        
        @return:  self.        
        """                
        self._x = float( xval )
        
        
    # property for x
    x = property( g_x, s_x )


    def g_y( self ):
        """
        Get point y value.
        
        @return:  y coordinate of Point instance - float.
        """
        return self._y
    

    def s_y( self, yval ):
        """
        Set point y value.
        
        @param  yval:  y coordinate to be set to the point.
        @type  yval:  number or string convertible to float.
        
        @return:  self.        
        """
        self._y = yval


    # property for y
    y = property( g_y, s_y )


    # z value
    def g_z( self ):
        """
        Get point z value.
        
        @return:  z coordinate of Point instance - float.
        """
        return self._z


    def s_z( self, zval ):
        """
        Set point z value.
        
        @param  zval:  z coordinate to be set to the point.
        @type  zval:  number or string convertible to float.
        
        @return:  self. 
        """
        self._z = zval


    # property for z
    z = property( g_z, s_z )
    
  
    def distance( self, another ):
        """
        Calculate Euclidean distance between two points.

        @param  another:  the Point instance for which the distance should be calculated
        @type  another:  Point.
        
        @return:  distance between the two points - float.        
        """
        
        assert self.z != np.nan
        assert another.z != np.nan
        return sqrt( ( self.x - another.x ) ** 2 + ( self.y - another.y ) ** 2 + ( self.z - another.z ) ** 2 )
    
    
    def hor_distance( self, another ):
        """
        Calculate Euclidean horizontal distance between two points.

        @param  another:  the Point instance for which the distance should be calculated
        @type  another:  Point.
        
        @return:  horizontal distance between the two points - float.        
        """
        return sqrt( ( self.x - another.x ) ** 2 + ( self.y - another.y ) ** 2 )
    
        
    def movedby( self, sx = 0.0 , sy = 0.0 , sz = 0.0 ):
        """
        Create a new point shifted by given amount from the self instance.

        @param  sx:  the shift to be applied along the x axis.
        @type  sx:  float.
        @param  sy:  the shift to be applied along the y axis.
        @type  sy:  float.
        @param  sz:  the shift to be applied along the z axis.
        @type  sz:  float.
                
        @return:  a new Point instance shifted by the given amounts with respect to the original one.        
        """
        return Point( self.x + sx , self.y + sy, self.z + sz )
    

    def displaced_by_vector(self, displacement_vector ):
        
        return Point( self.x + displacement_vector.x , self.y + displacement_vector.y, self.z + displacement_vector.z )
        

    def to_vector(self, end_pt):
        
        dx = end_pt.x - self.x
        dy = end_pt.y - self.y
        if self.z == np.nan or end_pt == np.nan: dz = 0.0            
        else: dz = end_pt.z - self.z 
        
        return Vector( dx, dy, dz )        
        
        
    def distance_uvector(self, end_pt):
        
        d_x = end_pt.x - self.x
        d_y = end_pt.y - self.y
        d_z = end_pt.z - self.z
         
        length_2d_squared = d_x*d_x +  d_y*d_y  
        length_2d = sqrt( length_2d_squared )         
        lenght_3d = sqrt( length_2d_squared +  d_z*d_z )
        
        if length_2d < MINIMUM_SEPARATION_THRESHOLD or lenght_3d < MINIMUM_SEPARATION_THRESHOLD:
            return 0.0, 0.0, Vector( 0.0, 0.0, 0.0 ), Vector( 0.0, 0.0, 0.0 )
        else:
            return length_2d, lenght_3d, Vector( d_x/length_2d, d_y/length_2d, 0.0 ), Vector( d_x/lenght_3d, d_y/lenght_3d, d_z/lenght_3d )


    def as_vector( self ):
        
        return Vector(self.x,
                      self.y,
                      self.z )
        
        
class Point4D( Point ):
    
    def __init__(self, x = np.nan, y = np.nan, z = 0.0, time = np.nan ):
        
        super(Point4D, self).__init__(x, y, z)
        self._t = time
            
    
    def g_t( self ):
        """
        Get point time value.
        
        @return:  time coordinate of Point instance - float.
        """        
        return self._t
    

    def s_t( self, time_val ):
        """
        Set point time value.
        
        @param  time_val:  time coordinate to be set to the point.
        @type  time_val:  number or string convertible to float.
        
        @return:  self.        
        """                
        self._t = float( time_val )
        
        
    # property for x
    t = property( g_t, s_t )
    
    
    def delta_time(self, anotherPt4D): 
        
        return anotherPt4D._t - self._t
    
    
    def speed(self, anotherPt4D):
        
        try:
            return self.distance( anotherPt4D ) / self.delta_time( anotherPt4D )
        except:
            return np.nan

    
class Vector( Point ):
    
    def __init__(self, x = np.nan, y = np.nan, z = 0.0):
        
        super(Vector, self).__init__(x, y, z)
          
         
    def lenght_hor(self):
        
        return sqrt( self.x * self.x + self.y * self.y )
    
       
    def lenght_3d(self):
        
        return sqrt( self.x * self.x + self.y * self.y + self.z * self.z )
        
        
    def scale(self, scale_factor):
                
        return Vector(self.x*scale_factor, 
                      self.y*scale_factor, 
                      self.z*scale_factor)


    def to_down_vector(self):
        
        if self.z > 0.0:
            return self.scale(-1.0)
        else:
            return self
        
  
    def to_unit_vector(self):
        
        scale_factor = self.lenght_3d()
        
        assert scale_factor is not None
        assert scale_factor > 1.0e-4
        
        return self.scale( 1.0 / scale_factor )
    
        
    def add(self, another):
        
        return Vector(self.x + another.x, 
                      self.y + another.y, 
                      self.z + another.z)
        

    def slope_radians(self):
        
        try:
            return atan( self.z / self.lenght_hor() )
        except:
            return np.nan
              

    def to_axis(self):
        
        if self.lenght_3d() < MINIMUM_VECTOR_MAGNITUDE:
            return None
        
        unit_vect = self.to_unit_vector()
        
        plunge = - degrees( asin( unit_vect.z ) ) # upward negative, downward positive
        
        trend = 90.0 - degrees( atan2( unit_vect.y, unit_vect.x ) )
        if trend < 0.0: trend += 360.0
        elif trend > 360.0: trend -= 360.0
        
        assert 0.0 <= trend < 360.0
        assert -90.0 <= plunge <= 90.0
         
        return Axis( trend, plunge )
    
    
    def scalar_product( self, another ):
        
        return self.x * another.x + self.y * another.y + self.z * another.z


    def vectors_cos_angle( self, another ):
 
        try:
            return self.scalar_product( another ) / ( self.lenght_3d() * another.lenght_3d() )
        except ZeroDivisionError:
            return None
        

    def vector_product( self, another ):
        
        x = (self.y*another.z) - (self.z*another.y)
        y = (self.z*another.x) - (self.x*another.z)
        z = (self.x*another.y) - (self.y*another.x)

        return Vector( x, y, z )


    def by_matrix( self, matrix3x3 ):
        
        vx = matrix3x3[0,0]*self.x + matrix3x3[0,1]*self.y + matrix3x3[0,2]*self.z 
        vy = matrix3x3[1,0]*self.x + matrix3x3[1,1]*self.y + matrix3x3[1,2]*self.z 
        vz = matrix3x3[2,0]*self.x + matrix3x3[2,1]*self.y + matrix3x3[2,2]*self.z 
        
        return Vector( vx, vy, vz )
        

    def as_point( self ):
        
        return Point( self.x, self.y, self.z )
        

class Axis( object ):
    """
    Structural axis,
    defined by trend and plunge (both in degrees)
    Trend range: [0.0, 360.0[  clockwise, from 0 (North) 
    Plunge: [-90.0, 90.0], negative value: upward axis, positive values: downward axis
    """
    
    def __init__(self, srcTrend, srcPlunge ):
        
        assert 0.0 <= float( srcTrend ) < 360.0
        assert -90.0 <= float( srcPlunge ) <= 90.0
        
        self._trend = float( srcTrend )
        self._plunge = float( srcPlunge )        
        

    def to_versor( self ):
        
        north_coord = cos( radians( self._plunge ) ) * cos( radians( self._trend ) )
        east_coord = cos( radians( self._plunge ) ) * sin( radians( self._trend ) )
        down_coord = sin( radians( self._plunge ) )
        
        return Vector( east_coord, north_coord, -down_coord )
        
             
    def to_down_axis( self ): 
        
        trend, plunge = self._trend, self._plunge
        if plunge < 0.0:
            trend += 180.0
            if trend > 360.0: trend -= 360.0
            plunge = - plunge
        
        return Axis( trend, plunge )


    def to_normal_geolplane( self ):
        
        down_axis = self.to_down_axis()
        
        dipdir = down_axis._trend + 180.0
        if dipdir >= 360.0: dipdir -= 360.0
        dipangle = 90.0 - down_axis._plunge
        
        return GeolPlane( dipdir, dipangle ) 



class CartesianPlane(object):
    """
    Cartesian plane, expressed by equation:
    ax + by + cz + d = 0

    """

    def __init__(self, a=None, b=None, c=None, d = 0.0 ):

        self._a = a
        self._b = b        
        self._c = c        
        self._d = d
                

    def from_points(self, pt1, pt2, pt3):
        
        matr_a = np.array([[pt1.y, pt1.z, 1], 
                           [pt2.y, pt2.z, 1],
                           [pt3.y, pt3.z, 1] ])

        matr_b = - np.array([[pt1.x, pt1.z, 1], 
                           [pt2.x, pt2.z, 1],
                           [pt3.x, pt3.z, 1] ])  
        
        matr_c = np.array([[pt1.x, pt1.y, 1], 
                           [pt2.x, pt2.y, 1],
                           [pt3.x, pt3.y, 1] ])
        
        matr_d = - np.array([[pt1.x, pt1.y, pt1.z], 
                           [pt2.x, pt2.y, pt2.z],
                           [pt3.x, pt3.y, pt3.z] ])
         
        self._a = np.linalg.det(matr_a)
        self._b = np.linalg.det(matr_b)        
        self._c = np.linalg.det(matr_c)        
        self._d = np.linalg.det(matr_d)
        
        return self
        
        
    def to_normal_versor(self):
        """
        return the normal versor to the cartesian plane
        """
        
        assert self._a is not None
        assert self._b is not None        
        assert self._c is not None
                
        return Vector(self._a,self._b,self._c).to_unit_vector()
        
        
    def to_geolplane_and_point(self):
        """
        converts a cartesian plane into a geological plane
        and a point lying in the plane (non-unique solution)
        """

        geol_plane = self.to_normal_versor().to_axis().to_normal_geolplane()
               
        point = Point( point_solution(np.array([[self._a,self._b,self._c]]),
                                      np.array([-self._d])))
        
        return geol_plane, point
    
        
    def intersect(self, another):
        """
        return intersection versor and point on intersection (non-unique solution)
        for two intersecting planes
        """

        # find the intersection vector 
        
        self_plane_normal = self.to_normal_versor()
        another_plane_normal = another.to_normal_versor()
    
        intersection_versor = self_plane_normal.vector_product(another_plane_normal).to_unit_vector()


        # find a point lying on the intersection line (this is a non-unique solution)    
        a = np.array([[self._a,self._b,self._c],[another._a,another._b,another._c]])
        b = np.array([-self._d,-another._d]) 
        x, y, z = point_solution( a, b ) 
              
        inters_point = Point( x, y, z )
        
        return intersection_versor, inters_point
        
        
class GeolPlane(object):
    """
    Structural plane, following geological conventions:
    dip direction and dip angle.
    
    """
    
    def __init__( self, srcDipDir, srcDipAngle ):
        """
        Class constructor
        
        @param  srcDipDir:  Dip direction of the plane (0-360�).
        @type  srcDipDir:  number or string convertible to float.
        @param  srcDipAngle:  Dip angle of the plane (0-90�).
        @type  srcDipAngle:  number or string convertible to float.
           
        @return:  GeolPlane.
    
        """
        
        assert 0.0 <= float( srcDipDir ) < 360.0
        assert 0.0 <= float( srcDipAngle ) <= 90.0
  
        self._dipdir = float( srcDipDir )
        self._dipangle = float( srcDipAngle )
        
        
    def to_normal_axis( self ):
        
        trend = self._dipdir + 180.0
        if trend >= 360.0:
            trend -= 360.0            
        plunge = 90.0 - self._dipangle
        
        return Axis( trend, plunge )
        
        
    def plane_x_coeff( self ):
        """
        Calculate the slope of a given plane along the x direction.
        The plane orientation  is expressed following the geological convention. 
               
        @return:  slope - float.    
        """ 
        return - sin( radians( self._dipdir ) ) * tan( radians( self._dipangle ) )


    def plane_y_coeff( self ):
        """
        Calculate the slope of a given plane along the y direction.
        The plane orientation  is expressed following the geological convention. 
               
        @return:  slope - float.     
        """ 
        return - cos( radians( self._dipdir ) ) * tan( radians( self._dipangle ) )

       
    def plane_from_geo( self, or_Pt ):
        """
        Closure that embodies the analytical formula for a given, non-vertical plane.
        This closure is used to calculate the z value from given horizontal coordinates (x, y).
    
        @param  or_Pt:  Point instance expressing a location point contained by the plane.
        @type  or_Pt:  Point.    
        
        @return:  lambda (closure) expressing an analytical formula for deriving z given x and y values.    
        """
    
        x0 =  or_Pt.x     
        y0 =  or_Pt.y
        z0 =  or_Pt.z
    
        # slope of the line parallel to the x axis and contained by the plane
        a = self.plane_x_coeff(  ) 
               
        # slope of the line parallel to the y axis and contained by the plane    
        b = self.plane_y_coeff(  )
                        
        return lambda x, y : a * ( x - x0 )  +  b * ( y - y0 )  +  z0
 
 
    def to_cartes_plane(self, point):
        
        normal_versor = self.to_normal_axis().to_down_axis().to_versor()
        
        a, b, c = normal_versor.x, normal_versor.y, normal_versor.z
        
        d = - ( a * point.x + b * point.y + c * point.z )
        
        return CartesianPlane( a, b, c, d )
    
      
class ArrCoord(object):
    """
    2D Array coordinates.
    Manages coordinates in the raster (array) space. 
    
    """
    
    def __init__( self, ival = 0.0, jval = 0.0 ):
        """
        Class constructor.

        @param  ival:  the i (-y) array coordinate of the point.
        @type  ival:  number or string convertible to float.
        @param  jval:  the j (x) array coordinate of the point.
        @type  jval:  number or string convertible to float.
               
        @return:  self.                
        """
        self._i = float( ival )
        self._j = float( jval )


    def g_i( self ):
        """
        Get i (row) coordinate value.
        
        @return:  the i (-y) array coordinate of the point - float.
        """
        return self._i
    

    def s_i( self, ival ):
        """
        Set i (row) coordinate value.
        
        @param  ival:  the i (-y) array coordinate of the point.
        @type  ival:  number or string convertible to float.
        
        @return:  self.        
        """
        self._i = float( ival )


    # set property for i
    i = property( g_i, s_i )
     
     
    def g_j( self ):
        """
        Get j (column) coordinate value.
        
        @return:  the j (x) array coordinate of the point - float.        
        """
        return self._j


    def s_j( self, jval ):
        """
        Set j (column) coordinate value.
        
        @param  jval:  the j (x) array coordinate of the point.
        @type  jval:  number or string convertible to float.
        
        @return:  self.         
        """
        self._j = jval


    # set property for j
    j = property( g_j, s_j )            
           
    
    def grid2geogcoord( self, currGeoGrid ):
        
        currPt_geogr_y = currGeoGrid.domain.g_trcorner().y - self.i*currGeoGrid.cellsize_y()
        currPt_geogr_x = currGeoGrid.domain.g_llcorner().x + self.j*currGeoGrid.cellsize_x()
        return Point( currPt_geogr_x, currPt_geogr_y )



class RectangularDomain(object):
    """
    Rectangular spatial domain class.
    
    """
    
    def __init__( self, pt_llc = None, pt_trc = None ): 
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

 
    def g_llcorner( self ):
        """
        Get lower-left corner of the spatial domain.
        
        @return:  lower-left corner of the spatial domain - Point.        
        """
        return self._llcorner


    def g_trcorner( self ):
        """
        Get top-right corner of the spatial domain.
        
        @return:  top-right corner of the spatial domain - Point.         
        """
        return self._trcorner
    
 
    def g_xrange( self ):
        """
        Get x range of spatial domain.
        
        @return:  x range - float.                 
        """
        return self._trcorner.x - self._llcorner.x

 
    def g_yrange( self ):
        """
        Get y range of spatial domain.
        
        @return:  y range - float.
        """
        return self._trcorner.y - self._llcorner.y


    def g_zrange( self ):
        """
        Get z range of spatial domain.
        
        @return:  z range - float.
        """
        return self._trcorner.z - self._llcorner.z


    def g_horiz_area( self ):
        """
        Get horizontal area of spatial domain.
        
        @return:  area - float.
        """
        return self.g_xrange() * self.g_yrange()


        
class Grid(object):
    """
    Grid class.
    Stores and manages the most of data and processing.
    
    """

    def __init__(self, source_filename= None, grid_params = None, grid_data = None):
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
            pt_llc = grid_params.llcorner()
            pt_trc = grid_params.trcorner()       
        else:
            pt_llc = None
            pt_trc = None
           
        self._grid_domain = RectangularDomain(pt_llc, pt_trc)          

        if grid_data is not None: 
            self._grid_data = grid_data.copy() 
        else:
            self._grid_data = None 
        
        
    def s_domain( self, domain ):
        """
        Set spatial domain.
        
        @param  domain:  Spatial domain to be attributed to the current Grid instance.
        @type  domain:  class RectangularDomain.
        
        @return: self
        """        
        del self._grid_domain
        self._grid_domain = copy.deepcopy(domain)
        
        
    def g_domain( self ):
        """
        Get spatial domain.
        
        @return: the spatial domain of the current Grid instance - class RectangularDomain.
        """        
        return self._grid_domain
    
        
    def d_domain( self ):
        """
        Delete current spatial domain of the Grid instance.
        
        @return: self
        """        
        del self._grid_domain


    # set property for spatial domain        
    domain = property( g_domain, s_domain, d_domain )

    
    def s_grid_data( self, data_array ):
        """
        Set grid data array.
        
        @param data_array: numpy.array of data values.
        @param type: 2D numpy.array.
        
        @return: self.
        """        
        if self._grid_data is not None:
            del self._grid_data
            
        self._grid_data = data_array.copy()  
        
            
    def g_grid_data( self ):
        """
        Get grid data array.
        
        @return: 2D numpy.array.
        """        
        return self._grid_data


    def d_grid_data( self ):
        """
        Delete grid data array.
        
        @return: self.
        """        
        del self._grid_data
            
    data = property( g_grid_data, s_grid_data, d_grid_data )
      

    def grid_extent( self ):
        """
        Return the xmin, xmax and ymin, ymax values as a dictionary
        """
        
        return dict(xmin=self.domain.g_llcorner().x, 
                    xmax=self.domain.g_trcorner().x,
                    ymin=self.domain.g_llcorner().y,
                    ymax=self.domain.g_trcorner().y )
        
        
    def _xmin(self):
        
        return self.grid_extent()['xmin']
    
    xmin = property( _xmin )
    
    
    def _xmax(self):
        
        return self.grid_extent()['xmax']
    
    xmax = property( _xmax )
        
        
    def _ymin(self):
        
        return self.grid_extent()['ymin']
    
    ymin = property( _ymin )
        
        
    def _ymax(self):
        
        return self.grid_extent()['ymax']   
            
    ymax = property( _ymax )
    
                
    def row_num( self ):
        """
        Get row number of the grid domain.     
        
        @return: number of rows of data array - int. 
        """        
        return np.shape( self.data )[0]        


    def col_num( self ):
        """
        Get column number of the grid domain.
        
        @return: number of columns of data array - int. 
        """        
        return np.shape( self.data )[1]        
            
 
    def cellsize_x( self ):
        """
        Get the cell size of the grid in the x direction.
        
        @return: cell size in the x (j) direction - float.
        """        
        return self.domain.g_xrange() / float( self.col_num() )

  
    def cellsize_y( self ):
        """
        Get the cell size of the grid in the y direction.
        
        @return: cell size in the y (-i) direction - float.
        """        
        return self.domain.g_yrange() / float( self.row_num() )
            
 
    def cellsize_h( self ):
        """
        Get the mean horizontal cell size.
        
        @return: mean horizontal cell size - float.
        """        
        return ( self.cellsize_x() + self.cellsize_y() ) / 2.0
              
   
    def geog2array_coord( self, curr_Pt ):
        """
        Converts from geographic to raster (array) coordinates.
        
        @param curr_Pt: point whose geographical coordinates will be converted to raster (array) ones.
        @type curr_Pt: Point.
        
        @return: point coordinates in raster (array) frame - class ArrCoord.
        """       
        currArrCoord_grid_i = ( self.domain.g_trcorner().y - curr_Pt.y ) / self.cellsize_y()
        currArrCoord_grid_j = ( curr_Pt.x - self.domain.g_llcorner().x ) / self.cellsize_x()
        
        return ArrCoord(currArrCoord_grid_i, currArrCoord_grid_j)


    def x( self ):
        """
        Creates an array storing the geographical coordinates of the cell centers along the x axis.
        Direction is from left to right.
        
        @return: numpy.array, shape: 1 x col_num.
        """        
        x_values = self.domain.g_llcorner().x + self.cellsize_x() * ( 0.5 + np.arange( self.col_num() ) ) 
        return x_values[ np.newaxis, : ]
        

    def y( self ):
        """
        Creates an array storing the geographical coordinates of the cell centers along the y axis.
        Direction is from top to bottom.
        
        @return: numpy.array, shape: row_num x 1.
        """        
        y_values = self.domain.g_trcorner().y - self.cellsize_y() * ( 0.5 + np.arange( self.row_num() ) )
        return y_values[ : , np.newaxis ]
        
              
    def grad_forward_y( self ):
        """
        Return an array representing the forward gradient in the y direction (top-wards), with values scaled by cell size.
        
        @return: numpy.array, same shape as current Grid instance        
        """        
        gf = np.zeros( np.shape( self.data ) ) * np.NaN
        gf[ 1:,:] = self.data[:-1,:] - self.data[1:,:] 
        
        return gf / float( self.cellsize_y() )

              
    def grad_forward_x( self ):
        """
        Return an array representing the forward gradient in the x direction (right-wards), with values scaled by cell size.

        @return: numpy.array, same shape as current Grid instance
        """
        gf = np.zeros( np.shape( self.data ),  ) * np.NaN
        gf[:, :-1] = self.data[:, 1:] - self.data[: , :-1] 
        
        return gf / float( self.cellsize_x() )
           
         
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

        assert currPt_cellcenter_i > 0,  currPt_cellcenter_j > 0
              
        grid_val_00 = self.data[int(floor(currPt_cellcenter_i)), int(floor(currPt_cellcenter_j))]
        grid_val_01 = self.data[int(floor(currPt_cellcenter_i)), int(ceil(currPt_cellcenter_j))]
        grid_val_10 = self.data[int(ceil(currPt_cellcenter_i)), int(floor(currPt_cellcenter_j))]    
        grid_val_11 = self.data[int(ceil(currPt_cellcenter_i)), int(ceil(currPt_cellcenter_j))]
                        
        delta_i = currPt_cellcenter_i - floor(currPt_cellcenter_i)
        delta_j = currPt_cellcenter_j - floor(currPt_cellcenter_j)                            
        
        grid_val_y0 = grid_val_00 + (grid_val_10-grid_val_00)*delta_i
        grid_val_y1 = grid_val_01 + (grid_val_11-grid_val_01)*delta_i        
                   
        grid_val_interp = grid_val_y0 + (grid_val_y1-grid_val_y0)*delta_j
                
        return grid_val_interp  
             

    def intersection_with_surface( self, surf_type, srcPt, srcPlaneAttitude ):
        """
        Calculates the intersections (as points) between DEM (self) and analytical surface.
        Currently it works only with planes as analytical surface cases.
        
        @param surf_type: type of considered surface (e.g., plane).
        @type surf_type: String.
        @param srcPt: point, expressed in geographical coordinates, that the plane must contain.
        @type srcPt: Point.
        @param srcPlaneAttitude: orientation of the surface (currently only planes).
        @type srcPlaneAttitude: class GeolPlane.
        
        @return: tuple of six arrays
        """
        
        if surf_type == 'plane': 
                        
            # lambdas to compute the geographic coordinates (in x- and y-) of a cell center 
            coord_grid2geog_x = lambda j : self.domain.g_llcorner().x + self.cellsize_x() * ( 0.5 + j )
            coord_grid2geog_y = lambda i : self.domain.g_trcorner().y - self.cellsize_y() * ( 0.5 + i )
             
            # arrays storing the geographical coordinates of the cell centers along the x- and y- axes    
            x_values = self.x()
            y_values = self.y()      

            ycoords_x, xcoords_y  = np.broadcast_arrays( x_values, y_values )
                        
            #### x-axis direction intersections
            
            # 2D array of DEM segment parameters                         
            x_dem_m = self.grad_forward_x()            
            x_dem_q = self.data - x_values * x_dem_m            
            
            # equation for the planar surface that, given (x,y), will be used to derive z  
            plane_z = srcPlaneAttitude.plane_from_geo( srcPt  )
            
            # 2D array of plane segment parameters
            x_plane_m = srcPlaneAttitude.plane_x_coeff(  )            
            x_plane_q = array_from_function( self.row_num(), 1, lambda j: 0, coord_grid2geog_y, plane_z )

            # 2D array that defines denominator for intersections between local segments
            x_inters_denomin =  np.where( x_dem_m != x_plane_m, x_dem_m-x_plane_m, np.NaN )            
            coincident_x = np.where( x_dem_q != x_plane_q, np.NaN, ycoords_x )            
            xcoords_x = np.where( x_dem_m != x_plane_m , (x_plane_q - x_dem_q ) / x_inters_denomin, coincident_x )
            
            xcoords_x = np.where( xcoords_x < ycoords_x , np.NaN, xcoords_x )           
            xcoords_x = np.where( xcoords_x >= ycoords_x + self.cellsize_x() , np.NaN, xcoords_x )  
                        
            
            #### y-axis direction intersections

            # 2D array of DEM segment parameters  
            y_dem_m = self.grad_forward_y()            
            y_dem_q = self.data - y_values * y_dem_m
 
            # 2D array of plane segment parameters
            y_plane_m = srcPlaneAttitude.plane_y_coeff(  )            
            y_plane_q = array_from_function( 1, self.col_num(), coord_grid2geog_x , lambda i: 0, plane_z )

            # 2D array that defines denominator for intersections between local segments
            y_inters_denomin =  np.where( y_dem_m != y_plane_m, y_dem_m - y_plane_m, np.NaN )
            coincident_y = np.where( y_dem_q != y_plane_q, np.NaN, xcoords_y )
                        
            ycoords_y = np.where( y_dem_m != y_plane_m, (y_plane_q - y_dem_q ) / y_inters_denomin, coincident_y )            

            # filter out cases where intersection is outside cell range
            ycoords_y = np.where( ycoords_y < xcoords_y , np.NaN, ycoords_y )           
            ycoords_y = np.where( ycoords_y >= xcoords_y + self.cellsize_y() , np.NaN, ycoords_y )            

            for i in xrange(xcoords_x.shape[0]):
                for j in xrange(xcoords_x.shape[1]):
                    if abs(xcoords_x[i,j]-ycoords_x[i,j])<1.0e-5 and abs(ycoords_y[i,j]-xcoords_y[i,j])<1.0e-5:
                        ycoords_y[i,j] = np.NaN
                                                         
            return xcoords_x, xcoords_y, ycoords_x, ycoords_y

    """
    def intersection_within_searchradius( self, surf_type, srcPt, srcPlaneAttitude ):

        
        if surf_type == 'plane': 
                        
            # lambdas to compute the geographic coordinates (in x- and y-) of a cell center 
            coord_grid2geog_x = lambda j : self.domain.g_llcorner().x + self.cellsize_x() * ( 0.5 + j )
            coord_grid2geog_y = lambda i : self.domain.g_trcorner().y - self.cellsize_y() * ( 0.5 + i )
             
            # arrays storing the geographical coordinates of the cell centers along the x- and y- axes    
            x_values = self.x()
            y_values = self.y()      

            ycoords_x, xcoords_y  = np.broadcast_arrays( x_values, y_values )
                        
            #### x-axis direction intersections
            
            # 2D array of DEM segment parameters                         
            x_dem_m = self.grad_forward_x()            
            x_dem_q = self.data - x_values * x_dem_m            
            
            # equation for the planar surface that, given (x,y), will be used to derive z  
            plane_z = srcPlaneAttitude.plane_from_geo( srcPt )
            
            # 2D array of plane segment parameters
            x_plane_m = srcPlaneAttitude.plane_x_coeff(  )            
            x_plane_q = array_from_function( self.row_num(), 1, lambda j: 0, coord_grid2geog_y, plane_z )

            # 2D array that defines denominator for intersections between local segments
            x_inters_denomin =  np.where( x_dem_m != x_plane_m, x_dem_m-x_plane_m, np.NaN )            
            coincident_x = np.where( x_dem_q != x_plane_q, np.NaN, ycoords_x )            
            xcoords_x = np.where( x_dem_m != x_plane_m , (x_plane_q - x_dem_q ) / x_inters_denomin, coincident_x )
            
            xcoords_x = np.where( xcoords_x < ycoords_x , np.NaN, xcoords_x )           
            xcoords_x = np.where( xcoords_x >= ycoords_x + self.cellsize_x() , np.NaN, xcoords_x )  
                        
            
            #### y-axis direction intersections

            # 2D array of DEM segment parameters  
            y_dem_m = self.grad_forward_y()            
            y_dem_q = self.data - y_values * y_dem_m
 
            # 2D array of plane segment parameters
            y_plane_m = srcPlaneAttitude.plane_y_coeff(  )            
            y_plane_q = array_from_function( 1, self.col_num(), coord_grid2geog_x , lambda i: 0, plane_z )

            # 2D array that defines denominator for intersections between local segments
            y_inters_denomin =  np.where( y_dem_m != y_plane_m, y_dem_m - y_plane_m, np.NaN )
            coincident_y = np.where( y_dem_q != y_plane_q, np.NaN, xcoords_y )
                        
            ycoords_y = np.where( y_dem_m != y_plane_m, (y_plane_q - y_dem_q ) / y_inters_denomin, coincident_y )            

            # filter out cases where intersection is outside cell range
            ycoords_y = np.where( ycoords_y < xcoords_y , np.NaN, ycoords_y )           
            ycoords_y = np.where( ycoords_y >= xcoords_y + self.cellsize_y() , np.NaN, ycoords_y )            

            for i in xrange(xcoords_x.shape[0]):
                for j in xrange(xcoords_x.shape[1]):
                    if abs(xcoords_x[i,j]-ycoords_x[i,j])<1.0e-5 and abs(ycoords_y[i,j]-xcoords_y[i,j])<1.0e-5:
                        ycoords_y[i,j] = np.NaN
                                                         
            return xcoords_x, xcoords_y, ycoords_x, ycoords_y
        """

