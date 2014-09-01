# -*- coding: utf-8 -*-

from __future__  import division

from math import sqrt, floor, ceil, sin, cos, tan, radians, asin, acos, atan, atan2, degrees

import numpy as np

import copy

from qgis.core import QgsCoordinateTransform, QgsPoint


from .utils import array_from_function
from .algebr_utils import point_solution




MINIMUM_SEPARATION_THRESHOLD = 1e-10
MINIMUM_VECTOR_MAGNITUDE = 1e-10


def qgs_point( x, y):
    
    return QgsPoint(x, y)
        

def project_qgs_point( qgsPt, srcCrs, destCrs ):
    
    return QgsCoordinateTransform( srcCrs, destCrs ).transform( qgsPt )


class Point_2D( object ):
    
    def __init__(self, x = np.nan, y = np.nan ):
        
        self._x = x
        self._y = y
    

    def copy(self):
        
        return Point_2D( self._x, self._y )
    
    
    def distance(self, another ):
        
        return sqrt( (self._x - another._x)**2 + (self._y - another._y)**2 )
    
    
    def to_3D( self, z = 0.0 ):
        
        return Point_3D( self._x, self._y, z )
    

    def displaced_by_vector(self, displacement_vector ):
        
        return Point_2D( self._x + displacement_vector._x , self._y + displacement_vector._y )
        

    def is_coincident_with( self, another, tolerance = 1.0e-7 ):
        
        if self.distance(another) > tolerance:
            return False
        else:
            return True
        

    def project_crs( self, srcCrs, destCrs ):
        
        qgis_pt = qgs_point( self._x, self._y )
        destCrs_qgis_pt = project_qgs_point( qgis_pt, srcCrs, destCrs )
        
        return Point_2D( destCrs_qgis_pt.x(), destCrs_qgis_pt.y() )
  
                        
class Segment_2D( object ):
    
    def __init__(self, start_pt_2d, end_pt_2d ): 
        
        self._start_pt = start_pt_2d.copy()
        self._end_pt = end_pt_2d.copy()
                   
    
    def copy( self ):
        
        return Segment_2D( self._start_pt, self._end_pt )

              
    def to_vector( self ):
        
        return Vector_2D( self._end_pt._x - self._start_pt._x,
                          self._end_pt._y - self._start_pt._y ) 
        
            
    def length( self ):
        
        return self._start_pt.distance( self._end_pt )
        

    def delta_x( self ):
        
        return self._end_pt._x - self._start_pt._x
    
    
    def delta_y( self ):

        return self._end_pt._y - self._start_pt._y
    
                    
    def scale( self, scale_factor ):
        
        delta_x = self.delta_x() * scale_factor
        delta_y = self.delta_y() * scale_factor
        
        return Segment_2D( self._start_pt, Point_2D( self._start_pt._x + delta_x, self._start_pt._y + delta_y ) )        
        
        
    def to_3D( self ):
        
        return Segment_3D( self._start_pt.to_3D(), self._end_pt.to_3D() ) 


    def densify_nodes( self, densify_distance ):
        
        assert densify_distance > 0.0
        segment_length = self.length()
        assert segment_length > 0.0
  
        generator_vector = self.to_vector().to_versor().scale( densify_distance )    
        interpolated_line = Line_2D( [ self._start_pt ] )       
        n = 0        
        while ( True ):        
            n += 1 
            new_pt = self._start_pt.displaced_by_vector( generator_vector.scale( n ) )
            if self._start_pt.distance(new_pt) >= segment_length:
                break        
            interpolated_line = interpolated_line.add_pt( new_pt )            
        interpolated_line = interpolated_line.add_pt( self._end_pt )
            
        return interpolated_line
            

class Vector_2D( object ):
    
    def __init__(self, x = np.nan, y = np.nan ):
        
        self._x = x
        self._y = y
       
    
    def copy(self):
        
        return Vector_2D( self._x, self._y )
    
       
    def length( self ):
       
        return sqrt( self._x * self._x + self._y * self._y )
    

    def scale(self, scale_factor ):
        
        return Vector_2D( self._x * scale_factor, self._y * scale_factor )
    

    def to_versor( self ):
        
        return self.scale( 1.0 / self.length() )
    

    def add(self, another ):
        
        return Vector_2D( self._x + another._x, self._y + another._y )
    
            
    def minus(self, another):
        
        return self.add( another.scale(-1) )
    
    
    def to_3D( self, z = 0.0 ):
        
        return Vector_3D( self._x, self._y, z )


class Line_2D( object ):
    
    def __init__( self, pt_2d_list = [] ):
        
        self._pts = [ pt_2d.copy() for pt_2d in pt_2d_list ]
        

    def copy( self ):
        
        return Line_2D( self._pts )
    
            
    def add_pt(self, pt_2d ):
        
        return Line_2D( self._pts + [ pt_2d ] )
        
        
    def add_pts(self, pt_2d_list ):
        
        return Line_2D( self._pts + pt_2d_list )
    

    def num_points( self ):
        
        return len( self._pts )
    

    def x_list(self):
        
        return [pt._x for pt in self._pts ]
    
    
    def y_list(self):
        
        return [pt._y for pt in self._pts ]
    
    
    def x_y_lists( self ):

        return self.x_list(), self.y_list()
    
       
    def x_min( self ):
        
        return min( self.x_list() )
    
    
    def x_max(self):
        
        return max( self.x_list() )
    
    
    def y_min(self):
        
        return min( self.y_list() )
    
    
    def y_max(self):
        
        return max( self.y_list() )
    
            
    def remove_coincident_successive_points( self ):
        
        assert self.num_points() > 0        
        new_line = Line_2D( [ self._pts[ 0 ] ] )
        for ndx in range(1, self.num_points() ):
            if not self._pts[ndx].is_coincident_with( new_line._pts[-1] ):
                new_line = new_line.add_pt( self._pts[ndx] )
        return new_line
    

    def to_segments( self ):

        pts_pairs = zip( self._pts[:-1], self._pts[1:] )
        return [ Segment_2D( pt_a, pt_b ) for ( pt_a, pt_b ) in pts_pairs ]
        
        
    def densify_nodes( self, sample_distance ):

        assert sample_distance > 0.0
        densified_line_list = [ segment.densify_nodes( sample_distance ) for segment in self.to_segments() ]
        assert len( densified_line_list ) > 0
        return MultiLine_2D( densified_line_list ).to_line().remove_coincident_successive_points()
        

    def length( self ):
        
        length = 0.0
        for ndx in range( self.num_points()-1 ):
            length += self._pts[ndx].distance( self._pts[ndx+1] )
        return length            


    def incremental_length( self ):
        
        incremental_length_list = []
        length = 0.0
        incremental_length_list.append( length )
        for ndx in range( self.num_points()-1 ):            
            length += self._pts[ndx].distance( self._pts[ndx+1] )
            incremental_length_list.append( length )        
            
        return incremental_length_list         
        
        
    def project_crs( self, srcCrs, destCrs ):

            points = []            
            for point in self._pts:
                destCrs_point = point.project_crs( srcCrs, destCrs )
                points.append( destCrs_point )
                
            return Line_2D( points )
                
                                          
class MultiLine_2D(object):
    # MultiLine_2D is a list of Line_2D objects
    
    def __init__( self, lines_list = [] ):
     
        self._lines = [ line_2d.copy() for line_2d in lines_list ]    
    
    
    def add( self, line ):
        
        return MultiLine_2D( self._lines + [ line ] )
    
    
    def copy( self ):
        
        return MultiLine_2D( self._lines )
    
    
    def num_parts( self ):
        
        return len( self._lines )
    
    
    def num_points( self ):

        num_elements = map( lambda x: len( x._pts ), self._lines )
        return reduce( lambda x, y: x + y, num_elements )
    

    def is_continuous( self ):
        
        for line_ndx in range( len(self._lines) - 1 ):
            if not self._lines[line_ndx]._pts[-1].is_coincident_with( self._lines[line_ndx+1]._pts[0] ) or \
               not self._lines[line_ndx]._pts[-1].is_coincident_with( self._lines[line_ndx+1]._pts[-1] ):
                return False            
        return True


    def x_min(self):
        
        return min( [ line.x_min() for line in self._lines ] )
    
    
    def x_max(self):
        
        return max( [ line.x_max() for line in self._lines ] )
     

    def y_min(self):
        
        return min( [ line.y_min() for line in self._lines ] )
    
    
    def y_max(self):
        
        return max( [ line.y_max() for line in self._lines ] )
     
            
    def is_unidirectional( self ):        
        
        for line_ndx in range( len(self._lines) - 1):
            if not self._lines[line_ndx]._pts[-1].is_coincident_with( self._lines[line_ndx+1]._pts[0] ):
                return False            
        return True


    def to_line( self ):

        return Line_2D( [ point for line in self._lines for point in line._pts ] )


    def project_crs( self, srcCrs, destCrs ):

        lines = []        
        for line_2d in self._lines:
            lines.append( line_2d.project_crs( srcCrs, destCrs ) )
                        
        return MultiLine_2D( lines ) 
        

    def densify( self, sample_distance ):
       
        densified_multiline_2d_list = []
        for line_2d in self._lines:
            densified_multiline_2d_list.append( line_2d.densify_nodes( sample_distance ) )
   
        return MultiLine_2D( densified_multiline_2d_list )
     

    def remove_coincident_points( self ):
        
        cleaned_lines = []
        for line_2d in self._lines:
            cleaned_lines.append( line_2d.remove_coincident_successive_points() )
        
        return MultiLine_2D( cleaned_lines )
        
        
class Point_3D( object ):

        
    # class constructor
    def __init__( self, x = np.nan, y = np.nan, z = np.nan ):
     
        self._x = x
        self._y = y
        self._z = z
    

    def copy( self ):
        
        return Point_3D( self._x, self._y, self._z )
    
    
    def to_2D( self ):
        
        return Point_2D( self._x, self._y ), self._z
            
                  
    def distance( self, another ):
        """
        Calculate Euclidean distance between two points.

        @param  another:  the Point_3D instance for which the distance should be calculated
        @type  another:  Point_3D.
        
        @return:  distance between the two points - float.        
        """
        
        return sqrt( ( self._x - another._x ) ** 2 + ( self._y - another._y ) ** 2 + ( self._z - another._z ) ** 2 )
        

    def distance_2D( self, another ):

        return sqrt( ( self._x - another._x ) ** 2 + ( self._y - another._y ) ** 2 )
    
    
    def is_coincident_with( self, another, tolerance = 1.0e-6 ):
        
        if self.distance(another) > tolerance:
            return False
        else:
            return True
        
                
    def movedby( self, sx = 0.0 , sy = 0.0 , sz = 0.0 ):
        """
        Create a new point shifted by given amount from the self instance.

        @param  sx:  the shift to be applied along the x axis.
        @type  sx:  float.
        @param  sy:  the shift to be applied along the y axis.
        @type  sy:  float.
        @param  sz:  the shift to be applied along the z axis.
        @type  sz:  float.
                
        @return:  a new Point_3D instance shifted by the given amounts with respect to the original one.        
        """
        
        return Point_3D( self._x + sx , self._y + sy, self._z + sz )
    

    def displaced_by_vector(self, displacement_vector ):
        
        return Point_3D( self._x + displacement_vector._x , self._y + displacement_vector._y, self._z + displacement_vector._z )


    def as_vector( self ):
        
        return Vector_3D( self._x, self._y, self._z )
                

class Segment_3D( object ):
    
    def __init__(self, start_point, end_point ):
        
        self._start_pt = start_point.copy()
        self._end_pt = end_point.copy()
          

    def copy( self ):
        
        return Segment_3D( self._start_pt, self._end_pt )
    
              
    def to_vector( self ):
        
        return Vector_3D( self._end_pt._x - self._start_pt._x,
                          self._end_pt._y - self._start_pt._y,
                          self._end_pt._z - self._start_pt._z ) 
        

    def length( self ):
        
        return self._start_pt.distance( self._end_pt )
    
            
    def trend_and_plunge( self ):
        
        as_geol_axis = self.to_vector().to_geol_axis()
        return as_geol_axis._trend, as_geol_axis._plunge
           

    def to_vertical_cartes_plane( self ):
        """
        Creates a vertical Cartesian plane passing through the self Segment_3D
        """
          
        trend, _ = self.trend_and_plunge()
        dip_dir = trend + 90.0
        if dip_dir >= 360.0:
            dip_dir -= 360.0        
        return GeolPlane( dip_dir, 90.0 ).to_cartes_plane( self._start_pt )    


    def densify_nodes( self, densify_distance ):
        
        length = self.length()
        assert length > 0.0
  
        generator_vector = self.to_vector().to_versor().scale( densify_distance )
    
        interpolated_line = Line_3D( [ self._start_pt ] )       
        n = 0        
        while ( True ):        
            n += 1 
            new_pt = self._start_pt.displaced_by_vector( generator_vector.scale( n ) )
            if self._start_pt.distance(new_pt) >= length:
                break        
            interpolated_line = interpolated_line.add_pt( new_pt )            
        interpolated_line = interpolated_line.add_pt( self._end_pt )
            
        return interpolated_line
            
                
class Vector_3D( object ):
    
    def __init__(self, x = np.nan, y = np.nan, z = np.nan ):
        
        self._x = x
        self._y = y
        self._z = z


    def copy( self ):
        
        return Vector_3D( self._x, self._y, self._z ) 
    
          
    def length(self):
        
        return sqrt( self._x * self._x + self._y * self._y + self._z * self._z )
        

    def length_horiz(self):
        
        return sqrt( self._x * self._x + self._y * self._y )
    
            
    def scale( self, scale_factor ):
                
        return Vector_3D( self._x * scale_factor, 
                          self._y * scale_factor, 
                          self._z * scale_factor )
        
  
    def to_versor(self):
        
        return self.scale( 1.0 / self.length() )
    

    def to_down_vector(self):
        
        if self._z > 0.0:
            return self.scale(-1.0)
        else:
            return self.copy()
    
        
    def add(self, another):
        
        return Vector_3D( self._x + another._x, 
                          self._y + another._y, 
                          self._z + another._z )
        

    def slope_radians(self):
        
        return atan( self._z / self.length_horiz() )
              

    def to_geol_axis(self):
        
        if self.length() < MINIMUM_VECTOR_MAGNITUDE:
            return None
        
        unit_vect = self.to_versor()
        
        plunge = - degrees( asin( unit_vect._z ) ) # upward negative, downward positive
        
        trend = 90.0 - degrees( atan2( unit_vect._y, unit_vect._x ) )
        if trend < 0.0: trend += 360.0
        elif trend > 360.0: trend -= 360.0
        
        assert 0.0 <= trend < 360.0
        assert -90.0 <= plunge <= 90.0
         
        return GeolAxis( trend, plunge )
    
    
    def scalar_product( self, another ):
        
        return self._x * another._x + self._y * another._y + self._z * another._z


    def vectors_cos_angle( self, another ):
 
        try:
            return self.scalar_product( another ) / ( self.length() * another.length() )
        except ZeroDivisionError:
            return np.nan
        
        
    def angle_degr( self, another ):
        
        return degrees( acos ( self.vectors_cos_angle( another ) ) )
    

    def vector_product( self, another ):
        
        x = (self._y * another._z) - (self._z * another._y)
        y = (self._z * another._x) - (self._x * another._z)
        z = (self._x * another._y) - (self._y * another._x)

        return Vector_3D( x, y, z )


    def by_matrix( self, matrix3x3 ):
        
        vx = matrix3x3[0,0] * self._x + matrix3x3[0,1] * self._y + matrix3x3[0,2] * self._z 
        vy = matrix3x3[1,0] * self._x + matrix3x3[1,1] * self._y + matrix3x3[1,2] * self._z 
        vz = matrix3x3[2,0] * self._x + matrix3x3[2,1] * self._y + matrix3x3[2,2] * self._z 
        
        return Vector_3D( vx, vy, vz )
        

    def as_point( self ):
        
        return Point_3D( self._x, self._y, self._z )


class Line_3D(object):
    # Line_3D is a list of Point_3D objects

    def __init__( self, pt_3d_list = [] ):
     
        self._pts = [ pt_3d.copy() for pt_3d in pt_3d_list ]
        
        
    def copy( self ):
        
        return Line_3D( self._pts )
        

    def add_pt( self, pt ):
        
        return Line_3D( self._pts + [ pt ] )
        
 
    def add_pts( self, pt_list ):

        return Line_3D( self._pts + pt_list )
                        

    def remove_coincident_successive_points( self ):
        
        new_line = Line_3D( self._pts[ : 1 ] )
        for ndx in range(1, self.num_points() ):
            if not self._pts[ndx].is_coincident_with( new_line._pts[-1] ):
                new_line = new_line.add_point( self._pts[ndx] )
        return new_line
    

    def join( self, another ):
        """
        Joins together two lines and returns the join as a new line without point changes, 
        with possible overlapping points 
        and orientation mismatches between the two original lines
        """
        
        return Line_3D( self._pts + another._pts )
   
                    
    def num_points( self ):
        
        return len( self._pts )
    
        
    def length( self ):
        
        length = 0.0
        for ndx in range( self.num_points()-1 ):
            length += self._pts[ndx].distance( self._pts[ndx+1] )
        return length            


    def length_2d( self ):
        
        length = 0.0
        for ndx in range( self.num_points()-1 ):
            length += self._pts[ndx].distance_2D( self._pts[ndx+1] )
        return length 
    
    
    def z_min(self):
        
        return min( [ pt._z for pt in self._pts ] )
    
    
    def z_max(self):
        
        return max( [ pt._z for pt in self._pts ] )
    
    
    def incremental_length( self ):
        
        incremental_length_list = []
        length = 0.0
        incremental_length_list.append( length )
        for ndx in range( self.num_points()-1 ):            
            length += self._pts[ndx].distance( self._pts[ndx+1] )
            incremental_length_list.append( length )        
            
        return incremental_length_list         
        

    def incremental_length_2d( self ):
        
        incremental_length_list = []
        length = 0.0
        incremental_length_list.append( length )
        for ndx in range( self.num_points()-1 ):            
            length += self._pts[ndx].distance_2D( self._pts[ndx+1] )
            incremental_length_list.append( length )        
            
        return incremental_length_list
    
    
    def reverse_direction( self ):
        
        new_line = self.copy()
        new_line._pts.reverse() # in-place operation on new_line
        
        return new_line
    
                
    def slopes_list( self ):
  
        slopes_list = []
        for ndx in range( self.num_points() - 1 ):            
            vector = Segment_3D( self._pts[ndx], self._pts[ndx+1] ).to_vector()
            slopes_list.append( degrees( vector.slope_radians() ) ) 
        slopes_list.append( np.nan ) # slope value for last point is unknown         
        return slopes_list        
 
  
class MultiLine_3D(object):
    # MultiLine_3D is a list of Line_3D objects
    
    def __init__( self, lines_list ):
     
        self._lines = lines_list    
    
    
    def num_parts( self ):
        
        return len( self._lines )
    
    
    def num_points( self ):
        
        num_points = 0
        for line in self._lines:
            num_points += line.num_points()            
        return num_points
    

    def is_continuous( self ):
        
        for line_ndx in range( len(self._lines) - 1 ):
            if not self._lines[line_ndx]._pts[-1].is_coincident_with( self._lines[line_ndx+1]._pts[0] ) or \
               not self._lines[line_ndx]._pts[-1].is_coincident_with( self._lines[line_ndx+1]._pts[-1] ):
                return False
            
        return True
        
            
    def is_unidirectional( self ):        
        
        for line_ndx in range( len(self._lines) - 1):
            if not self._lines[line_ndx]._pts[-1].is_coincident_with( self._lines[line_ndx+1]._pts[0] ):
                return False
            
        return True


    def to_line( self ):

        return Line_3D( [ point for line in self._lines for point in line._pts ] )
      
                        
class Point_4D( Point_3D ):
    
    def __init__(self, x = np.nan, y = np.nan, z = 0.0, time = np.nan ):
        
        super(Point_4D, self).__init__(x, y, z)
        self._t = time
                
    
    def delta_time(self, anotherPt4D): 
        
        return anotherPt4D._t - self._t
    
    
    def speed(self, anotherPt4D):
        
        try:
            return self.distance( anotherPt4D ) / self.delta_time( anotherPt4D )
        except:
            return np.nan


class ParamLine( object ):
    # parametric line
    # srcPt: source Point_3D
    # l, m, n: .....
    
    def __init__(self, srcPt, l, m, n ):
        
        assert  -1.0 <= l <= 1.0
        assert  -1.0 <= m <= 1.0        
        assert  -1.0 <= n <= 1.0        
               
        self._srcPt = srcPt
        self._l = l     
        self._m = m    
        self._n = n    
    

    def intersect_cartes_plane( self, cartes_plane ):
        """
        Return intersection point between parameteric line and Cartesian plane
        """
         
        # line parameters
        x1, y1, z1 = self._srcPt._x, self._srcPt._y, self._srcPt._z
        l, m, n = self._l, self._m, self._n
        
        # Cartesian plane parameters
        a, b, c, d = cartes_plane._a, cartes_plane._b, cartes_plane._c, cartes_plane._d
        
        try:
            k = ( a * x1 + b * y1 + c * z1 + d ) / ( a * l + b * m + c * n )
        except ZeroDivisionError:
            return None
        
        return Point_3D( x1 - l * k,
                      y1 - m * k,
                      z1 - n * k )
        
        
class GeolAxis( object ):
    """
    Structural axis,
    defined by trend and plunge (both in degrees)
    Trend range: [0.0, 360.0[  clockwise, from 0 (North) 
    Plunge: [-90.0, 90.0], negative value: upward axis, positive values: downward axis
    """
    
    def __init__(self, srcTrend, srcPlunge ):

        assert 0.0 <= srcTrend < 360.0
        assert -90.0 <= srcPlunge <= 90.0
        
        self._trend = srcTrend 
        self._plunge = srcPlunge         
        

    def to_versor( self ):
        
        north_coord = cos( radians( self._plunge ) ) * cos( radians( self._trend ) )
        east_coord = cos( radians( self._plunge ) ) * sin( radians( self._trend ) )
        down_coord = sin( radians( self._plunge ) )
        
        return Vector_3D( east_coord, north_coord, -down_coord )
        
             
    def to_down_axis( self ): 
        
        trend, plunge = self._trend, self._plunge
        if plunge < 0.0:
            trend += 180.0
            if trend > 360.0: trend -= 360.0
            plunge = - plunge
        
        return GeolAxis( trend, plunge )


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
                
        
    def to_normal_versor(self):
        """
        return the normal versor to the cartesian plane
        """
                
        return Vector_3D( self._a, self._b, self._c ).to_versor()
        
        
    def to_geolplane_and_point(self):
        """
        converts a cartesian plane into a geological plane
        and a point lying in the plane (non-unique solution)
        """

        geol_plane = self.to_normal_versor().to_geol_axis().to_normal_geolplane()               
        point = Point_3D( point_solution( np.array([[self._a,self._b,self._c]]),
                                          np.array([-self._d])))        
        return geol_plane, point
    
        
    def intersection_versor(self, another):
        """
        return intersection versor for two intersecting planes
        """

        return self.to_normal_versor().vector_product( another.to_normal_versor() ).to_versor()


    def intersection_point(self, another):
        """
        return point on intersection line (obviously non-unique solution)
        for two planes
        """
        
        # find a point lying on the intersection line (this is a non-unique solution)    
        a = np.array([[self._a, self._b, self._c], [another._a, another._b, another._c]])
        b = np.array([-self._d, -another._d]) 
        x, y, z = point_solution( a, b ) 
              
        return Point_3D( x, y, z )


    def set_point_inside( self, pt ):
        
        return self._a * pt._x + self._b * pt._y + self._c * pt._z + self._d
    

    def angle_degr( self, another ):
        
        angle_degr = self.to_normal_versor().angle_degr( another.to_normal_versor() )
        assert angle_degr > 0.0
        if angle_degr > 90.0:
            angle_degr = 180.0 - angle_degr
        return angle_degr   
                
        
def cartes_plane_from_points( pt1, pt2, pt3 ):
    
    matr_a = np.array([[pt1._y, pt1._z, 1], 
                       [pt2._y, pt2._z, 1],
                       [pt3._y, pt3._z, 1] ])

    matr_b = - np.array([[pt1._x, pt1._z, 1], 
                       [pt2._x, pt2._z, 1],
                       [pt3._x, pt3._z, 1] ])  
    
    matr_c = np.array([[pt1._x, pt1._y, 1], 
                       [pt2._x, pt2._y, 1],
                       [pt3._x, pt3._y, 1] ])
    
    matr_d = - np.array([[pt1._x, pt1._y, pt1._z], 
                       [pt2._x, pt2._y, pt2._z],
                       [pt3._x, pt3._y, pt3._z] ])
    
    return CartesianPlane( np.linalg.det(matr_a),
                           np.linalg.det(matr_b),
                           np.linalg.det(matr_c),
                           np.linalg.det(matr_d) ) 
     
        
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
        
        return GeolAxis( trend, plunge )
        
        
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
    
        @param  or_Pt:  Point_3D instance expressing a location point contained by the plane.
        @type  or_Pt:  Point_3D.    
        
        @return:  lambda (closure) expressing an analytical formula for deriving z given x and y values.    
        """
    
        x0 =  or_Pt._x     
        y0 =  or_Pt._y
        z0 =  or_Pt._z
    
        # slope of the line parallel to the x axis and contained by the plane
        a = self.plane_x_coeff(  ) 
               
        # slope of the line parallel to the y axis and contained by the plane    
        b = self.plane_y_coeff(  )
                        
        return lambda x, y : a * ( x - x0 )  +  b * ( y - y0 )  +  z0
 
 
    def to_cartes_plane(self, point):
        
        normal_versor = self.to_normal_axis().to_down_axis().to_versor()        
        a, b, c = normal_versor._x, normal_versor._y, normal_versor._z        
        d = - ( a * point._x + b * point._y + c * point._z )        
        return CartesianPlane( a, b, c, d )
    

def eq_xy_pair( xy_pair_1, xy_pair_2 ):

    if xy_pair_1[0] == xy_pair_2[0] and xy_pair_1[1] == xy_pair_2[1]:
        return True
    
    return False

 
def remove_equal_consecutive_xypairs( xy_list ):
    
    out_xy_list = [ xy_list[0] ]
    
    for n in range( 1, len( xy_list ) ):
        if not eq_xy_pair( xy_list[n], out_xy_list[-1] ):
            out_xy_list.append( xy_list[n] )
            
    return out_xy_list
    
           
def xytuple_list_to_Line_2D( xytuple_list ):

    return Line_2D( [ Point_2D(x,y) for (x,y) in xytuple_list ] )

    
def xytuple_list2_to_MultiLine_2D( xytuple_list2 ):
    # input is a list of list of (x,y) values
    
    assert len( xytuple_list2 ) > 0
    lines_list = []
    for xy_list in xytuple_list2:
        assert len( xy_list ) > 0
        lines_list.append( xytuple_list_to_Line_2D( xy_list ) )
        
    return MultiLine_2D( lines_list )


def list2_to_list( list2 ):
    """
    input: a list of list of (x,y) tuples
    output: a list of (x,y) tuples
    """

    out_list = []
    for list1 in list2:
        for el in list1:
            out_list.append( el )
        
    return out_list


def list3_to_list( list3 ):
    """
    input: a list of list of (x,y) tuples
    output: a list of (x,y) tuples
    """

    out_list = []
    for list2 in list3:
        for list1 in list2:
            out_list += list1
        
    return out_list

     
def merge_lines( lines, progress_ids ):
    """
    lines: a list of list of (x,y,z) tuples for multilines
    """

    sorted_line_list = [line for (_, line) in sorted(zip( progress_ids, lines))]

    line_list = []
    for line in sorted_line_list:
       
        line_type, line_geometry = line 
     
        if line_type == 'multiline': 
            path_line = xytuple_list2_to_MultiLine_2D( line_geometry ).to_line()
        elif line_type == 'line':            
            path_line = xytuple_list_to_Line_2D( line_geometry ) 
        line_list.append( path_line )  # now a list of Lines     
                
    # now the list of Lines is transformed into a single Line_2D 
    return MultiLine_2D( line_list ).to_line().remove_coincident_successive_points()
               
                   
class ArrCoord(object):
    """
    2D Array coordinates.
    Manages coordinates in the raster (array) space. 
    
    """
    
    def __init__( self, ival = 0.0, jval = 0.0 ):
        """
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
        
        currPt_geogr_y = currGeoGrid.domain.g_trcorner()._y - self.i * currGeoGrid.cellsize_y()
        currPt_geogr_x = currGeoGrid.domain.g_llcorner()._x + self.j * currGeoGrid.cellsize_x()
        return Point_3D( currPt_geogr_x, currPt_geogr_y )



class RectangularDomain(object):
    """
    Rectangular spatial domain class.
    
    """
    
    def __init__( self, pt_llc = None, pt_trc = None ): 
        """
        Class constructor.
        
        @param  pt_llc:  lower-left corner of the domain.
        @type  pt_llc:  Point_3D.
        @param  pt_trc:  top-right corner of the domain.
        @type  pt_trc:  Point_3D.
                        
        @return:  RectangularDomain instance.
        """     
        self._llcorner = pt_llc  
        self._trcorner = pt_trc 

 
    def g_llcorner( self ):
        """
        Get lower-left corner of the spatial domain.
        
        @return:  lower-left corner of the spatial domain - Point_3D.        
        """
        return self._llcorner


    def g_trcorner( self ):
        """
        Get top-right corner of the spatial domain.
        
        @return:  top-right corner of the spatial domain - Point_3D.         
        """
        return self._trcorner
    
 
    def g_xrange( self ):
        """
        Get x range of spatial domain.
        
        @return:  x range - float.                 
        """
        return self._trcorner._x - self._llcorner._x

 
    def g_yrange( self ):
        """
        Get y range of spatial domain.
        
        @return:  y range - float.
        """
        return self._trcorner._y - self._llcorner._y


    def g_zrange( self ):
        """
        Get z range of spatial domain.
        
        @return:  z range - float.
        """
        return self._trcorner._z - self._llcorner._z


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
        
        return dict(xmin=self.domain.g_llcorner()._x, 
                    xmax=self.domain.g_trcorner()._x,
                    ymin=self.domain.g_llcorner()._y,
                    ymax=self.domain.g_trcorner()._y )
        
        
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
        @type curr_Pt: Point_3D.
        
        @return: point coordinates in raster (array) frame - class ArrCoord.
        """       
        currArrCoord_grid_i = ( self.domain.g_trcorner()._y - curr_Pt._y ) / self.cellsize_y()
        currArrCoord_grid_j = ( curr_Pt._x - self.domain.g_llcorner()._x ) / self.cellsize_x()
        
        return ArrCoord(currArrCoord_grid_i, currArrCoord_grid_j)


    def x( self ):
        """
        Creates an array storing the geographical coordinates of the cell centers along the x axis.
        Direction is from left to right.
        
        @return: numpy.array, shape: 1 x col_num.
        """        
        x_values = self.domain.g_llcorner()._x + self.cellsize_x() * ( 0.5 + np.arange( self.col_num() ) ) 
        return x_values[ np.newaxis, : ]
        

    def y( self ):
        """
        Creates an array storing the geographical coordinates of the cell centers along the y axis.
        Direction is from top to bottom.
        
        @return: numpy.array, shape: row_num x 1.
        """        
        y_values = self.domain.g_trcorner()._y - self.cellsize_y() * ( 0.5 + np.arange( self.row_num() ) )
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
        
        grid_val_y0 = grid_val_00 + (grid_val_10-grid_val_00) * delta_i
        grid_val_y1 = grid_val_01 + (grid_val_11-grid_val_01) * delta_i        
                   
        grid_val_interp = grid_val_y0 + (grid_val_y1-grid_val_y0) * delta_j
                
        return grid_val_interp  
             

    def intersection_with_surface( self, surf_type, srcPt, srcPlaneAttitude ):
        """
        Calculates the intersections (as points) between DEM (the self object) and an analytical surface.
        Currently it works only with planes.
        
        @param surf_type: type of considered surface (i.e., plane, the only case implemented at present).
        @type surf_type: String.
        @param srcPt: point, expressed in geographical coordinates, that the plane must contain.
        @type srcPt: Point_3D.
        @param srcPlaneAttitude: orientation of the surface (currently only planes).
        @type srcPlaneAttitude: class GeolPlane.
        
        @return: tuple of four arrays
        """
        
        if surf_type == 'plane': 
                        
            # closures to compute the geographic coordinates (in x- and y-) of a cell center
            # the grid coordinates of the cell center are expressed by i and j 
            grid_coord_to_geogr_coord_x_closure = lambda j : self.domain.g_llcorner()._x + self.cellsize_x() * ( 0.5 + j )
            grid_coord_to_geogr_coord_y_closure = lambda i : self.domain.g_trcorner()._y - self.cellsize_y() * ( 0.5 + i )
             
            # arrays storing the geographical coordinates of the cell centers along the x- and y- axes    
            cell_center_x_array = self.x()
            cell_center_y_array = self.y()      

            ycoords_x, xcoords_y  = np.broadcast_arrays( cell_center_x_array, cell_center_y_array )
                        
            #### x-axis direction intersections
            
            # 2D array of DEM segment parameters                         
            x_dem_m = self.grad_forward_x()            
            x_dem_q = self.data - cell_center_x_array * x_dem_m            
            
            # closure for the planar surface that, given (x,y), will be used to derive z  
            plane_z_closure = srcPlaneAttitude.plane_from_geo( srcPt  )
            
            # 2D array of plane segment parameters
            x_plane_m = srcPlaneAttitude.plane_x_coeff(  )            
            x_plane_q = array_from_function( self.row_num(), 1, lambda j: 0, grid_coord_to_geogr_coord_y_closure, plane_z_closure )

            # 2D array that defines denominator for intersections between local segments
            x_inters_denomin =  np.where( x_dem_m != x_plane_m, x_dem_m-x_plane_m, np.NaN )  
                      
            coincident_x = np.where( x_dem_q != x_plane_q, np.NaN, ycoords_x )
                        
            xcoords_x = np.where( x_dem_m != x_plane_m , (x_plane_q - x_dem_q ) / x_inters_denomin, coincident_x )            
            xcoords_x = np.where( xcoords_x < ycoords_x , np.NaN, xcoords_x )           
            xcoords_x = np.where( xcoords_x >= ycoords_x + self.cellsize_x() , np.NaN, xcoords_x )  
                        
            
            #### y-axis direction intersections

            # 2D array of DEM segment parameters  
            y_dem_m = self.grad_forward_y()            
            y_dem_q = self.data - cell_center_y_array * y_dem_m
 
            # 2D array of plane segment parameters
            y_plane_m = srcPlaneAttitude.plane_y_coeff(  )            
            y_plane_q = array_from_function( 1, self.col_num(), grid_coord_to_geogr_coord_x_closure , lambda i: 0, plane_z_closure )

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
            grid_coord_ to_geogr_coord_x = lambda j : self.domain.g_llcorner()._x + self.cellsize_x() * ( 0.5 + j )
            grid_coord_to_geogr_coord_y_closure = lambda i : self.domain.g_trcorner()._y - self.cellsize_y() * ( 0.5 + i )
             
            # arrays storing the geographical coordinates of the cell centers along the x- and y- axes    
            cell_center_x_array = self._x()
            cell_center_y_array = self._y()      

            ycoords_x, xcoords_y  = np.broadcast_arrays( cell_center_x_array, cell_center_y_array )
                        
            #### x-axis direction intersections
            
            # 2D array of DEM segment parameters                         
            x_dem_m = self.grad_forward_x()            
            x_dem_q = self.data - cell_center_x_array * x_dem_m            
            
            # equation for the planar surface that, given (x,y), will be used to derive z  
            plane_z_closure = srcPlaneAttitude.plane_from_geo( srcPt )
            
            # 2D array of plane segment parameters
            x_plane_m = srcPlaneAttitude.plane_x_coeff(  )            
            x_plane_q = array_from_function( self.row_num(), 1, lambda j: 0, grid_coord_to_geogr_coord_y_closure, plane_z_closure )

            # 2D array that defines denominator for intersections between local segments
            x_inters_denomin =  np.where( x_dem_m != x_plane_m, x_dem_m-x_plane_m, np.NaN )            
            coincident_x = np.where( x_dem_q != x_plane_q, np.NaN, ycoords_x )            
            xcoords_x = np.where( x_dem_m != x_plane_m , (x_plane_q - x_dem_q ) / x_inters_denomin, coincident_x )
            
            xcoords_x = np.where( xcoords_x < ycoords_x , np.NaN, xcoords_x )           
            xcoords_x = np.where( xcoords_x >= ycoords_x + self.cellsize_x() , np.NaN, xcoords_x )  
                        
            
            #### y-axis direction intersections

            # 2D array of DEM segment parameters  
            y_dem_m = self.grad_forward_y()            
            y_dem_q = self.data - cell_center_y_array * y_dem_m
 
            # 2D array of plane segment parameters
            y_plane_m = srcPlaneAttitude.plane_y_coeff(  )            
            y_plane_q = array_from_function( 1, self.col_num(), grid_coord_ to_geogr_coord_x , lambda i: 0, plane_z_closure )

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

