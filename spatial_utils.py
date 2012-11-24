# -*- coding: utf-8 -*-

from __future__  import division

from math import *
import numpy as np

import copy

from function_utils import *


from PyQt4.QtCore import *
from PyQt4.QtGui import *


MINIMUM_DISTANCE = 1e-10

        
def plane_x_coeff( or_plane_attitude ):
    """
    Calculate the slope of a given plane along the x direction.
    The plane orientation  is expressed following the geological convention. 
 
    @param  or_plane_attitude:  orientation of a plane.
    @type  or_plane_attitude:  StructPlane.
           
    @return:  slope - float.    
    """ 
    return - sin( radians( or_plane_attitude._dipdir ) ) * tan( radians( or_plane_attitude._dipangle ) )


def plane_y_coeff( or_plane_attitude ):
    """
    Calculate the slope of a given plane along the y direction.
    The plane orientation  is expressed following the geological convention. 

    @param  or_plane_attitude:  orientation of a plane.
    @type  or_plane_attitude:  StructPlane.
           
    @return:  slope - float.     
    """ 
    return - cos( radians( or_plane_attitude._dipdir ) ) * tan( radians( or_plane_attitude._dipangle ) )

        
def plane_from_geo( or_Pt, or_plane_attitude ):
    """
    Closure that embodies the analytical formula for the given plane.
    This closure is used to calculate the z value from given horizontal coordinates (x, y).

    @param  or_Pt:  Point instance expressing a location point contained by the plane.
    @type  or_Pt:  Point.    
    @param  or_plane_attitude:  orientation of a plane.
    @type  or_plane_attitude:  StructPlane.
    
    @return:  lambda (closure) expressing an analytical formula for deriving z given x and y values.    
    """

    x0 =  or_Pt.x     
    y0 =  or_Pt.y
    z0 =  or_Pt.z

    # slope of the line parallel to the x axis and contained by the plane
    a = plane_x_coeff( or_plane_attitude ) 
           
    # slope of the line parallel to the y axis and contained by the plane    
    b = plane_y_coeff( or_plane_attitude )
                    
    return lambda x, y : a * ( x - x0 )  +  b * ( y - y0 )  +  z0


def plane_from_pts( p1, p2, p3 ):
    """
    Closure.
    Creates plane equation given three points (p1-p3).
    """
    
    a_arr = np.array([[p1.y, p1.z, 1.0], [p2.y, p2.z, 1.0], [p3.y, p3.z, 1.0]])
    b_arr = np.array([[p1.x, p1.z, 1.0], [p2.x, p2.z, 1.0], [p3.x, p3.z, 1.0]])
    c_arr = np.array([[p1.x, p1.y, 1.0], [p2.x, p2.y, 1.0], [p3.x, p3.y, 1.0]])
    d_arr = np.array([[p1.x, p1.y, p1.z], [p2.x, p2.y, p2.z], [p3.x, p3.y, p3.z]]) 
            
    a = np.linalg.det(a_arr)
    b = np.linalg.det(b_arr)*(-1.0)
    c = np.linalg.det(c_arr)
    d = np.linalg.det(d_arr)*(-1.0)

    return lambda x, y, z : a * x  +  b * y  + c * z + d    


def point_in_plane( pt, plane, tolerance = 1.0e-1 ):
    """
    Determines if a point is within a plane.
    """

    if abs( plane (pt.x, pt.y, pt.z) ) < tolerance:
        return True
    return False


class Point(object):
    """
    3D points.
    
    """
        
    # class constructor
    def __init__( self, x = np.nan, y = np.nan, z = 0.0 ):
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
    
 
    def distance( self, pt ):
        """
        Calculate Euclidean distance between two points.

        @param  pt:  the Point instance for which the distance should be calculated
        @type  pt:  Point.
        
        @return:  distance between the two points - float.        
        """
        d_x = pt.x - self.x
        d_y = pt.y - self.y
        d_z = pt.z - self.z        
        
        return sqrt( d_x*d_x +  d_y*d_y +  d_z*d_z )


    def distance_uvector(self, end_pt):
        
        d_x = end_pt.x - self.x
        d_y = end_pt.y - self.y
        d_z = end_pt.z - self.z
                    
        length = sqrt( d_x*d_x +  d_y*d_y +  d_z*d_z )
        
        if length < MINIMUM_DISTANCE:
            return 0.0, Vector( 0.0, 0.0, 0.0 )
        else:
            return length, Vector( d_x/length, d_y/length, d_z/length )
        
    
    def displaced_by_deltas( self, sx = 0.0 , sy = 0.0 , sz = 0.0 ):
        """
        Create a new point shifted by a given amount from the self Point instance.

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
        

class Vector(Point):
    
    def __init__(self, x = np.nan, y = np.nan, z = 0.0):
        
        super(Vector, self).__init__(x, y, z)
        
    
    def length(self):
        
        return sqrt( self.x * self.x + self.y * self.y + self.z * self.z )
        
        
    def scale(self, scale_factor):
                
        return Vector(self.x*scale_factor, 
                      self.y*scale_factor, 
                      self.z*scale_factor)
  
  
    def to_unit_vector(self):
        
        scale_factor = self.lenght
        
        return self.scale(1.0/scale_factor)
    
        
    def add(self, another):
        
        return Vector(self.x + another.x, 
                      self.y + another.y, 
                      self.z + another.z)
        
        
class StructPlane(object):
    """
    Structural plane, following geological conventions.
    
    """
    
    def __init__( self, srcDipDir, srcDipAngle ):
        """
        Class constructor
        
        @param  srcDipDir:  Dip direction of the plane (0-360�).
        @type  srcDipDir:  number or string convertible to float.
        @param  srcDipAngle:  Dip angle of the plane (0-90�).
        @type  srcDipAngle:  number or string convertible to float.
           
        @return:  StructPlane.
    
        """
        
        self._dipdir = float( srcDipDir )
        self._dipangle = float( srcDipAngle )
  
   
  
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


class SpatialDomain(object):
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
                        
        @return:  SpatialDomain instance.
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

    """

    def __init__(self, grid_params = None, grid_data = None):
        """
        Grid class constructor.
        
        @param  grid_params:  the geo-parameters of the grid.
        @type  grid_params:  class GDALParameters.
        @param  grid_data:  the array storing the data.
        @type  grid_data:  2D np.array.
               
        @return:  self.
        """ 
                  
        if grid_params is not None:
            pt_llc = grid_params.llcorner()
            pt_trc = grid_params.trcorner()       
        else:
            pt_llc = None
            pt_trc = None
           
        self._grid_domain = SpatialDomain(pt_llc, pt_trc)          

        if grid_data is not None: 
            self._grid_data = grid_data.copy() 
        else:
            self._grid_data = None 
        
        
    def s_domain( self, domain ):
        """
        Set spatial domain.
        
        @param  domain:  Spatial domain to be attributed to the current Grid instance.
        @type  domain:  class SpatialDomain.
        
        @return: self
        """        
        del self._grid_domain
        self._grid_domain = copy.deepcopy(domain)
        
        
    def g_domain( self ):
        """
        Get spatial domain.
        
        @return: the spatial domain of the current Grid instance - class SpatialDomain.
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
        @type srcPlaneAttitude: class StructPlane.
        
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
            plane_z = plane_from_geo( srcPt, srcPlaneAttitude )
            
            # 2D array of plane segment parameters
            x_plane_m = plane_x_coeff( srcPlaneAttitude )            
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
            y_plane_m = plane_y_coeff( srcPlaneAttitude )            
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
        
        
    def calculate_profile_from_path( self, path, profile_resolution ):
        
        # result: list of [ Point, height, incremental_distance ]
        
        assert len(path) == 1
        
        profile = []
        
        single_path = path[0]
        
        cumulated_distance = 0.0

        for point_ndx in xrange( len( single_path) - 1 ):

            point_start = single_path[ point_ndx ]
            point_end = single_path[ point_ndx+1 ]
            distance, unit_vector = point_start.distance_uvector(point_end)
            
            if distance < MINIMUM_DISTANCE:
                continue                

            point_height = self.interpolate_bilinear( self.geog2array_coord( point_start ) )
            profile.append( [point_start, cumulated_distance, point_height ] )
            
            num_iter = 0
            potential_displacement = (num_iter+1)*profile_resolution
                        
            while potential_displacement < distance:

                displacement_vector = unit_vector.scale(potential_displacement)
                intermediate_point = point_start.displaced_by_vector(displacement_vector) 
                intermediate_cumulated_distance = cumulated_distance + potential_displacement
                point_height = self.interpolate_bilinear( self.geog2array_coord( intermediate_point ) )
                profile.append( [intermediate_point, intermediate_cumulated_distance, point_height ] )
                            
                num_iter += 1
                potential_displacement = (num_iter+1)*profile_resolution

            cumulated_distance += distance
            
        point_height = self.interpolate_bilinear( self.geog2array_coord( single_path[-1] ) )                                      
        profile.append( [single_path[-1], cumulated_distance, point_height ] )
    
        return profile
            
            
            
            
            
            
            
            