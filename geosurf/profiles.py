
from .spatial import Line_2D


class Profiles( object ):
    
    def __init__(self, max_spacing ):

        self.max_spacing = max_spacing # max spacing along profile; float
        self.init_elements()
        

    def init_elements(self):
        
        self.profile_dems = []
        self.profile_line = None                
        self.topo_profiles = []
        self.plane_attitudes = []
        self.curves = []
        self.curves_ids = []
        self.topo_pts = []
        self.topo_lines = []
        
        
    def add_topo_profile(self, topo_profile ):
        
        self.topo_profiles.append( topo_profile )
 

    def get_current_dem_names( self ):
        
        return [ topo_profile.dem_name for topo_profile in self.topo_profiles ]
    

    def get_min_s( self ):
        
        return min( [  topo_profile.length_2d() for topo_profile in self.topo_profiles ] )


    def get_max_s( self ):
        
        return max( [  topo_profile.length_2d() for topo_profile in self.topo_profiles ] )


    def min_z_topo_profiles( self ):
        
        return min( [  topo_profile.get_min_z() for topo_profile in self.topo_profiles ] )
        

    def max_z_topo_profiles( self ):
        
        return max( [  topo_profile.get_max_z() for topo_profile in self.topo_profiles ] )


    def min_z_plane_attitudes( self ):
           
        # TODO:  manage case for possible nan z values
        return min( [ plane_attitude.pt_3d._z for plane_attitude_set in self.plane_attitudes for plane_attitude in plane_attitude_set if  0.0 <= plane_attitude.sign_hor_dist <= self.get_max_s() ] )
        

    def max_z_plane_attitudes( self ):
            
        # TODO:  manage case for possible nan z values
        return max( [ plane_attitude.pt_3d._z for plane_attitude_set in self.plane_attitudes for plane_attitude in plane_attitude_set if  0.0 <= plane_attitude.sign_hor_dist <= self.get_max_s() ] )


    def min_z_curves( self ):
        
        return min( [ pt_2d._y for multiline_2d_list in self.curves for multiline_2d in multiline_2d_list for line_2d in multiline_2d._lines for pt_2d in line_2d._pts if 0.0 <= pt_2d._x <= self.get_max_s() ] )
                     

    def max_z_curves( self ):
        
        return max( [ pt_2d._y for multiline_2d_list in self.curves for multiline_2d in multiline_2d_list for line_2d in multiline_2d._lines for pt_2d in line_2d._pts if 0.0 <= pt_2d._x <= self.get_max_s() ] )
        
                                        
    def get_min_z( self ):
        
        min_z = self.min_z_topo_profiles()
        
        if len( self.plane_attitudes ) > 0:            
            min_z = min( [ min_z, self.min_z_plane_attitudes() ])
            
        if len( self.curves ) > 0:
            min_z = min( [ min_z, self.min_z_curves() ])
            
        return min_z


    def get_max_z( self ):
        
        max_z = self.max_z_topo_profiles()

        if len( self.plane_attitudes ) > 0:            
            max_z = max( [max_z, self.max_z_plane_attitudes() ])

        if len( self.curves ) > 0:
            max_z = max( [ max_z, self.max_z_curves() ])
                        
        return max_z
    
                    
    def add_plane_attitudes(self, plane_attitudes ):
        
        self.plane_attitudes.append( plane_attitudes )


    def add_curves(self, multiline_2d_list, ids_list ):
        
        self.curves.append( multiline_2d_list )
        self.curves_ids.append( ids_list )
 

    def add_topo_pts(self, topo_pts ):
        
        self.topo_pts.append( topo_pts )


    def add_topo_lines(self, topo_lines ):
        
        self.topo_lines.append( topo_lines )
   

class ProfileDEM( object ):
    
    def __init__( self, layer, params ):

        self.layer = layer
        self.params = params
        
                     
class TopoProfile( object ):
    
    def __init__( self, dem_name, line_3d ):
        
        self.dem_name = dem_name     
        self.profile_3d = line_3d
        

    def get_min_z( self ):
        
        return self.profile_3d.z_min()
    
    
    def get_max_z( self ):
        
        return self.profile_3d.z_max()
        

    def x_list( self ):
        
        return [ pt_3d._x for pt_3d in self.profile_3d._pts ]
    
    
    def y_list( self ):
        
        return [ pt_3d._y for pt_3d in self.profile_3d._pts ]
    
    
    def z_list( self ):
        
        return [ pt_3d._z for pt_3d in self.profile_3d._pts ]
    
    
    def slope_list( self ):
        
        return self.profile_3d.slopes_list()
    
            
    def length_2d( self ):
        
        return self.profile_3d.length_2d()
        
                   
    def get_increm_dist_3d( self ):
        
        return self.profile_3d.incremental_length()
        
       
    def get_increm_dist_2d( self ):
        
        return self.profile_3d.incremental_length_2d()

                   
class PlaneAttitude( object ):
    
    def __init__(self, rec_id, source_point_3d, source_geol_plane, point_3d, slope_rad, dwnwrd_sense, sign_hor_dist ):
        
        self.id = rec_id
        self.src_pt_3d = source_point_3d
        self.src_geol_plane = source_geol_plane
        self.pt_3d = point_3d
        self.slope_rad = slope_rad
        self.dwnwrd_sense = dwnwrd_sense
        self.sign_hor_dist = sign_hor_dist


class TopoPoints( object ):
    
    def __init__(self, topo_pt_list, attrs ):
        
        self.topo_points = topo_pt_list
        self.attrs = attrs


class TopoLines( object ):
    
    def __init__( self, line_2d_list, cats ):
        
        self.topo_lines_2d = line_2d_list
        self.cats = cats
        
        
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
