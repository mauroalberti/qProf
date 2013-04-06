# -*- coding: utf-8 -*-


import os
import math
import numpy as np

from osgeo import ogr

from PyQt4.QtCore import *
from PyQt4.QtGui import *

import qgis.core as QgsCore

from gis_utils import read_line_shapefile, read_raster_band, \
                      Vector_Input_Errors, Raster_Parameters_Errors
                      
from spatial_utils import *

from qProf_mplwidget import MplWidget

import unicodedata

        
class qProfDialog( QDialog ):
    """
    Constructor
    
    """
    
    def __init__( self ):

        super( qProfDialog, self ).__init__()        
        self.current_directory = os.path.dirname(__file__)
        self.initialize_parameters()                 
        self.setup_gui()        
           
    
    def get_layers_from_qgis(self):        
        
        # filter raster and vector layers
        curr_map_layers = QgsCore.QgsMapLayerRegistry.instance().mapLayers()
        mapLayers = zip(unicode(curr_map_layers.keys()), curr_map_layers.values())       
        rasterLayers = filter( lambda layer: layer[1].type() == QgsCore.QgsMapLayer.RasterLayer, mapLayers )
        vectorLayers = filter( lambda layer: layer[1].type() == QgsCore.QgsMapLayer.VectorLayer, mapLayers )
        linevectLayers = filter( lambda layer: layer[1].geometryType() == QgsCore.QGis.Line, vectorLayers )

        return rasterLayers, linevectLayers


    def initialize_parameters(self):
        
        self.rasterLayers, self.linevectLayers = self.get_layers_from_qgis()
        
        self.profile_windows = []
        
        self.current_profiles = None
        

    def setup_gui( self ): 

        self.dialog_layout = QVBoxLayout()

        self.main_widget = QTabWidget()
        
        self.setup_input_tab()
        self.setup_about_tab()
        
        self.dialog_layout.addWidget(self.main_widget)                             
        self.setLayout(self.dialog_layout)            
        self.adjustSize()               
        self.setWindowTitle( 'qProf' )        
                
  
    def setup_input_tab( self ):  

        inputWidget = QWidget() 
        inputLayout = QGridLayout()
 
        inputLayout.addWidget( QLabel( self.tr("Path shapefile:") ), 0, 0, 1, 1)        

        self.Trace_comboBox = QComboBox()
        inputLayout.addWidget(self.Trace_comboBox, 0, 1, 1, 3) 
        for ( name, layer ) in self.linevectLayers:
            self.Trace_comboBox.addItem( layer.name() )
                     
        inputLayout.addWidget(QLabel( "Use DEMs:" ), 1, 0, 1, 1)

        # input data tree view        
        self.listDEMs_treeWidget = QTreeWidget()
        self.listDEMs_treeWidget.setHorizontalScrollBarPolicy(Qt.ScrollBarAsNeeded)
        self.listDEMs_treeWidget.setDragEnabled(False)
        self.listDEMs_treeWidget.setDragDropMode(QAbstractItemView.NoDragDrop)
        self.listDEMs_treeWidget.setAlternatingRowColors(True)
        self.listDEMs_treeWidget.setSelectionMode(QAbstractItemView.SingleSelection)
        self.listDEMs_treeWidget.setTextElideMode(Qt.ElideNone)
        self.listDEMs_treeWidget.headerItem().setText( 0, self.tr( "Files" ) )
        self.listDEMs_treeWidget.header().setVisible(False)
                                     
        self.listDEMs_treeWidget.clear()                             
        for (name, layer) in self.rasterLayers:
            tree_item = QTreeWidgetItem( self.listDEMs_treeWidget, [ layer.name() ] )       
            tree_item.setFlags( tree_item.flags() | Qt.ItemIsUserCheckable )
            tree_item.setCheckState( 0, 0 )
            
        inputLayout.addWidget(self.listDEMs_treeWidget, 1, 1, 1, 3 )
            
        inputLayout.addWidget( QLabel( self.tr("Max sampling distance") ), 2, 0, 1, 1 ) 
        
        self.sample_distance_lineedit = QLineEdit()
        inputLayout.addWidget( self.sample_distance_lineedit, 2, 1, 1, 3 )
        
        self.CalcProf_pushbutton = QPushButton(self.tr("Calculate profile")) 
        QObject.connect(self.CalcProf_pushbutton, SIGNAL("clicked()"), self.calculate_profiles )               
        inputLayout.addWidget( self.CalcProf_pushbutton, 3, 0, 1, 4 )

        self.PlotProf_pushbutton = QPushButton(self.tr("Plot")) 
        QObject.connect(self.PlotProf_pushbutton, SIGNAL("clicked()"), self.plot_profiles )               
        inputLayout.addWidget( self.PlotProf_pushbutton, 4, 0, 1, 2 )        
     
        self.PlotHeight_checkbox = QCheckBox( self.tr( "height"))
        inputLayout.addWidget( self.PlotHeight_checkbox, 4, 2, 1, 1 )  
        self.PlotHeight_checkbox.setChecked( True ) 

        self.PlotSlope_checkbox = QCheckBox( self.tr( "slope"))
        inputLayout.addWidget( self.PlotSlope_checkbox, 4, 3, 1, 1 )   
        self.PlotSlope_checkbox.setChecked( True ) 
                               
        self.SaveData_pushbutton = QPushButton(self.tr("Save full results as")) 
        QObject.connect(self.SaveData_pushbutton, SIGNAL("clicked()"), self.save_full_data )               
        inputLayout.addWidget( self.SaveData_pushbutton, 5, 0, 1, 2 )        

        self.SaveCSV_checkbox = QCheckBox( self.tr( "csv"))
        inputLayout.addWidget( self.SaveCSV_checkbox, 5, 2, 1, 1 )  

        self.SavePtShp_checkbox = QCheckBox( self.tr( "2D point shp"))
        inputLayout.addWidget( self.SavePtShp_checkbox, 5, 3, 1, 1 ) 

        self.Export3DLine_pushbutton = QPushButton(self.tr("Export 3D line from DEM ")) 
        QObject.connect(self.Export3DLine_pushbutton, SIGNAL("clicked()"), self.export3dline )               
        inputLayout.addWidget( self.Export3DLine_pushbutton, 6, 0, 1, 2 )
         
        self.DEM_3D_Export_comboBox = QComboBox()
        inputLayout.addWidget(self.DEM_3D_Export_comboBox, 6, 2, 1, 2)
             
        inputWidget.setLayout(inputLayout)
                                
        self.main_widget.addTab( inputWidget, "Input layers" ) 
  
  
    def setup_about_tab(self):
        
        aboutWidget = QWidget()  
        aboutLayout = QVBoxLayout( )
        
        htmlText = """
        <h3>qProf release 0.1.3 (2013-04-02, experimental)</h3>
        Created by M. Alberti and M. Zanieri.
        <br />Concept: M. Zanieri, implementation: M. Alberti.
        <br />We thank S. Peduzzi for his vigorous testing.
        <br />Plugin for creating profiles, with slope calculation and 2D-3D output.
        <br /><br />For info see www.malg.eu.        
        """
        
        aboutQTextBrowser = QTextBrowser( aboutWidget )        
        aboutQTextBrowser.insertHtml( htmlText ) 
        aboutLayout.addWidget( aboutQTextBrowser )  
        aboutWidget.setLayout(aboutLayout)              
        self.main_widget.addTab( aboutWidget, "About" ) 
             
                          
    def calculate_profiles(self):

        # checking sample distance
        if self.sample_distance_lineedit.text() == '':
            QMessageBox.critical( self, "Error while reading sample distance", "Check max sampling distance value" )
            return            
        try:
            self.sample_distance = float( self.sample_distance_lineedit.text() )
        except:
            QMessageBox.critical( self, "Error while reading sample distance", "Check sample distance value" )
            return 
        
        # checking input path shapefile
        path_shape_qgis_ndx = self.Trace_comboBox.currentIndex()
        if path_shape_qgis_ndx < 0:
            QMessageBox.critical( self, "Error while reading path shapefile", "No selected shapefile" )             
            return       
        # trying to read input path shapefile
        try:
            line_shape_fpath = self.linevectLayers[ path_shape_qgis_ndx ][ 1 ].source()  
            lines_points, layer_extent_x, layer_extent_y = read_line_shapefile( line_shape_fpath )
        except ( IOError, TypeError, Vector_Input_Errors ), e:                    
            QMessageBox.critical( self, "Error while reading input shapefile", unicode(e) )
            return
                
        # checking if any selected DEMs
        selected_dems = False
        for dem_qgis_ndx in range( self.listDEMs_treeWidget.topLevelItemCount () ): 
            if self.listDEMs_treeWidget.topLevelItem ( dem_qgis_ndx ).checkState ( 0 ) == 2: 
                selected_dems = True
                break  
        if selected_dems == False:
            QMessageBox.critical( self, "Selected DEMs", "No DEM is selected" )
            return            
 
        profile_list = []                         
        for dem_qgis_ndx in range( self.listDEMs_treeWidget.topLevelItemCount () ): 
            if self.listDEMs_treeWidget.topLevelItem ( dem_qgis_ndx ).checkState ( 0 ) != 2: continue           
            dem_fpath = self.rasterLayers[ dem_qgis_ndx ][ 1 ].source()
            try:
                dem_params, dem_array = read_raster_band( dem_fpath )          
                dem_params.check_params()
            except ( IOError, TypeError, Raster_Parameters_Errors ), e:                    
                QMessageBox.critical( self, "Error while reading input DEM", unicode(e) )
                return
            
            dem = Grid( dem_params, dem_array ) 
            dem_name = self.listDEMs_treeWidget.topLevelItem ( dem_qgis_ndx ).text(0)   
            profile_result = dem.calculate_profile_from_2d_path( lines_points, self.sample_distance )
            
            if not profile_result[0]:
                QMessageBox.critical( self, "Error while reading path", unicode( profile_result[1] ) )
                return                

            profile_list.append( [ profile_result[1], dem_name ] )

        self.current_profiles = profile_list
        
        self.current_results = self.process_results()
        
        self.update_3d_dem_export_combobx()
        
        QMessageBox.information( self, "Profile calculated", "Now you can plot/save/export" )


    def process_results(self):
        
        # definition of output results
 
        profiles = [ profile for profile, _ in self.current_profiles ]        
                   
        x_list = [ vertex_3d.x for vertex_3d, _, _, _ in self.current_profiles[0][0] ]
        y_list = [ vertex_3d.y for vertex_3d, _, _, _ in self.current_profiles[0][0] ]
        cumul_2d_distance_list = [ cumulated_2D_distance for _, cumulated_2D_distance, _, _ in self.current_profiles[0][0] ]

        profiles_z_list = []
        profiles_cumul3d_dist_list = []
        profiles_slope_list = []
        for profile in profiles:
            profile_z = [ vertex_3d.z for vertex_3d, _, _, _ in profile ]
            profile_cumul3d_dist = [ cumulated_3D_distance for _, _, cumulated_3D_distance, _ in profile ]
            profile_slopes = [ slope_degr for _, _, _, slope_degr in profile ]
            profiles_z_list.append( profile_z ) 
            profiles_cumul3d_dist_list.append( profile_cumul3d_dist ) 
            profiles_slope_list.append( profile_slopes )          

        zs_list = zip( *profiles_z_list )
        cum3d_dists_list = zip( *profiles_cumul3d_dist_list ) 
        slopes_list = zip( *profiles_slope_list )  

        result_data = []
        for x, y, cum_2d_dist, zs, cum3d_dists, slopes in zip( x_list, y_list, cumul_2d_distance_list, zs_list, cum3d_dists_list, slopes_list ):
            record = [ x, y, cum_2d_dist]
            for z, cum3d_dist, slope in zip( zs, cum3d_dists, slopes ):
                if math.isnan(z): z = ''
                if math.isnan(cum3d_dist): cum3d_dist = ''
                if math.isnan(slope): slope = ''
                record += [z, cum3d_dist, slope ]
            result_data.append( record )
         
        return result_data
     
            
    def update_3d_dem_export_combobx(self):
        
        self.DEM_3D_Export_comboBox.clear()        
        self.DEM_3D_Export_comboBox.addItems( [dem_name for profile, dem_name in self.current_profiles])
        
                    
    def plot_profiles( self ):
 
        if self.current_profiles is None:            
            QMessageBox.critical( self, "Plotting profile", "No profile has still been calculated" )
            return 
 
        if not ( self.PlotHeight_checkbox.isChecked() or self.PlotSlope_checkbox.isChecked() ):
            QMessageBox.critical( self, "Plotting profile", "Neither height or slope options are selected" )
            return         
                   
        profile_list = self.current_profiles

        # defines the input values
        dem_names_list = []; z_values_list = []; slope_values_list = []      
        for profile, dem_name in profile_list:
            dem_names_list.append( dem_name )
            z_values_list.append( [ profile_point_3d.z for profile_point_3d, _, _, _ in profile ] )
            slope_values_list.append( [ slope for _, _, _, slope in profile ] )

        # defines the extent for the plot window: s min and max 
        s_2d_values = [ cum_dist_2d for _, cum_dist_2d, _, _ in profile_list[0][0] ]      
        profiles_s_min, profiles_s_max = 0, s_2d_values[-1]        

        # defines z min and max values
        profiles_z_list = [ z for z_values in z_values_list for z in z_values if not math.isnan( z ) ]
        profiles_z_min, profiles_z_max = min( profiles_z_list ), max( profiles_z_list )
        delta_z = profiles_z_max - profiles_z_min 
        profiles_z_min, profiles_z_max = profiles_z_min - delta_z * 0.05, profiles_z_max + delta_z * 0.05

        # defines slope min and max values
        profiles_slope_list = [ slope for slope_values in slope_values_list for slope in slope_values if not math.isnan( slope ) ]
        profiles_slope_min, profiles_slope_max = min( profiles_slope_list ), max( profiles_slope_list )
        delta_slope = profiles_slope_max - profiles_slope_min 
        profiles_slope_min, profiles_slope_max = profiles_slope_min - delta_slope*0.2, profiles_slope_max + delta_slope*0.2  

        # map
        profile_window = MplWidget()  

        if self.PlotHeight_checkbox.isChecked() and self.PlotSlope_checkbox.isChecked():
            mpl_code_list = [ 211, 212 ]
        else:
            mpl_code_list = [ 111 ]            

        # for mpl_code in mpl_code_list:
        colors = ['r', 'b', 'g', 'y', 'c', 'k', 'm', 'orange']  

        s_2d_values_array = np.array( s_2d_values )
        
        subplot_code = mpl_code_list[0]   
        if self.PlotHeight_checkbox.isChecked():                            
            axes_height = profile_window.canvas.fig.add_subplot( subplot_code )
            axes_height.set_xlim( profiles_s_min, profiles_s_max )
            axes_height.set_ylim( profiles_z_min, profiles_z_max ) 
            axes_height.set_color_cycle( colors )
            for dem_name, z_values, color in zip( dem_names_list, z_values_list, colors ):              
                z_values_array = np.array( z_values )
                valid_s_2d_values = s_2d_values_array[ np.isfinite( z_values_array ) ]
                valid_z_values = z_values_array[ np.isfinite( z_values_array ) ]                
                axes_height.fill_between( valid_s_2d_values, profiles_z_min, valid_z_values, facecolor=color, alpha=0.1 )                       
                axes_height.plot(s_2d_values, z_values,'-', label = unicode(dem_name) )
                
            axes_height.grid(True)
            axes_height.legend(loc = 'upper left', shadow=True)              
            
        if self.PlotSlope_checkbox.isChecked():            
            if len(mpl_code_list) == 2: subplot_code = mpl_code_list[1]
            axes_slope = profile_window.canvas.fig.add_subplot( subplot_code )
            axes_slope.set_xlim( profiles_s_min, profiles_s_max )
            axes_slope.set_ylim( profiles_slope_min, profiles_slope_max )
            axes_slope.set_color_cycle( colors )
            for dem_name, slope_values, color in zip( dem_names_list, slope_values_list, colors):
                
                slope_values_array = np.array( slope_values )
                valid_s_2d_values = s_2d_values_array[ np.isfinite( slope_values_array ) ]
                valid_slope_values = slope_values_array[ np.isfinite( slope_values_array ) ]  
                
                axes_slope.fill_between( valid_s_2d_values, 0, valid_slope_values, facecolor=color, alpha=0.1 )
                axes_slope.plot(s_2d_values, slope_values,'-', label = unicode(dem_name) )
                
            axes_slope.grid(True)
            axes_slope.legend(loc = 'upper left', shadow=True)              
            
                    
        profile_window.canvas.draw() 
        
        self.profile_windows.append(profile_window) 
          
    
    def save_full_data(self):
        
        if self.current_profiles is None:            
            QMessageBox.critical( self, "Saving results", "No profile has still been calculated" )
            return 
        
        if not ( self.SaveCSV_checkbox.isChecked() or self.SavePtShp_checkbox.isChecked() ):
            QMessageBox.critical( self, "Saving results", "No output format is selected" )
            return   
                          
        # definition of field names        
        dem_names = [ dem_name for _, dem_name in self.current_profiles ]  
        cum3ddist_headers = ['cumulated_3d_distance']*len(dem_names)
        slopes_headers = ['slopes (degr)']*len(dem_names)
        header_list = ['x', 'y', 'cumulated_2d_distance'] + [ name for sublist in zip(dem_names, cum3ddist_headers, slopes_headers) for name in sublist ]              

        # output for csv file
        if self.SaveCSV_checkbox.isChecked():            
            self.write_csv( header_list )
        if self.SavePtShp_checkbox.isChecked():
            self.write_ptshp( dem_names )
        
        QMessageBox.information( self, "Saving profile results", "Results saved" )
        
                      
    def write_ptshp( self, dem_names ):
        
        fileName = QFileDialog.getSaveFileName(self, 
                                               self.tr("Save results as 2D point shapefile"),
                                                "*.shp",
                                                self.tr("shapefile (*.shp)"))

        if fileName is None or fileName == '':
            QMessageBox.critical( self, "Saving results", "No file has been defined" )
            return           

        shape_driver_name = "ESRI Shapefile"
        shape_driver = ogr.GetDriverByName( shape_driver_name )
        if shape_driver is None:
            QMessageBox.critical( self, "Saving results", "%s driver is not available" % shape_driver_name )
            return             
            
        ptshp_datasource = shape_driver.CreateDataSource( unicode( fileName ) )
        if ptshp_datasource is None:
            QMessageBox.critical( self, "Saving results", "Creation of %s shapefile failed" % os.path.split( fileName )[1] )
            return         
        
        ptshp_layer = ptshp_datasource.CreateLayer( 'profile', geom_type=ogr.wkbPoint )
        if ptshp_layer is None:
            QMessageBox.critical( self, "Saving results", "Layer creation failed" )
            return 
                  
        # creates required fields         
        ptshp_layer.CreateField( ogr.FieldDefn( 'id', ogr.OFTInteger ) )
        ptshp_layer.CreateField( ogr.FieldDefn( 'x', ogr.OFTReal ) )
        ptshp_layer.CreateField( ogr.FieldDefn( 'y', ogr.OFTReal ) )       
        ptshp_layer.CreateField( ogr.FieldDefn( 'cum2dis', ogr.OFTReal ) )
        for dem_ndx, dem_name in enumerate( dem_names ):       
            ptshp_layer.CreateField( ogr.FieldDefn( unicodedata.normalize('NFKD', unicode( dem_name[:10] )).encode('ascii', 'ignore'), ogr.OFTReal ) )       
            ptshp_layer.CreateField( ogr.FieldDefn( 'cum3dis'+str(dem_ndx+1), ogr.OFTReal ) )        
            ptshp_layer.CreateField( ogr.FieldDefn( 'slope'+str(dem_ndx+1), ogr.OFTReal ) )        
        
        ptshp_featureDefn = ptshp_layer.GetLayerDefn()
        
        # loops through output records                       
        for ndx, rec in enumerate( self.current_results ):
            
            pt_feature = ogr.Feature( ptshp_featureDefn )
            
            pt = ogr.Geometry(ogr.wkbPoint)
            pt.SetPoint_2D( 0, rec[0], rec[1] )        
            pt_feature.SetGeometry( pt )
            
            pt_feature.SetField('id', ndx+1 )
            pt_feature.SetField('x', rec[0] )   
            pt_feature.SetField('y', rec[1] ) 
            pt_feature.SetField('cum2dis', rec[2] )  
            for dem_ndx, dem_name in enumerate( dem_names ):
                dem_height = rec[3+dem_ndx*3+0]
                if dem_height != '': pt_feature.SetField( unicodedata.normalize('NFKD', unicode( dem_name[:10] )).encode('ascii', 'ignore'), dem_height )
                cum3ddist = rec[3+dem_ndx*3+1]
                if cum3ddist != '': pt_feature.SetField( 'cum3dis'+str(dem_ndx+1), cum3ddist )
                slope = rec[3+dem_ndx*3+2]
                if slope != '': pt_feature.SetField( 'slope'+str(dem_ndx+1), slope )                  
  
            ptshp_layer.CreateFeature( pt_feature )
            
            pt_feature.Destroy()
            
        ptshp_datasource.Destroy()
                                
        
    def write_csv( self, header_list ):
               
        fileName = QFileDialog.getSaveFileName(self, 
                                               self.tr("Save results"),
                                                "*.csv",
                                                self.tr("csv (*.csv)"))

        if fileName is None or fileName == '':
            QMessageBox.critical( self, "Saving results", "No file has been defined" )
            return   

        header_list = [ unicodedata.normalize('NFKD', unicode(header)).encode('ascii', 'ignore') for header in header_list]
        with open( unicode( fileName ), 'w') as f:
            f.write( ','.join( header_list )+'\n' )
            for rec in self.current_results:
                out_rec_string = ''
                for val in rec:
                    out_rec_string += str( val ) + ','
                f.write( out_rec_string[:-1]+'\n' )
                
        
    def export3dline(self):
           
        fileName = QFileDialog.getSaveFileName(self, 
                                               self.tr("Save results as line shapefile"),
                                                "*.shp",
                                                self.tr("shapefile (*.shp)"))

        if fileName is None or fileName == '':
            QMessageBox.critical( self, "Saving results", "No file has been defined" )
            return       
    
        shape_driver_name = "ESRI Shapefile"
        shape_driver = ogr.GetDriverByName( shape_driver_name )
        if shape_driver is None:
            QMessageBox.critical( self, "Saving results", "%s driver is not available" % shape_driver_name )
            return             
            
        lnshp_datasource = shape_driver.CreateDataSource( unicode( fileName ) )
        if lnshp_datasource is None:
            QMessageBox.critical( self, "Saving results", "Creation of %s shapefile failed" % os.path.split( fileName )[1] )
            return         
        
        lnshp_layer = lnshp_datasource.CreateLayer( 'profile', geom_type=ogr.wkbLineString25D )
        if lnshp_layer is None:
            QMessageBox.critical( self, "Saving results", "Layer creation failed" )
            return     

        current_dem_ndx = self.DEM_3D_Export_comboBox.currentIndex()
            
        # creates required fields         
        lnshp_layer.CreateField( ogr.FieldDefn( 'id', ogr.OFTInteger ) )      
        lnshp_layer.CreateField( ogr.FieldDefn( 'cum2dis', ogr.OFTReal ) )
        lnshp_layer.CreateField( ogr.FieldDefn( 'cum3dis', ogr.OFTReal ) )
        lnshp_layer.CreateField( ogr.FieldDefn( 'slope', ogr.OFTReal ) )     
    
        lnshp_featureDefn = lnshp_layer.GetLayerDefn()
        
        # loops through output records                     
        for ndx in range( len( self.current_results )-1 ):
                
            rec_a = self.current_results[ndx]
            rec_b = self.current_results[ndx+1]
                            
            x0, y0, z0 = rec_a[0], rec_a[1], rec_a[3+current_dem_ndx*3]
            x1, y1, z1 = rec_b[0], rec_b[1], rec_b[3+current_dem_ndx*3]
            
            if z0 == '' or z1 == '': continue

            ln_feature = ogr.Feature( lnshp_featureDefn )            
            segment_3d = ogr.CreateGeometryFromWkt('LINESTRING(%f %f %f, %f %f %f)' % (x0, y0, z0, x1, y1, z1))       
            ln_feature.SetGeometry( segment_3d ) 
            
            ln_feature.SetField('id', ndx+1 )
            ln_feature.SetField('cum2dis', rec_b[2] )  
            cum3ddist = rec_b[3+current_dem_ndx*3+1]
            if cum3ddist != '': ln_feature.SetField( 'cum3dis', cum3ddist )
            slope = rec_b[3+current_dem_ndx*3+2]
            if slope != '': ln_feature.SetField( 'slope', slope )  
                            
            lnshp_layer.CreateFeature( ln_feature )
            
            ln_feature.Destroy()
            
        lnshp_datasource.Destroy()
        
        QMessageBox.information( self, "Saving 3d line", "Shapefile saved" )
 
