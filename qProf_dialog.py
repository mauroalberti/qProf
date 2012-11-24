# -*- coding: utf-8 -*-


import os
import sys

from PyQt4.QtCore import *
from PyQt4.QtGui import *

from matplotlib.offsetbox import AnchoredOffsetbox, TextArea 

import qgis.core as QgsCore

import guiqt_utils as guiqt
from gis_utils import read_line_shapefile, read_raster_band, \
                      Vector_Input_Errors, Raster_Parameters_Errors
                      
from spatial_utils import *

from qProf_mplwidget import MplWidget

        
class qProfDialog( QDialog ):
    """
    Constructor
    
    """

    # static attributes
    colormaps = ["jet", "gray", "bone", "hot","autumn", "cool","copper", "hsv", "pink", "spring", "summer", "winter", "spectral", "flag", ] 
    trace_colors = [ "white", "red", "blue", "yellow", "orange", "brown",]
    original_extent = [0,100]
    
    def __init__( self ):

        super( qProfDialog, self ).__init__()
        
        self.current_directory = os.path.dirname(__file__)

        self.initialize_parameters() 
                
        self.setup_gui()
        
        self.set_layers_from_qgis()
        
        #self.setConnections() 


    def initialize_parameters(self):

        # dem and map extent settings      
        self.dem_extent_x = self.dem_extent_y = qProfDialog.original_extent
        self.map_extent_x = self.map_extent_y = qProfDialog.original_extent
        
        # colormaps
        self.dem_colormap = qProfDialog.colormaps[0]
        self.intersection_color = qProfDialog.trace_colors[0]        
        
        # initialize intersection drawing
        self.valid_intersections = False
             

    def setup_gui( self ):
        
        self.setWindowTitle( 'qProf' )

        self.button_size = (30,30)
        
        self.icon_directory = os.path.join( self.current_directory, 'icons' )
                
        self.openfolder_icon_path = os.path.join( self.icon_directory, 'folder_orange_open.png' )
        #self.calcprof_icon_path = os.path.join( self.icon_directory, 'folder_orange_open.png' )
               
        # set layout for dialog
        dialog_splitter = QSplitter( Qt.Horizontal, self )
        
        # define left widget and layout
        left_widget = QWidget()
        left_layout = QVBoxLayout()
        left_widget.setLayout( left_layout )

        # define right widget and layout        
        right_widget = QWidget()
        right_layout = QVBoxLayout()
        right_widget.setLayout( right_layout )
        
        # setup input data section
        self.setup_input_section( left_layout )

        # setup look section
        self.setup_look_section( left_layout )      
        
        # setup profile parameters
        self.setup_profileinfo_section( left_layout )

        # setup output settings
        self.setup_output_section( left_layout )
        
        # map
        self.mplwidget = MplWidget( self.map_extent_x, self.map_extent_y )
        right_layout.addWidget( self.mplwidget )

        # add left and right widgets to dialog layout
        dialog_splitter.addWidget( left_widget )
        dialog_splitter.addWidget( right_widget )        
        
        self.resize( 1100, 700 )


    def setup_input_section( self, layout ):        

                
        inputGroupBox = QGroupBox( "Input layers" )       
        inputLayout = QGridLayout()
        inputGroupBox.setLayout(inputLayout) 
 
        inputLayout.addWidget( QLabel( self.tr("Path") ), 0, 0, 1, 1)        

        self.Trace_comboBox = QComboBox()
        inputLayout.addWidget(self.Trace_comboBox, 0, 1, 1, 2) 

        self.TraceBrowse_pushbutton = guiqt.push_button_with_icon(self.openfolder_icon_path, self.button_size, self.find_line_layer )       
        inputLayout.addWidget( self.TraceBrowse_pushbutton, 0, 3, 1, 1 )
                       
        inputLayout.addWidget(QLabel( "DEM" ), 1, 0, 1, 1)
        
        self.DEM_comboBox = QComboBox()
        inputLayout.addWidget(self.DEM_comboBox, 1, 1, 1, 2)

        self.DEMBrowse_pushbutton = guiqt.push_button_with_icon(self.openfolder_icon_path, self.button_size, self.find_DEM_layer )        
        inputLayout.addWidget( self.DEMBrowse_pushbutton, 1, 3, 1, 1)
        
        inputLayout.addWidget( QLabel( self.tr("Profile resolution") ), 2, 0, 1, 1) 
        
        self.profile_resolution_lineedit = QLineEdit()
        inputLayout.addWidget( self.profile_resolution_lineedit, 2, 1, 1, 1)        
        
        self.CalcProf_pushbutton = QPushButton(self.tr("Get profile")) 
        QObject.connect(self.CalcProf_pushbutton, SIGNAL("clicked()"), self.calculate_profile)               
        inputLayout.addWidget( self.CalcProf_pushbutton, 2, 2, 1, 2)
                        
        layout.addWidget( inputGroupBox )      


            
    def set_layers_from_qgis(self):        
        
        # filter raster and vector layers
        curr_map_layers = QgsCore.QgsMapLayerRegistry.instance().mapLayers()
        self.mapLayers = zip(unicode(curr_map_layers.keys()), curr_map_layers.values())       
        self.rasterLayers = filter( lambda layer: layer[1].type() == QgsCore.QgsMapLayer.RasterLayer, self.mapLayers )
        self.vectorLayers = filter( lambda layer: layer[1].type() == QgsCore.QgsMapLayer.VectorLayer, self.mapLayers )
        self.linevectLayers = filter( lambda layer: layer[1].geometryType() == QgsCore.QGis.Line, self.vectorLayers )
        
        for (name,layer) in self.rasterLayers:
            self.DEM_comboBox.addItem(layer.name())
            
        for (name,layer) in self.linevectLayers:
            self.Trace_comboBox.addItem(layer.name())


    def find_DEM_layer(self):
        
        pass
    

    def find_line_layer(self):
        
        pass
    
        
    def calculate_profile(self):

        try:
            self.profile_resolution = float( self.profile_resolution_lineedit.text() )
        except:
            QMessageBox.critical( self, "Error while reading profile resolution", "Found non-numeric value" )
            return            
        
        # reading input line shapefile
        path_shape_qgis_ndx = self.Trace_comboBox.currentIndex()
        if path_shape_qgis_ndx < 0: 
            return

        try:
            line_shape_fpath = self.linevectLayers[ path_shape_qgis_ndx ][ 1 ].source()  
            lines_points, layer_extent_x, layer_extent_y = read_line_shapefile( line_shape_fpath )
        except ( IOError, TypeError, Vector_Input_Errors ), e:                    
            QMessageBox.critical( self, "Error while reading input shapefile", str(e) )
            return
           
        # reading input DEM            
        dem_qgis_ndx = self.DEM_comboBox.currentIndex()
        if dem_qgis_ndx < 0: 
            return

        dem_fpath = self.rasterLayers[ dem_qgis_ndx ][ 1 ].source()
        try:
            dem_params, dem_array = read_raster_band( dem_fpath )
            dem_params.check_params()
        except ( IOError, TypeError, Raster_Parameters_Errors ), e:                    
            QMessageBox.critical( self, "Error while reading input DEM", str(e) )
            return
        
        dem = Grid( dem_params, dem_array )

        self.profile = dem.calculate_profile_from_path( lines_points, self.profile_resolution )
        
        self.plot_profile()
        
    
    
    def get_profile_range(self):
        
        x_min_init = self.profile[0][1]
        x_max_init = self.profile[0][1]
        y_min_init = self.profile[0][2]
        y_max_init = self.profile[0][2]
        
        for profile_point in self.profile[1:]:
            if profile_point[1] < x_min_init:
                x_min_init = profile_point[1]
            if profile_point[1] > x_max_init:
                x_max_init = profile_point[1]    
            if profile_point[2] < y_min_init:
                y_min_init = profile_point[2]
            if profile_point[2] > y_max_init:
                y_max_init = profile_point[2]
                
        return x_min_init,x_max_init, y_min_init,y_max_init


    def extract_profile_values(self):
        
        x_values = []
        y_values = []
        for profile_point in self.profile:
            
            x_values.append(profile_point[1])
            y_values.append(profile_point[2])
            
        return x_values, y_values            
      
        
    def plot_profile(self):

        self.mplwidget.canvas.ax.cla()
        
        x_values, y_values = self.extract_profile_values()
        
        x_min_init,x_max_init, y_min_init,y_max_init = self.get_profile_range()
        
        delta_y = y_max_init - y_min_init
        
        y_min_init = y_min_init - delta_y*0.05
        y_max_init = y_max_init + delta_y*0.05        

        self.mplwidget.canvas.ax.set_xlim(x_min_init,x_max_init)
        self.mplwidget.canvas.ax.set_ylim(y_min_init,y_max_init) 
                
        self.mplwidget.canvas.ax.plot(x_values, y_values,'-')
            
        self.mplwidget.canvas.draw()                                    
            
            
    def setup_look_section( self, layout ): 

        lookGroupBox = QGroupBox( "Rendering" )       
        lookLayout = QGridLayout()
        lookGroupBox.setLayout(lookLayout) 

        lookLayout.addWidget( QLabel( self.tr("X axis characteristics") ), 0, 1) 
        self.vest_axis_x_icon_path = os.path.join( self.icon_directory, 'vest_axis_x.png' )                 
        self.vest_axis_x_pushbutton = guiqt.push_button_with_icon(self.vest_axis_x_icon_path, self.button_size, self.do_vest_axis_x )       
        lookLayout.addWidget( self.vest_axis_x_pushbutton, 0, 0 )

        lookLayout.addWidget( QLabel( self.tr("Y axis characteristics") ), 1, 1) 
        self.vest_axis_y_icon_path = os.path.join( self.icon_directory, 'vest_axis_y.png' )                 
        self.vest_axis_y_pushbutton = guiqt.push_button_with_icon(self.vest_axis_y_icon_path, self.button_size, self.do_vest_axis_y )       
        lookLayout.addWidget( self.vest_axis_y_pushbutton, 1, 0 )        

        lookLayout.addWidget( QLabel( self.tr("Vertical exageration") ), 2, 1) 
        self.vest_vertex_icon_path = os.path.join( self.icon_directory, 'vest_exag_v.png' )                 
        self.vest_vertex_pushbutton = guiqt.push_button_with_icon(self.vest_vertex_icon_path, self.button_size, self.do_vert_exag )       
        lookLayout.addWidget( self.vest_vertex_pushbutton, 2, 0 )

        lookLayout.addWidget( QLabel( self.tr("Text") ), 3, 1) 
        self.vest_font_icon_path = os.path.join( self.icon_directory, 'vest_font.png' )                 
        self.vest_font_pushbutton = guiqt.push_button_with_icon(self.vest_font_icon_path, self.button_size, self.do_font )       
        lookLayout.addWidget( self.vest_font_pushbutton, 3, 0 )

        lookLayout.addWidget( QLabel( self.tr("Measure unit") ), 4, 1) 
        self.vest_measure_icon_path = os.path.join( self.icon_directory, 'vest_measure.png' )                 
        self.vest_measure_pushbutton = guiqt.push_button_with_icon(self.vest_measure_icon_path, self.button_size, self.do_mesure_unit )       
        lookLayout.addWidget( self.vest_measure_pushbutton, 4, 0 )        

        lookLayout.addWidget( QLabel( self.tr("Line style") ), 0, 3) 
        self.vest_linestyle_icon_path = os.path.join( self.icon_directory, 'vest_line_style.png' )                 
        self.vest_linestyle_pushbutton = guiqt.push_button_with_icon(self.vest_linestyle_icon_path, self.button_size, self.do_line_style )       
        lookLayout.addWidget( self.vest_linestyle_pushbutton, 0, 2 )                

        lookLayout.addWidget( QLabel( self.tr("Axis style") ), 1, 3) 
        self.vest_axisstyle_icon_path = os.path.join( self.icon_directory, 'vest_axis_style.png' )                 
        self.vest_axisstyle_pushbutton = guiqt.push_button_with_icon(self.vest_axisstyle_icon_path, self.button_size, self.do_axis_style )       
        lookLayout.addWidget( self.vest_axisstyle_pushbutton, 1, 2 )
        
        lookLayout.addWidget( QLabel( self.tr("Retino complessivo") ), 2, 3) 
        self.vest_retincompl_icon_path = os.path.join( self.icon_directory, 'vest_retin_compl.png' )                 
        self.vest_retincompl_pushbutton = guiqt.push_button_with_icon(self.vest_retincompl_icon_path, self.button_size, self.do_retin_compl )       
        lookLayout.addWidget( self.vest_retincompl_pushbutton, 2, 2 )
        
        lookLayout.addWidget( QLabel( self.tr("Retino classificato") ), 3, 3) 
        self.vest_retinclassif_icon_path = os.path.join( self.icon_directory, 'vest_retin_class.png' )                 
        self.vest_retinclassif_pushbutton = guiqt.push_button_with_icon(self.vest_retinclassif_icon_path, self.button_size, self.do_retin_classif )       
        lookLayout.addWidget( self.vest_retinclassif_pushbutton, 3, 2 )
        
        layout.addWidget( lookGroupBox )
                
        
    def do_vest_axis_x(self):
        
        pass
        
        
    def do_vest_axis_y(self):
        
        pass
    
    
    def do_vert_exag(self):
        
        pass
                
 
    def do_font(self):
        
        pass
    

    def do_mesure_unit(self):
        
        pass
    
 
    def do_line_style(self):
        
        pass
    

    def do_axis_style(self):
        
        pass
    
 
    def do_retin_compl(self):
        
        pass
    
 
    def do_retin_classif(self):
        
        pass
    

    def setup_profileinfo_section(self, layout ):
        
        profinfoGroupBox = QGroupBox( "Show profile parameters" )       
        profinfoLayout = QGridLayout()
        profinfoGroupBox.setLayout(profinfoLayout) 

        profinfoLayout.addWidget( QLabel( self.tr("Slope") ), 0, 1) 

        self.prof_show_slope_min = QCheckBox("min")        
        profinfoLayout.addWidget( self.prof_show_slope_min, 1, 1) 
        self.prof_show_slope_mean = QCheckBox("mean")        
        profinfoLayout.addWidget( self.prof_show_slope_mean, 2, 1)
        self.prof_show_slope_max = QCheckBox("max")        
        profinfoLayout.addWidget( self.prof_show_slope_max, 3, 1)

        profinfoLayout.addWidget( QLabel( self.tr("Elevation") ), 0, 0) 

        self.prof_show_elev_start = QCheckBox("start")        
        profinfoLayout.addWidget( self.prof_show_elev_start, 1, 0) 
        self.prof_show_elev_stop = QCheckBox("stop")        
        profinfoLayout.addWidget( self.prof_show_elev_stop, 2, 0)
        self.prof_show_elev_min = QCheckBox("min")        
        profinfoLayout.addWidget( self.prof_show_elev_min, 3, 0) 
        self.prof_show_elev_mean = QCheckBox("mean")        
        profinfoLayout.addWidget( self.prof_show_elev_mean, 4, 0)
        self.prof_show_elev_max = QCheckBox("max")        
        profinfoLayout.addWidget( self.prof_show_elev_max, 5, 0)       
  
        layout.addWidget( profinfoGroupBox )
              
        
    def setup_output_section( self, layout ):
 
        outputGroupBox = QGroupBox( "Output" )
        outputLayout = QGridLayout()
        outputGroupBox.setLayout(outputLayout)         

        outputLayout.addWidget( QLabel( self.tr("Save plot as:") ), 0, 0) 
        
        self.OutputPlotFilename_lineedit = QLineEdit()
        outputLayout.addWidget( self.OutputPlotFilename_lineedit, 0, 1 )

        self.OutputPlotBrowse_pushbutton = guiqt.push_button_with_icon(self.openfolder_icon_path, self.button_size, self.set_plot_export )       
        outputLayout.addWidget( self.OutputPlotBrowse_pushbutton, 0, 2 )
        
        self.OutputPlotSave_pushbutton = QPushButton(self.tr("Save"))
        QObject.connect(self.OutputPlotSave_pushbutton, SIGNAL("clicked()"), self.export_plot)               
        outputLayout.addWidget( self.OutputPlotSave_pushbutton, 0, 3 )        

        outputLayout.addWidget( QLabel( self.tr("Save data as:") ), 1, 0) 
        
        self.OutputDataFileName_lineedit = QLineEdit()
        outputLayout.addWidget( self.OutputDataFileName_lineedit, 1, 1 )

        self.OutputDataBrowse_pushbutton = guiqt.push_button_with_icon(self.openfolder_icon_path, 
                                                                       self.button_size, 
                                                                       self.set_data_export)       
        outputLayout.addWidget( self.OutputDataBrowse_pushbutton, 1, 2 )

        self.OutputDataSave_pushbutton = QPushButton(self.tr("Save"))
        QObject.connect(self.OutputDataSave_pushbutton, SIGNAL("clicked()"), self.export_data)               
        outputLayout.addWidget( self.OutputDataSave_pushbutton, 1, 3 ) 
                                              
        layout.addWidget( outputGroupBox )                               


    def set_plot_export(self):
        
        pass


    def export_plot(self):
        
        pass
    
    
    def set_data_export(self):
        
        fileName = QFileDialog.getSaveFileName(self, 
                                               self.tr("Save data"),
                                                "*.txt",
                                                self.tr("Text (*.txt)"))
        
        self.OutputDataFileName_lineedit.setText(fileName)
         

    def export_data(self):
        
        fileName = self.OutputDataFileName_lineedit.text()
                          
        with open(fileName, 'w') as f:
            f.write('cumulated_distance, point_height\n')
            for result in self.profile:
                f.write('%.1f, %.2f\n' % (result[1], result[2]))

    
            
                
         
class AnchoredText(AnchoredOffsetbox):
    """
    Creation of an info box in the plot
    
    """
    def __init__(self, s, loc, pad=0.4, borderpad=0.5, prop=None, frameon=True):

        self.txt = TextArea( s, minimumdescent=False )

        super(AnchoredText, self).__init__(loc, pad=pad, borderpad=borderpad,
                                           child=self.txt,
                                           prop=prop,
                                           frameon=frameon)


            
 