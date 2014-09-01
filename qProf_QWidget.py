# -*- coding: utf-8 -*-


import os
from math import isnan, sin, cos, asin, radians, degrees, floor, ceil
import numpy as np

import copy

import xml.dom.minidom

import unicodedata

from osgeo import ogr

from PyQt4.QtCore import *
from PyQt4.QtGui import *

from qgis.core import QgsPoint, QgsRaster, QgsMapLayerRegistry, QgsMapLayer, QGis
from qgis.gui import QgsRubberBand

from geosurf.qt_utils import lastUsedDir, setLastUsedDir

from geosurf.spatial import Point_2D, Segment_2D, Vector_2D, MultiLine_2D, Line_2D
from geosurf.spatial import Point_3D, Segment_3D, Vector_3D, MultiLine_3D, Line_3D                           
from geosurf.spatial import merge_lines, cartes_plane_from_points, GeolPlane, GeolAxis, CartesianPlane, ParamLine
from geosurf.spatial import xytuple_list2_to_MultiLine_2D
                            
from geosurf.geoio import QGisRasterParameters
from geosurf.profiles import Profiles, TopoProfile, ProfileDEM
from geosurf.geodetic import TrackPointGPX
from geosurf.intersections import map_struct_pts_on_section, calculate_distance_with_sign
from geosurf.errors import Vector_Input_Errors, GPXIOException, VectorIOException
           
from geosurf.qgs_tools import loaded_line_layers, loaded_point_layers, pt_geoms_attrs, \
                            raster_qgis_params, loaded_monoband_raster_layers, \
                            line_geoms_with_id, qgs_point, project_qgs_point, vect_attrs, \
                            line_geoms_attrs, field_values, MapDigitizeTool
            
from mpl.mpl_widget import MplWidget, plot_line, plot_filled_line
from mpl.utils import valid_intervals
from projections import project_line_2d

        
class qprof_QWidget( QWidget ):
    
    colors = ['orange', 'green', 'red', 'grey', 'brown', 'yellow', 'magenta', 'black', 'blue', 'white', 'cyan', 'chartreuse' ]

    def __init__( self, canvas ):

        super( qprof_QWidget, self ).__init__() 
        self.mapcanvas = canvas
        self.initialize_parameters()      
        self.setup_gui()
           

    def initialize_parameters(self):
 
        self.profile_windows = []
        self.cross_section_windows = [] 
        self.profile_line = None   
        self.current_directory = os.path.dirname(__file__) 
        self.profiles = None
        self.DEM_data_export = None 
        self.profile_GPX = None
        self.selected_dem_colors = []
        self.plane_attitudes_colors = []
        self.curve_colors = [] 
        
        
    def setup_gui( self ): 

        self.dialog_layout = QVBoxLayout()
        self.main_widget = QTabWidget()        
        self.main_widget.addTab( self.setup_topoprofile_tab(), "Profile" )         
        self.main_widget.addTab( self.setup_project_section_tab(), "Project" )
        self.main_widget.addTab( self.setup_about_tab(), "About" )

        QgsMapLayerRegistry.instance().layerWasAdded.connect( self.refresh_struct_point_lyr_combobox )        
        QgsMapLayerRegistry.instance().layerWasAdded.connect( self.refresh_struct_line_lyr_combobox )

        QgsMapLayerRegistry.instance().layerRemoved.connect( self.refresh_struct_point_lyr_combobox )                
        QgsMapLayerRegistry.instance().layerRemoved.connect( self.refresh_struct_line_lyr_combobox )
                
        self.dialog_layout.addWidget(self.main_widget)                             
        self.setLayout(self.dialog_layout)            
        self.adjustSize()               
        self.setWindowTitle( 'qProf' )        
                
  
    def setup_topoprofile_tab( self ):  

        profile_widget = QWidget() 
        profile_layout = QVBoxLayout()
 
        profile_toolbox = QToolBox()        
        
        ## input section
        
        profileDEM_QWidget = QWidget()
        profileDEM_Layout = QVBoxLayout()        

        ## Input from DEM
        
        inputDEM_QGroupBox = QGroupBox( profileDEM_QWidget )
        inputDEM_QGroupBox.setTitle("Input DEMs")
        
        inputDEM_Layout = QGridLayout()

        self.DefineSourceDEMs_pushbutton = QPushButton(self.tr("Define source DEMs")) 
        self.DefineSourceDEMs_pushbutton.clicked.connect( self.define_source_DEMs )
            
        inputDEM_Layout.addWidget(self.DefineSourceDEMs_pushbutton, 0, 0, 1, 3 ) 

        inputDEM_QGroupBox.setLayout( inputDEM_Layout )
        
        profileDEM_Layout.addWidget( inputDEM_QGroupBox )
                
        ## Line layer input
        
        inputLine_QGroupBox = QGroupBox( profileDEM_QWidget )
        inputLine_QGroupBox.setTitle("Input line")
        
        inputLine_Layout = QGridLayout()
        
        self.DefineLine_pushbutton = QPushButton(self.tr("Define"))  
        self.DefineLine_pushbutton.clicked.connect( self.define_line ) 
        inputLine_Layout.addWidget(self.DefineLine_pushbutton, 0, 0, 1, 2 )
 
        self.LoadLineLayer_checkbox = QRadioButton( self.tr("from layer") )
        self.LoadLineLayer_checkbox.setChecked( True )
        inputLine_Layout.addWidget( self.LoadLineLayer_checkbox, 0, 2, 1, 1 )   
               
        self.DigitizeLine_checkbox = QRadioButton( self.tr("by digitization") )
        inputLine_Layout.addWidget( self.DigitizeLine_checkbox, 0, 3, 1, 1 )
                        
        inputLine_QGroupBox.setLayout( inputLine_Layout )
        
        profileDEM_Layout.addWidget( inputLine_QGroupBox )
                  
         
        ## create profile section
        
        plotDEM_QGroupBox = QGroupBox( profileDEM_QWidget )
        plotDEM_QGroupBox.setTitle( 'Create profile')
        
        plotDEM_Layout = QGridLayout()                

        # profile options
        
        # trace sampling distance                 
        plotDEM_Layout.addWidget( QLabel( self.tr("Line densify distance") ), 0, 0, 1, 1 )         
        self.profile_densify_distance_lineedit = QLineEdit()
        plotDEM_Layout.addWidget( self.profile_densify_distance_lineedit, 0, 1, 1, 3 )
         
        self.DEM_plot_height_checkbox = QCheckBox( self.tr( "height"))
        self.DEM_plot_height_checkbox.setChecked( True ) 
        plotDEM_Layout.addWidget( self.DEM_plot_height_checkbox, 1, 0, 1, 1 )  

        self.DEM_plot_height_filled_checkbox = QCheckBox( self.tr( "(filled)"))
        plotDEM_Layout.addWidget( self.DEM_plot_height_filled_checkbox, 1, 1, 1, 1 ) 

        self.DEM_exageration_1_1_checkbox = QCheckBox( self.tr( "scale ratio 1:1"))
        plotDEM_Layout.addWidget( self.DEM_exageration_1_1_checkbox, 1, 2, 1, 2 ) 
                
        self.DEM_plot_slope_checkbox = QCheckBox( self.tr( "slope (degrees)"))
        plotDEM_Layout.addWidget( self.DEM_plot_slope_checkbox, 2, 0, 1, 1 ) 
       
        self.DEM_plot_slope_filled_checkbox = QCheckBox( self.tr( "(filled)"))
        plotDEM_Layout.addWidget( self.DEM_plot_slope_filled_checkbox, 2, 1, 1, 1 ) 
        
        self.CreateProfDEM_pushbutton = QPushButton(self.tr("Create profile")) 
        self.CreateProfDEM_pushbutton.clicked.connect( self.create_topo_profiles_from_DEMs )
                       
        plotDEM_Layout.addWidget( self.CreateProfDEM_pushbutton, 2, 2, 1, 3 )

     
        plotDEM_QGroupBox.setLayout( plotDEM_Layout )

        profileDEM_Layout.addWidget( plotDEM_QGroupBox )
         
                 
        ## Export section 
        
        exportDEM_QGroupBox = QGroupBox( profileDEM_QWidget )
        exportDEM_QGroupBox.setTitle( 'Export')
        
        exportDEM_Layout = QGridLayout()         

        self.Export_fromDEMData_pushbutton = QPushButton(self.tr("Export profiles as")) 
        self.Export_fromDEMData_pushbutton.clicked.connect( self.export_from_DEM_data )
        exportDEM_Layout.addWidget( self.Export_fromDEMData_pushbutton, 0, 0, 1, 2 )        

        self.ExportfromDEM_asCSV_checkbox = QCheckBox( self.tr( "csv"))
        exportDEM_Layout.addWidget( self.ExportfromDEM_asCSV_checkbox, 0, 2, 1, 1 )  

        self.ExportfromDEM_asPtShp_checkbox = QCheckBox( self.tr( "2D point shp"))
        exportDEM_Layout.addWidget( self.ExportfromDEM_asPtShp_checkbox, 0, 3, 1, 1 ) 

        self.Export3DLine_pushbutton = QPushButton(self.tr("Export 3D line from DEM ")) 
        self.Export3DLine_pushbutton.clicked.connect( self.write_DEM_3D_lnshp )
        exportDEM_Layout.addWidget( self.Export3DLine_pushbutton, 1, 0, 1, 2 )
         
        self.DEM_3D_Export_comboBox = QComboBox()
        exportDEM_Layout.addWidget(self.DEM_3D_Export_comboBox, 1, 2, 1, 2)   

        exportDEM_QGroupBox.setLayout( exportDEM_Layout )
        
        profileDEM_Layout.addWidget( exportDEM_QGroupBox )
         
        ## 
               
        profileDEM_QWidget.setLayout( profileDEM_Layout )
                
        profile_toolbox.addItem ( profileDEM_QWidget, "Topographic profiles from DEMs" ) 
        
                        
        ## Input from GPX
        
        profileGPX_QWidget = QWidget()
        profileGPX_Layout = QVBoxLayout()        

        ## Input from DEM
        
        inputGPX_QGroupBox = QGroupBox( profileDEM_QWidget )
        inputGPX_QGroupBox.setTitle('Input')
        
        inputGPX_Layout = QGridLayout()
                        
        inputGPX_Layout.addWidget( QLabel( self.tr("Input GPX file with track points:") ), 0, 0, 1, 3)       

        self.input_gpx_lineEdit = QLineEdit()
        self.input_gpx_lineEdit.setPlaceholderText("my_track.gpx")
        inputGPX_Layout.addWidget(self.input_gpx_lineEdit, 1, 0, 1, 2)
        
        self.input_gpx_QPButt = QPushButton( "..." )
        self.input_gpx_QPButt.clicked.connect( self.select_input_gpxFile )
        inputGPX_Layout.addWidget(self.input_gpx_QPButt, 1, 2, 1, 1)
               
        inputGPX_QGroupBox.setLayout( inputGPX_Layout )
        
        profileGPX_Layout.addWidget( inputGPX_QGroupBox )
        
        # Plot section

        plotGPX_QGroupBox = QGroupBox( profileGPX_QWidget )
        plotGPX_QGroupBox.setTitle('Plot')
        
        plotGPX_Layout = QGridLayout()
                
        plotGPX_Layout.addWidget( QLabel( self.tr("Plot:") ), 2, 0, 1, 1 )  
        
        self.GPX_plot_height_checkbox = QCheckBox( self.tr( "height"))
        self.GPX_plot_height_checkbox.setChecked( True ) 
        plotGPX_Layout.addWidget( self.GPX_plot_height_checkbox, 2, 1, 1, 1 )  

        self.GPX_plot_slope_checkbox = QCheckBox( self.tr( "slope"))
        self.GPX_plot_slope_checkbox.setChecked( False )
        plotGPX_Layout.addWidget( self.GPX_plot_slope_checkbox, 2, 2, 1, 1 ) 
        
        self.CalcProf3D_pushbutton = QPushButton(self.tr("Create profile")) 
        self.CalcProf3D_pushbutton.clicked.connect( self.create_topo_profile_from_GPX )
             
        plotGPX_Layout.addWidget( self.CalcProf3D_pushbutton, 3, 0, 1, 3 )                       

        plotGPX_QGroupBox.setLayout( plotGPX_Layout )
        
        profileGPX_Layout.addWidget( plotGPX_QGroupBox )
        
        ## Export Section

        exportGPX_QGroupBox = QGroupBox( profileGPX_QWidget )
        exportGPX_QGroupBox.setTitle('Export')
        
        exportGPX_Layout = QGridLayout()
        
        self.Export_fromGPXData_pushbutton = QPushButton(self.tr("Export profiles as:")) 
        self.Export_fromGPXData_pushbutton.clicked.connect( self.export_from_GPX_data )
        exportGPX_Layout.addWidget( self.Export_fromGPXData_pushbutton, 4, 0, 1, 3 )         
        
        self.ExportfromGPX_asCSV_checkbox = QCheckBox( self.tr( "csv"))
        exportGPX_Layout.addWidget( self.ExportfromGPX_asCSV_checkbox, 5, 0, 1, 1 ) 
        
        self.ExportfromGPX_asPtShp_checkbox = QCheckBox( self.tr( "2D point shp"))
        exportGPX_Layout.addWidget( self.ExportfromGPX_asPtShp_checkbox, 5, 1, 1, 1 )         
        
        self.ExportfromGPX_asLnShp_checkbox = QCheckBox( self.tr( "3D line shp"))
        exportGPX_Layout.addWidget( self.ExportfromGPX_asLnShp_checkbox, 5, 2, 1, 1 )         
           
        exportGPX_QGroupBox.setLayout( exportGPX_Layout )
        
        profileGPX_Layout.addWidget( exportGPX_QGroupBox )
         
        ## 
               
        profileGPX_QWidget.setLayout( profileGPX_Layout )                
        profile_toolbox.addItem ( profileGPX_QWidget, "Topographic profiles from GPXs" ) 
                 
        # widget final setup 
                       
        profile_layout.addWidget(profile_toolbox)
        profile_widget.setLayout(profile_layout) 
        
        return profile_widget     
        

    def setup_project_section_tab( self ):
        
        section_project_QWidget = QWidget()  
        section_project_layout = QVBoxLayout() 

        project_toolbox = QToolBox()
        
        ### Point project toolbox
        
        xs_point_proj_QWidget = QWidget()
        xs_point_proj_Layout = QVBoxLayout()         

        ## input section
        
        xs_input_point_proj_QGroupBox = QGroupBox( xs_point_proj_QWidget )
        xs_input_point_proj_QGroupBox.setTitle('Input')
        
        xs_input_point_proj_Layout = QGridLayout()
                
        # input point geological layer
                
        xs_input_point_proj_Layout.addWidget( QLabel("Layer"), 0, 0, 1, 1 )
        self.prj_struct_point_comboBox = QComboBox()
        self.prj_struct_point_comboBox.currentIndexChanged[int].connect( self.update_point_layers_boxes )
               
        xs_input_point_proj_Layout.addWidget( self.prj_struct_point_comboBox, 0, 1, 1, 3 )        
        self.refresh_struct_point_lyr_combobox()
                   
        xs_input_point_proj_Layout.addWidget( QLabel("Id field"), 1, 0, 1, 1 )
        self.proj_point_id_fld_comboBox = QComboBox()
        xs_input_point_proj_Layout.addWidget( self.proj_point_id_fld_comboBox, 1, 1, 1, 3 )
                     
        xs_input_point_proj_Layout.addWidget( QLabel("Dip dir. field"), 2, 0, 1, 1 )
        self.proj_point_dipdir_fld_comboBox = QComboBox()
        xs_input_point_proj_Layout.addWidget( self.proj_point_dipdir_fld_comboBox, 2, 1, 1, 1 ) 
                
        xs_input_point_proj_Layout.addWidget( QLabel("Dip ang. field"), 2, 2, 1, 1 )
        self.proj_point_dipang_fld_comboBox = QComboBox()
        xs_input_point_proj_Layout.addWidget( self.proj_point_dipang_fld_comboBox, 2, 3, 1, 1 )        
        
        xs_input_point_proj_QGroupBox.setLayout( xs_input_point_proj_Layout )        
        xs_point_proj_Layout.addWidget( xs_input_point_proj_QGroupBox )        


        ## interpolation method
        
        xs_method_point_proj_QGroupBox = QGroupBox( xs_point_proj_QWidget )
        xs_method_point_proj_QGroupBox.setTitle('Projection method')
        
        xs_method_point_proj_Layout = QGridLayout()
        
        self.nearest_point_proj_choice = QRadioButton( "Nearest intersection")
        xs_method_point_proj_Layout.addWidget( self.nearest_point_proj_choice, 0, 0, 1, 3 )        

        self.axis_common_point_proj_choice = QRadioButton( "Along axis projection - common axis")
        xs_method_point_proj_Layout.addWidget( self.axis_common_point_proj_choice, 1, 0, 1, 3 )
        
        xs_method_point_proj_Layout.addWidget( QLabel("Trend"), 2, 0, 1, 1 )
        
        self.common_axis_point_trend_SpinBox = QDoubleSpinBox()
        self.common_axis_point_trend_SpinBox.setMinimum( 0.0 )
        self.common_axis_point_trend_SpinBox.setMaximum( 359.9 ) 
        self.common_axis_point_trend_SpinBox.setDecimals( 1 )
        xs_method_point_proj_Layout.addWidget( self.common_axis_point_trend_SpinBox, 2, 1, 1, 1 )

        xs_method_point_proj_Layout.addWidget( QLabel("Plunge"), 2, 2, 1, 1 )
                
        self.common_axis_point_plunge_SpinBox = QDoubleSpinBox()
        self.common_axis_point_plunge_SpinBox.setMinimum( 0.0 )
        self.common_axis_point_plunge_SpinBox.setMaximum( 89.9 ) 
        self.common_axis_point_plunge_SpinBox.setDecimals( 1 )
        xs_method_point_proj_Layout.addWidget( self.common_axis_point_plunge_SpinBox, 2, 3, 1, 1 )                       

        self.axis_individual_point_proj_choice = QRadioButton( "Along axis projection - individual axes")
        xs_method_point_proj_Layout.addWidget( self.axis_individual_point_proj_choice, 3, 0, 1, 4 )

        xs_method_point_proj_Layout.addWidget( QLabel("Trend field"), 4, 0, 1, 1 )
        self.proj_point_indivax_trend_fld_comboBox = QComboBox()
        xs_method_point_proj_Layout.addWidget( self.proj_point_indivax_trend_fld_comboBox, 4, 1, 1, 1 ) 

        xs_method_point_proj_Layout.addWidget( QLabel("Plunge field"), 4, 2, 1, 1 )
        self.proj_point_indivax_plunge_fld_comboBox = QComboBox()
        xs_method_point_proj_Layout.addWidget( self.proj_point_indivax_plunge_fld_comboBox, 4, 3, 1, 1 )
       
        xs_method_point_proj_QGroupBox.setLayout( xs_method_point_proj_Layout )        
        xs_point_proj_Layout.addWidget( xs_method_point_proj_QGroupBox )        

        
        ## Plot groupbox
                
        xs_plot_proj_QGroupBox = QGroupBox( xs_point_proj_QWidget )
        xs_plot_proj_QGroupBox.setTitle('Plot')
        
        xs_plot_proj_Layout = QGridLayout() 

        xs_plot_proj_Layout.addWidget( QLabel("Color"), 0, 0, 1, 1 )
        self.proj_point_color_comboBox = QComboBox()
        self.proj_point_color_comboBox.addItems( qprof_QWidget.colors )
        xs_plot_proj_Layout.addWidget( self.proj_point_color_comboBox, 0, 1, 1, 1 )        
        
        xs_plot_proj_Layout.addWidget( QLabel("Add labels"), 1, 0, 1, 1 )
        self.plot_prj_add_pt_id_label = QCheckBox( "Id")
        xs_plot_proj_Layout.addWidget( self.plot_prj_add_pt_id_label, 1, 1, 1, 1 ) 
        self.plot_prj_add_trendplunge_label = QCheckBox( "Dip dir/plunge")
        xs_plot_proj_Layout.addWidget( self.plot_prj_add_trendplunge_label, 1, 2, 1, 1 )       

        self.project_point_pushbutton = QPushButton(self.tr("Project geological attitudes"))
        self.project_point_pushbutton.clicked.connect( self.create_struct_point_projection )
        
        xs_plot_proj_Layout.addWidget( self.project_point_pushbutton, 2, 0, 1, 4 )

        self.reset_point_pushbutton = QPushButton(self.tr("Reset geological attitudes"))
        self.reset_point_pushbutton.clicked.connect( self.reset_struct_point_projection )

        xs_plot_proj_Layout.addWidget( self.reset_point_pushbutton, 3, 0, 1, 4 )
                                                        
        xs_plot_proj_QGroupBox.setLayout( xs_plot_proj_Layout )        
        xs_point_proj_Layout.addWidget( xs_plot_proj_QGroupBox )                

        self.flds_prj_point_comboBoxes = [self.proj_point_id_fld_comboBox,
                                          self.proj_point_dipdir_fld_comboBox,
                                          self.proj_point_dipang_fld_comboBox,
                                          self.proj_point_indivax_trend_fld_comboBox,
                                          self.proj_point_indivax_plunge_fld_comboBox ]
                                     
        
        ## output section
        
        xs_output_point_proj_QGroupBox = QGroupBox( xs_point_proj_QWidget )
        xs_output_point_proj_QGroupBox.setTitle( 'Output' )
        
        xs_output_point_proj_Layout = QGridLayout()

        self.save_proj_point_results_asCSV_checkbox = QCheckBox( self.tr( "csv"))
        xs_output_point_proj_Layout.addWidget( self.save_proj_point_results_asCSV_checkbox, 0, 0, 1, 1 ) 
       
        self.save_proj_point_results_asPtShp_checkbox = QCheckBox( self.tr( "3D point shp"))
        xs_output_point_proj_Layout.addWidget( self.save_proj_point_results_asPtShp_checkbox, 0, 1, 1, 1 )  
                        
        self.save_proj_point_results_pushbutton = QPushButton(self.tr("Export projected attitudes")) 
        self.save_proj_point_results_pushbutton.clicked.connect( self.save_proj_points_results )
        xs_output_point_proj_Layout.addWidget( self.save_proj_point_results_pushbutton, 1, 0, 1, 4 )

        xs_output_point_proj_QGroupBox.setLayout(xs_output_point_proj_Layout)        
        xs_point_proj_Layout.addWidget( xs_output_point_proj_QGroupBox ) 
        
        xs_point_proj_QWidget.setLayout( xs_point_proj_Layout )                
        project_toolbox.addItem ( xs_point_proj_QWidget, "Geological attitudes" ) 
                
        ## END Point project toolbox
        
        ### Line project toolbox
        
        xs_line_proj_QWidget = QWidget()
        xs_line_proj_Layout = QVBoxLayout()         

        ## input section
        
        xs_input_line_proj_QGroupBox = QGroupBox( xs_line_proj_QWidget )
        xs_input_line_proj_QGroupBox.setTitle('Input')
        
        xs_input_line_proj_Layout = QGridLayout()
                
        # input geological layer
                
        xs_input_line_proj_Layout.addWidget( QLabel("Layer"), 0, 0, 1, 1 )
        self.prj_input_line_comboBox = QComboBox()
        self.prj_input_line_comboBox.currentIndexChanged[int].connect( self.update_line_layers_boxes )
               
        xs_input_line_proj_Layout.addWidget( self.prj_input_line_comboBox, 0, 1, 1, 3 )        
        self.refresh_struct_line_lyr_combobox()
                   
        xs_input_line_proj_Layout.addWidget( QLabel("Id field"), 1, 0, 1, 1 )
        self.id_fld_line_prj_comboBox = QComboBox()
        xs_input_line_proj_Layout.addWidget( self.id_fld_line_prj_comboBox, 1, 1, 1, 3 )

        xs_input_line_proj_Layout.addWidget( QLabel( "Line densify distance" ), 2, 0, 1, 1 )         
        self.project_line_densify_distance_lineedit = QLineEdit()
        xs_input_line_proj_Layout.addWidget( self.project_line_densify_distance_lineedit, 2, 1, 1, 3 )

        xs_input_line_proj_QGroupBox.setLayout( xs_input_line_proj_Layout )        
        xs_line_proj_Layout.addWidget( xs_input_line_proj_QGroupBox )       

        ## interpolation method
        
        xs_method_line_proj_QGroupBox = QGroupBox( xs_line_proj_QWidget )
        xs_method_line_proj_QGroupBox.setTitle('Project')
        
        xs_method_line_proj_Layout = QGridLayout()

        xs_method_line_proj_Layout.addWidget( QLabel("Projection axis"), 0, 0, 1, 4 )
                
        xs_method_line_proj_Layout.addWidget( QLabel("Trend"), 1, 0, 1, 1 )
        
        self.common_axis_line_trend_SpinBox = QDoubleSpinBox()
        self.common_axis_line_trend_SpinBox.setMinimum( 0.0 )
        self.common_axis_line_trend_SpinBox.setMaximum( 359.9 ) 
        self.common_axis_line_trend_SpinBox.setDecimals( 1 )
        xs_method_line_proj_Layout.addWidget( self.common_axis_line_trend_SpinBox, 1, 1, 1, 1 )

        xs_method_line_proj_Layout.addWidget( QLabel("Plunge"), 1, 2, 1, 1 )
                
        self.common_axis_line_plunge_SpinBox = QDoubleSpinBox()
        self.common_axis_line_plunge_SpinBox.setMinimum( 0.0 )
        self.common_axis_line_plunge_SpinBox.setMaximum( 89.9 ) 
        self.common_axis_line_plunge_SpinBox.setDecimals( 1 )
        xs_method_line_proj_Layout.addWidget( self.common_axis_line_plunge_SpinBox, 1, 3, 1, 1 )                       

        # calculate profile

        xs_method_line_proj_Layout.addWidget( QLabel("Color"), 2, 0, 1, 1 )
         
        self.project_line_color_comboBox = QComboBox()
        self.project_line_color_comboBox.addItems( qprof_QWidget.colors )
        xs_method_line_proj_Layout.addWidget( self.project_line_color_comboBox, 2, 1, 1, 3 )
                         
        self.project_line_pushbutton = QPushButton(self.tr("Plot traces"))
        self.project_line_pushbutton.clicked.connect( self.create_struct_line_projection )
        
        xs_method_line_proj_Layout.addWidget( self.project_line_pushbutton, 3, 0, 1, 4 )

        self.reset_curves_pushbutton = QPushButton(self.tr("Reset traces"))
        self.reset_curves_pushbutton.clicked.connect( self.reset_structural_lines_projection )

        xs_method_line_proj_Layout.addWidget( self.reset_curves_pushbutton, 4, 0, 1, 4 )
                                                
        xs_method_line_proj_QGroupBox.setLayout( xs_method_line_proj_Layout )        
        xs_line_proj_Layout.addWidget( xs_method_line_proj_QGroupBox )                

        self.flds_prj_line_comboBoxes = [ self.id_fld_line_prj_comboBox ]                                     

        ## output section
        
        xs_output_line_proj_QGroupBox = QGroupBox( xs_line_proj_QWidget )
        xs_output_line_proj_QGroupBox.setTitle( 'Output' )
        
        xs_output_line_proj_Layout = QGridLayout()
              
        self.save_proj_line_results_pushbutton = QPushButton(self.tr("Export projected lines as csv")) 
        self.save_proj_line_results_pushbutton.clicked.connect( self.save_proj_lines_results )
        xs_output_line_proj_Layout.addWidget( self.save_proj_line_results_pushbutton, 1, 0, 1, 4 )

        xs_output_line_proj_QGroupBox.setLayout(xs_output_line_proj_Layout)        
        xs_line_proj_Layout.addWidget( xs_output_line_proj_QGroupBox ) 
                
        ## 
        
        xs_line_proj_QWidget.setLayout( xs_line_proj_Layout )                
        project_toolbox.addItem ( xs_line_proj_QWidget, "Geological traces" ) 
    
        ## END Line project toolbox
                
        # widget final setup
                
        section_project_layout.addWidget( project_toolbox ) 
        # section_project_layout.addWidget( xs_line_proj_toolbox )     
        section_project_QWidget.setLayout( section_project_layout )
        
        return section_project_QWidget


    def setup_about_tab(self):
        
        about_widget = QWidget()  
        about_layout = QVBoxLayout( )
        
        htmlText = """
        <h3>qProf</h3>
        Created by M. Alberti (www.malg.eu) and M. Zanieri.
        <br />Concept: M. Zanieri and M. Alberti, implementation: M. Alberti.
        <br />We thank S. Peduzzi for his vigorous testing.
        <br />Plugin for creating profiles from DEM and GPX files, and as an aid in creating geological cross-sections.       
        <br />
        <h4>Profile</h4>
        It is possible to create a topographic profile from one or more DEMs and a line digitized or stored in a layer, or directly 
        from a GPX file storing track points. Data sources can be projected in different CRS,
        and the created profiles will be in the project CRS.
        When directly digitizing a line in the map, you add points wtih left clicks, and stops a line with a right click. 
        <br /><br />When using a line layer, multiple lines will be merged into a single line, based on the 
        chosen order field when available, or otherwise on the internal line order.
        Some artifacts in derived profiles can be due to erroneous line ordering, not corrected by 
        defining an order in an integer field (order values start from 1).
        <br /><br />When calculating a profile from DEM, you must define the 'Line densify distance',
        i.e. the distance between consecutive sampling point automatically added
        when the original vertices of the used path are distanced more than this value.
        It is suggested to use a value comparable to the resolution of the used DEM,
        for instance 30 m for Aster DEMs. 
        <br /><br />After having calculated the profile, you can plot its elevations and slopes (as degrees),
        and save the results as a csv file, a 2D point shapefile or a 3D line shapefile.
        <br />
        <h4>Project</h4>
        Having defined and calculated a profile as previously described,
        it is also possible to project geological attitudes or traces on the section.
         <br /> <br /><u>Projection of geological attitudes</u>
        <br /><br />The geological attitudes source is a point layer: just selected points will be projected, 
        unless, in case of no selection, all points will be projected.
        <br /> <br />Required fields are the geological point <i>id</i> and its surface orientation, expressed by <i>dip direction</i> and 
        <i>dip angle</i> values.
        <br />The geological attitudes can be projected on the section plane according to three
        methods: 1) nearest point; 2) projection along a common axis; 3) projection along individual 
        axes, for each geological record.
        <br /> <br />When choosing the <i>Along axis projection - individual axes</i> option, <i>trend</i> and <i>plunge</i> fields,
        storing the fold axis values along which to project each observation, are required in the source point layer. 
        <br /><br />The results can be exported both as a csv file or as a 3D point shapefile.
                
        <br /> <br /><u>Projection of trace lines</u>
        <br /><br />Geological traces can be projected on the section plane, based on a fold axis for which 
        trend and plunge values have to be defined.
        
        
        

        <br />
        """
        
        aboutQTextBrowser = QTextBrowser( about_widget )        
        aboutQTextBrowser.insertHtml( htmlText )
         
        about_layout.addWidget( aboutQTextBrowser )  
        about_widget.setLayout(about_layout) 
        
        return about_widget


    def refresh_struct_point_lyr_combobox(self):
        
        self.pointLayers = loaded_point_layers()
        self.prj_struct_point_comboBox.clear()        
        message = "choose"
        self.prj_struct_point_comboBox.addItem( message )
        self.prj_struct_point_comboBox.addItems( [ layer.name() for layer in self.pointLayers ] )              


    def refresh_struct_line_lyr_combobox( self ):
        
        self.current_line_layers = loaded_line_layers()
        self.prj_input_line_comboBox.clear()        
        message = "choose"
        self.prj_input_line_comboBox.addItem( message )
        self.prj_input_line_comboBox.addItems( [ layer.name() for layer in self.current_line_layers ] )              
   

    def define_source_DEMs(self):  
        
        current_raster_layers = loaded_monoband_raster_layers()  
        if len( current_raster_layers ) == 0:
            QMessageBox.critical( self, 
                                      "DEMs sources", 
                                      "No available DEM" )
            return            

        dialog = SourceDEMsDialog( current_raster_layers )

        if dialog.exec_():
            selected_dems, selected_dem_colors = self.get_selected_dems_params( dialog )
        else:
            QMessageBox.critical( self, 
                                      "DEMs", 
                                      "No DEM chosen" )
            return
            
        if  len( selected_dems ) == 0:       
            QMessageBox.critical( self, 
                                      "DEMs sources", 
                                      "No selected DEM" )
        else:
            self.selected_dems = selected_dems
            self.selected_dem_colors = selected_dem_colors 
  
        
    def get_selected_dems_params( self, dialog ):   

        selected_dems = []
        selected_dem_colors = [] 
        for dem_qgis_ndx in range( dialog.listDEMs_treeWidget.topLevelItemCount () ):
            curr_DEM_item = dialog.listDEMs_treeWidget.topLevelItem ( dem_qgis_ndx ) 
            if curr_DEM_item.checkState ( 0 ) == 2:
                selected_dems.append( dialog.singleband_raster_layers_in_project[ dem_qgis_ndx ] )
                selected_dem_colors.append( dialog.listDEMs_treeWidget.itemWidget( curr_DEM_item, 1 ).currentText() )  
         
        return selected_dems, selected_dem_colors
        
 
    def define_line( self ):
        
        if self.DigitizeLine_checkbox.isChecked():
            self.digitize_line()
        else:
            self.load_line_layer()
        
        
    def get_line_layer_params( self, dialog ):
        
        line_layer = dialog.line_shape
        order_field_ndx = dialog.Trace2D_order_field_comboBox.currentIndex() 
        
        return line_layer, order_field_ndx
        
   
    def digitize_line( self ):
        
        try:
            self.rubberband.reset( QGis.Line )
        except:
            pass

        self.previous_maptool = self.mapcanvas.mapTool( )            # Save the standard map tool for restoring it at the end
        self.digitize_maptool = MapDigitizeTool( self.mapcanvas )        #  mouse listener
        self.mapcanvas.setMapTool( self.digitize_maptool )
        self.connect_digitize_maptool()
        
        self.polygon = False
        self.rubberband = QgsRubberBand( self.mapcanvas, self.polygon )
        self.rubberband.setWidth( 1 )
        self.rubberband.setColor( QColor( Qt.red ) )

        self.profile_canvas_points = []


    def connect_digitize_maptool( self ):
        
        QObject.connect( self.digitize_maptool, SIGNAL( "moved"), self.canvas_refresh_profile_line )
        QObject.connect( self.digitize_maptool, SIGNAL( "leftClicked"), self.canvas_add_point_to_profile )        
        QObject.connect( self.digitize_maptool, SIGNAL( "rightClicked"), self.canvas_end_profile_line )


    def disconnect_digitize_maptool( self):
        
        QObject.disconnect( self.digitize_maptool, SIGNAL( "moved"), self.canvas_refresh_profile_line )        
        QObject.disconnect( self.digitize_maptool, SIGNAL( "leftClicked"), self.canvas_add_point_to_profile  )
        QObject.disconnect( self.digitize_maptool, SIGNAL( "rightClicked"), self.canvas_end_profile_line )
         

    def xy_from_canvas( self, position ):
        
        mapPos = self.mapcanvas.getCoordinateTransform().toMapCoordinates( position["x"], position["y"] )
        return  mapPos.x(), mapPos.y()


    def refresh_rubberband(self, xy_list ):
        
        self.rubberband.reset( QGis.Line ) 
        for x,y in xy_list:
            self.rubberband.addPoint( QgsPoint( x, y ) )
                
        
    def canvas_refresh_profile_line(  self, position ):  
   
        if len( self.profile_canvas_points) == 0:
            return     

        x, y = self.xy_from_canvas( position ) 
        self.refresh_rubberband( self.profile_canvas_points + [[x,y]] )           


    def profile_add_point(self, position ):

        x, y = self.xy_from_canvas( position )        
        self.profile_canvas_points.append( [x,y] )        


    def canvas_add_point_to_profile(  self, position ):
        
        if len( self.profile_canvas_points) == 0:
            self.rubberband.reset( self.polygon)
            
        self.profile_add_point( position )
                         

    def canvas_end_profile_line(  self, position ):
        
        self.refresh_rubberband( self.profile_canvas_points )

        self.profile_line = Line_2D( [ Point_2D(x,y) for x,y in self.profile_canvas_points ] )
        #self.profile_points_from_canvas = self.profile_canvas_points
        self.profile_canvas_points = []
       

    def restore_previous_map_tool( self ):
        
        self.mapcanvas.unsetMapTool( self.digitize_maptool )
        self.mapcanvas.setMapTool( self.previous_maptool )


    def load_line_layer( self ):
        
        try:
            self.disconnect_digitize_maptool()
        except:
            pass
        
        try:
            self.rubberband.reset( QGis.Line )
        except:
            pass
        
        
        current_line_layers = loaded_line_layers()   

        if len( current_line_layers ) == 0:
            QMessageBox.critical( self, 
                                      "Line sources", 
                                      "No available layers" )
            return            

        dialog = SourceLineLayerDialog( current_line_layers )

        if dialog.exec_():
            line_layer, order_field_ndx = self.get_line_layer_params( dialog )
        else:
            QMessageBox.critical( self, 
                                      "Line source", 
                                      "No choice made" )
            return

        line_fld_ndx = int( order_field_ndx ) - 1
        # get profile path from input line layer
        success, result = self.get_line_trace( line_layer, line_fld_ndx )
        if not success:
            raise VectorIOException, result
        
        profile_orig_lines, mergeorder_ids = result

        profile_processed_line_2d = merge_lines( profile_orig_lines, mergeorder_ids )   
        
        on_the_fly_projection, project_crs = self.get_on_the_fly_projection()         

        # process input line layer
        profile_projected_line_2d = self.create_projected_line( profile_processed_line_2d, line_layer.crs(), on_the_fly_projection, project_crs )
        
        self.profile_line = profile_projected_line_2d.remove_coincident_successive_points()
        
           
    def get_on_the_fly_projection(self):
        
        on_the_fly_projection = True if self.mapcanvas.hasCrsTransformEnabled() else False
        if on_the_fly_projection:
            project_crs = self.mapcanvas.mapRenderer().destinationCrs()
        else:
            project_crs = None
            
        return on_the_fly_projection, project_crs                   
 

    def get_dem_parameters(self, dem ):
    
        return QGisRasterParameters( *raster_qgis_params( dem ) )

    
    def get_z(self, dem_layer, point ):
        
        identification = dem_layer.dataProvider().identify( QgsPoint( point._x, point._y ), QgsRaster.IdentifyFormatValue )
        if not identification.isValid(): 
            return np.nan
        else:
            try: 
                result_map = identification.results()
                return float( result_map[1] )
            except:
                return np.nan
                          

    def interpolate_bilinear(self, dem, dem_params, point):  
      
        array_coords_dict = dem_params.geogr2raster( point )
        floor_x_raster = floor( array_coords_dict["x"] )
        ceil_x_raster = ceil( array_coords_dict["x"] )
        floor_y_raster = floor( array_coords_dict["y"] )
        ceil_y_raster = ceil( array_coords_dict["y"] )                
         
        # bottom-left center
        p1 = dem_params.raster2geogr( dict(x=floor_x_raster,
                                           y=floor_y_raster))
        # bottom-right center       
        p2 = dem_params.raster2geogr( dict(x=ceil_x_raster,
                                           y=floor_y_raster))
        # top-left center
        p3 = dem_params.raster2geogr( dict(x=floor_x_raster,
                                           y=ceil_y_raster)) 
        # top-right center       
        p4 = dem_params.raster2geogr( dict(x=ceil_x_raster,
                                           y=ceil_y_raster))
         
        z1 = self.get_z( dem, p1 )
        z2 = self.get_z( dem, p2 )         
        z3 = self.get_z( dem, p3 )
        z4 = self.get_z( dem, p4 ) 
        
        delta_x = point._x - p1._x
        delta_y = point._y - p1._y 

        z_x_a = z1 + (z2-z1)*delta_x/dem_params.cellsizeEW
        z_x_b = z3 + (z4-z3)*delta_x/dem_params.cellsizeEW        
        
        return z_x_a + (z_x_b-z_x_a)*delta_y/dem_params.cellsizeNS
            
        
    def interpolate_point_z( self, dem, dem_params, point ):
        
        if dem_params.point_in_interpolation_area( point ):
            return self.interpolate_bilinear( dem, dem_params, point  )
        elif dem_params.point_in_dem_area( point ):
            return self.get_z( dem, point )
        else:
            return np.nan            
                        

    def check_trace_sampling_distance(self, sampling_distance ):
        
        if sampling_distance == '':
            return False, "Line densify distance is not set" 
        try:
            sample_distance = float( sampling_distance )
        except:
            return False, "Input cannot be converted to a float value"        
        if not sample_distance > 0.0:
            return False, "Input value is negative or zero"
        return True, ''
            

    def get_selected_dem_indices( self ):   

        selected_dem_indices = []
        for dem_qgis_ndx in range( self.listDEMs_treeWidget.topLevelItemCount () ):
            curr_DEM_item = self.listDEMs_treeWidget.topLevelItem ( dem_qgis_ndx ) 
            if curr_DEM_item.checkState ( 0 ) == 2:
                selected_dem_indices.append( dem_qgis_ndx )
        return selected_dem_indices
           

    def get_line_trace( self, line_shape, order_field_ndx ):

        try:
            profile_orig_lines, mergeorder_ids = line_geoms_with_id( line_shape, order_field_ndx )
        except Vector_Input_Errors as error_msg:
            return False, error_msg
        return True, ( profile_orig_lines, mergeorder_ids )        
        

    def create_projected_line( self, profile_processed_line, line_layer_crs, on_the_fly_projection, project_crs ):
        
        if not on_the_fly_projection:
            return profile_processed_line
        else: 
            return project_line_2d( profile_processed_line, line_layer_crs, project_crs )            


    def get_DEM_3Dlines( self, dems, dems_params, resampled_line_2d, on_the_fly_projection, project_crs ):
    
        dem_3Dlines = []
        for dem, dem_params in zip( dems, dems_params ):
            if on_the_fly_projection and dem.crs() != project_crs:
                line_with_dem_crs = project_line_2d( resampled_line_2d, project_crs, dem.crs() ) 
            else:
                line_with_dem_crs = resampled_line_2d
                
            profile_3d = Line_3D()
            for pt_dem_crs, pt_project_crs in zip(line_with_dem_crs._pts, resampled_line_2d._pts):
                profile_3d = profile_3d.add_pt( Point_3D( pt_project_crs._x, 
                                                         pt_project_crs._y, 
                                                         self.interpolate_point_z( dem, dem_params, pt_dem_crs ) ) )
            dem_3Dlines.append( profile_3d )
        return dem_3Dlines


    def check_profile_creation_parameters( self, selected_dems, profile_line, sample_distance, plot_height_choice, plot_slope_choice ):

        if len( selected_dems ) == 0:
            return False, "No DEM layer available/selected"

        try:
            if profile_line.num_points() < 2:
                return False, "No line defined"
        except:
                return False, "No line defined"            
        
        sampling_distance_state = self.check_trace_sampling_distance( sample_distance )
        if not sampling_distance_state[0]:
            return False, sampling_distance_state[1]
        
        if not ( plot_height_choice or plot_slope_choice ):
            return False, "Neither height or slope options are selected"

        return True, ''
    
              
    def create_topo_profiles_from_DEMs( self ):

        try:
            selected_dems = self.selected_dems
        except:
            QMessageBox.critical( self, "Error", "DEM(s) not yet defined" )
            return            
        
        try:
            profile_line = self.profile_line
        except:
            QMessageBox.critical( self, "Error", "Profile line not yet defined" )
            return       
        
        # get profile creation parameters                        
        sample_distance = self.profile_densify_distance_lineedit.text()                                
        plot_height_choice = self.DEM_plot_height_checkbox.isChecked()
        plot_slope_choice = self.DEM_plot_slope_checkbox.isChecked() 
                        
        # check profile creation parameters 
        input_parameters_state_ok, msg = self.check_profile_creation_parameters(selected_dems,
                                                                                profile_line, 
                                                                                sample_distance,
                                                                                plot_height_choice,
                                                                                plot_slope_choice )
        
        if not ( input_parameters_state_ok ):
            QMessageBox.critical( self, "Error", msg )
            return
                                 
        try:               
            self.profiles = self.calculate_elevations_from_DEMs( self.selected_dems,
                                                                 float( sample_distance )  )        
        except VectorIOException, msg:
            self.profiles = None
            QMessageBox.critical( self, "Error", msg )
            return             
                
        # update combobox for 3D line export       
        self.DEM_3D_Export_comboBox.clear()        
        self.DEM_3D_Export_comboBox.addItems( self.profiles.get_current_dem_names() )
        
        self.plot_profile_elements( )        

                                                                        
    def calculate_elevations_from_DEMs(self, selected_dems, sample_distance ):

        # get project CRS information
        on_the_fly_projection, project_crs = self.get_on_the_fly_projection() 
                
        # get geodata
        selected_dem_parameters = [ self.get_dem_parameters(dem) for dem in selected_dems ]        

        resampled_line_2d = self.profile_line.densify_nodes( sample_distance ) # line resampled by sample distance
       
        # calculate 3D profiles from DEMs
        dem_3Dlines = self.get_DEM_3Dlines( selected_dems, selected_dem_parameters, resampled_line_2d, on_the_fly_projection, project_crs )
        
        return self.parse_topo_profiles_from_DEM( sample_distance, selected_dems, selected_dem_parameters, resampled_line_2d, dem_3Dlines )


    def parse_topo_profiles_from_DEM(self, sample_distance, dems, dems_parameters, resampled_line_2d, dem_3Dlines ):
        
        profiles = Profiles( sample_distance )

        # extract parameters for final processings       
        profiles.profile_dems = [ ProfileDEM( dem, params ) for (dem, params) in zip(dems, dems_parameters ) ]
        
        profiles.profile_line = resampled_line_2d
                
        for dem_name, line_3d in zip( [ dem.name() for dem in dems ], dem_3Dlines ):
            profiles.add_topo_profile( TopoProfile( dem_name, line_3d ) )
        
        return profiles


    def parse_DEM_results_for_export( self, profiles ):
                
        # definition of output results         
        x_list = profiles.topo_profiles[0].x_list()
        y_list = profiles.topo_profiles[0].y_list() 
        elev_list = [ topo_profile.z_list() for topo_profile in profiles.topo_profiles ]              
        cumdist_2D_list = profiles.topo_profiles[0].get_increm_dist_2d()
        cumdist_3d_list = [ topo_profile.get_increm_dist_3d() for topo_profile in profiles.topo_profiles ]
        slope_list = [ topo_profile.slope_list() for topo_profile in profiles.topo_profiles ]

        elev_list_zip = zip( *elev_list )
        cumdist_3d_list_zip = zip( *cumdist_3d_list ) 
        slope_list_zip = zip( *slope_list )  

        result_data = []
        for x, y, cum_2d_dist, zs, cum3d_dists, slopes in zip( x_list, y_list, cumdist_2D_list, elev_list_zip, cumdist_3d_list_zip, slope_list_zip ):
            record = [ x, y, cum_2d_dist]
            for z, cum3d_dist, slope in zip( zs, cum3d_dists, slopes ):
                if isnan(z): z = ''
                if isnan(cum3d_dist): cum3d_dist = ''
                if isnan(slope): slope = ''
                record += [z, cum3d_dist, slope ]
            result_data.append( record )
         
        return profiles.get_current_dem_names(), result_data

     
    def export_from_DEM_data( self ):
        
        if self.profiles is None:            
            QMessageBox.critical( self, "Saving results", "No DEM-derived profile is available for export" )
            return 
        
        if not ( self.ExportfromDEM_asCSV_checkbox.isChecked() or self.ExportfromDEM_asPtShp_checkbox.isChecked() ):
            QMessageBox.critical( self, "Saving results", "No output format is selected" )
            return   
        
        # process results for data export         
        dem_names, export_data = self.parse_DEM_results_for_export( self.profiles )
                                  
        # definition of field names         
        cum3ddist_headers = ['cumulated_3d_distance']*len(dem_names)
        slopes_headers = ['slopes (degr)']*len(dem_names)
        header_list = ['x', 'y', 'cumulated_2d_distance'] + [ name for sublist in zip(dem_names, cum3ddist_headers, slopes_headers) for name in sublist ]              

        # output for csv file
        if self.ExportfromDEM_asCSV_checkbox.isChecked():            
            self.write_DEM_2D_csv( header_list, export_data )

        # output for 2D pt shapefile            
        if self.ExportfromDEM_asPtShp_checkbox.isChecked():
            self.write_DEM_2D_ptshp( dem_names, export_data )
        
        QMessageBox.information( self, "Profile export", "Finished" )
        

    def write_DEM_2D_csv( self, header_list, export_data ):
               
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
            for rec in export_data:
                out_rec_string = ''
                for val in rec:
                    out_rec_string += str( val ) + ','
                f.write( out_rec_string[:-1]+'\n' )
                
                      
    def write_DEM_2D_ptshp( self, dem_names, export_data ):
        
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
        
        try:    
            ptshp_datasource = shape_driver.CreateDataSource( unicode( fileName ) )
        except TypeError:
            ptshp_datasource = shape_driver.CreateDataSource( str( fileName ) )
            
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
        for ndx, rec in enumerate( export_data ):
            
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
                                
                
    def write_DEM_3D_lnshp( self ):

        if self.profiles is None:            
            QMessageBox.critical( self, "Saving results", "No DEM-derived profile is available for export" )
            return        

        # process results for data export         
        _, export_data = self.parse_DEM_results_for_export( self.profiles ) 
                   
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
            
        try:
            lnshp_datasource = shape_driver.CreateDataSource( unicode( fileName ) )
        except TypeError:
            lnshp_datasource = shape_driver.CreateDataSource( str( fileName ) )
                        
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
        for ndx in range( len( export_data )-1 ):
                
            rec_a = export_data[ndx]
            rec_b = export_data[ndx+1]
                            
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
 
            
    def select_input_gpxFile( self ):
            
        fileName = QFileDialog.getOpenFileName( self, 
                                                self.tr( "Open GPX file" ), 
                                                lastUsedDir(), 
                                                "GPX (*.gpx *.GPX)" )
        if not fileName:
            return
        setLastUsedDir( fileName )    
        self.input_gpx_lineEdit.setText( fileName )


    def check_GPX_profile_parameters(self):
        
        source_gpx_path = unicode( self.input_gpx_lineEdit.text() )
        if source_gpx_path == '':
            return False, 'Source GPX file is not set'

        plot_height_choice = self.GPX_plot_height_checkbox.isChecked()
        plot_slope_choice = self.GPX_plot_slope_checkbox.isChecked()
        
        if not ( plot_height_choice or plot_slope_choice ):
            return False, 'One of height or slope plot options are to be chosen'        
        
        return True, 'OK'
    
        
    def create_topo_profile_from_GPX( self ):

        preliminar_check = self.check_GPX_profile_parameters()
        
        if not preliminar_check[0]:
            QMessageBox.critical( self, "GPX input error", preliminar_check[1] )
            return   
        
        source_gpx_path = unicode( self.input_gpx_lineEdit.text() )
         
        try:       
            self.profile_GPX = self.calculate_profile_from_GPX( source_gpx_path )
        except GPXIOException, msg:
            self.profile_GPX = None
            QMessageBox.critical( self, "GPX input error", msg )
            return 
        
        self.plot_GPX_profile( self.profile_GPX )  
        
 
    def calculate_profile_from_GPX( self, source_gpx_path ):
                                 
        doc = xml.dom.minidom.parse( source_gpx_path )

        # get track point values (lat, lon, elev, time)
        trk_measures = []
        for trk_node in doc.getElementsByTagName( 'trk'):
            trkname = trk_node.getElementsByTagName( 'name' )[0].firstChild.data
            trksegments = trk_node.getElementsByTagName( 'trkseg' )
            for trksegment in trksegments:
                trk_pts = trksegment.getElementsByTagName( 'trkpt' )
                for tkr_pt in trk_pts:               
                    gpx_trackpoint = TrackPointGPX( tkr_pt.getAttribute("lat"),
                                                    tkr_pt.getAttribute("lon"),
                                                    tkr_pt.getElementsByTagName("ele")[0].childNodes[0].data,
                                                    tkr_pt.getElementsByTagName("time")[0].childNodes[0].data )
                    trk_measures.append( gpx_trackpoint )
        
        # check for the presence of track points
        if len( trk_measures ) == 0:
            raise GPXIOException, "No track point found in this file"            
        
        # calculate delta elevations between consecutive points
        delta_elev_values = [ np.nan ]
        for ndx in range( 1, len ( trk_measures ) ):
            delta_elev_values.append( trk_measures[ndx].elev - trk_measures[ndx-1].elev )        
        
        # covert values into ECEF values (x, y, z in ECEF global coordinate system)        
        trk_ECEFpoints = [ trk_value.toPoint4D() for trk_value in  trk_measures ]
        
        # calculate 3D distances between consecutive points
        dist_3D_values = [ np.nan ]
        for ndx in range( 1, len ( trk_ECEFpoints ) ):
            dist_3D_values.append( trk_ECEFpoints[ndx].distance( trk_ECEFpoints[ndx-1] ) ) 
                    
        # calculate slope along track
        slopes = [ ]
        for delta_elev, dist_3D in zip( delta_elev_values, dist_3D_values ):
            try:
                slopes.append( degrees( asin( delta_elev / dist_3D ) ) )
            except:
                slopes.append( np.nan ) 
        
        # calculate horizontal distance along track
        horiz_dist_values = []
        for slope, dist_3D in zip( slopes, dist_3D_values ):
            try:
                horiz_dist_values.append( dist_3D * cos( radians( slope ) ) )
            except: 
                horiz_dist_values.append( np.nan )
                        
        # defines the cumulative 2D distance values
        cum_distances_2D = [ 0.0 ]
        for ndx in range( 1, len( horiz_dist_values ) ):
            cum_distances_2D.append( cum_distances_2D[-1] + horiz_dist_values[ndx] )          

        # defines the cumulative 3D distance values
        cum_distances_3D = [ 0.0 ]
        for ndx in range( 1, len( dist_3D_values ) ):
            cum_distances_3D.append( cum_distances_3D[-1] + dist_3D_values[ndx] ) 
            
        # define GPX names, elevations and slopes
        dataset_name = [ trkname ]
        lat_values = [ track.lat for track in trk_measures ]
        lon_values = [ track.lon for track in trk_measures ]
        time_values = [ track.time for track in trk_measures ]        
        elevations = [ track.elev for track in trk_measures ]           
        
        # define variable for plotting                
        profiles = dict(  lats=lat_values,
                          lons=lon_values,
                          times=time_values,
                          dataset_names=dataset_name,
                          cum_distances_2D=cum_distances_2D,
                          cum_distances_3D=[cum_distances_3D], # [] required for compatibility with DEM plotting
                          elevations=[elevations], # [] required for compatibility with DEM plotting
                          slopes=[slopes]) # [] required for compatibility with DEM plotting

        return profiles


    def plot_GPX_profile( self, profiles ):
 
        dataset_names = profiles['dataset_names']
        cum_distances_2D = profiles['cum_distances_2D']
        elevations = profiles['elevations']
        slopes = profiles['slopes']
                               
        # defines the extent for the plot window: s min and max     
        profiles_s_min, profiles_s_max = cum_distances_2D[0], cum_distances_2D[-1] 

        # defines z min and max values
        elev_list = [ z for z_values in elevations for z in z_values if not isnan( z ) ]
        plot_z_min, plot_z_max = min( elev_list ), max( elev_list )
        delta_z = plot_z_max - plot_z_min 
        plot_z_min, plot_z_max = plot_z_min - delta_z * 0.05, plot_z_max + delta_z * 0.05

        # defines slope min and max values
        slope_list = [ slope for profile_slopes in slopes for slope in profile_slopes if not isnan( slope ) ]
        profiles_slope_min, profiles_slope_max = min( slope_list ), max( slope_list )
        delta_slope = profiles_slope_max - profiles_slope_min 
        profiles_slope_min, profiles_slope_max = profiles_slope_min - delta_slope*0.2, profiles_slope_max + delta_slope*0.2 

        # map
        profile_window = MplWidget()  

        plot_height_choice = self.GPX_plot_height_checkbox.isChecked()
        plot_slope_choice = self.GPX_plot_slope_checkbox.isChecked() 
            
        if plot_height_choice and plot_slope_choice:
            mpl_code_list = [ 211, 212 ]
        else:
            mpl_code_list = [ 111 ]          

        s_2d_values_array = np.array( cum_distances_2D )
        
        subplot_code = mpl_code_list[0]   
        if plot_height_choice:                            
            axes_height = profile_window.canvas.fig.add_subplot( subplot_code )
            axes_height.set_xlim( profiles_s_min, profiles_s_max )
            axes_height.set_ylim( plot_z_min, plot_z_max ) 
            axes_height.set_color_cycle( qprof_QWidget.colors )
            for dem_name, z_values, color in zip( dataset_names, elevations, qprof_QWidget.colors ):              
                z_values_array = np.array( z_values )
                for val_int in valid_intervals( z_values_array ):               
                    axes_height.fill_between( s_2d_values_array[ val_int['start'] : val_int['end']+1 ], 
                                              plot_z_min, 
                                              z_values_array[ val_int['start'] : val_int['end']+1 ], 
                                              facecolor=color, 
                                              alpha=0.1 )                       
                axes_height.plot(cum_distances_2D, z_values,'-', label = unicode(dem_name) )
                
            axes_height.grid(True)
            axes_height.legend(loc = 'upper left', shadow=True)              
  
        if plot_slope_choice:            
            if len(mpl_code_list) == 2: subplot_code = mpl_code_list[1]
            axes_slope = profile_window.canvas.fig.add_subplot( subplot_code )
            axes_slope.set_xlim( profiles_s_min, profiles_s_max )
            axes_slope.set_ylim( profiles_slope_min, profiles_slope_max )
            axes_slope.set_color_cycle( qprof_QWidget.colors )
            for dem_name, profile_slopes, color in zip( dataset_names, slopes, qprof_QWidget.colors):
                
                slope_values_array = np.array( profile_slopes ) 
                for val_int in valid_intervals( slope_values_array ):              
                    axes_slope.fill_between( s_2d_values_array[ val_int['start'] : val_int['end']+1 ], 
                                             0, 
                                             slope_values_array[ val_int['start'] : val_int['end']+1 ], 
                                             facecolor=color, 
                                             alpha=0.1 )                
                axes_slope.plot(cum_distances_2D, profile_slopes,'-', label = unicode(dem_name) )
                
            axes_slope.grid(True)
            axes_slope.legend(loc = 'upper left', shadow=True)  
                    
        profile_window.canvas.draw() 
        
        self.profile_windows.append( profile_window )
    
            
    def parse_GPX_results_for_export( self, GPXprofile ):        

        # definition of output results        
        lat_list = GPXprofile['lats']
        lon_list = GPXprofile['lons'] 
        time_list = GPXprofile['times']          
        cumdist_2D_list = GPXprofile['cum_distances_2D']
        elev_list = GPXprofile['elevations'][0] # [0] required for compatibility with DEM processing                  
        cumdist_3d_list = GPXprofile['cum_distances_3D'] [0] # [0] required for compatibility with DEM processing
        slope_list = GPXprofile['slopes'][0] # [0] required for compatibility with DEM processing

        result_data = []
        for lat, lon, time, elev, cumdist_2D, cumdist_3D, slope in zip( lat_list, lon_list, time_list, elev_list, cumdist_2D_list, cumdist_3d_list, slope_list ):
            if isnan( elev ): elev = ''
            if isnan( cumdist_3D ): cumdist_3D = ''
            if isnan( slope ): slope = ''
            record = [ lat, lon, time, elev, cumdist_2D, cumdist_3D, slope ]
            result_data.append( record )
         
        return result_data

 
    def export_from_GPX_data( self ):
        
        if self.profile_GPX is None:            
            QMessageBox.critical( self, "Saving results", "No GPX-derived profile is available for export" )
            return 
    
        if not ( self.ExportfromGPX_asCSV_checkbox.isChecked() or \
                 self.ExportfromGPX_asPtShp_checkbox.isChecked() or \
                 self.ExportfromGPX_asLnShp_checkbox.isChecked()):
            QMessageBox.critical( self, "Saving results", "No output format is selected" )
            return  

        # process results from export         
        gpx_parsed_results = self.parse_GPX_results_for_export( self.profile_GPX )
                
        # definition of field names        
        header_list = ['lat', 'lon', 'time', 'elev', 'cumulated_2d_distance', 'cumulated_3d_distance', 'slopes (degr)' ]              
        
        # output for csv file
        if self.ExportfromGPX_asCSV_checkbox.isChecked():         
            self.write_results_as_csv( header_list, gpx_parsed_results )

        # output for 2D pt shapefile            
        if self.ExportfromGPX_asPtShp_checkbox.isChecked():
            self.write_GPX_2D_ptshp( gpx_parsed_results )

        if self.ExportfromGPX_asLnShp_checkbox.isChecked():
            self.write_GPX_3D_lnshp( gpx_parsed_results )
                    
        QMessageBox.information( self, "Profile export", "Finished" )        
        
                      
    def write_GPX_2D_ptshp( self, gpx_parsed_results ):
                
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
            
        try:
            ptshp_datasource = shape_driver.CreateDataSource( unicode( fileName ) )
        except TypeError:
            ptshp_datasource = shape_driver.CreateDataSource( str( fileName ) )
                        
        if ptshp_datasource is None:
            QMessageBox.critical( self, "Saving results", "Creation of %s shapefile failed" % os.path.split( fileName )[1] )
            return         
        
        ptshp_layer = ptshp_datasource.CreateLayer( 'profile', geom_type=ogr.wkbPoint )
        if ptshp_layer is None:
            QMessageBox.critical( self, "Saving results", "Layer creation failed" )
            return 
                  
        # creates required fields         
        ptshp_layer.CreateField( ogr.FieldDefn( 'id', ogr.OFTInteger ) )
        ptshp_layer.CreateField( ogr.FieldDefn( 'lat', ogr.OFTReal ) )
        ptshp_layer.CreateField( ogr.FieldDefn( 'lon', ogr.OFTReal ) )
        time_field = ogr.FieldDefn( 'time', ogr.OFTString )
        time_field.SetWidth(20)  
        ptshp_layer.CreateField( time_field ) 
        ptshp_layer.CreateField( ogr.FieldDefn( 'elev', ogr.OFTReal ) )      
        ptshp_layer.CreateField( ogr.FieldDefn( 'cum2dist', ogr.OFTReal ) )
        ptshp_layer.CreateField( ogr.FieldDefn( 'cum3dist', ogr.OFTReal ) )        
        ptshp_layer.CreateField( ogr.FieldDefn( 'slope', ogr.OFTReal ) )        
        
        ptshp_featureDefn = ptshp_layer.GetLayerDefn()
        
        # loops through output records                       
        for ndx, rec in enumerate( gpx_parsed_results ):
            
            pt_feature = ogr.Feature( ptshp_featureDefn )
            
            pt = ogr.Geometry(ogr.wkbPoint)
            pt.SetPoint_2D( 0, rec[1], rec[0] )        
            pt_feature.SetGeometry( pt )
            
            pt_feature.SetField('id', ndx+1 )
            pt_feature.SetField('lon', rec[1] )   
            pt_feature.SetField('lat', rec[0] )
                    
            pt_feature.SetField('time', str( rec[2] ) )
            if rec[3] != '': pt_feature.SetField('elev', rec[3] ) 
            pt_feature.SetField('cum2dist', rec[4] )
            if rec[5] != '': pt_feature.SetField('cum3dist', rec[5] )
            if rec[6] != '': pt_feature.SetField('slope', rec[6] )
                                                              
            ptshp_layer.CreateFeature( pt_feature )
            
            pt_feature.Destroy()
            
        ptshp_datasource.Destroy()
                                
       
    def write_GPX_3D_lnshp( self, gpx_parsed_results ):
     
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
        
        try:    
            lnshp_datasource = shape_driver.CreateDataSource( unicode( fileName ) )
        except TypeError:
            lnshp_datasource = shape_driver.CreateDataSource( str( fileName ) )
                        
        if lnshp_datasource is None:
            QMessageBox.critical( self, "Saving results", "Creation of %s shapefile failed" % os.path.split( fileName )[1] )
            return         
        
        lnshp_layer = lnshp_datasource.CreateLayer( 'profile', geom_type=ogr.wkbLineString25D )
        if lnshp_layer is None:
            QMessageBox.critical( self, "Saving results", "Layer creation failed" )
            return     
 
        # creates required fields         
        lnshp_layer.CreateField( ogr.FieldDefn( 'id', ogr.OFTInteger ) )
        time_beg_field = ogr.FieldDefn( 'time_beg', ogr.OFTString )
        time_beg_field.SetWidth(20)  
        lnshp_layer.CreateField( time_beg_field )
        time_end_field = ogr.FieldDefn( 'time_end', ogr.OFTString )
        time_end_field.SetWidth(20)  
        lnshp_layer.CreateField( time_end_field )     
        lnshp_layer.CreateField( ogr.FieldDefn( 'cum2dist', ogr.OFTReal ) )
        lnshp_layer.CreateField( ogr.FieldDefn( 'cum3dist', ogr.OFTReal ) )
        lnshp_layer.CreateField( ogr.FieldDefn( 'slope', ogr.OFTReal ) )     
    
        lnshp_featureDefn = lnshp_layer.GetLayerDefn()
        
        # loops through output records                     
        for ndx in range( len( gpx_parsed_results )-1 ):
                
            rec_a = gpx_parsed_results[ndx]
            rec_b = gpx_parsed_results[ndx+1]
                            
            lon0, lat0, z0 = rec_a[1], rec_a[0], rec_a[3]
            lon1, lat1, z1 = rec_b[1], rec_b[0], rec_b[3]
            
            if z0 == '' or z1 == '': continue

            ln_feature = ogr.Feature( lnshp_featureDefn )            
            segment_3d = ogr.CreateGeometryFromWkt('LINESTRING(%f %f %f, %f %f %f)' % (lon0, lat0, z0, lon1, lat1, z1))       
            ln_feature.SetGeometry( segment_3d ) 
            
            ln_feature.SetField('id', ndx+1 )
            ln_feature.SetField('time_beg', str( rec_a[2] ) ) 
            ln_feature.SetField('time_end', str( rec_b[2] ) )           
            ln_feature.SetField('cum2dist', rec_b[4] )  
            if rec_b[5] != '': ln_feature.SetField( 'cum3dist', rec_b[5] )
            if rec_b[6] != '': ln_feature.SetField( 'slope', rec_b[6] )  
                            
            lnshp_layer.CreateFeature( ln_feature )
            
            ln_feature.Destroy()
            
        lnshp_datasource.Destroy()        


    def update_point_layers_boxes( self ):
        
        if len(self.pointLayers) == 0:
            return
        
        shape_qgis_ndx = self.prj_struct_point_comboBox.currentIndex() - 1 # minus 1 to account for initial text in combo box
        if shape_qgis_ndx < 0: 
            return
        
        layer = self.pointLayers[ shape_qgis_ndx ]
        fields = layer.dataProvider().fields()     
        field_names = [ field.name() for field in fields.toList()]
                
        for ndx, combobox in enumerate( self.flds_prj_point_comboBoxes ):
            combobox.clear()
            if ndx == 0:
                combobox.addItems( ["none"])
            combobox.addItems( field_names )


    def update_line_layers_boxes( self ):
        
        if len(self.current_line_layers) == 0:
            return
        
        shape_qgis_ndx = self.prj_input_line_comboBox.currentIndex() - 1 # minus 1 to account for initial text in combo box
        if shape_qgis_ndx < 0: 
            return
        
        layer = self.current_line_layers[ shape_qgis_ndx ]
        fields = layer.dataProvider().fields()     
        field_names = [ field.name() for field in fields.toList()]
                
        for ndx, combobox in enumerate( self.flds_prj_line_comboBoxes ):
            combobox.clear()
            if ndx == 0:
                combobox.addItems( ["none"])
            combobox.addItems( field_names )
            
            
    def get_current_combobox_values( self, combobox_list ):
        
        return [ combobox.currentText() for combobox in combobox_list ]


    def define_plot_structural_segment( self, structural_attitude, profile_length ):
        
        slope_radians = structural_attitude.slope_rad
        intersection_downward_sense = structural_attitude.dwnwrd_sense
        intersection_point = structural_attitude.pt_3d
        horiz_distance = structural_attitude.sign_hor_dist        
        
        segment_horiz_scale_factor = 50.0   
        segment_emilength = profile_length / segment_horiz_scale_factor
        
        delta_height = segment_emilength * sin( float( slope_radians ) )
        delta_distance = segment_emilength * cos( float( slope_radians ) )
        
        z0 = intersection_point._z

        structural_segment_s = [ horiz_distance - delta_distance, horiz_distance + delta_distance ]
        structural_segment_z = [ z0 + delta_height, z0 - delta_height ]
                    
        if intersection_downward_sense == "left":
            structural_segment_z = [ z0 - delta_height, z0 + delta_height ]       
        
        return structural_segment_s, structural_segment_z
        

    def get_z_from_dem( self, struct_pts_2d, demObj ):      
        
        z_list = []
        for point_2d in struct_pts_2d:
            interp_z = self.interpolate_point_z( demObj.layer, demObj.params, point_2d )
            z_list.append( interp_z )
            
        return z_list
            

    def calculate_pts_in_projection( self, pts_in_orig_crs, srcCrs, destCrs ):

        pts_in_prj_crs = []
        for pt in pts_in_orig_crs:
            qgs_pt = qgs_point(pt._x,pt._y)
            qgs_pt_prj_crs = project_qgs_point( qgs_pt, srcCrs, destCrs )
            pts_in_prj_crs.append(  Point_3D( qgs_pt_prj_crs.x(), qgs_pt_prj_crs.y() ) )        
        return pts_in_prj_crs
        

    def calculate_projected_3d_pts( self, struct_pts, structural_pts_crs, demObj ):

        demCrs = demObj.params.crs
                        
        # check if on-the-fly-projection is set on
        on_the_fly_projection, project_crs = self.get_on_the_fly_projection()

        # set points in the project crs                    
        if on_the_fly_projection and structural_pts_crs != project_crs:
            struct_pts_in_prj_crs = self.calculate_pts_in_projection( struct_pts, structural_pts_crs, project_crs )
        else:
            struct_pts_in_prj_crs = copy.deepcopy(struct_pts)    
        
        # project the source points from point layer crs to DEM crs
        # if the two crs are different       
        if structural_pts_crs != demCrs:
            struct_pts_in_dem_crs = self.calculate_pts_in_projection( struct_pts, structural_pts_crs, demCrs )
        else:
            struct_pts_in_dem_crs = copy.deepcopy(struct_pts)    
            
        # - 3D structural points, with x, y, and z extracted from the current DEM
        struct_pts_z = self.get_z_from_dem( struct_pts_in_dem_crs, demObj )
        
        assert len(struct_pts_in_prj_crs) == len(struct_pts_z)
        
        return [ Point_3D(pt._x,pt._y,z) for (pt,z) in zip(struct_pts_in_prj_crs, struct_pts_z)] 
       

    def calculate_section_data( self ):
        
        sect_pt_1, sect_pt_2 = self.profile_line._pts
        
        section_init_pt = Point_3D( sect_pt_1._x, sect_pt_1._y, 0.0 )
        section_final_pt = Point_3D( sect_pt_2._x, sect_pt_2._y, 0.0 )

        section_final_pt_up = Point_3D( section_final_pt._x, section_final_pt._y, 1000.0 ) # arbitrary point on the same vertical as sect_pt_2    
        section_cartes_plane = cartes_plane_from_points(section_init_pt, section_final_pt, section_final_pt_up)    
        section_vector = Segment_3D( section_init_pt, section_final_pt ).to_vector( )
        
        return { 'init_pt': section_init_pt, 'cartes_plane': section_cartes_plane, 'vector': section_vector }
                             

    def get_mapping_method(self):
        
        if self.nearest_point_proj_choice.isChecked ():
            return { 'method': 'nearest' }
        
        if self.axis_common_point_proj_choice.isChecked ():
            return { 'method': 'common axis', 
                     'trend': float( self.common_axis_point_trend_SpinBox.value() ),
                     'plunge': float( self.common_axis_point_plunge_SpinBox.value() )}
        
        if self.axis_individual_point_proj_choice.isChecked ():
            return { 'method': 'individual axes', 
                     'trend field': unicode( self.proj_point_indivax_trend_fld_comboBox.currentText() ),
                     'plunge field': unicode( self.proj_point_indivax_plunge_fld_comboBox.currentText()  )}          



    def check_struct_point_proj_parameters(self):
        
        # check if profile exists
        if self.profile_line is None:                         
            return False, "Profile not calculated"
        
        # check that section is made up of only two points
        if self.profile_line.num_points() != 2:                   
            return False, "Profile not made up by only two points"
                        
        # dem number
        if len( self.profiles.topo_profiles ) > 1:           
            return False, "One DEM (and only one DEM) has to be in the profile section" 

        # get point structural layer with parameter fields
        prj_struct_point_qgis_ndx = self.prj_struct_point_comboBox.currentIndex() - 1 # minus 1 to account for initial text in combo box
        if prj_struct_point_qgis_ndx < 0:            
            return False, "No defined point layer for structural data"
        
        return True, "OK"
    
                                                                         
    def create_struct_point_projection( self ):

        parameters_check_ok, parameters_check_msg = self.check_struct_point_proj_parameters()
        if not parameters_check_ok:
            QMessageBox.critical( self, 
                                  "Error", 
                                  parameters_check_msg )            
            return            

        # get color for projected points
        color = self.proj_point_color_comboBox.currentText()

        # define structural layer 
        prj_struct_point_qgis_ndx = self.prj_struct_point_comboBox.currentIndex() - 1 # minus 1 to account for initial text in combo box       
        structural_layer = self.pointLayers[ prj_struct_point_qgis_ndx ]
        structural_layer_crs = structural_layer.crs()
        structural_field_list = self.get_current_combobox_values( self.flds_prj_point_comboBoxes ) 
                  
        # retrieve selected structural points with their attributes        
        structural_pts_attrs = pt_geoms_attrs( structural_layer, structural_field_list ) 
                             
        # list of structural points with original crs
        struct_pts_in_orig_crs = [ Point_3D( rec[0], rec[1] ) for rec in structural_pts_attrs ]
        
        # IDs of structural points
        struct_pts_ids = [ rec[2] for rec in structural_pts_attrs ]
        
        # - geological planes (3D), as geological planes
        try:
            structural_planes = [ GeolPlane( rec[3], rec[4] ) for rec in structural_pts_attrs ]
        except:
            QMessageBox.critical( self, "Error", "Check chosen fields for possible errors" )            
            return
        
        struct_pts_3d = self.calculate_projected_3d_pts(struct_pts_in_orig_crs, 
                                                        structural_layer_crs, 
                                                        self.profiles.profile_dems[0] )

        # - zip together the point value data sets                     
        assert len( struct_pts_3d ) == len( structural_planes )
        structural_data = zip( struct_pts_3d, structural_planes, struct_pts_ids )   
               
        ### map points onto section ###
        
        # calculation of Cartesian plane expressing section plane        
        self.section_data = self.calculate_section_data( )
        
        # calculation of projected structural points
        
        # get chosen mapping method
        mapping_method = self.get_mapping_method()
        if mapping_method['method'] == 'individual axes':
            trend_field_name, plunge_field_name = mapping_method['trend field'], mapping_method['plunge field']
            # retrieve structural points mapping axes        
            mapping_method['individual_axes_values'] = vect_attrs( structural_layer, [trend_field_name, plunge_field_name]) 

        self.profiles.add_plane_attitudes( map_struct_pts_on_section( structural_data, self.section_data, mapping_method  ) )
        self.plane_attitudes_colors.append( color )
        ### plot structural points in section ###
        self.plot_profile_elements()


    def reset_struct_point_projection(self):
        
        try:
            self.profiles.plane_attitudes = []
            self.plane_attitudes_colors = []
        except:
            pass
         


    def check_structural_line_projection_inputs( self ):

        # dem parameters
        try:
            num_dems_in_profile = len( self.profiles.topo_profiles )
        except:
            return False, "Profile has not been calculated"
        else:
            if num_dems_in_profile == 0:
                return False, "Profile has not been calculated"
            elif num_dems_in_profile > 1:
                return False, "One DEM (and only one DEM) has to be used in the profile section"
                        
        # check if profile exists
        if self.profile_line is None:                          
            return False, "Profile has not been calculated"
        
        # check that section is made up of only two points
        if self.profile_line.num_points() != 2:                   
            return False, "Current profile is not made up by only two points"

        # line structural layer with parameter fields
        prj_struct_line_qgis_ndx = self.prj_input_line_comboBox.currentIndex() - 1 # minus 1 to account for initial text in combo box
        if prj_struct_line_qgis_ndx < 0:            
            return False, "No defined point layer for structural data"

        try:
            densify_distance = float( self.project_line_densify_distance_lineedit.text() )
        except:
            return False, "No numeric value for densify line distance"
        else:
            if densify_distance == 0.0:
                return False, "Densify line distance cannot be zero"                
                        
        return True, "OK"

                                                     
    def create_struct_line_projection(self):

        input_params_valid = self.check_structural_line_projection_inputs()
        if not input_params_valid[0]:         
            QMessageBox.critical( self, 
                                  "Error", 
                                  input_params_valid[1] )            
            return
               
        assert len( self.profiles.profile_dems ) == 1
        demLayer = self.profiles.profile_dems[0].layer      
        demParams = self.profiles.profile_dems[0].params

        # get line structural layer
        prj_struct_line_qgis_ndx = self.prj_input_line_comboBox.currentIndex() - 1 # minus 1 to account for initial text in combo box
        
        # get id field
        prj_struct_line_id_field_ndx = self.id_fld_line_prj_comboBox.currentIndex() - 1 # minus 1 to account for initial text in combo box
        
        # get color 
        color = self.project_line_color_comboBox.currentText()
        
        # define structural layer        
        structural_line_layer = self.current_line_layers[ prj_struct_line_qgis_ndx ]
        structural_line_layer_crs = structural_line_layer.crs()
                   
        # read structural line values
        id_list = field_values( structural_line_layer, prj_struct_line_id_field_ndx ) 
        line_orig_crs_geom_data = line_geoms_attrs( structural_line_layer)
        assert len( id_list ) == len( line_orig_crs_geom_data )
        line_orig_geom_list3 = [ geom_data[0] for geom_data in line_orig_crs_geom_data ]
        line_orig_crs_MultiLine_2D_list = [ xytuple_list2_to_MultiLine_2D( xy_list2 ) for xy_list2 in line_orig_geom_list3 ]
        line_orig_crs_clean_MultiLine_2D_list = [ multiline_2d.remove_coincident_points() for multiline_2d in line_orig_crs_MultiLine_2D_list ]

        # get project CRS information
        on_the_fly_projection, project_crs = self.get_on_the_fly_projection() 
                
        # project input line layer to project CRS
        if on_the_fly_projection:
            line_proj_crs_MultiLine_2D_list = [ multiline2d.project_crs( structural_line_layer_crs, project_crs ) for multiline2d in line_orig_crs_clean_MultiLine_2D_list ]
        else:
            line_proj_crs_MultiLine_2D_list = line_orig_crs_clean_MultiLine_2D_list
               
        # densify with provided distance
        densify_proj_crs_distance = float( self.project_line_densify_distance_lineedit.text() )
        densified_proj_crs_MultiLine_2D_list = [ multiline_2d.densify( densify_proj_crs_distance ) for multiline_2d in line_proj_crs_MultiLine_2D_list ]
                    
        # project to Dem CRS
        if on_the_fly_projection and demParams.crs != project_crs:
            densified_dem_crs_MultiLine_2D_list = [ multiline_2d.project_crs( project_crs, demParams.crs ) for multiline_2d in densified_proj_crs_MultiLine_2D_list ]
        else:
            densified_dem_crs_MultiLine_2D_list = densified_proj_crs_MultiLine_2D_list
        
        # interpolate z values from Dem
        z_list = [ self.interpolate_point_z( demLayer, demParams, pt_2d ) for multiline_2d in densified_dem_crs_MultiLine_2D_list for line_2d in multiline_2d._lines for pt_2d in line_2d._pts ]

        # extract x-y pairs for creation of 3D points 
        xy_list = [ ( pt_2d._x, pt_2d._y ) for multiline_2d in densified_proj_crs_MultiLine_2D_list for line_2d in multiline_2d._lines for pt_2d in line_2d._pts ]
        
        # debug: verify length of two lists             
        assert len( z_list ) == len ( xy_list )  
                    
        # replicate MultiLine list structure with 3D points with project CRS
        ndx = -1
        multiline_3d_proj_crs_list = []
        for multiline_2d in densified_proj_crs_MultiLine_2D_list:
            multiline_3d_list = []
            for line_2d in multiline_2d._lines:
                line_3d_pts_list = []
                for pt_2d in line_2d._pts:
                    ndx += 1
                    line_3d_pts_list.append( Point_3D( xy_list[ndx][0], xy_list[ndx][1], z_list[ndx]) )
                multiline_3d_list.append( Line_3D( line_3d_pts_list ) )
            multiline_3d_proj_crs_list.append( MultiLine_3D( multiline_3d_list ) )

        # create projection vector        
        trend = float( self.common_axis_line_trend_SpinBox.value() )
        plunge = float( self.common_axis_line_plunge_SpinBox.value() )
        axis_versor = GeolAxis( trend, plunge ).to_versor()
        l, m, n = axis_versor._x, axis_versor._y, axis_versor._z
        
        # calculation of Cartesian plane expressing section plane        
        self.section_data = self.calculate_section_data( )
                
        # project MultiLine_3D points to section
        intersection_point_list = []
        for multiline_3d in multiline_3d_proj_crs_list:
            for line_3d in multiline_3d._lines:
                for pt_3d in line_3d._pts:
                    srcPt = pt_3d
                    param_line = ParamLine( srcPt, l, m, n )
                    intersection_point_list.append( param_line.intersect_cartes_plane( self.section_data['cartes_plane'] )  )
                                     
        # replicate MultiLine list structure with 3D points with project CRS
        ndx = -1
        multiline_3d_proj_crs_section_list = []
        for multiline_3d in multiline_3d_proj_crs_list:
            multiline_3d_list = []
            for line_3d in multiline_3d._lines:
                line_3d_pts_list = []
                for pt_3d in line_3d._pts:
                    ndx += 1
                    line_3d_pts_list.append( intersection_point_list[ndx] ) 
                multiline_3d_list.append( Line_3D( line_3d_pts_list ) )
            multiline_3d_proj_crs_section_list.append( MultiLine_3D( multiline_3d_list ) )
        

        section_start_point, section_vector = self.section_data['init_pt'], self.section_data['vector']
        curves_2d_list = []
        for multiline_3d in multiline_3d_proj_crs_section_list:
            multiline_2d_list = []
            for line_3d in multiline_3d._lines:
                line_2d_pts_list = []
                for pt_3d in line_3d._pts:
                    s = calculate_distance_with_sign( pt_3d, section_start_point, section_vector )
                    z = pt_3d._z
                    line_2d_pts_list.append( Point_2D( s, z ) )
                multiline_2d_list.append( Line_2D( line_2d_pts_list ))
            curves_2d_list.append( MultiLine_2D( multiline_2d_list ) )
         
        assert len( curves_2d_list ) == len( line_orig_geom_list3 )
                       
        self.profiles.add_curves( curves_2d_list, id_list )
        self.curve_colors.append( color )
                                       
        # plot new cross section
        self.plot_profile_elements()


    def reset_structural_lines_projection(self):
        
        try:
            self.profiles.curves = []
            self.profiles.curves_ids = []
            self.curve_colors = []
        except:
            pass
        
                                         
    def save_proj_points_results(self):
        
        try:
            num_plane_attitudes_sets = len( self.profiles.plane_attitudes )
        except:
            QMessageBox.critical( self, "Saving results", "No available geological attitudes" )
            return 
        else:
            if num_plane_attitudes_sets == 0:
                QMessageBox.critical( self, "Saving results", "No available geological attitudes" )
                return                              
    
        if not ( self.save_proj_point_results_asCSV_checkbox.isChecked() or \
                 self.save_proj_point_results_asPtShp_checkbox.isChecked()):
            QMessageBox.critical( self, "Saving results", "No output format is selected" )
            return  
        
        # definition of field names
        header_list = ['id', 
                       'or_strpt_x',
                       'or_strpt_y',
                       'or_strpt_z', 
                       'prj_strpt_x',
                       'prj_strpt_y',
                       'prj_strpt_z',
                       's',
                       'or_dipdir',
                       'or_dipangle',
                       'trc_dipangle',
                       'trc_dipdir']
                                               
        parsed_crosssect_results = self.parse_crosssect_results_for_export( self.profiles.plane_attitudes )

        # output for csv file
        if self.save_proj_point_results_asCSV_checkbox.isChecked():
            self.write_results_as_csv( header_list, parsed_crosssect_results )

        # output for 3D pt shapefile            
        if self.save_proj_point_results_asPtShp_checkbox.isChecked():
            self.write_crosssect_result_ptshp( header_list, parsed_crosssect_results )
                    
        QMessageBox.information( self, "Saving projected attitudes", "Completed" )
    
    
    def parse_crosssect_results_for_export( self, plane_attitudes_datasets ):

        result_data = []  
              
        for dataset in plane_attitudes_datasets:  
                      
            for plane_attitude_rec in dataset:
                
                pt_id = plane_attitude_rec.id
                or_pt_x = plane_attitude_rec.src_pt_3d._x
                or_pt_y = plane_attitude_rec.src_pt_3d._y            
                or_pt_z = plane_attitude_rec.src_pt_3d._z            
                pr_pt_x = plane_attitude_rec.pt_3d._x
                pr_pt_y = plane_attitude_rec.pt_3d._y            
                pr_pt_z = plane_attitude_rec.pt_3d._z            
                s = plane_attitude_rec.sign_hor_dist
                or_dipdir = plane_attitude_rec.src_geol_plane._dipdir
                or_dipangle = plane_attitude_rec.src_geol_plane._dipangle
                tr_dipangle = degrees( plane_attitude_rec.slope_rad )
                tr_dipdir = plane_attitude_rec.dwnwrd_sense
                
                record = [ pt_id, or_pt_x, or_pt_y, or_pt_z, pr_pt_x, pr_pt_y, pr_pt_z, s, or_dipdir, or_dipangle, tr_dipangle, tr_dipdir ]

                result_data.append( record )
         
        return result_data


    def parse_curves_for_export(self):
        
        data_list = []
        for curve_set, id_set in zip( self.profiles.curves, self.profiles.curves_ids ):
            for curve, id in zip( curve_set, id_set):
                for line in curve._lines:
                    for pt in line._pts:
                        data_list.append( [ id, pt._x, pt._y ] )
        return data_list
     
        
    def save_proj_lines_results(self):
        
        try:
            num_proj_lines_sets = len( self.profiles.curves )
        except:
            QMessageBox.critical( self, "Saving results", "No available geological traces" )
            return 
        else:
            if num_proj_lines_sets == 0:
                QMessageBox.critical( self, "Saving results", "No available geological traces to save" )
                return 
    
        parsed_curves_for_export = self.parse_curves_for_export()
        header_list = ['id', 's', 'z']
        
        self.write_results_as_csv( header_list, parsed_curves_for_export )
               
               
    def write_results_as_csv( self, header_list, parsed_crosssect_results ):
        
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
            for rec in parsed_crosssect_results:
                out_rec_string = ''
                for val in rec:
                    out_rec_string += str( val ) + ','
                f.write( out_rec_string[:-1]+'\n' )
    

    def write_crosssect_result_ptshp( self, header_list, parsed_crosssect_results ):

        fileName = QFileDialog.getSaveFileName(self,
                                               self.tr("Save results as a 3D point shapefile"),
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

        try:
            ptshp_datasource = shape_driver.CreateDataSource( unicode( fileName ) )
        except TypeError:
            ptshp_datasource = shape_driver.CreateDataSource( str( fileName ) )
            
        if ptshp_datasource is None:
            QMessageBox.critical( self, "Saving results", "Creation of %s shapefile failed" % os.path.split( fileName )[1] )
            return

        ptshp_layer = ptshp_datasource.CreateLayer( 'profile', geom_type=ogr.wkbPoint25D )
        if ptshp_layer is None:
            QMessageBox.critical( self, "Saving results", "Layer creation failed" )
            return

        # creates required fields
        ptshp_layer.CreateField( ogr.FieldDefn( 'id', ogr.OFTString ) )
        ptshp_layer.CreateField( ogr.FieldDefn( 'or_pt_x', ogr.OFTReal ) )
        ptshp_layer.CreateField( ogr.FieldDefn( 'or_pt_y', ogr.OFTReal ) )
        ptshp_layer.CreateField( ogr.FieldDefn( 'or_pt_z', ogr.OFTReal ) )
        ptshp_layer.CreateField( ogr.FieldDefn( 'prj_pt_x', ogr.OFTReal ) )
        ptshp_layer.CreateField( ogr.FieldDefn( 'prj_pt_y', ogr.OFTReal ) )
        ptshp_layer.CreateField( ogr.FieldDefn( 'prj_pt_z', ogr.OFTReal ) )
        ptshp_layer.CreateField( ogr.FieldDefn( 's', ogr.OFTReal ) )
        ptshp_layer.CreateField( ogr.FieldDefn( 'or_dpdir', ogr.OFTReal ) )
        ptshp_layer.CreateField( ogr.FieldDefn( 'or_dpang', ogr.OFTReal ) )
        ptshp_layer.CreateField( ogr.FieldDefn( 'tr_dpang', ogr.OFTReal ) )
        ptshp_layer.CreateField( ogr.FieldDefn( 'tr_dpdir', ogr.OFTString ) )

        ptshp_featureDefn = ptshp_layer.GetLayerDefn()

        # loops through output records
        for rec in parsed_crosssect_results:

            pt_id, or_pt_x, or_pt_y, or_pt_z, pr_pt_x, pr_pt_y, pr_pt_z, s, or_dipdir, or_dipangle, tr_dipangle, tr_dipdir = rec

            pt_feature = ogr.Feature( ptshp_featureDefn )

            pt = ogr.Geometry( ogr.wkbPoint25D )
            pt.SetPoint( 0, pr_pt_x, pr_pt_y, pr_pt_z )
            pt_feature.SetGeometry( pt )

            pt_feature.SetField( 'id', str( pt_id ) )
            pt_feature.SetField( 'or_pt_x', or_pt_x )
            pt_feature.SetField( 'or_pt_y', or_pt_y )
            pt_feature.SetField( 'or_pt_z', or_pt_z )
            pt_feature.SetField( 'prj_pt_x', pr_pt_x )
            pt_feature.SetField( 'prj_pt_y', pr_pt_y )
            pt_feature.SetField( 'prj_pt_z', pr_pt_z )
            pt_feature.SetField( 's', s )
            pt_feature.SetField( 'or_dpdir',or_dipdir )
            pt_feature.SetField( 'or_dpang', or_dipangle )
            pt_feature.SetField( 'tr_dpang', tr_dipangle )
            pt_feature.SetField( 'tr_dpdir', str( tr_dipdir ) )

            ptshp_layer.CreateFeature( pt_feature )

            pt_feature.Destroy()

        ptshp_datasource.Destroy()

    
    def plot_profile_elements( self, elev_type="DEM" ):
                               
        # defines the extent for the plot window: s min and max     
        plot_s_min, plot_s_max = 0, self.profiles.get_max_s() 

        # defines z min and max values
        profile_z_min, profile_z_max = self.profiles.get_min_z(), self.profiles.get_max_z()
        delta_z = profile_z_max - profile_z_min 
        plot_z_min, plot_z_max = profile_z_min - delta_z * 0.05, profile_z_max + delta_z * 0.05

        # defines slope min and max values
        slope_list = [ topo_profile.profile_3d.slopes_list() for topo_profile in self.profiles.topo_profiles ]
        profiles_slope_min, profiles_slope_max = min( [ min(slist) for slist in slope_list ] ), max( [ max(slist) for slist in slope_list ] )
        delta_slope = profiles_slope_max - profiles_slope_min 
        plot_slope_min, plot_slope_max = profiles_slope_min - delta_slope*0.2, profiles_slope_max + delta_slope*0.2 

        # map
        profile_window = MplWidget()  

        if elev_type == 'DEM':
            plot_height_choice = self.DEM_plot_height_checkbox.isChecked()
            plot_slope_choice = self.DEM_plot_slope_checkbox.isChecked() 
        elif elev_type == 'GPX':
            plot_height_choice = self.GPX_plot_height_checkbox.isChecked()
            plot_slope_choice = self.GPX_plot_slope_checkbox.isChecked() 

        if plot_height_choice and plot_slope_choice:
            mpl_code_list = [ 211, 212 ]
        else:
            mpl_code_list = [ 111 ]            
        subplot_code = mpl_code_list[0] 
  
        if plot_height_choice:
            
            self.axes_elevation = self.plot_topo_profile_lines( subplot_code, 
                                                                  profile_window, 
                                                                  self.profiles.topo_profiles, 
                                                                  'elevation', 
                                                                  (plot_s_min, plot_s_max), 
                                                                  (plot_z_min, plot_z_max),
                                                                  self.selected_dem_colors,
                                                                  self.DEM_plot_height_filled_checkbox.isChecked() )
            
            if self.DEM_exageration_1_1_checkbox.isChecked():
                self.axes_elevation.set_aspect('equal')
            
        if plot_slope_choice:
                        
            if len(mpl_code_list) == 2: subplot_code = mpl_code_list[1]            
            self.axes_slopes = self.plot_topo_profile_lines( subplot_code, 
                                                              profile_window, 
                                                              self.profiles.topo_profiles, 
                                                              'slope', 
                                                              (plot_s_min, plot_s_max), 
                                                              (plot_slope_min, plot_slope_max), 
                                                              self.selected_dem_colors,
                                                              self.DEM_plot_slope_filled_checkbox.isChecked() )
                        
        if len( self.profiles.plane_attitudes ) > 0: 
                       
            for plane_attitude_set, color in zip( self.profiles.plane_attitudes, self.plane_attitudes_colors):                
                self.plot_structural_attitude( self.axes_elevation, plot_s_max, plane_attitude_set, color )                   
                   
        if len( self.profiles.curves ) > 0: 

            for curve_set, color in zip( self.profiles.curves, self.curve_colors) :                
                self.plot_curve_set( self.axes_elevation, curve_set, color )          
                     
        profile_window.canvas.draw() 
        
        self.profile_windows.append( profile_window )


    def plot_topo_profile_lines(self, subplot_code, profile_window, topo_profiles, topo_type, plot_x_range, plot_y_range, dem_colors, filled_choice ):
        
        axes = self.create_axes( subplot_code,
                                  profile_window, 
                                  plot_x_range, 
                                  plot_y_range  )

        # label = unicode(dem_name)
        for topo_profile, dem_color in zip( topo_profiles, dem_colors ): 
            
            if topo_type == 'elevation':
                y_list = topo_profile.z_list()
                plot_y_min = plot_y_range[0]
            elif topo_type == 'slope':
                y_list = topo_profile.slope_list()
                plot_y_min = 0.0
            
            if filled_choice:    
                plot_filled_line( axes,
                                topo_profile.get_increm_dist_2d(), 
                                y_list, 
                                plot_y_min, 
                                dem_color  )
                
            plot_line( axes,
                        topo_profile.get_increm_dist_2d(), 
                        y_list, 
                        dem_color ) 

        return axes
    
    
    def create_axes(self, subplot_code, profile_window, plot_x_range, plot_y_range ):

            x_min, x_max = plot_x_range
            y_min, y_max = plot_y_range
            axes = profile_window.canvas.fig.add_subplot( subplot_code )
            axes.set_xlim( x_min, x_max )
            axes.set_ylim( y_min, y_max )

            axes.grid(True)
                       
            return axes
                
        
    def plot_structural_attitude( self, axes, section_length, structural_attitude_list, color ):
        
        # TODO:  manage case for possible nan z values
        projected_z = [ structural_attitude.pt_3d._z for structural_attitude in structural_attitude_list if  0.0 <= structural_attitude.sign_hor_dist <= section_length ]
                 
        # TODO:  manage case for possible nan z values
        projected_s = [ structural_attitude.sign_hor_dist for structural_attitude in structural_attitude_list if 0.0 <= structural_attitude.sign_hor_dist <= section_length ]
        
        projected_ids = [ structural_attitude.id for structural_attitude in structural_attitude_list if 0.0 <= structural_attitude.sign_hor_dist <= section_length ]
        
                
        axes.plot( projected_s, projected_z,'o', color=color ) 
        
        # plot segments representing structural data       
        for structural_attitude in structural_attitude_list:
            if 0.0 <= structural_attitude.sign_hor_dist <= section_length:
                structural_segment_s, structural_segment_z = self.define_plot_structural_segment( structural_attitude, section_length )            
                
                axes.plot( structural_segment_s, structural_segment_z,'-', color=color )
        
        if self.plot_prj_add_trendplunge_label.isChecked() or self.plot_prj_add_pt_id_label.isChecked():
            
            src_dip_dirs = [ structural_attitude.src_geol_plane._dipdir for structural_attitude in structural_attitude_list if 0.0 <= structural_attitude.sign_hor_dist <= section_length ]
            src_dip_angs = [ structural_attitude.src_geol_plane._dipangle for structural_attitude in structural_attitude_list if 0.0 <= structural_attitude.sign_hor_dist <= section_length ]
        
            for rec_id, src_dip_dir, src_dip_ang, s, z in zip( projected_ids, src_dip_dirs, src_dip_angs, projected_s, projected_z):
                
                if self.plot_prj_add_trendplunge_label.isChecked() and self.plot_prj_add_pt_id_label.isChecked():
                    label = "%s-%03d/%02d" % ( rec_id, src_dip_dir, src_dip_ang )
                elif self.plot_prj_add_pt_id_label.isChecked():
                    label = "%s" % rec_id
                elif self.plot_prj_add_trendplunge_label.isChecked():
                    label = "%03d/%02d" % ( src_dip_dir, src_dip_ang )
                    
                axes.annotate( label, (s+15,z+15) )                    
            

    def plot_curve_set(self, axes, curve_set, color ):        
        
        for multiline_2d in curve_set:
                for line_2d in multiline_2d._lines:
                    plot_line( axes, line_2d.x_list(), line_2d.y_list(), color )
                    
                    
    def closeEvent( self, event ):
        
        try:
            self.rubberband.reset( QGis.Line )
        except:
            pass
        
        try:
            self.disconnect_digitize_maptool()
        except:
            pass
 
        try:       
            QgsMapLayerRegistry.instance().layerWasAdded.disconnect( self.refresh_struct_point_lyr_combobox ) 
        except:
            pass
                    
        try:       
            QgsMapLayerRegistry.instance().layerWasAdded.disconnect( self.refresh_struct_line_lyr_combobox )
        except:
            pass
           
        try:
            QgsMapLayerRegistry.instance().layerRemoved.disconnect( self.refresh_struct_point_lyr_combobox ) 
        except:
            pass               
            
        try:
            QgsMapLayerRegistry.instance().layerRemoved.disconnect( self.refresh_struct_line_lyr_combobox )
        except:
            pass
     

        
class SourceDEMsDialog( QDialog ):
    
    def __init__(self, raster_layers, parent=None):
        
        super( SourceDEMsDialog, self ).__init__(parent)        
        
        self.singleband_raster_layers_in_project = raster_layers

                                       
        self.listDEMs_treeWidget = QTreeWidget()
        self.listDEMs_treeWidget.setColumnCount( 2 )
        self.listDEMs_treeWidget.setColumnWidth ( 0, 200 )
        self.listDEMs_treeWidget.headerItem().setText( 0, "Name" )
        self.listDEMs_treeWidget.headerItem().setText( 1, "Plot color" )
        self.listDEMs_treeWidget.setHorizontalScrollBarPolicy(Qt.ScrollBarAsNeeded)
        self.listDEMs_treeWidget.setDragEnabled(False)
        self.listDEMs_treeWidget.setDragDropMode(QAbstractItemView.NoDragDrop)
        self.listDEMs_treeWidget.setAlternatingRowColors(True)
        self.listDEMs_treeWidget.setSelectionMode(QAbstractItemView.SingleSelection)
        self.listDEMs_treeWidget.setTextElideMode( Qt.ElideLeft )
         
        self.refresh_raster_layer_treewidget()
        
        okButton = QPushButton("&OK")
        cancelButton = QPushButton("Cancel")

        buttonLayout = QHBoxLayout()
        buttonLayout.addStretch()
        buttonLayout.addWidget(okButton)
        buttonLayout.addWidget(cancelButton)
        
        layout = QGridLayout()

        layout.addWidget( self.listDEMs_treeWidget, 0, 0, 1, 3 )                 
        layout.addLayout( buttonLayout, 1, 0, 1, 3 )
        
        self.setLayout(layout)

        self.connect(okButton, SIGNAL("clicked()"),
                     self,  SLOT("accept()") )
        self.connect(cancelButton, SIGNAL("clicked()"),
                     self, SLOT("reject()"))
        
        self.setWindowTitle("Define source DEMs")
         

    def refresh_raster_layer_treewidget( self ):

        self.listDEMs_treeWidget.clear() 
                                    
        for raster_layer in self.singleband_raster_layers_in_project:
            
            tree_item = QTreeWidgetItem( self.listDEMs_treeWidget )
            tree_item.setText(0, raster_layer.name() )
            combo_box = QComboBox()
            combo_box.setSizeAdjustPolicy ( 0 )
            combo_box.addItems( qprof_QWidget.colors)
            self.listDEMs_treeWidget.setItemWidget( tree_item, 1, combo_box )       
            tree_item.setFlags( tree_item.flags() | Qt.ItemIsUserCheckable )
            tree_item.setCheckState( 0, 0 )


        
class SourceLineLayerDialog( QDialog ):
    
    def __init__(self, current_line_layers, parent=None):
                
        super( SourceLineLayerDialog, self ).__init__(parent)
        
        self.current_line_layers = current_line_layers
 
        layout = QGridLayout()
                                              
        layout.addWidget( QLabel( self.tr("Line layer:") ), 0, 0, 1, 1) 
        self.LineLayers_comboBox = QComboBox()                         
        layout.addWidget(self.LineLayers_comboBox, 0, 1, 1, 3)         
        self.refresh_input_profile_layer_combobox( )

        layout.addWidget( QLabel( self.tr("Line order field:") ), 1, 0, 1, 1) 
        
        self.Trace2D_order_field_comboBox = QComboBox()                        
        layout.addWidget(self.Trace2D_order_field_comboBox, 1, 1, 1, 3) 
                
        self.refresh_order_field_combobox( )
        
        self.LineLayers_comboBox.currentIndexChanged[int].connect (self.refresh_order_field_combobox )
        
        okButton = QPushButton("&OK")
        cancelButton = QPushButton("Cancel")

        buttonLayout = QHBoxLayout()
        buttonLayout.addStretch()
        buttonLayout.addWidget(okButton)
        buttonLayout.addWidget(cancelButton)
              
        layout.addLayout( buttonLayout, 2, 0, 1, 3 )
        
        self.setLayout( layout )

        self.connect(okButton, SIGNAL("clicked()"),
                     self,  SLOT("accept()") )
        self.connect(cancelButton, SIGNAL("clicked()"),
                     self, SLOT("reject()"))
        
        self.setWindowTitle("Define source line layer")


    def refresh_input_profile_layer_combobox( self ):
        
        self.LineLayers_comboBox.clear()
      
        for layer in self.current_line_layers:
            self.LineLayers_comboBox.addItem( layer.name() )        

        shape_qgis_ndx = self.LineLayers_comboBox.currentIndex()
        self.line_shape = self.current_line_layers[ shape_qgis_ndx ]            
         
         
    def refresh_order_field_combobox( self ):       
        
        self.Trace2D_order_field_comboBox.clear()
        self.Trace2D_order_field_comboBox.addItem('--optional--')  

        shape_qgis_ndx = self.LineLayers_comboBox.currentIndex()
        self.line_shape = self.current_line_layers[ shape_qgis_ndx ]
       
        line_layer_field_list = self.line_shape.dataProvider().fields().toList( )        
        for field in line_layer_field_list:
            self.Trace2D_order_field_comboBox.addItem( field.name() )      
            
              

        



        
        
        
        
