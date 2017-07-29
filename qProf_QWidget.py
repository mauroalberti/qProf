# -*- coding: utf-8 -*-

import os
import unicodedata

from qgis.core import QgsGeometry, QgsVectorLayer

from .gsf.geometry import Plane, GPlane, GAxis
from .gsf.array_utils import to_float

from .gis_utils.features import Segment, MultiLine, Line, \
    merge_lines, ParamLine3D, xytuple_list_to_Line
from .gis_utils.intersections import map_struct_pts_on_section, calculate_distance_with_sign
from .gis_utils.profile import Profile_Elements, topoprofiles_from_dems, topoprofiles_from_gpxfile, \
    intersect_with_dem, calculate_profile_lines_intersection, intersection_distances_by_profile_start_list, \
    extract_multiline2d_list, profile_polygon_intersection, calculate_projected_3d_pts
from .gis_utils.qgs_tools import *
from .gis_utils.statistics import get_statistics
from .gis_utils.errors import VectorInputException, VectorIOException

from .qt_utils.filesystem import update_directory_key, new_file_path, old_file_path
from .qt_utils.tools import info, warn, error, update_ComboBox

from .string_utils.utils_string import clean_string

from qProf_plotting import plot_profile_elements
from qProf_export import write_intersection_polygon_lnshp, write_topography_allDEMs_csv, write_topography_singleDEM_csv, \
          write_generic_csv, write_line_csv, write_topography_allDEMs_ptshp, write_topography_allDEMs_lnshp, \
          write_geological_attitudes_ptshp, write_profile_lnshp, write_topography_GPX_lnshp, write_topography_singleDEM_lnshp, \
          write_topography_singleDEM_ptshp, write_topography_GPX_ptshp, write_intersection_line_ptshp


class qprof_QWidget(QWidget):

    colors = ['orange', 'green', 'red', 'grey', 'brown', 'yellow', 'magenta', 'black', 'blue', 'white', 'cyan',
              'chartreuse']

    map_digitations = 0

    def __init__(self, plugin_name, canvas):

        super(qprof_QWidget, self).__init__()

        self.plugin_name = plugin_name
        self.canvas = canvas

        self.current_directory = os.path.dirname(__file__)

        self.settings = QSettings("alberese", self.plugin_name)
        self.settings_gpxdir_key = "gpx/last_used_dir"

        self.choose_message = "choose"

        self.demline_source = "demline"
        self.gpxfile_source = "gpxfile"
        self.digitized_profile_line2dt = None
        self.polygon_classification_colors = None

        self.profile_elements = None

        self.profile_windows = []

        self.plane_attitudes_colors = []
        self.curve_colors = []

        self.setup_gui()

    def setup_gui(self):

        self.dialog_layout = QVBoxLayout()
        self.main_widget = QTabWidget()
        self.main_widget.addTab(self.setup_topoprofile_tab(), "Topography")
        self.main_widget.addTab(self.setup_project_section_tab(), "Projections")
        self.main_widget.addTab(self.setup_intersect_section_tab(), "Intersections")
        self.main_widget.addTab(self.setup_export_section_tab(), "Export")

        self.prj_input_line_comboBox.currentIndexChanged[int].connect(self.update_linepoly_layers_boxes)
        self.inters_input_line_comboBox.currentIndexChanged[int].connect(self.update_linepoly_layers_boxes)
        self.inters_input_polygon_comboBox.currentIndexChanged[int].connect(self.update_linepoly_layers_boxes)

        self.struct_line_refresh_lyr_combobox()
        self.struct_polygon_refresh_lyr_combobox()

        QgsMapLayerRegistry.instance().layerWasAdded.connect(self.struct_point_refresh_lyr_combobox)
        QgsMapLayerRegistry.instance().layerWasAdded.connect(self.struct_line_refresh_lyr_combobox)
        QgsMapLayerRegistry.instance().layerWasAdded.connect(self.struct_polygon_refresh_lyr_combobox)

        QgsMapLayerRegistry.instance().layerRemoved.connect(self.struct_point_refresh_lyr_combobox)
        QgsMapLayerRegistry.instance().layerRemoved.connect(self.struct_line_refresh_lyr_combobox)
        QgsMapLayerRegistry.instance().layerRemoved.connect(self.struct_polygon_refresh_lyr_combobox)

        self.dialog_layout.addWidget(self.main_widget)
        self.setLayout(self.dialog_layout)
        self.adjustSize()
        self.setWindowTitle(self.plugin_name)

    def select_input_gpxFile(self):

        gpx_last_used_dir = self.settings.value(self.settings_gpxdir_key,
                                                "")
        fileName = QFileDialog.getOpenFileName(self,
                                               self.tr("Open GPX file"),
                                               gpx_last_used_dir,
                                               "GPX (*.gpx *.GPX)")
        if not fileName:
            return
        else:
            update_directory_key(self.settings,
                                 self.settings_gpxdir_key,
                                 fileName)
            self.qlneInputGPXFile.setText(fileName)


    def setup_topoprofile_tab(self):

        def read_input_gpxFile(qlneInputGPXFile, qcbxInvertProfile):

            try:
                source_gpx_path = unicode(qlneInputGPXFile.text())
                if source_gpx_path == '':
                    warn(self,
                         self.plugin_name,
                         "Source GPX file is not set")
                    return
            except Exception as e:
                warn(self,
                     self.plugin_name,
                     "Source GPX file not correctly set: {}".format(e.message))
                return

            invert_profile_choice = qcbxInvertProfile.isChecked()

            try:
                topo_profiles = topoprofiles_from_gpxfile(source_gpx_path,
                                                          invert_profile_choice,
                                                          self.gpxfile_source)
            except Exception as e:
                warn(self,
                     self.plugin_name,
                     "Error with profile calculation from GPX file: {}".format(e.message))
                return

            if topo_profiles is None:
                warn(self,
                     self.plugin_name,
                     "Debug: profile not created")
                return

            topo_source_type = self.gpxfile_source


        qwdgTopoProfile = QWidget()
        qlytTopoProfile = QVBoxLayout()

        ## Input data

        qgbxTopoSources = QGroupBox(qwdgTopoProfile)
        qgbxTopoSources.setTitle("Topographic profile sources")

        qlyTopoSources = QVBoxLayout()

        ####

        qtbxGPXInput = QToolBox()
        qwdgGPXInput = QWidget()
        qlytGPXInput = QGridLayout()

        qlytGPXInput.addWidget(QLabel(self.tr("Input GPX file with track points:")), 0, 0, 1, 1)

        self.qlneInputGPXFile = QLineEdit()
        self.qlneInputGPXFile.setPlaceholderText("my_track.gpx")
        qlytGPXInput.addWidget(self.qlneInputGPXFile, 0, 1, 1, 1)

        self.qphbInputGPXFile = QPushButton("...")
        self.qphbInputGPXFile.clicked.connect(self.select_input_gpxFile)
        qlytGPXInput.addWidget(self.qphbInputGPXFile, 0, 0, 1, 2)

        self.qcbxInvertProfile = QCheckBox("Invert profile orientation")
        qlytGPXInput.addWidget(self.qcbxInvertProfile, 1, 0, 1, 1)

        self.qphbGetInputGPXFile = QPushButton("Get")
        self.qphbGetInputGPXFile.clicked.connect(self.read_input_gpxFile)
        qlytGPXInput.addWidget(self.qphbGetInputGPXFile, 1, 2, 1, 1)

        qwdgGPXInput.setLayout(qlytGPXInput)
        qtbxGPXInput.addItem(qwdgGPXInput, "GPX input")
        qlyTopoSources.addWidget(qtbxGPXInput)

        #####

        qtbxDEMLineInput = QToolBox()
        qwdgDEMLineInput = QWidget()
        qlytDEMLineInput = QHBoxLayout()


        qwdgDEMLineInput.setLayout(qlytDEMLineInput)
        qtbxDEMLineInput.addItem(qwdgDEMLineInput, "DEM and line input")
        qlyTopoSources.addWidget(qtbxDEMLineInput)

        """
        self.prof_toposources_fromdems_checkbox = QRadioButton(self.tr("DEM-line"))
        self.prof_toposources_fromdems_checkbox.setChecked(True)
        qlyTopoSources.addWidget(self.prof_toposources_fromdems_checkbox, 1, 0, 1, 1)

        self.prof_digitizeline_pushbutton = QPushButton(self.tr("Digitize line"))
        self.prof_digitizeline_pushbutton.clicked.connect(self.digitize_line)
        self.prof_digitizeline_pushbutton.setToolTip("Digitize a line on the map.\n"
                                                     "Left click: add point\n"
                                                     "Right click: end adding point\n"
                                                     "From: Define topographic sources (below)\n"
                                                     "you can use also an existing line\n"
                                                     "or a point list")
        qlyTopoSources.addWidget(self.prof_digitizeline_pushbutton, 2, 1, 1, 1)

        self.prof_clearline_pushbutton = QPushButton(self.tr("Clear"))
        self.prof_clearline_pushbutton.clicked.connect(self.clear_rubberband)
        qlyTopoSources.addWidget(self.prof_clearline_pushbutton, 2, 2, 1, 1)

        self.prof_clearline_pushbutton = QPushButton(self.tr("Save"))
        self.prof_clearline_pushbutton.clicked.connect(self.save_rubberband)
        qlyTopoSources.addWidget(self.prof_clearline_pushbutton, 2, 3, 1, 1)

        self.prof_toposources_reverse_direction_checkbox = QCheckBox(self.tr("Invert source line orientation"))
        qlyTopoSources.addWidget(self.prof_toposources_reverse_direction_checkbox, 3, 1, 1, 2)

        self.prof_toposources_pushbutton = QPushButton(self.tr("Define topographic sources"))
        self.prof_toposources_pushbutton.clicked.connect(self.define_topographic_sources)
        qlyTopoSources.addWidget(self.prof_toposources_pushbutton, 4, 0, 1, 4)
        """

        qgbxTopoSources.setLayout(qlyTopoSources)

        qlytTopoProfile.addWidget(qgbxTopoSources)

        ## Profile statistics

        prof_stats_QGroupBox = QGroupBox(qwdgTopoProfile)
        prof_stats_QGroupBox.setTitle("Profile statistics")

        prof_stats_Layout = QGridLayout()

        self.profile_stats_pushbutton = QPushButton(self.tr("Calculate profile statistics"))
        self.profile_stats_pushbutton.clicked.connect(self.calculate_profile_statistics)

        prof_stats_Layout.addWidget(self.profile_stats_pushbutton, 0, 0, 1, 3)

        prof_stats_QGroupBox.setLayout(prof_stats_Layout)

        qlytTopoProfile.addWidget(prof_stats_QGroupBox)

        ## Create profile section

        plotProfile_QGroupBox = QGroupBox(qwdgTopoProfile)
        plotProfile_QGroupBox.setTitle('Profile plot')

        plotProfile_Layout = QGridLayout()

        self.plotProfile_create_pushbutton = QPushButton(self.tr("Create profile"))
        self.plotProfile_create_pushbutton.clicked.connect(self.plot_topo_profiles)

        plotProfile_Layout.addWidget(self.plotProfile_create_pushbutton, 0, 0, 1, 4)

        plotProfile_QGroupBox.setLayout(plotProfile_Layout)

        qlytTopoProfile.addWidget(plotProfile_QGroupBox)

        ###################

        qwdgTopoProfile.setLayout(qlytTopoProfile)

        return qwdgTopoProfile

    def setup_project_section_tab(self):

        section_project_QWidget = QWidget()
        section_project_layout = QVBoxLayout()

        project_toolbox = QToolBox()

        ### Point project toolbox

        xs_point_proj_QWidget = QWidget()
        xs_point_proj_Layout = QVBoxLayout()

        ## input section

        xs_input_point_proj_QGroupBox = QGroupBox(xs_point_proj_QWidget)
        xs_input_point_proj_QGroupBox.setTitle('Input')

        xs_input_point_proj_Layout = QGridLayout()

        # input point geological layer

        xs_input_point_proj_Layout.addWidget(QLabel("Layer "), 0, 0, 1, 1)
        self.prj_struct_point_comboBox = QComboBox()
        self.prj_struct_point_comboBox.currentIndexChanged[int].connect(self.update_point_layers_boxes)

        xs_input_point_proj_Layout.addWidget(self.prj_struct_point_comboBox, 0, 1, 1, 6)
        self.struct_point_refresh_lyr_combobox()

        xs_input_point_proj_Layout.addWidget(QLabel("Fields:"), 1, 0, 1, 1)

        xs_input_point_proj_Layout.addWidget(QLabel("Id"), 1, 1, 1, 1)

        self.proj_point_id_fld_comboBox = QComboBox()
        xs_input_point_proj_Layout.addWidget(self.proj_point_id_fld_comboBox, 1, 2, 1, 1)

        self.plot_prj_use_rhrstrike_QRadioButt = QRadioButton("RHR str.")
        self.plot_prj_use_rhrstrike_QRadioButt.setChecked(True)
        self.plot_prj_use_dipdire_QRadioButt = QRadioButton("Dip dir.")
        xs_input_point_proj_Layout.addWidget(self.plot_prj_use_rhrstrike_QRadioButt, 1, 3, 1, 1)
        xs_input_point_proj_Layout.addWidget(self.plot_prj_use_dipdire_QRadioButt, 2, 3, 1, 1)
        self.proj_point_orient_fld_comboBox = QComboBox()
        xs_input_point_proj_Layout.addWidget(self.proj_point_orient_fld_comboBox, 1, 4, 1, 1)

        xs_input_point_proj_Layout.addWidget(QLabel("Dip"), 1, 5, 1, 1)
        self.proj_point_dipang_fld_comboBox = QComboBox()
        xs_input_point_proj_Layout.addWidget(self.proj_point_dipang_fld_comboBox, 1, 6, 1, 1)

        xs_input_point_proj_QGroupBox.setLayout(xs_input_point_proj_Layout)
        xs_point_proj_Layout.addWidget(xs_input_point_proj_QGroupBox)

        ## interpolation method

        xs_method_point_proj_QGroupBox = QGroupBox(xs_point_proj_QWidget)
        xs_method_point_proj_QGroupBox.setTitle('Project along')

        xs_method_point_proj_Layout = QGridLayout()

        self.nearest_point_proj_choice = QRadioButton("nearest intersection")
        self.nearest_point_proj_choice.setChecked(True)
        xs_method_point_proj_Layout.addWidget(self.nearest_point_proj_choice, 0, 0, 1, 3)

        self.axis_common_point_proj_choice = QRadioButton("axis with trend")
        xs_method_point_proj_Layout.addWidget(self.axis_common_point_proj_choice, 1, 0, 1, 1)

        self.common_axis_point_trend_SpinBox = QDoubleSpinBox()
        self.common_axis_point_trend_SpinBox.setMinimum(0.0)
        self.common_axis_point_trend_SpinBox.setMaximum(359.9)
        self.common_axis_point_trend_SpinBox.setDecimals(1)
        xs_method_point_proj_Layout.addWidget(self.common_axis_point_trend_SpinBox, 1, 1, 1, 1)

        xs_method_point_proj_Layout.addWidget(QLabel("and plunge"), 1, 2, 1, 1)

        self.common_axis_point_plunge_SpinBox = QDoubleSpinBox()
        self.common_axis_point_plunge_SpinBox.setMinimum(0.0)
        self.common_axis_point_plunge_SpinBox.setMaximum(89.9)
        self.common_axis_point_plunge_SpinBox.setDecimals(1)
        xs_method_point_proj_Layout.addWidget(self.common_axis_point_plunge_SpinBox, 1, 3, 1, 1)

        self.axis_individual_point_proj_choice = QRadioButton("axes from trend field")
        xs_method_point_proj_Layout.addWidget(self.axis_individual_point_proj_choice, 2, 0, 1, 1)

        self.proj_point_indivax_trend_fld_comboBox = QComboBox()
        xs_method_point_proj_Layout.addWidget(self.proj_point_indivax_trend_fld_comboBox, 2, 1, 1, 1)

        xs_method_point_proj_Layout.addWidget(QLabel("and plunge field"), 2, 2, 1, 1)
        self.proj_point_indivax_plunge_fld_comboBox = QComboBox()
        xs_method_point_proj_Layout.addWidget(self.proj_point_indivax_plunge_fld_comboBox, 2, 3, 1, 1)

        xs_method_point_proj_QGroupBox.setLayout(xs_method_point_proj_Layout)
        xs_point_proj_Layout.addWidget(xs_method_point_proj_QGroupBox)

        ## Plot groupbox

        xs_plot_proj_QGroupBox = QGroupBox(xs_point_proj_QWidget)
        xs_plot_proj_QGroupBox.setTitle('Plot geological attitudes')

        xs_plot_proj_Layout = QGridLayout()

        xs_plot_proj_Layout.addWidget(QLabel("Labels"), 0, 0, 1, 1)

        self.plot_prj_add_trendplunge_label = QCheckBox("or./dip")
        xs_plot_proj_Layout.addWidget(self.plot_prj_add_trendplunge_label, 0, 1, 1, 1)

        self.plot_prj_add_pt_id_label = QCheckBox("id")
        xs_plot_proj_Layout.addWidget(self.plot_prj_add_pt_id_label, 0, 2, 1, 1)

        xs_plot_proj_Layout.addWidget(QLabel("Color"), 0, 3, 1, 1)

        self.proj_point_color_QgsColorButtonV2 = QgsColorButtonV2()
        self.proj_point_color_QgsColorButtonV2.setColor(QColor('orange'))
        xs_plot_proj_Layout.addWidget(self.proj_point_color_QgsColorButtonV2, 0, 4, 1, 1)

        self.project_point_pushbutton = QPushButton(self.tr("Plot"))
        self.project_point_pushbutton.clicked.connect(self.create_struct_point_projection)
        xs_plot_proj_Layout.addWidget(self.project_point_pushbutton, 1, 0, 1, 3)

        self.reset_point_pushbutton = QPushButton(self.tr("Reset plot"))
        self.reset_point_pushbutton.clicked.connect(self.reset_struct_point_projection)

        xs_plot_proj_Layout.addWidget(self.reset_point_pushbutton, 1, 3, 1, 2)

        xs_plot_proj_QGroupBox.setLayout(xs_plot_proj_Layout)
        xs_point_proj_Layout.addWidget(xs_plot_proj_QGroupBox)

        self.flds_prj_point_comboBoxes = [self.proj_point_id_fld_comboBox,
                                          self.proj_point_orient_fld_comboBox,
                                          self.proj_point_dipang_fld_comboBox,
                                          self.proj_point_indivax_trend_fld_comboBox,
                                          self.proj_point_indivax_plunge_fld_comboBox]

        ##

        xs_point_proj_QWidget.setLayout(xs_point_proj_Layout)
        project_toolbox.addItem(xs_point_proj_QWidget, "Geological attitudes")

        ## END Point project toolbox

        ### Line project toolbox

        xs_line_proj_QWidget = QWidget()
        xs_line_proj_Layout = QVBoxLayout()

        ## input section

        xs_input_line_proj_QGroupBox = QGroupBox(xs_line_proj_QWidget)
        xs_input_line_proj_QGroupBox.setTitle('Input')

        xs_input_line_proj_Layout = QGridLayout()

        # input geological layer

        xs_input_line_proj_Layout.addWidget(QLabel("Layer"), 0, 0, 1, 1)
        self.prj_input_line_comboBox = QComboBox()

        xs_input_line_proj_Layout.addWidget(self.prj_input_line_comboBox, 0, 1, 1, 3)

        xs_input_line_proj_Layout.addWidget(QLabel("Id field"), 1, 0, 1, 1)
        self.id_fld_line_prj_comboBox = QComboBox()
        xs_input_line_proj_Layout.addWidget(self.id_fld_line_prj_comboBox, 1, 1, 1, 3)

        xs_input_line_proj_Layout.addWidget(QLabel("Line densify distance"), 2, 0, 1, 1)
        self.project_line_densify_distance_lineedit = QLineEdit()
        xs_input_line_proj_Layout.addWidget(self.project_line_densify_distance_lineedit, 2, 1, 1, 3)

        self.flds_prj_line_comboBoxes = [self.id_fld_line_prj_comboBox]

        xs_input_line_proj_QGroupBox.setLayout(xs_input_line_proj_Layout)
        xs_line_proj_Layout.addWidget(xs_input_line_proj_QGroupBox)

        ## interpolation method

        xs_method_line_proj_QGroupBox = QGroupBox(xs_line_proj_QWidget)
        xs_method_line_proj_QGroupBox.setTitle('Project')

        xs_method_line_proj_Layout = QGridLayout()

        xs_method_line_proj_Layout.addWidget(QLabel("Projection axis:"), 0, 0, 1, 1)

        xs_method_line_proj_Layout.addWidget(QLabel("trend"), 0, 1, 1, 1)

        self.common_axis_line_trend_SpinBox = QDoubleSpinBox()
        self.common_axis_line_trend_SpinBox.setMinimum(0.0)
        self.common_axis_line_trend_SpinBox.setMaximum(359.9)
        self.common_axis_line_trend_SpinBox.setDecimals(1)
        xs_method_line_proj_Layout.addWidget(self.common_axis_line_trend_SpinBox, 0, 2, 1, 1)

        xs_method_line_proj_Layout.addWidget(QLabel("plunge"), 0, 3, 1, 1)

        self.common_axis_line_plunge_SpinBox = QDoubleSpinBox()
        self.common_axis_line_plunge_SpinBox.setMinimum(0.0)
        self.common_axis_line_plunge_SpinBox.setMaximum(89.9)
        self.common_axis_line_plunge_SpinBox.setDecimals(1)
        xs_method_line_proj_Layout.addWidget(self.common_axis_line_plunge_SpinBox, 0, 4, 1, 1)

        # calculate profile

        self.project_line_pushbutton = QPushButton(self.tr("Plot traces"))
        self.project_line_pushbutton.clicked.connect(self.create_struct_line_projection)
        xs_method_line_proj_Layout.addWidget(self.project_line_pushbutton, 1, 0, 1, 5)

        self.reset_curves_pushbutton = QPushButton(self.tr("Reset traces"))
        self.reset_curves_pushbutton.clicked.connect(self.reset_structural_lines_projection)

        xs_method_line_proj_Layout.addWidget(self.reset_curves_pushbutton, 2, 0, 1, 5)

        xs_method_line_proj_QGroupBox.setLayout(xs_method_line_proj_Layout)
        xs_line_proj_Layout.addWidget(xs_method_line_proj_QGroupBox)

        ## 

        xs_line_proj_QWidget.setLayout(xs_line_proj_Layout)
        project_toolbox.addItem(xs_line_proj_QWidget, "Geological traces")

        ## END Line project toolbox

        # widget final setup

        section_project_layout.addWidget(project_toolbox)

        section_project_QWidget.setLayout(section_project_layout)

        return section_project_QWidget

    def setup_intersect_section_tab(self):

        intersect_widget = QWidget()
        intersect_layout = QVBoxLayout()

        intersect_toolbox = QToolBox()

        ### Line intersection section

        line_intersect_QWidget = QWidget()
        line_intersect_Layout = QVBoxLayout()

        ## input section

        inters_line_input_QGroupBox = QGroupBox(line_intersect_QWidget)
        inters_line_input_QGroupBox.setTitle('Input')

        inters_line_input_Layout = QGridLayout()

        # input traces layer

        inters_line_input_Layout.addWidget(QLabel("Line layer"), 0, 0, 1, 1)

        self.inters_input_line_comboBox = QComboBox()

        inters_line_input_Layout.addWidget(self.inters_input_line_comboBox, 0, 1, 1, 3)
        self.struct_line_refresh_lyr_combobox()

        inters_line_input_Layout.addWidget(QLabel("Id field"), 1, 0, 1, 1)
        self.inters_input_id_fld_line_comboBox = QComboBox()
        inters_line_input_Layout.addWidget(self.inters_input_id_fld_line_comboBox, 1, 1, 1, 3)

        self.flds_inters_line_comboBoxes = [self.inters_input_id_fld_line_comboBox]

        inters_line_input_Layout.addWidget(QLabel("Color"), 2, 0, 1, 1)
        self.inters_line_point_color_QgsColorButtonV2 = QgsColorButtonV2()
        self.inters_line_point_color_QgsColorButtonV2.setColor(QColor('blue'))
        inters_line_input_Layout.addWidget(self.inters_line_point_color_QgsColorButtonV2, 2, 1, 1, 1)

        inters_line_input_QGroupBox.setLayout(inters_line_input_Layout)
        line_intersect_Layout.addWidget(inters_line_input_QGroupBox)

        ## do section

        inters_line_do_QGroupBox = QGroupBox(line_intersect_QWidget)
        inters_line_do_QGroupBox.setTitle('Intersect')

        inters_line_do_Layout = QGridLayout()

        self.inters_line_do_pushbutton = QPushButton(self.tr("Intersect"))
        self.inters_line_do_pushbutton.clicked.connect(self.do_line_intersection)
        inters_line_do_Layout.addWidget(self.inters_line_do_pushbutton, 1, 0, 1, 4)

        self.line_inters_reset_pushbutton = QPushButton(self.tr("Reset intersections"))
        self.line_inters_reset_pushbutton.clicked.connect(self.line_intersection_reset)
        inters_line_do_Layout.addWidget(self.line_inters_reset_pushbutton, 2, 0, 1, 4)

        inters_line_do_QGroupBox.setLayout(inters_line_do_Layout)
        line_intersect_Layout.addWidget(inters_line_do_QGroupBox)

        # END do section

        line_intersect_QWidget.setLayout(line_intersect_Layout)
        intersect_toolbox.addItem(line_intersect_QWidget, "Intersect line layer")

        # END Line intersection section

        ### Polygon intersection section

        polygon_intersect_QWidget = QWidget()
        polygon_intersect_Layout = QVBoxLayout()

        ## input section

        inters_polygon_input_QGroupBox = QGroupBox(polygon_intersect_QWidget)
        inters_polygon_input_QGroupBox.setTitle('Input')

        inters_polygon_input_Layout = QGridLayout()

        # input traces layer

        inters_polygon_input_Layout.addWidget(QLabel("Polygon layer"), 0, 0, 1, 1)

        self.inters_input_polygon_comboBox = QComboBox()

        inters_polygon_input_Layout.addWidget(self.inters_input_polygon_comboBox, 0, 1, 1, 3)
        self.struct_polygon_refresh_lyr_combobox()

        inters_polygon_input_Layout.addWidget(QLabel("Classification field"), 1, 0, 1, 1)
        self.inters_polygon_classifaction_field_comboBox = QComboBox()
        inters_polygon_input_Layout.addWidget(self.inters_polygon_classifaction_field_comboBox, 1, 1, 1, 3)

        self.flds_inters_polygon_comboBoxes = [self.inters_polygon_classifaction_field_comboBox]

        inters_polygon_input_QGroupBox.setLayout(inters_polygon_input_Layout)
        polygon_intersect_Layout.addWidget(inters_polygon_input_QGroupBox)

        ## do section

        inters_polygon_do_QGroupBox = QGroupBox(polygon_intersect_QWidget)
        inters_polygon_do_QGroupBox.setTitle('Intersect')

        inters_polygon_do_Layout = QGridLayout()

        self.inters_polygon_do_pushbutton = QPushButton(self.tr("Intersect"))
        self.inters_polygon_do_pushbutton.clicked.connect(self.do_polygon_intersection)
        inters_polygon_do_Layout.addWidget(self.inters_polygon_do_pushbutton, 1, 0, 1, 4)

        self.polygon_inters_reset_pushbutton = QPushButton(self.tr("Reset intersections"))
        self.polygon_inters_reset_pushbutton.clicked.connect(self.polygon_intersection_reset)
        inters_polygon_do_Layout.addWidget(self.polygon_inters_reset_pushbutton, 2, 0, 1, 4)

        inters_polygon_do_QGroupBox.setLayout(inters_polygon_do_Layout)
        polygon_intersect_Layout.addWidget(inters_polygon_do_QGroupBox)

        # END do section

        polygon_intersect_QWidget.setLayout(polygon_intersect_Layout)
        intersect_toolbox.addItem(polygon_intersect_QWidget, "Intersect polygon layer")

        # END Polygon intersection section


        intersect_layout.addWidget(intersect_toolbox)

        intersect_widget.setLayout(intersect_layout)

        return intersect_widget

    def setup_export_section_tab(self):

        impexp_widget = QWidget()
        impexp_layout = QVBoxLayout()

        # Export section        

        export_QGroupBox = QGroupBox(impexp_widget)
        export_QGroupBox.setTitle('Export')

        export_inner_Layout = QGridLayout()

        self.export_image_QPushButton = QPushButton("Figure")
        export_inner_Layout.addWidget(self.export_image_QPushButton, 1, 0, 1, 4)
        self.export_image_QPushButton.clicked.connect(self.do_export_image)

        self.export_topographic_profile_QPushButton = QPushButton("Topographic profile data")
        export_inner_Layout.addWidget(self.export_topographic_profile_QPushButton, 2, 0, 1, 4)
        self.export_topographic_profile_QPushButton.clicked.connect(self.do_export_topo_profiles)

        self.export_project_geol_attitudes_QPushButton = QPushButton("Projected geological attitude data")
        export_inner_Layout.addWidget(self.export_project_geol_attitudes_QPushButton, 3, 0, 1, 4)
        self.export_project_geol_attitudes_QPushButton.clicked.connect(self.do_export_project_geol_attitudes)

        self.export_project_geol_lines__QPushButton = QPushButton("Projected geological line data")
        export_inner_Layout.addWidget(self.export_project_geol_lines__QPushButton, 4, 0, 1, 4)
        self.export_project_geol_lines__QPushButton.clicked.connect(self.do_export_project_geol_lines)

        self.export_line_intersections_QPushButton = QPushButton("Line intersection data")
        export_inner_Layout.addWidget(self.export_line_intersections_QPushButton, 5, 0, 1, 4)
        self.export_line_intersections_QPushButton.clicked.connect(self.do_export_line_intersections)

        self.export_polygon_intersections_QPushButton = QPushButton("Polygon intersection data")
        export_inner_Layout.addWidget(self.export_polygon_intersections_QPushButton, 6, 0, 1, 4)
        self.export_polygon_intersections_QPushButton.clicked.connect(self.do_export_polygon_intersections)

        export_QGroupBox.setLayout(export_inner_Layout)
        impexp_layout.addWidget(export_QGroupBox)

        impexp_widget.setLayout(impexp_layout)

        return impexp_widget

    def stop_rubberband(self):

        try:
            self.canvas_end_profile_line()
        except:
            pass

        try:
            self.clear_rubberband()
        except:
            pass

    def define_topographic_sources(self):

        def create_topo_profiles():

            selected_dems = None
            selected_dem_parameters = None

            sample_distance = None
            source_profile_line2dt = None

            if topo_source_type == self.demline_source:

                try:
                    selected_dems = dialog.selected_dems
                    selected_dem_parameters = dialog.selected_dem_parameters
                    topoline_colors = dialog.selected_dem_colors
                except Exception as e:
                    warn(self,
                         self.plugin_name,
                         "Input DEMs definition not correct: {}".format(e.message))
                    return

                try:
                    sample_distance = float(dialog.profile_densify_distance_lineedit.text())
                    assert sample_distance > 0.0
                except Exception as e:
                    warn(self,
                         self.plugin_name,
                         "Sample distance value not correct: {}".format(e.message))
                    return

                if dialog.DigitizeLine_checkbox.isChecked():
                    if self.digitized_profile_line2dt is None or \
                       self.digitized_profile_line2dt.num_pts < 2:
                        warn(self,
                             self.plugin_name,
                             "No digitized line available")
                        return
                    else:
                        source_profile_line2dt = self.digitized_profile_line2dt
                else:
                    self.stop_rubberband()
                    try:
                        source_profile_line2dt = dialog.dem_source_profile_line2dt
                    except:
                        warn(self,
                             self.plugin_name,
                             "DEM-line profile source not correctly created [1]")
                        return
                    if source_profile_line2dt is None:
                        warn(self,
                             self.plugin_name,
                             "DEM-line profile source not correctly created [2]")
                        return

            elif topo_source_type == self.gpxfile_source:
                self.stop_rubberband()
                try:
                    source_gpx_path = unicode(dialog.input_gpx_lineEdit.text())
                    if source_gpx_path == '':
                        warn(self,
                             self.plugin_name,
                             "Source GPX file is not set")
                        return
                except Exception as e:
                    warn(self,
                         self.plugin_name,
                         "Source GPX file not correctly set: {}".format(e.message))
                    return
                topoline_colors = [qcolor2rgbmpl(dialog.inputGPX_color_button.color())]

            else:
                warn(self,
                     self.plugin_name,
                     "Debug: uncorrect type source for topo sources def")
                return

            # calculates profiles
            invert_profile = self.prof_toposources_reverse_direction_checkbox.isChecked()
            if topo_source_type == self.demline_source:  # sources are DEM(s) and line
                try:
                    topo_profiles = topoprofiles_from_dems(self.canvas,
                                                            source_profile_line2dt,
                                                            sample_distance,
                                                            selected_dems,
                                                            selected_dem_parameters,
                                                            invert_profile)
                except Exception as e:
                     warn(self,
                         self.plugin_name,
                         e.message)
                     return
            elif topo_source_type == self.gpxfile_source:  # source is GPX file
                try:
                    topo_profiles = topoprofiles_from_gpxfile(source_gpx_path,
                                                           topoline_colors,
                                                           invert_profile)
                except Exception as e:
                    warn(self,
                         self.plugin_name,
                         "Error with profile calculation from GPX file: {}".format(e.message))
                    return
            else:  # source error
                error(self,
                      self.plugin_name,
                     "Algorithm error: profile calculation not defined")
                return

            if topo_profiles is None:
                warn(self,
                     self.plugin_name,
                     "Debug: profile not created")
                return

            profile_elements = Profile_Elements()
            profile_elements.profile_source_type = topo_source_type
            profile_elements.topoline_colors = topoline_colors
            profile_elements.source_profile_line = source_profile_line2dt
            profile_elements.sample_distance = sample_distance
            profile_elements.set_topo_profiles(topo_profiles)

            return profile_elements

        if self.prof_toposources_fromdems_checkbox.isChecked():
            dialog = TopoSourceFromDEMAndLineDialog(self.plugin_name, self.canvas)
            topo_source_type = self.demline_source
        elif self.prof_toposources_fromgpxfile_checkbox.isChecked():
            dialog = TopoSourceFromGPXFileDialog(self.plugin_name,
                                                 self.settings,
                                                 self.settings_gpxdir_key)
            topo_source_type = self.gpxfile_source
        else:
            warn(self,
                 self.plugin_name,
                 "Debug error: incorrect source definition")
            return

        if dialog.exec_():
            self.profile_elements = create_topo_profiles()
            if self.profile_elements is None:
                warn(self,
                     self.plugin_name,
                     "Debug error: self.profile_elements is None")
                return
            else:
                info(self,
                     self.plugin_name,
                     "Topographic sources defined")
        else:
            warn(self,
                 self.plugin_name,
                 "No topographic source defined")
            return


    def clear_rubberband(self):

        self.profile_canvas_points = []
        self.digitized_profile_line2dt = None
        try:
            self.rubberband.reset()
        except:
            pass


    def digitize_line(self):

        qprof_QWidget.map_digitations += 1

        self.clear_rubberband()
        self.profile_canvas_points = []

        if qprof_QWidget.map_digitations == 1:
            info(self,
                 self.plugin_name,
                 "Now you can digitize a line on the map.\nLeft click: add point\nRight click: end adding point")

        self.previous_maptool = self.canvas.mapTool()  # Save the standard map tool for restoring it at the end
        self.digitize_maptool = MapDigitizeTool(self.canvas)  # mouse listener
        self.canvas.setMapTool(self.digitize_maptool)
        self.connect_digitize_maptool()

        self.rubberband = QgsRubberBand(self.canvas)
        self.rubberband.setWidth(2)
        self.rubberband.setColor(QColor(Qt.red))

    def connect_digitize_maptool(self):

        QObject.connect(self.digitize_maptool, SIGNAL("moved"), self.canvas_refresh_profile_line)
        QObject.connect(self.digitize_maptool, SIGNAL("leftClicked"), self.profile_add_point)
        QObject.connect(self.digitize_maptool, SIGNAL("rightClicked"), self.canvas_end_profile_line)

    def disconnect_digitize_maptool(self):

        QObject.disconnect(self.digitize_maptool, SIGNAL("moved"), self.canvas_refresh_profile_line)
        QObject.disconnect(self.digitize_maptool, SIGNAL("leftClicked"), self.profile_add_point)
        QObject.disconnect(self.digitize_maptool, SIGNAL("rightClicked"), self.canvas_end_profile_line)


    def refresh_rubberband(self, xy_list):

        self.rubberband.reset(QGis.Line)
        for x, y in xy_list:
            self.rubberband.addPoint(QgsPoint(x, y))

    def canvas_refresh_profile_line(self, position):

        if len(self.profile_canvas_points) == 0:
            return

        x, y = xy_from_canvas(self.canvas, position)
        self.refresh_rubberband(self.profile_canvas_points + [[x, y]])

    def profile_add_point(self, position):

        x, y = xy_from_canvas(self.canvas, position)
        self.profile_canvas_points.append([x, y])

    def canvas_end_profile_line(self):

        self.refresh_rubberband(self.profile_canvas_points)

        if len(self.profile_canvas_points) > 1:
            self.digitized_profile_line2dt = Line(
                [Point(x, y) for x, y in self.profile_canvas_points])
        else:
            self.digitized_profile_line2dt = None

        self.profile_canvas_points = []

        self.restore_previous_map_tool()

    def restore_previous_map_tool(self):

        self.canvas.unsetMapTool(self.digitize_maptool)
        self.canvas.setMapTool(self.previous_maptool)

    def do_export_topo_profiles(self):

        def get_source_type():

            if dialog.src_allselecteddems_QRadioButton.isChecked():
                return ["all_dems"]
            elif dialog.src_singledem_QRadioButton.isChecked():
                return ["single_dem", dialog.src_singledemlist_QComboBox.currentIndex()]
            elif dialog.src_singlegpx_QRadioButton.isChecked():
                return ["gpx_file"]
            else:
                return []

        def get_format_type():

            if dialog.outtype_shapefile_line_QRadioButton.isChecked():
                return "shapefile - line"
            elif dialog.outtype_shapefile_point_QRadioButton.isChecked():
                return "shapefile - point"
            elif dialog.outtype_csv_QRadioButton.isChecked():
                return "csv"
            else:
                return ""

        try:
            self.profile_elements.topo_profiles.s
        except:
            warn(self,
                 self.plugin_name,
                 "Profile not yet calculated")
            return

        selected_dems_params = self.profile_elements.topo_profiles.dem_params
        dialog = TopographicProfileExportDialog(self.plugin_name,
                                                selected_dems_params)

        if dialog.exec_():

            output_source = get_source_type()
            if not output_source:
                warn(self,
                     self.plugin_name,
                     "Error in output source")
                return

            output_format = get_format_type()
            if output_format == "":
                warn(self,
                     self.plugin_name,
                     "Error in output format")
                return

            output_filepath = dialog.outpath_QLineEdit.text()
            if len(output_filepath) == 0:
                warn(self,
                     self.plugin_name,
                     "Error in output path")
                return
            add_to_project = dialog.load_output_checkBox.isChecked()
        else:
            warn(self,
                 self.plugin_name,
                 "No export defined")
            return

        # get project CRS information
        project_crs_osr = get_prjcrs_as_proj4str(self.canvas)

        self.output_topography(output_source, output_format, output_filepath, project_crs_osr)

        # add theme to QGis project
        if 'shapefile' in output_format and add_to_project:
            try:
                layer = QgsVectorLayer(output_filepath,
                                      QFileInfo(output_filepath).baseName(),
                                      "ogr")
                QgsMapLayerRegistry.instance().addMapLayer(layer)
            except:
                QMessageBox.critical(self, "Result", "Unable to load layer in project")
                return

    def output_topography(self, output_source, output_format, output_filepath, project_crs_osr):

        if output_source[0] == "all_dems":
            self.export_topography_all_dems(output_format, output_filepath, project_crs_osr)
        elif output_source[0] == "single_dem":
            ndx_dem_to_export = output_source[1]
            self.export_topography_single_dem(output_format, ndx_dem_to_export, output_filepath, project_crs_osr)
        elif output_source[0] == "gpx_file":
            self.export_topography_gpx_data(output_format, output_filepath, project_crs_osr)
        else:
            error(self,
                  self.plugin_name,
                  "Debug: output choice not correctly defined")
            return

    def do_export_project_geol_attitudes(self):

        def get_format_type():

            if dialog.outtype_shapefile_point_QRadioButton.isChecked():
                return "shapefile - point"
            elif dialog.outtype_csv_QRadioButton.isChecked():
                return "csv"
            else:
                return ""

        try:
            num_plane_attitudes_sets = len(self.profile_elements.plane_attitudes)
        except:
            warn(self,
                 self.plugin_name,
                 "No available geological attitudes")
            return
        else:
            if num_plane_attitudes_sets == 0:
                warn(self,
                     self.plugin_name,
                     "No available geological attitudes")
                return

        dialog = PointDataExportDialog(self.plugin_name)

        if dialog.exec_():

            output_format = get_format_type()
            if output_format == "":
                warn(self,
                     self.plugin_name,
                     "Error in output format")
                return
            output_filepath = dialog.outpath_QLineEdit.text()
            if len(output_filepath) == 0:
                warn(self,
                     self.plugin_name,
                     "Error in output path")
                return
            add_to_project = dialog.load_output_checkBox.isChecked()
        else:
            warn(self,
                 self.plugin_name,
                 "No export defined")
            return

        # get project CRS information
        project_crs_osr = get_prjcrs_as_proj4str(self.canvas)

        self.output_geological_attitudes(output_format, output_filepath, project_crs_osr)

        # add theme to QGis project
        if 'shapefile' in output_format and add_to_project:
            try:
                layer = QgsVectorLayer(output_filepath,
                                       QFileInfo(output_filepath).baseName(),
                                       "ogr")
                QgsMapLayerRegistry.instance().addMapLayer(layer)
            except:
                QMessageBox.critical(self, "Result", "Unable to load layer in project")
                return

    def output_geological_attitudes(self, output_format, output_filepath, project_crs_osr):

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

        parsed_geologicalattitudes_results = self.export_parse_geologicalattitudes_results(
            self.profile_elements.plane_attitudes)

        # output for csv file
        if output_format == "csv":
            success, msg = write_generic_csv(output_filepath, header_list, parsed_geologicalattitudes_results)
            if not success:
                warn(self,
                     self.plugin_name,
                     msg)
        elif output_format == "shapefile - point":
            success, msg = write_geological_attitudes_ptshp(output_filepath, parsed_geologicalattitudes_results, project_crs_osr)
            if not success:
                warn(self,
                     self.plugin_name,
                     msg)
        else:
            error(self,
                  self.plugin_name,
                  "Debug: error in export format")
            return

        if success:
            info(self,
                 self.plugin_name,
                 "Projected attitudes saved")

    def do_export_project_geol_lines(self):

        try:
            num_proj_lines_sets = len(self.profile_elements.curves)
        except:
            warn(self,
                 self.plugin_name,
                 "No available geological traces")
            return
        else:
            if num_proj_lines_sets == 0:
                warn(self,
                     self.plugin_name,
                     "No available geological traces to save")
                return

        fileName = QFileDialog.getSaveFileName(self,
                                               self.tr("Save results"),
                                               "*.csv",
                                               self.tr("csv (*.csv)"))

        if fileName is None or fileName == '':
            warn(self,
                 self.plugin_name,
                 "No output file has been defined")
            return

        parsed_curves_for_export = self.export_parse_geologicalcurves()
        header_list = ['id', 's', 'z']

        write_generic_csv(fileName, header_list, parsed_curves_for_export)

        info(self,
             self.plugin_name,
             "Projected lines saved")

    def do_export_line_intersections(self):

        def get_format_type():

            if dialog.outtype_shapefile_point_QRadioButton.isChecked():
                return "shapefile - point"
            elif dialog.outtype_csv_QRadioButton.isChecked():
                return "csv"
            else:
                return ""

        try:
            num_intersection_pts = len(self.profile_elements.intersection_pts)
        except:
            warn(self,
                 self.plugin_name,
                 "No available profile-line intersections")
            return
        else:
            if num_intersection_pts == 0:
                warn(self,
                     self.plugin_name,
                     "No available profile-line intersections")
                return

        dialog = PointDataExportDialog(self.plugin_name)

        if dialog.exec_():
            output_format = get_format_type()
            if output_format == "":
                warn(self,
                     self.plugin_name,
                     "Error in output format")
                return
            output_filepath = dialog.outpath_QLineEdit.text()
            if len(output_filepath) == 0:
                warn(self,
                     self.plugin_name,
                     "Error in output path")
                return
            add_to_project = dialog.load_output_checkBox.isChecked()
        else:
            warn(self,
                     self.plugin_name,
                     "No export defined")
            return

        # get project CRS information
        project_crs_osr = get_prjcrs_as_proj4str(self.canvas)

        self.output_profile_lines_intersections(output_format, output_filepath, project_crs_osr)

        # add theme to QGis project
        if 'shapefile' in output_format and add_to_project:
            try:
                layer = QgsVectorLayer(output_filepath,
                                       QFileInfo(output_filepath).baseName(),
                                       "ogr")
                QgsMapLayerRegistry.instance().addMapLayer(layer)
            except:
                QMessageBox.critical(self, "Result", "Unable to load layer in project")
                return

    def output_profile_lines_intersections(self, output_format, output_filepath, project_crs_osr):

        # definition of field names
        header_list = ['id',
                       's',
                       'x',
                       'y',
                       'z']

        parsed_profilelineintersections = self.export_parse_lineintersections(self.profile_elements.intersection_pts)

        # output for csv file
        if output_format == "csv":
            success, msg = write_generic_csv(output_filepath, header_list, parsed_profilelineintersections)
            if not success:
                warn(self,
                     self.plugin_name,
                     msg)
        elif output_format == "shapefile - point":
            success, msg = write_intersection_line_ptshp(output_filepath, header_list, parsed_profilelineintersections, project_crs_osr)
            if not success:
                warn(self,
                     self.plugin_name,
                     msg)
        else:
            error(self,
                  self.plugin_name,
                  "Debug: error in export format")
            return

        if success:
            info(self,
                 self.plugin_name,
                 "Line intersections saved")

    def save_rubberband(self):

        def get_format_type():

            if dialog.outtype_shapefile_line_QRadioButton.isChecked():
                return "shapefile - line"
            elif dialog.outtype_csv_QRadioButton.isChecked():
                return "csv"
            else:
                return ""

        if self.digitized_profile_line2dt is None:
            warn(self,
                     self.plugin_name,
                     "No available line to save [1]")
            return
        elif self.digitized_profile_line2dt.num_pts < 2:
            warn(self,
                 self.plugin_name,
                 "No available line to save [2]")
            return

        dialog = LineDataExportDialog(self.plugin_name)
        if dialog.exec_():
            output_format = get_format_type()
            if output_format == "":
                warn(self,
                     self.plugin_name,
                     "Error in output format")
                return
            output_filepath = dialog.outpath_QLineEdit.text()
            if len(output_filepath) == 0:
                warn(self,
                     self.plugin_name,
                     "Error in output path")
                return
            add_to_project = dialog.load_output_checkBox.isChecked()
        else:
            warn(self,
                     self.plugin_name,
                     "No export defined")
            return

        # get project CRS information
        project_crs_osr = get_prjcrs_as_proj4str(self.canvas)

        self.output_profile_line(output_format, output_filepath, self.digitized_profile_line2dt.pts, project_crs_osr)

        # add theme to QGis project
        if output_format == "shapefile - line" and add_to_project:
            try:
                digitized_line_layer = QgsVectorLayer(output_filepath,
                                                      QFileInfo(output_filepath).baseName(),
                                                      "ogr")
                QgsMapLayerRegistry.instance().addMapLayer(digitized_line_layer)
            except:
                QMessageBox.critical(self, "Result", "Unable to load layer in project")
                return

    def output_profile_line(self, output_format, output_filepath, pts2dt, proj_sr):

        points = [[n, pt2dt.x, pt2dt.y] for n, pt2dt in enumerate(pts2dt)]
        if output_format == "csv":
            success, msg = write_generic_csv(output_filepath,
                                             ['id', 'x', 'y'],
                                             points)
            if not success:
                warn(self,
                     self.plugin_name,
                     msg)
        elif output_format == "shapefile - line":
            success, msg = write_profile_lnshp(output_filepath,
                                               ['id'],
                                               points,
                                               proj_sr)
            if not success:
                warn(self,
                     self.plugin_name,
                     msg)
        else:
            error(self,
                  self.plugin_name,
                  "Debug: error in export format")
            return

        if success:
            info(self,
                 self.plugin_name,
                 "Line saved")

    def do_export_polygon_intersections(self):

        def get_format_type():

            if dialog.outtype_shapefile_line_QRadioButton.isChecked():
                return "shapefile - line"
            elif dialog.outtype_csv_QRadioButton.isChecked():
                return "csv"
            else:
                return ""

        try:
            num_intersection_lines = len(self.profile_elements.intersection_lines)
        except:
            warn(self,
                     self.plugin_name,
                     "No available profile-polygon intersections")
            return
        else:
            if num_intersection_lines == 0:
                warn(self,
                     self.plugin_name,
                     "No available profile-polygon intersections")
                return

        dialog = LineDataExportDialog(self.plugin_name)
        if dialog.exec_():
            output_format = get_format_type()
            if output_format == "":
                warn(self,
                     self.plugin_name,
                     "Error in output format")
                return
            output_filepath = dialog.outpath_QLineEdit.text()
            if len(output_filepath) == 0:
                warn(self,
                     self.plugin_name,
                     "Error in output path")
                return
            add_to_project = dialog.load_output_checkBox.isChecked()
        else:
            warn(self,
                     self.plugin_name,
                     "No export defined")
            return

        # get project CRS information
        project_crs_osr = get_prjcrs_as_proj4str(self.canvas)

        self.output_profile_polygons_intersections(output_format, output_filepath, project_crs_osr)

        # add theme to QGis project
        if 'shapefile' in output_format and add_to_project:
            try:
                layer = QgsVectorLayer(output_filepath,
                                       QFileInfo(output_filepath).baseName(),
                                       "ogr")
                QgsMapLayerRegistry.instance().addMapLayer(layer)
            except:
                QMessageBox.critical(self, "Result", "Unable to load layer in project")
                return

    def output_profile_polygons_intersections(self, output_format, output_filepath, sr):

        # definition of field names
        header_list = ['class_fld',
                       's',
                       'x',
                       'y',
                       'z']

        # output for csv file
        if output_format == "csv":
            success, msg = write_line_csv(output_filepath, header_list, self.profile_elements.intersection_lines)
            if not success:
                warn(self,
                     self.plugin_name,
                     msg)
        elif output_format == "shapefile - line":
            success, msg = write_intersection_polygon_lnshp(output_filepath, header_list, self.profile_elements.intersection_lines, sr)
            if not success:
                warn(self,
                     self.plugin_name,
                     msg)
        else:
            error("Debug: error in export format")
            return

        if success:
            info(self,
                 self.plugin_name,
                 "Polygon intersections saved")

    def struct_point_refresh_lyr_combobox(self):

        self.pointLayers = loaded_point_layers()

        update_ComboBox(self.prj_struct_point_comboBox,
                        self.choose_message,
                        [layer.name() for layer in self.pointLayers])

    def struct_polygon_refresh_lyr_combobox(self):

        self.current_polygon_layers = loaded_polygon_layers()
        update_ComboBox(self.inters_input_polygon_comboBox,
                        self.choose_message,
                        [layer.name() for layer in self.current_polygon_layers])

    def struct_line_refresh_lyr_combobox(self):

        self.current_line_layers = loaded_line_layers()
        update_ComboBox(self.prj_input_line_comboBox,
                        self.choose_message,
                        [layer.name() for layer in self.current_line_layers])
        update_ComboBox(self.inters_input_line_comboBox,
                        self.choose_message,
                        [layer.name() for layer in self.current_line_layers])


    def calculate_profile_statistics(self):

        if not self.check_pre_statistics():
            return

        topo_profiles = self.profile_elements.topo_profiles
        self.profile_elements.topo_profiles.statistics_elev = map(lambda p: get_statistics(p), topo_profiles.elevs)
        self.profile_elements.topo_profiles.statistics_dirslopes = map(lambda p: get_statistics(p), topo_profiles.dir_slopes)
        self.profile_elements.topo_profiles.statistics_slopes = map(lambda p: get_statistics(p), np.absolute(topo_profiles.dir_slopes))

        self.profile_elements.topo_profiles.profile_length = topo_profiles.s[-1] - topo_profiles.s[0]
        statistics_elev = self.profile_elements.topo_profiles.statistics_elev
        self.profile_elements.topo_profiles.natural_elev_range = (np.nanmin(np.array(map(lambda ds_stats: ds_stats["min"], statistics_elev))),
                                                                  np.nanmax(np.array(map(lambda ds_stats: ds_stats["max"], statistics_elev))))

        dialog = StatisticsDialog(self.plugin_name,
                                  self.profile_elements.topo_profiles)
        dialog.exec_()

        self.profile_elements.topo_profiles.statistics_defined = True


    def check_pre_statistics(self):

        if self.profile_elements is None or \
           self.profile_elements.profile_source_type is None or \
           self.profile_elements.topo_profiles is None:
            warn(self,
                     self.plugin_name,
                     "Source profile not yet defined")
            return False

        return True

    def check_pre_profile(self):

        if not self.check_pre_statistics():
            return False

        if not self.profile_elements.topo_profiles.statistics_defined:
            warn(self,
                     self.plugin_name,
                     "Statistics not yet calculated")
            return False

        return True

    def check_post_profile(self):

        if not self.check_pre_profile():
            return False

        if not self.profile_elements.topo_profiles.profile_defined:
            warn(self,
                     self.plugin_name,
                     "Topographic profile not yet created")
            return False

        return True

    def plot_topo_profiles(self):

        if not self.check_pre_profile():
            return

        natural_elev_min, natural_elev_max = self.profile_elements.topo_profiles.natural_elev_range
        profile_length = self.profile_elements.topo_profiles.profile_length
        dialog = PlotTopoProfileDialog(self.plugin_name,
                                       profile_length,
                                       natural_elev_min,
                                       natural_elev_max)

        if dialog.exec_():
            self.profile_elements.plot_params = get_profile_plot_params(dialog)
        else:
            return

        self.profile_elements.topo_profiles.profile_defined = True

        # plot profiles

        plot_addit_params = dict()
        plot_addit_params["add_trendplunge_label"] = self.plot_prj_add_trendplunge_label.isChecked()
        plot_addit_params["add_ptid_label"] = self.plot_prj_add_pt_id_label.isChecked()
        plot_addit_params["polygon_class_colors"] = self.polygon_classification_colors
        plot_addit_params["plane_attitudes_colors"] = self.plane_attitudes_colors

        profile_window = plot_profile_elements(self.profile_elements,
                                               plot_addit_params)
        self.profile_windows.append(profile_window)

    def export_parse_DEM_results(self, profiles_elements):

        # definition of output results         
        x_list = profiles_elements.topo_profiles.xs
        y_list = profiles_elements.topo_profiles.ys
        elev_list = profiles_elements.topo_profiles.elevs
        cumdist_2D_list = profiles_elements.topo_profiles.s
        cumdist_3d_list = profiles_elements.topo_profiles.s3d
        slopes = profiles_elements.topo_profiles.dir_slopes

        elev_list_zip = zip(*elev_list)
        cumdist_3d_list_zip = zip(*cumdist_3d_list)
        slope_list_zip = zip(*slopes)

        result_data = []
        rec_id = 0
        for x, y, cum_2d_dist, zs, cum3d_dists, slopes in zip(x_list, y_list, cumdist_2D_list, elev_list_zip,
                                                              cumdist_3d_list_zip, slope_list_zip):
            rec_id += 1
            record = [rec_id, x, y, cum_2d_dist]
            for z, cum3d_dist, slope in zip(zs, cum3d_dists, slopes):
                if isnan(z): z = ''
                if isnan(cum3d_dist): cum3d_dist = ''
                if isnan(slope): slope = ''
                record += [z, cum3d_dist, slope]
            result_data.append(record)

        return profiles_elements.get_current_dem_names(), result_data

    def export_topography_all_dems(self, out_format, outfile_path, proj_sr):

        if self.profile_elements.profile_source_type != self.demline_source:
            warn(self,
                     self.plugin_name,
                     "No DEM-derived profile defined")
            return

            # process results for data export
        dem_names, export_data = self.export_parse_DEM_results(self.profile_elements)

        # definition of field names
        dem_headers = []
        cum3ddist_headers = []
        slopes_headers = []
        for ndx in range(len(dem_names)):
            dem_headers.append(unicodedata.normalize('NFKD', unicode(dem_names[ndx][:10])).encode('ascii', 'ignore'))
            cum3ddist_headers.append("cds3d_" + str(ndx + 1))
            slopes_headers.append("slopd_" + str(ndx + 1))

        header_list = ["id", "x", "y", "cds2d"] + [name for sublist in
                                                   zip(dem_headers, cum3ddist_headers, slopes_headers) for name in
                                                   sublist]

        if out_format == "csv":
            success, msg = write_topography_allDEMs_csv(outfile_path, header_list, export_data)
            if not success:
                warn(self,
                     self.plugin_name,
                     msg)
        elif out_format == "shapefile - point":
            success, msg = write_topography_allDEMs_ptshp(outfile_path, header_list, dem_names, export_data, proj_sr)
            if not success:
                warn(self,
                     self.plugin_name,
                     msg)
        elif out_format == "shapefile - line":
            success, msg = write_topography_allDEMs_lnshp(outfile_path, header_list, dem_names, export_data, proj_sr)
            if not success:
                warn(self,
                     self.plugin_name,
                     msg)
        else:
            error("Debug: error in export all DEMs")
            return

        if success:
            info(self,
                 self.plugin_name,
                 "Profile export completed")

    def export_topography_single_dem(self, out_format, ndx_dem_to_export, outfile_path, prj_srs):

        if self.profile_elements.profile_source_type != self.demline_source:
            warn(self,
                 self.plugin_name,
                 "No DEM-derived profile defined")
            return

        # process results for data export
        _, export_data = self.export_parse_DEM_results(self.profile_elements)

        # definition of field names         
        header_list = ["id", "x", "y", "cds2d", "z", "cds3d", "dirslop"]

        if out_format == "csv":
            success, msg = write_topography_singleDEM_csv(outfile_path, header_list, export_data, ndx_dem_to_export)
            if not success:
                warn(self,
                     self.plugin_name,
                     msg)
        elif out_format == "shapefile - point":
            success, msg = write_topography_singleDEM_ptshp(outfile_path, header_list, export_data, ndx_dem_to_export, prj_srs)
            if not success:
                warn(self,
                     self.plugin_name,
                     msg)
        elif out_format == "shapefile - line":
            success, msg = write_topography_singleDEM_lnshp(outfile_path, header_list, export_data, ndx_dem_to_export, prj_srs)
            if not success:
                warn(self,
                     self.plugin_name,
                     msg)
        else:
            error(self,
                  self.plugin_name,
                 "Debug: error in export single DEM")
            return

        if success:
            info(self,
                 self.plugin_name,
                 "Profile export completed")

    def export_topography_gpx_data(self, out_format, output_filepath, prj_srs):

        if self.profile_elements.profile_source_type != self.gpxfile_source:
            warn(self,
                     self.plugin_name,
                     "No GPX-derived profile defined")
            return

        # process results from export
        gpx_parsed_results = self.export_parse_gpx_results()

        # definition of field names        
        header_list = ["id", "lat", "lon", "time", "elev", "cds2d", "cds3d", "dirslop"]

        if out_format == "csv":
            success, msg = write_generic_csv(output_filepath, header_list, gpx_parsed_results)
            if not success:
                warn(self,
                     self.plugin_name,
                     msg)
        elif out_format == "shapefile - point":
            success, msg = write_topography_GPX_ptshp(output_filepath, header_list, gpx_parsed_results, prj_srs)
            if not success:
                warn(self,
                     self.plugin_name,
                     msg)
        elif out_format == "shapefile - line":
            success, msg = write_topography_GPX_lnshp(output_filepath, header_list, gpx_parsed_results, prj_srs)
            if not success:
                warn(self,
                     self.plugin_name,
                     msg)
        else:
            error(self,
                  self.plugin_name,
                 "Debug: error in export single DEM")
            return

        if success:
            info(self,
                 self.plugin_name,
                 "Profile export completed")


    def gpx_profile_check_parameters(self):

        source_gpx_path = unicode(self.qlneInputGPXFile.text())
        if source_gpx_path == '':
            return False, 'Source GPX file is not set'

        plot_height_choice = self.GPX_plot_height_checkbox.isChecked()
        plot_slope_choice = self.GPX_plot_slope_checkbox.isChecked()

        if not (plot_height_choice or plot_slope_choice):
            return False, 'One of height or slope plot options are to be chosen'

        return True, 'OK'

    def export_parse_gpx_results(self):

        # definition of output results
        topo_profile = self.profile_elements.topo_profiles
        lat_list = topo_profile.lats
        lon_list = topo_profile.lons
        time_list = topo_profile.times
        cumdist_2D_list = topo_profile.s
        elev_list = topo_profile.elevs[0]  # [0] required for compatibility with DEM processing
        cumdist_3d_list = topo_profile.s3d[0]  # [0] required for compatibility with DEM processing
        dirslope_list = topo_profile.dir_slopes[0]  # [0] required for compatibility with DEM processing

        result_data = []
        rec_id = 0
        for lat, lon, time, elev, cumdist_2D, cumdist_3D, slope in zip(lat_list, lon_list, time_list, elev_list,
                                                                       cumdist_2D_list, cumdist_3d_list, dirslope_list):
            rec_id += 1
            if isnan(elev):
                elev = ''
            if isnan(cumdist_3D):
                cumdist_3D = ''
            if isnan(slope):
                slope = ''
            record = [rec_id, lat, lon, time, elev, cumdist_2D, cumdist_3D, slope]
            result_data.append(record)

        return result_data

    def update_point_layers_boxes(self):

        if len(self.pointLayers) == 0:
            return

        shape_qgis_ndx = self.prj_struct_point_comboBox.currentIndex() - 1  # minus 1 to account for initial text in combo box
        if shape_qgis_ndx < 0:
            return

        layer = self.pointLayers[shape_qgis_ndx]
        fields = layer.dataProvider().fields()
        field_names = [field.name() for field in fields.toList()]

        for ndx, combobox in enumerate(self.flds_prj_point_comboBoxes):
            combobox.clear()
            if ndx == 0:
                combobox.addItems(["none"])
            combobox.addItems(field_names)

    def update_linepoly_layers_boxes(self):

        def update_field_combo_boxes():

            for combobox in field_combobox_list:
                combobox.clear()

            if shape_qgis_ndx < 0 or len(layer_list) == 0:
                return

            fields = layer_list[shape_qgis_ndx].dataProvider().fields()
            field_names = [field.name() for field in fields.toList()]

            for combobox in field_combobox_list:
                combobox.addItems(["none"] + field_names)

        if self.sender() is self.prj_input_line_comboBox:
            shape_qgis_ndx = self.prj_input_line_comboBox.currentIndex() - 1  # minus 1 to account for initial text in combo box
            field_combobox_list = self.flds_prj_line_comboBoxes
            layer_list = self.current_line_layers
        elif self.sender() is self.inters_input_line_comboBox:
            shape_qgis_ndx = self.inters_input_line_comboBox.currentIndex() - 1  # minus 1 to account for initial text in combo box
            field_combobox_list = self.flds_inters_line_comboBoxes
            layer_list = self.current_line_layers
        elif self.sender() is self.inters_input_polygon_comboBox:
            shape_qgis_ndx = self.inters_input_polygon_comboBox.currentIndex() - 1  # minus 1 to account for initial text in combo box
            field_combobox_list = self.flds_inters_polygon_comboBoxes
            layer_list = self.current_polygon_layers

        update_field_combo_boxes()

    def get_current_combobox_values(self, combobox_list):

        return [combobox.currentText() for combobox in combobox_list]

    def calculate_section_data(self):

        sect_pt_1, sect_pt_2 = self.profile_elements.source_profile_line2dt.pts

        section_init_pt = Point(sect_pt_1.x, sect_pt_1.y, 0.0)
        section_final_pt = Point(sect_pt_2.x, sect_pt_2.y, 0.0)

        section_final_pt_up = Point(section_final_pt.x, section_final_pt.y,
                                       1000.0)  # arbitrary point on the same vertical as sect_pt_2
        section_cartes_plane = Plane.from_points(section_init_pt, section_final_pt, section_final_pt_up)
        section_vector = Segment(section_init_pt, section_final_pt).vector()

        return {'init_pt': section_init_pt, 'cartes_plane': section_cartes_plane, 'vector': section_vector}

    def struct_prjct_get_mapping_method(self):

        if self.nearest_point_proj_choice.isChecked():
            return {'method': 'nearest'}

        if self.axis_common_point_proj_choice.isChecked():
            return {'method': 'common axis',
                    'trend': float(self.common_axis_point_trend_SpinBox.value()),
                    'plunge': float(self.common_axis_point_plunge_SpinBox.value())}

        if self.axis_individual_point_proj_choice.isChecked():
            return {'method': 'individual axes',
                    'trend field': unicode(self.proj_point_indivax_trend_fld_comboBox.currentText()),
                    'plunge field': unicode(self.proj_point_indivax_plunge_fld_comboBox.currentText())}

    def check_for_struc_process(self):

        if not self.check_post_profile():
            return False

        # check that section is made up of only two points
        if self.profile_elements.source_profile_line2dt.num_pts != 2:
            warn(self,
                 self.plugin_name,
                 "Profile not made up by only two points")
            return False

        # check dem number is 1
        if len(self.profile_elements.topo_profiles.s3d) > 1:
            warn(self,
                 self.plugin_name,
                 "One (and only) topographic surface has to be used in the profile section")
            return False

        return True

    def check_struct_point_proj_parameters(self):

        if not self.check_for_struc_process():
            return False

        # get point structural layer with parameter fields
        prj_struct_point_qgis_ndx = self.prj_struct_point_comboBox.currentIndex() - 1  # minus 1 to account for initial text in combo box
        if prj_struct_point_qgis_ndx < 0:
            warn(self,
                 self.plugin_name,
                 "No defined point layer for structural data")
            return False

        return True

    def create_struct_point_projection(self):

        if not self.check_struct_point_proj_parameters():
            return

        # get color for projected points
        color = qcolor2rgbmpl(self.proj_point_color_QgsColorButtonV2.color())

        # define structural layer 
        prj_struct_point_qgis_ndx = self.prj_struct_point_comboBox.currentIndex() - 1  # minus 1 to account for initial text in combo box
        structural_layer = self.pointLayers[prj_struct_point_qgis_ndx]
        structural_layer_crs = structural_layer.crs()
        structural_field_list = self.get_current_combobox_values(self.flds_prj_point_comboBoxes)
        isRHRStrike = self.plot_prj_use_rhrstrike_QRadioButt.isChecked()

        # retrieve selected structural points with their attributes
        structural_pts_attrs = pt_geoms_attrs(structural_layer, structural_field_list)

        # list of structural points with original crs
        struct_pts_in_orig_crs = [Point(float(rec[0]), float(rec[1])) for rec in structural_pts_attrs]

        # IDs of structural points
        struct_pts_ids = [rec[2] for rec in structural_pts_attrs]

        # - geological planes (3D), as geological planes
        try:
            structural_planes = [GPlane(float(rec[3]), float(rec[4]), isRHRStrike) for rec in structural_pts_attrs]
        except:
            warn(self,
                 self.plugin_name,
                 "Check defined fields for possible errors")
            return

        struct_pts_3d = calculate_projected_3d_pts(self.canvas,
                                                    struct_pts_in_orig_crs,
                                                    structural_layer_crs,
                                                    self.profile_elements.topo_profiles.dem_params[0])

        # - zip together the point value data sets                     
        assert len(struct_pts_3d) == len(structural_planes)
        structural_data = zip(struct_pts_3d, structural_planes, struct_pts_ids)

        ### map points onto section ###

        # calculation of Cartesian plane expressing section plane        
        self.section_data = self.calculate_section_data()

        # calculation of projected structural points

        # get chosen mapping method
        mapping_method = self.struct_prjct_get_mapping_method()
        if mapping_method['method'] == 'individual axes':
            trend_field_name, plunge_field_name = mapping_method['trend field'], mapping_method['plunge field']
            # retrieve structural points mapping axes        
            mapping_method['individual_axes_values'] = vect_attrs(structural_layer,
                                                                  [trend_field_name, plunge_field_name])

        self.profile_elements.add_plane_attitudes(map_struct_pts_on_section(structural_data, self.section_data, mapping_method))
        self.plane_attitudes_colors.append(color)

        # plot profiles

        plot_addit_params = dict()
        plot_addit_params["add_trendplunge_label"] = self.plot_prj_add_trendplunge_label.isChecked()
        plot_addit_params["add_ptid_label"] = self.plot_prj_add_pt_id_label.isChecked()
        plot_addit_params["polygon_class_colors"] = self.polygon_classification_colors
        plot_addit_params["plane_attitudes_colors"] = self.plane_attitudes_colors

        profile_window = plot_profile_elements(self.profile_elements,
                                               plot_addit_params)
        self.profile_windows.append(profile_window)


    def reset_struct_point_projection(self):

        try:
            self.profile_elements.plane_attitudes = []
            self.plane_attitudes_colors = []
        except:
            pass

    def check_structural_line_projection_inputs(self):

        if not self.check_for_struc_process():
            return False

        # line structural layer with parameter fields
        prj_struct_line_qgis_ndx = self.prj_input_line_comboBox.currentIndex() - 1  # minus 1 to account for initial text in combo box
        if prj_struct_line_qgis_ndx < 0:
            warn(self,
                 self.plugin_name,
                 "No defined structural line layer")
            return False

        try:
            densify_distance = float(self.project_line_densify_distance_lineedit.text())
        except:
            warn(self,
                 self.plugin_name,
                 "No valid numeric value for densify line distance")
            return False
        else:
            if densify_distance <= 0.0:
                warn(self,
                     self.plugin_name,
                     "Densify line distance must be larger than zero")
                return False

        return True

    def create_struct_line_projection(self):

        # check input values
        if not self.check_structural_line_projection_inputs():
            return

        # input dem parameters
        demLayer =self.profile_elements.topo_profiles.dem_params[0].layer
        demParams =self.profile_elements.topo_profiles.dem_params[0].params

        # get line structural layer
        prj_struct_line_qgis_ndx = self.prj_input_line_comboBox.currentIndex() - 1  # minus 1 to account for initial text in combo box

        # get id field
        prj_struct_line_id_field_ndx = self.id_fld_line_prj_comboBox.currentIndex() - 1  # minus 1 to account for initial text in combo box

        # define structural layer        
        structural_line_layer = self.current_line_layers[prj_struct_line_qgis_ndx]

        on_the_fly_projection, project_crs = get_on_the_fly_projection(self.canvas)

        # read structural line values
        id_list = field_values(structural_line_layer, prj_struct_line_id_field_ndx)
        line_proj_crs_MultiLine2D_list = extract_multiline2d_list(structural_line_layer, on_the_fly_projection,
                                                                       project_crs)

        # densify with provided spat_distance
        densify_proj_crs_distance = float(self.project_line_densify_distance_lineedit.text())
        densified_proj_crs_MultiLine2D_list = [multiline_2d.densify_2d_multiline(densify_proj_crs_distance) for multiline_2d in
                                               line_proj_crs_MultiLine2D_list]

        # project to Dem CRS
        if on_the_fly_projection and demParams.crs != project_crs:
            densified_dem_crs_MultiLine2D_list = [multiline_2d.crs_project(project_crs, demParams.crs) for
                                                  multiline_2d in densified_proj_crs_MultiLine2D_list]
        else:
            densified_dem_crs_MultiLine2D_list = densified_proj_crs_MultiLine2D_list

        # interpolate z values from Dem
        z_list = [interpolate_z(demLayer, demParams, pt_2d) for multiline_2d in densified_dem_crs_MultiLine2D_list
                  for line_2d in multiline_2d.lines for pt_2d in line_2d.pts]

        # extract x-y pairs for creation of 3D points
        xy_list = [(pt_2d.x, pt_2d.y) for multiline_2d in densified_proj_crs_MultiLine2D_list for line_2d in
                   multiline_2d.lines for pt_2d in line_2d.pts]

        # replicate MultiLine list structure with 3D points with project CRS
        ndx = -1
        multiline_3d_proj_crs_list = []
        for multiline_2d in densified_proj_crs_MultiLine2D_list:
            multiline_3d_list = []
            for line_2d in multiline_2d.lines:
                line_3d_pts_list = []
                for _ in line_2d.pts:
                    ndx += 1
                    line_3d_pts_list.append(Point(xy_list[ndx][0], xy_list[ndx][1], z_list[ndx]))
                multiline_3d_list.append(Line(line_3d_pts_list))
            multiline_3d_proj_crs_list.append(MultiLine(multiline_3d_list))

        # create projection vector        
        trend = float(self.common_axis_line_trend_SpinBox.value())
        plunge = float(self.common_axis_line_plunge_SpinBox.value())
        axis_versor = GAxis(trend, plunge).as_vect().versor
        l, m, n = axis_versor.x, axis_versor.y, axis_versor.z

        # calculation of Cartesian plane expressing section plane        
        self.section_data = self.calculate_section_data()

        # project CartesianMultiLine3DT points to section
        intersection_point_list = []
        for multiline_3d in multiline_3d_proj_crs_list:
            for line_3d in multiline_3d.lines:
                for pt_3d in line_3d.pts:
                    srcPt = pt_3d
                    param_line = ParamLine3D(srcPt, l, m, n)
                    intersection_point_list.append(param_line.intersect_cartes_plane(self.section_data['cartes_plane']))

        # replicate MultiLine list structure with 3D points with project CRS
        ndx = -1
        multiline_3d_proj_crs_section_list = []
        for multiline_3d in multiline_3d_proj_crs_list:
            multiline_3d_list = []
            for line_3d in multiline_3d.lines:
                line_3d_pts_list = []
                for _ in line_3d.pts:
                    ndx += 1
                    line_3d_pts_list.append(intersection_point_list[ndx])
                multiline_3d_list.append(Line(line_3d_pts_list))
            multiline_3d_proj_crs_section_list.append(MultiLine(multiline_3d_list))

        section_start_point, section_vector = self.section_data['init_pt'], self.section_data['vector']
        curves_2d_list = []
        for multiline_3d in multiline_3d_proj_crs_section_list:
            multiline_2d_list = []
            for line_3d in multiline_3d.lines:
                line_2d_pts_list = []
                for pt_3d in line_3d.pts:
                    s = calculate_distance_with_sign(pt_3d, section_start_point, section_vector)
                    z = pt_3d.z
                    line_2d_pts_list.append(Point(s, z))
                multiline_2d_list.append(Line(line_2d_pts_list))
            curves_2d_list.append(MultiLine(multiline_2d_list))

        self.profile_elements.add_curves(curves_2d_list, id_list)

        # plot profiles

        plot_addit_params = dict()
        plot_addit_params["add_trendplunge_label"] = self.plot_prj_add_trendplunge_label.isChecked()
        plot_addit_params["add_ptid_label"] = self.plot_prj_add_pt_id_label.isChecked()
        plot_addit_params["polygon_class_colors"] = self.polygon_classification_colors
        plot_addit_params["plane_attitudes_colors"] = self.plane_attitudes_colors

        profile_window = plot_profile_elements(self.profile_elements,
                                               plot_addit_params)
        self.profile_windows.append(profile_window)

    def reset_structural_lines_projection(self):

        try:
            self.profile_elements.curves = []
            self.profile_elements.curves_ids = []
            self.curve_colors = []
        except:
            pass

    def export_parse_geologicalattitudes_results(self, plane_attitudes_datasets):

        result_data = []

        for dataset in plane_attitudes_datasets:

            for plane_attitude_rec in dataset:
                pt_id = plane_attitude_rec.id
                or_pt_x = plane_attitude_rec.src_pt_3d.x
                or_pt_y = plane_attitude_rec.src_pt_3d.y
                or_pt_z = plane_attitude_rec.src_pt_3d.z
                pr_pt_x = plane_attitude_rec.pt_3d.x
                pr_pt_y = plane_attitude_rec.pt_3d.y
                pr_pt_z = plane_attitude_rec.pt_3d.z
                s = plane_attitude_rec.sign_hor_dist
                or_dipdir = plane_attitude_rec.src_geol_plane.dd
                or_dipangle = plane_attitude_rec.src_geol_plane.da
                tr_dipangle = degrees(plane_attitude_rec.slope_rad)
                tr_dipdir = plane_attitude_rec.dwnwrd_sense

                record = [pt_id, or_pt_x, or_pt_y, or_pt_z, pr_pt_x, pr_pt_y, pr_pt_z, s, or_dipdir, or_dipangle,
                          tr_dipangle, tr_dipdir]

                result_data.append(record)

        return result_data

    def export_parse_geologicalcurves(self):

        data_list = []
        for curve_set, id_set in zip(self.profile_elements.curves, self.profile_elements.curves_ids):
            for curve, rec_id in zip(curve_set, id_set):
                for line in curve.lines:
                    for pt in line.pts:
                        data_list.append([rec_id, pt.x, pt.y])
        return data_list

    def export_parse_lineintersections(self, profile_intersection_pts):

        result_data = []

        for distances_from_profile_start, intersection_point3d, intersection_id, _ in profile_intersection_pts:
            result_data.append(
                [intersection_id, distances_from_profile_start, intersection_point3d.x, intersection_point3d.y,
                 intersection_point3d.z])

        return result_data

    def line_intersection_reset(self):

        if self.profile_elements is not None:
            self.profile_elements.intersection_pts = []

    def polygon_intersection_reset(self):

        if self.profile_elements is not None:
            self.profile_elements.intersection_lines = []

    def check_intersection_polygon_inputs(self):

        if not self.check_for_struc_process():
            return False

        # polygon layer with parameter fields
        intersection_polygon_qgis_ndx = self.inters_input_polygon_comboBox.currentIndex() - 1  # minus 1 to account for initial text in combo box
        if intersection_polygon_qgis_ndx < 0:
            warn(self,
                 self.plugin_name,
                 "No defined polygon layer")
            return False

        return True

    def check_intersection_line_inputs(self):

        if not self.check_for_struc_process():
            return False

        # line structural layer with parameter fields
        intersection_line_qgis_ndx = self.inters_input_line_comboBox.currentIndex() - 1  # minus 1 in order to account for initial text in combo box
        if intersection_line_qgis_ndx < 0:
            warn(self,
                 self.plugin_name,
                 "No defined geological line layer")
            return False

        return True

    def do_polygon_intersection(self):

        # check input values
        if not self.check_intersection_polygon_inputs():
            return

        # get dem parameters
        demLayer =self.profile_elements.topo_profiles.dem_params[0].layer
        demParams =self.profile_elements.topo_profiles.dem_params[0].params

        # profile line2d, in project CRS and densified 
        profile_line2d_prjcrs_densif = self.profile_elements.source_profile_line2dt.densify_2d_line(self.profile_elements.sample_distance)

        # polygon layer
        intersection_polygon_qgis_ndx = self.inters_input_polygon_comboBox.currentIndex() - 1  # minus 1 to account for initial text in combo box
        inters_polygon_classifaction_field_ndx = self.inters_polygon_classifaction_field_comboBox.currentIndex() - 1  # minus 1 to account for initial text in combo box
        polygon_layer = self.current_polygon_layers[intersection_polygon_qgis_ndx]
        polygon_layer_crs = polygon_layer.crs()

        on_the_fly_projection, project_crs = get_on_the_fly_projection(self.canvas)

        if on_the_fly_projection and polygon_layer_crs != project_crs:
            profile_line2d_polycrs_densif = profile_line2d_prjcrs_densif.crs_project(project_crs,
                                                                                     polygon_layer_crs)
        else:
            profile_line2d_polycrs_densif = profile_line2d_prjcrs_densif

        profile_qgsgeometry = QgsGeometry.fromPolyline(
            [QgsPoint(pt2d.x, pt2d.y) for pt2d in profile_line2d_polycrs_densif.pts])

        success, return_data = profile_polygon_intersection(profile_qgsgeometry,
                                                            polygon_layer,
                                                            inters_polygon_classifaction_field_ndx)

        if not success:
            error(self,
                  self.plugin_name,
                  return_data)
            return

        lIntersectPolylinePolygonCrs = return_data

        if len(lIntersectPolylinePolygonCrs) == 0:
            warn(self,
                 self.plugin_name,
                 "No intersection found")
            return

        # transform polyline intersections into prj crs line2d & classification list
        lIntersLine2dPrjCrs = []
        for intersection_polyline_polygon_crs in lIntersectPolylinePolygonCrs:
            rec_classification, xy_tuple_list = intersection_polyline_polygon_crs
            intersection_polygon_crs_line2d = xytuple_list_to_Line(xy_tuple_list)
            if on_the_fly_projection and polygon_layer_crs != project_crs:
                intersection_prj_crs_line2d = intersection_polygon_crs_line2d.crs_project(polygon_layer_crs,
                                                                                          project_crs)
            else:
                intersection_prj_crs_line2d = intersection_polygon_crs_line2d
            lIntersLine2dPrjCrs.append([rec_classification, intersection_prj_crs_line2d])

        # create Point lists from intersection with source DEM

        polygon_classification_set = set()
        sect_pt_1, sect_pt_2 = self.profile_elements.source_profile_line2dt.pts
        formation_list = []
        intersection_line3d_list = []
        intersection_polygon_s_list2 = []
        for polygon_classification, line2d in lIntersLine2dPrjCrs:
            polygon_classification_set.add(polygon_classification)

            lptIntersPts3d = intersect_with_dem(demLayer, demParams, on_the_fly_projection, project_crs, line2d.pts)
            lineIntersectionLine3d = Line(lptIntersPts3d)

            s0_list = lineIntersectionLine3d.incremental_length_2d()
            s_start = sect_pt_1.dist_2d(lineIntersectionLine3d.pts[0])
            s_list = [s + s_start for s in s0_list]

            formation_list.append(polygon_classification)
            intersection_line3d_list.append(lineIntersectionLine3d)
            intersection_polygon_s_list2.append(s_list)

        if len(intersection_polygon_s_list2) == 0:
            warn(self,
                 self.plugin_name,
                 "No reprojected intersection")
            return

        # create windows for user_definition of intersection colors in profile
        if polygon_classification_set != set() and polygon_classification_set != set([None]):

            dialog = PolygonIntersectionRepresentationDialog(self.plugin_name,
                                                             polygon_classification_set)
            if dialog.exec_():
                polygon_classification_colors_dict = self.classification_colors(dialog)
            else:
                warn(self,
                     self.plugin_name,
                     "No color chosen")
                return
            if len(polygon_classification_colors_dict) == 0:
                warn(self,
                     self.plugin_name,
                     "No defined colors")
                return
            else:
                self.polygon_classification_colors = polygon_classification_colors_dict
        else:
            self.polygon_classification_colors = None

        self.profile_elements.add_intersections_lines(formation_list, intersection_line3d_list, intersection_polygon_s_list2)

        # plot profiles
        plot_addit_params = dict()
        plot_addit_params["add_trendplunge_label"] = self.plot_prj_add_trendplunge_label.isChecked()
        plot_addit_params["add_ptid_label"] = self.plot_prj_add_pt_id_label.isChecked()
        plot_addit_params["polygon_class_colors"] = self.polygon_classification_colors
        plot_addit_params["plane_attitudes_colors"] = self.plane_attitudes_colors

        profile_window = plot_profile_elements(self.profile_elements,
                                               plot_addit_params)
        self.profile_windows.append(profile_window)

    def classification_colors(self, dialog):

        polygon_classification_colors_dict = dict()
        for classification_ndx in range(dialog.polygon_classifications_treeWidget.topLevelItemCount()):
            class_itemwidget = dialog.polygon_classifications_treeWidget.topLevelItem(classification_ndx)
            classification = unicode(class_itemwidget.text(0))
            # get color
            color = qcolor2rgbmpl(dialog.polygon_classifications_treeWidget.itemWidget(class_itemwidget, 1).color())
            polygon_classification_colors_dict[classification] = color

        return polygon_classification_colors_dict


    def do_line_intersection(self):

        # check input values
        if not self.check_intersection_line_inputs():
            return

        # get color for projected points
        color = qcolor2rgbmpl(self.inters_line_point_color_QgsColorButtonV2.color())

        # get dem parameters
        demLayer = self.profile_elements.topo_profiles.dem_params[0].layer
        demParams = self.profile_elements.topo_profiles.dem_params[0].params

        # get line structural layer
        intersection_line_qgis_ndx = self.inters_input_line_comboBox.currentIndex() - 1  # minus 1 to account for initial text in combo box

        # get id field
        intersection_line_id_field_ndx = self.inters_input_id_fld_line_comboBox.currentIndex() - 1  # minus 1 in order to account for initial text in combo box

        # define structural layer        
        structural_line_layer = self.current_line_layers[intersection_line_qgis_ndx]

        on_the_fly_projection, project_crs = get_on_the_fly_projection(self.canvas)

        # read structural line values
        if intersection_line_id_field_ndx == -1:
            id_list = None
        else:
            id_list = field_values(structural_line_layer, intersection_line_id_field_ndx)

        line_proj_crs_MultiLine2D_list = extract_multiline2d_list(structural_line_layer, on_the_fly_projection,
                                                                       project_crs)

        # calculated Point intersection list
        intersection_point_id_list = calculate_profile_lines_intersection(line_proj_crs_MultiLine2D_list,
                                                                               id_list,
                                                                               self.profile_elements.source_profile_line2dt)

        # sort intersection points by spat_distance from profile start point
        distances_from_profile_start_list = intersection_distances_by_profile_start_list(self.profile_elements.source_profile_line2dt,
                                                                                              intersection_point_id_list)

        # create CartesianPoint from intersection with source DEM
        intersection_point_list = [pt2d for pt2d, _ in intersection_point_id_list]
        intersection_id_list = [id for _, id in intersection_point_id_list]
        intersection_point3d_list = intersect_with_dem(demLayer, demParams, on_the_fly_projection, project_crs,
                                                            intersection_point_list)
        intersection_colors = [color] * len(intersection_point_list)

        self.profile_elements.add_intersections_pts(
            zip(distances_from_profile_start_list, intersection_point3d_list, intersection_id_list, intersection_colors))

        # plot profiles

        plot_addit_params = dict()
        plot_addit_params["add_trendplunge_label"] = self.plot_prj_add_trendplunge_label.isChecked()
        plot_addit_params["add_ptid_label"] = self.plot_prj_add_pt_id_label.isChecked()
        plot_addit_params["polygon_class_colors"] = self.polygon_classification_colors
        plot_addit_params["plane_attitudes_colors"] = self.plane_attitudes_colors

        profile_window = plot_profile_elements(self.profile_elements,
                                               plot_addit_params)
        self.profile_windows.append(profile_window)


    def do_export_image(self):

        try:
            profile_window = self.profile_windows[-1]
        except:
            warn(self,
                 self.plugin_name,
                 "Profile not yet calculated")
            return

        dialog = FigureExportDialog(self.plugin_name)

        if dialog.exec_():

            try:
                fig_width_inches = float(dialog.figure_width_inches_QLineEdit.text())
            except:
                warn(self,
                     self.plugin_name,
                     "Error in figure width value")
                return

            try:
                fig_resolution_dpi = int(dialog.figure_resolution_dpi_QLineEdit.text())
            except:
                warn(self,
                     self.plugin_name,
                     "Error in figure resolution value")
                return

            try:
                fig_font_size_pts = float(dialog.figure_fontsize_pts_QLineEdit.text())
            except:
                warn(self,
                     self.plugin_name,
                     "Error in font size value")

            try:
                fig_outpath = unicode(dialog.figure_outpath_QLineEdit.text())
            except:
                warn(self,
                     self.plugin_name,
                     "Error in figure output path")
                return

            try:
                top_space_value = float(dialog.top_space_value_QDoubleSpinBox.value())
            except:
                warn(self,
                     self.plugin_name,
                     "Error in figure top space value")
                return

            try:
                left_space_value = float(dialog.left_space_value_QDoubleSpinBox.value())
            except:
                warn(self,
                     self.plugin_name,
                     "Error in figure left space value")
                return

            try:
                right_space_value = float(dialog.right_space_value_QDoubleSpinBox.value())
            except:
                warn(self,
                     self.plugin_name,
                     "Error in figure right space value")
                return

            try:
                bottom_space_value = float(dialog.bottom_space_value_QDoubleSpinBox.value())
            except:
                warn(self,
                     self.plugin_name,
                     "Error in figure bottom space value")
                return

            try:
                blank_width_space = float(dialog.blank_width_space_value_QDoubleSpinBox.value())
            except:
                warn(self,
                     self.plugin_name,
                     "Error in figure blank widht space value")
                return

            try:
                blank_height_space = float(dialog.blank_height_space_value_QDoubleSpinBox.value())
            except:
                warn(self,
                     self.plugin_name,
                     "Error in figure blank height space value")
                return

        else:

            warn(self,
                 self.plugin_name,
                 "No export figure defined")
            return

        figure = profile_window.canvas.fig

        fig_current_width, fig_current_height = figure.get_size_inches()
        fig_scale_factor = fig_width_inches / fig_current_width
        figure.set_size_inches(fig_width_inches, fig_scale_factor * fig_current_height)

        for axis in figure.axes:
            for label in (axis.get_xticklabels() + axis.get_yticklabels()):
                label.set_fontsize(fig_font_size_pts)

        figure.subplots_adjust(wspace=blank_width_space, hspace=blank_height_space, left=left_space_value,
                               right=right_space_value, top=top_space_value, bottom=bottom_space_value)

        try:
            figure.savefig(str(fig_outpath), dpi=fig_resolution_dpi)
        except:
            warn(self,
                 self.plugin_name,
                 "Error with image saving")
        else:
            info(self,
                 self.plugin_name,
                 "Image saved")

    def closeEvent(self, event):

        try:
            self.reset_profile_defs()
        except:
            pass

        try:
            self.clear_rubberband()
        except:
            pass

        try:
            QgsMapLayerRegistry.instance().layerWasAdded.disconnect(self.struct_polygon_refresh_lyr_combobox)
        except:
            pass

        try:
            QgsMapLayerRegistry.instance().layerWasAdded.disconnect(self.struct_line_refresh_lyr_combobox)
        except:
            pass

        try:
            QgsMapLayerRegistry.instance().layerWasAdded.disconnect(self.struct_point_refresh_lyr_combobox)
        except:
            pass

        try:
            QgsMapLayerRegistry.instance().layerRemoved.disconnect(self.struct_polygon_refresh_lyr_combobox)
        except:
            pass

        try:
            QgsMapLayerRegistry.instance().layerRemoved.disconnect(self.struct_line_refresh_lyr_combobox)
        except:
            pass

        try:
            QgsMapLayerRegistry.instance().layerRemoved.disconnect(self.struct_point_refresh_lyr_combobox)
        except:
            pass


class TopoSourceFromDEMAndLineDialog(QDialog):

    def __init__(self, tPluginName, canvas, parent=None):

        super(TopoSourceFromDEMAndLineDialog, self).__init__(parent)

        self.plugin_name = tPluginName
        self.canvas = canvas
        self.on_the_fly_projection, self.project_crs = get_on_the_fly_projection(self.canvas)

        profileDEM_QWidget = QWidget()
        profileDEM_Layout = QVBoxLayout()

        ## input DEM section

        inputDEM_QGroupBox = QGroupBox()
        inputDEM_QGroupBox.setTitle("Input DEMs")

        inputDEM_Layout = QVBoxLayout()
        self.DefineSourceDEMs_pushbutton = QPushButton(self.tr("Define source DEMs"))
        self.DefineSourceDEMs_pushbutton.clicked.connect(self.define_source_DEMs)
        inputDEM_Layout.addWidget(self.DefineSourceDEMs_pushbutton)
        inputDEM_QGroupBox.setLayout(inputDEM_Layout)

        profileDEM_Layout.addWidget(inputDEM_QGroupBox)

        ## input Line layer section

        inputLine_QGroupBox = QGroupBox()
        inputLine_QGroupBox.setTitle("Input line")
        inputLine_Layout = QGridLayout()

        self.DigitizeLine_checkbox = QRadioButton(self.tr("digitized line"))
        self.DigitizeLine_checkbox.setChecked(True)
        inputLine_Layout.addWidget(self.DigitizeLine_checkbox, 0, 0, 1, 1)

        self.LoadLineLayer_checkbox = QRadioButton(self.tr("line layer"))
        inputLine_Layout.addWidget(self.LoadLineLayer_checkbox, 1, 0, 1, 1)
        self.DefineLineLayer_pushbutton = QPushButton(self.tr("Choose layer"))
        self.DefineLineLayer_pushbutton.clicked.connect(self.load_line_layer)
        inputLine_Layout.addWidget(self.DefineLineLayer_pushbutton, 1, 1, 1, 2)

        self.PointListforLine_checkbox = QRadioButton(self.tr("point list"))
        inputLine_Layout.addWidget(self.PointListforLine_checkbox, 2, 0, 1, 1)
        self.DefinePointList_pushbutton = QPushButton(self.tr("Create"))
        self.DefinePointList_pushbutton.clicked.connect(self.load_point_list)
        inputLine_Layout.addWidget(self.DefinePointList_pushbutton, 2, 1, 1, 2)

        # trace sampling spat_distance
        inputLine_Layout.addWidget(QLabel(self.tr("line densify distance")), 3, 0, 1, 1)
        self.profile_densify_distance_lineedit = QLineEdit()
        inputLine_Layout.addWidget(self.profile_densify_distance_lineedit, 3, 1, 1, 2)

        inputLine_QGroupBox.setLayout(inputLine_Layout)

        profileDEM_Layout.addWidget(inputLine_QGroupBox)
        profileDEM_QWidget.setLayout(profileDEM_Layout)

        # accept/reject section

        okButton = QPushButton("&OK")
        cancelButton = QPushButton("Cancel")

        buttonLayout = QHBoxLayout()
        buttonLayout.addStretch()
        buttonLayout.addWidget(okButton)
        buttonLayout.addWidget(cancelButton)

        profileDEM_Layout.addLayout(buttonLayout)

        self.setLayout(profileDEM_Layout)

        self.connect(okButton, SIGNAL("clicked()"),
                     self, SLOT("accept()"))
        self.connect(cancelButton, SIGNAL("clicked()"),
                     self, SLOT("reject()"))

        self.setWindowTitle("DEM(s) and line sources")

    def define_source_DEMs(self):

        self.selected_dems = None
        self.selected_dem_colors = None
        self.selected_dem_parameters = []

        current_raster_layers = loaded_monoband_raster_layers()
        if len(current_raster_layers) == 0:
            warn(self,
                 self.plugin_name,
                 "No loaded DEM")
            return

        dialog = SourceDEMsDialog(self.plugin_name, current_raster_layers)

        if dialog.exec_():
            selected_dems, selected_dem_colors = self.get_selected_dems_params(dialog)
        else:
            warn(self,
                 self.plugin_name,
                 "No chosen DEM")
            return

        if len(selected_dems) == 0:
            warn(self,
                 self.plugin_name,
                 "No selected DEM")
            return
        else:
            self.selected_dems = selected_dems
            self.selected_dem_colors = selected_dem_colors

        # get geodata
        self.selected_dem_parameters = [self.get_dem_parameters(dem) for dem in selected_dems]

        # get DEMs resolutions in project CRS and choose the min value
        dem_resolutions_prj_crs_list = []
        for dem, dem_params in zip(self.selected_dems, self.selected_dem_parameters):
            dem_resolutions_prj_crs_list.append(
                self.get_dem_resolution_in_prj_crs(dem, dem_params, self.on_the_fly_projection, self.project_crs))

        min_dem_resolution = min(dem_resolutions_prj_crs_list)
        if min_dem_resolution > 1:
            min_dem_proposed_resolution = round(min_dem_resolution)
        else:
            min_dem_proposed_resolution = min_dem_resolution
        self.profile_densify_distance_lineedit.setText(str(min_dem_proposed_resolution))

    def get_dem_parameters(self, dem):

        return QGisRasterParameters(*raster_qgis_params(dem))

    def get_selected_dems_params(self, dialog):

        selected_dems = []
        selected_dem_colors = []
        for dem_qgis_ndx in range(dialog.listDEMs_treeWidget.topLevelItemCount()):
            curr_DEM_item = dialog.listDEMs_treeWidget.topLevelItem(dem_qgis_ndx)
            if curr_DEM_item.checkState(0) == 2:
                selected_dems.append(dialog.singleband_raster_layers_in_project[dem_qgis_ndx])
                qcolor = dialog.listDEMs_treeWidget.itemWidget(curr_DEM_item, 2).color()
                mpl_color = qcolor2rgbmpl(qcolor)
                selected_dem_colors.append(mpl_color)

        return selected_dems, selected_dem_colors


    def get_line_layer_params(self, dialog):

        line_layer = dialog.line_shape
        order_field_ndx = dialog.Trace2D_order_field_comboBox.currentIndex()

        return line_layer, order_field_ndx

    def get_point_list(self, dialog):
        raw_point_string = dialog.point_list_qtextedit.toPlainText()
        raw_point_list = raw_point_string.split("\n")
        raw_point_list = map(lambda unicode_txt: clean_string(str(unicode_txt)), raw_point_list)
        data_list = filter(lambda rp: rp != "", raw_point_list)

        point_list = [to_float(xy_pair.split(",")) for xy_pair in data_list]
        line2d = xytuple_list_to_Line(point_list)

        return line2d

    def load_point_list(self):

        dialog = LoadPointListDialog(self.plugin_name)

        if dialog.exec_():
            line2d = self.get_point_list(dialog)
        else:
            warn(self,
                 self.plugin_name,
                 "No defined line source")
            return
        try:
            npts = line2d.num_pts
            if npts < 2:
                warn(self,
                     self.plugin_name,
                     "Defined line source with less than two points")
                return
        except:
            warn(self,
                 self.plugin_name,
                 "No defined line source")
            return

        self.dem_source_profile_line2dt = line2d

    def get_dem_resolution_in_prj_crs(self, dem, dem_params, on_the_fly_projection, prj_crs):

        def distance_projected_pts(x, y, delta_x, delta_y, src_crs, dest_crs):

            qgspt_start_src_crs = qgs_pt(x, y)
            qgspt_end_src_crs = qgs_pt(x + delta_x, y + delta_y)

            qgspt_start_dest_crs = project_qgs_point(qgspt_start_src_crs, src_crs, dest_crs)
            qgspt_end_dest_crs = project_qgs_point(qgspt_end_src_crs, src_crs, dest_crs)

            pt2_start_dest_crs = Point(qgspt_start_dest_crs.x(), qgspt_start_dest_crs.y())
            pt2d_end_dest_crs = Point(qgspt_end_dest_crs.x(), qgspt_end_dest_crs.y())

            return pt2_start_dest_crs.dist_2d(pt2d_end_dest_crs)

        cellsizeEW, cellsizeNS = dem_params.cellsizeEW, dem_params.cellsizeNS
        xMin, yMin = dem_params.xMin, dem_params.yMin

        if on_the_fly_projection and dem.crs() != prj_crs:
            cellsizeEW_prj_crs = distance_projected_pts(xMin, yMin, cellsizeEW, 0, dem.crs(), prj_crs)
            cellsizeNS_prj_crs = distance_projected_pts(xMin, yMin, 0, cellsizeNS, dem.crs(), prj_crs)
        else:
            cellsizeEW_prj_crs = cellsizeEW
            cellsizeNS_prj_crs = cellsizeNS

        return 0.5 * (cellsizeEW_prj_crs + cellsizeNS_prj_crs)

    def stop_profile_digitize_tool(self):

        try:
            self.disconnect_digitize_maptool()
        except:
            pass

        try:
            self.canvas.setMapTool(self.previous_maptool)
        except:
            pass

    def reset_rubber_band(self):

        try:
            self.rubberband.reset(QGis.Line)
        except:
            pass

    def reset_profile_defs(self):

        self.dem_source_profile_line2dt = None
        self.reset_rubber_band()
        self.stop_profile_digitize_tool()

    def load_line_layer(self):

        current_line_layers = loaded_line_layers()

        if len(current_line_layers) == 0:
            warn(self,
                 self.plugin_name,
                 "No available line layers")
            return

        dialog = SourceLineLayerDialog(self.plugin_name,
                                       current_line_layers)

        if dialog.exec_():
            line_layer, order_field_ndx = self.get_line_layer_params(dialog)
        else:
            warn(self,
                 self.plugin_name,
                 "No defined line source")
            return

        line_fld_ndx = int(order_field_ndx) - 1
        # get profile path from input line layer
        success, result = self.get_line_trace(line_layer, line_fld_ndx)
        if not success:
            raise VectorIOException(result)

        profile_orig_lines, mergeorder_ids = result

        profile_processed_line2d = merge_lines(profile_orig_lines, mergeorder_ids)

        # process input line layer
        profile_projected_line_2d = self.create_line_in_project_crs(profile_processed_line2d,
                                                                    line_layer.crs(),
                                                                    self.on_the_fly_projection,
                                                                    self.project_crs)

        self.dem_source_profile_line2dt = profile_projected_line_2d.remove_coincident_points()

    def get_line_trace(self, line_shape, order_field_ndx):

        try:
            profile_orig_lines, mergeorder_ids = line_geoms_with_id(line_shape, order_field_ndx)
        except VectorInputException as error_msg:
            return False, error_msg
        return True, (profile_orig_lines, mergeorder_ids)

    def create_line_in_project_crs(self, profile_processed_line, line_layer_crs, on_the_fly_projection,
                                   project_crs):

        if not on_the_fly_projection:
            return profile_processed_line
        else:
            return profile_processed_line.crs_project(line_layer_crs, project_crs)


class TopoSourceFromGPXFileDialog(QDialog):

    def __init__(self, plugin_name, settings, settings_gpxdir_key, parent=None):

        super(TopoSourceFromGPXFileDialog, self).__init__(parent)

        self.plugin_name = plugin_name
        self.settings = settings
        self.settings_gpxdir_key = settings_gpxdir_key

        profileGPX_Qwidget = QWidget()
        profileGPX_layout = QVBoxLayout()

        inputGPX_QGroupBox = QGroupBox()
        inputGPX_QGroupBox.setTitle('Input GPX file')

        inputGPX_Layout = QGridLayout()
        inputGPX_Layout.addWidget(QLabel(self.tr("Input GPX file with track points:")), 0, 0, 1, 1)

        self.input_gpx_lineEdit = QLineEdit()
        self.input_gpx_lineEdit.setPlaceholderText("my_track.gpx")
        inputGPX_Layout.addWidget(self.input_gpx_lineEdit, 0, 1, 1, 1)

        self.input_gpx_QPButt = QPushButton("...")
        self.input_gpx_QPButt.clicked.connect(self.select_input_gpxFile)
        inputGPX_Layout.addWidget(self.input_gpx_QPButt, 0, 2, 1, 1)

        inputGPX_Layout.addWidget(QLabel(self.tr("Profile color")), 1, 0, 1, 1)
        self.inputGPX_color_button = QgsColorButtonV2()
        self.inputGPX_color_button.setColor(QColor('red'))
        inputGPX_Layout.addWidget(self.inputGPX_color_button, 1, 1, 1, 1)

        inputGPX_QGroupBox.setLayout(inputGPX_Layout)

        profileGPX_layout.addWidget(inputGPX_QGroupBox)

        profileGPX_Qwidget.setLayout(profileGPX_layout)

        # accept/reject section

        okButton = QPushButton("&OK")
        cancelButton = QPushButton("Cancel")

        buttonLayout = QHBoxLayout()
        buttonLayout.addStretch()
        buttonLayout.addWidget(okButton)
        buttonLayout.addWidget(cancelButton)

        profileGPX_layout.addLayout(buttonLayout)

        self.setLayout(profileGPX_layout)

        self.connect(okButton, SIGNAL("clicked()"),
                     self, SLOT("accept()"))
        self.connect(cancelButton, SIGNAL("clicked()"),
                     self, SLOT("reject()"))

        self.setWindowTitle("GPX file source")


    def select_input_gpxFile(self):

        gpx_last_used_dir = self.settings.value(self.settings_gpxdir_key,
                                                "")
        fileName = QFileDialog.getOpenFileName(self,
                                               self.tr("Open GPX file"),
                                               gpx_last_used_dir,
                                               "GPX (*.gpx *.GPX)")
        if not fileName:
            return
        else:
            update_directory_key(self.settings,
                                 self.settings_gpxdir_key,
                                 fileName)
            self.input_gpx_lineEdit.setText(fileName)


class SourceDEMsDialog(QDialog):

    def __init__(self, plugin_name, raster_layers, parent=None):

        super(SourceDEMsDialog, self).__init__(parent)

        self.plugin_name = plugin_name

        self.singleband_raster_layers_in_project = raster_layers

        self.listDEMs_treeWidget = QTreeWidget()
        self.listDEMs_treeWidget.setColumnCount(3)
        self.listDEMs_treeWidget.setColumnWidth(0, 43)
        self.listDEMs_treeWidget.setColumnWidth(1, 135)
        self.listDEMs_treeWidget.setColumnWidth(2, 58)
        self.listDEMs_treeWidget.headerItem().setText(0, "Select")
        self.listDEMs_treeWidget.headerItem().setText(1, "Name")
        self.listDEMs_treeWidget.headerItem().setText(2, "Color")
        self.listDEMs_treeWidget.setHorizontalScrollBarPolicy(Qt.ScrollBarAsNeeded)
        self.listDEMs_treeWidget.setDragEnabled(False)
        self.listDEMs_treeWidget.setDragDropMode(QAbstractItemView.NoDragDrop)
        self.listDEMs_treeWidget.setAlternatingRowColors(True)
        self.listDEMs_treeWidget.setSelectionMode(QAbstractItemView.SingleSelection)
        self.listDEMs_treeWidget.setTextElideMode(Qt.ElideLeft)

        self.refresh_raster_layer_treewidget()

        okButton = QPushButton("&OK")
        cancelButton = QPushButton("Cancel")

        buttonLayout = QHBoxLayout()
        buttonLayout.addStretch()
        buttonLayout.addWidget(okButton)
        buttonLayout.addWidget(cancelButton)

        layout = QGridLayout()

        layout.addWidget(self.listDEMs_treeWidget, 0, 0, 1, 3)
        layout.addLayout(buttonLayout, 1, 0, 1, 3)

        self.setLayout(layout)

        self.connect(okButton, SIGNAL("clicked()"),
                     self, SLOT("accept()"))
        self.connect(cancelButton, SIGNAL("clicked()"),
                     self, SLOT("reject()"))

        self.setWindowTitle("Define source DEMs")

    def refresh_raster_layer_treewidget(self):

        self.listDEMs_treeWidget.clear()

        for raster_layer in self.singleband_raster_layers_in_project:
            tree_item = QTreeWidgetItem(self.listDEMs_treeWidget)
            tree_item.setText(1, raster_layer.name())
            color_button = QgsColorButtonV2()
            color_button.setColor(QColor('red'))

            self.listDEMs_treeWidget.setItemWidget(tree_item, 2, color_button)
            tree_item.setFlags(tree_item.flags() | Qt.ItemIsUserCheckable)
            tree_item.setCheckState(0, 0)


class SourceLineLayerDialog(QDialog):

    def __init__(self, plugin_name, current_line_layers, parent=None):

        super(SourceLineLayerDialog, self).__init__(parent)

        self.plugin_name = plugin_name
        self.current_line_layers = current_line_layers

        layout = QGridLayout()

        layout.addWidget(QLabel(self.tr("Line layer:")), 0, 0, 1, 1)
        self.LineLayers_comboBox = QComboBox()
        layout.addWidget(self.LineLayers_comboBox, 0, 1, 1, 3)
        self.refresh_input_profile_layer_combobox()

        layout.addWidget(QLabel(self.tr("Line order field:")), 1, 0, 1, 1)

        self.Trace2D_order_field_comboBox = QComboBox()
        layout.addWidget(self.Trace2D_order_field_comboBox, 1, 1, 1, 3)

        self.refresh_order_field_combobox()

        self.LineLayers_comboBox.currentIndexChanged[int].connect(self.refresh_order_field_combobox)

        okButton = QPushButton("&OK")
        cancelButton = QPushButton("Cancel")

        buttonLayout = QHBoxLayout()
        buttonLayout.addStretch()
        buttonLayout.addWidget(okButton)
        buttonLayout.addWidget(cancelButton)

        layout.addLayout(buttonLayout, 2, 0, 1, 3)

        self.setLayout(layout)

        self.connect(okButton, SIGNAL("clicked()"),
                     self, SLOT("accept()"))
        self.connect(cancelButton, SIGNAL("clicked()"),
                     self, SLOT("reject()"))

        self.setWindowTitle("Define source line layer")

    def refresh_input_profile_layer_combobox(self):

        self.LineLayers_comboBox.clear()

        for layer in self.current_line_layers:
            self.LineLayers_comboBox.addItem(layer.name())

        shape_qgis_ndx = self.LineLayers_comboBox.currentIndex()
        self.line_shape = self.current_line_layers[shape_qgis_ndx]

    def refresh_order_field_combobox(self):

        self.Trace2D_order_field_comboBox.clear()
        self.Trace2D_order_field_comboBox.addItem('--optional--')

        shape_qgis_ndx = self.LineLayers_comboBox.currentIndex()
        self.line_shape = self.current_line_layers[shape_qgis_ndx]

        line_layer_field_list = self.line_shape.dataProvider().fields().toList()
        for field in line_layer_field_list:
            self.Trace2D_order_field_comboBox.addItem(field.name())


class LoadPointListDialog(QDialog):

    def __init__(self, plugin_name, parent=None):

        super(LoadPointListDialog, self).__init__(parent)

        self.plugin_name = plugin_name
        layout = QGridLayout()

        layout.addWidget(QLabel(self.tr("Point list, with at least two points.")), 0, 0, 1, 1)
        layout.addWidget(
            QLabel(self.tr("Each point is defined by a comma-separated, x-y coordinate pair, one for each row")), 1, 0,
            1, 1)
        layout.addWidget(QLabel(self.tr("Example:\n549242.7, 242942.2\n578370.3, 322634.5")), 2, 0, 1, 1)

        self.point_list_qtextedit = QTextEdit()
        layout.addWidget(self.point_list_qtextedit, 3, 0, 1, 1)

        okButton = QPushButton("&OK")
        cancelButton = QPushButton("Cancel")

        buttonLayout = QHBoxLayout()
        buttonLayout.addStretch()
        buttonLayout.addWidget(okButton)
        buttonLayout.addWidget(cancelButton)

        layout.addLayout(buttonLayout, 4, 0, 1, 3)

        self.setLayout(layout)

        self.connect(okButton, SIGNAL("clicked()"),
                     self, SLOT("accept()"))
        self.connect(cancelButton, SIGNAL("clicked()"),
                     self, SLOT("reject()"))

        self.setWindowTitle("Point list")


class PolygonIntersectionRepresentationDialog(QDialog):

    colors = ["darkseagreen", "darkgoldenrod", "darkviolet", "hotpink", "powderblue", "yellowgreen", "palevioletred",
              "seagreen", "darkturquoise", "beige", "darkkhaki", "red", "yellow", "magenta", "blue", "cyan",
              "chartreuse"]

    def __init__(self, plugin_name, polygon_classification_set, parent=None):

        super(PolygonIntersectionRepresentationDialog, self).__init__(parent)

        self.plugin_name = plugin_name
        self.polygon_classifications = list(polygon_classification_set)

        self.polygon_classifications_treeWidget = QTreeWidget()
        self.polygon_classifications_treeWidget.setColumnCount(2)
        self.polygon_classifications_treeWidget.setColumnWidth(0, 200)
        self.polygon_classifications_treeWidget.headerItem().setText(0, "Name")
        self.polygon_classifications_treeWidget.headerItem().setText(1, "Plot color")
        self.polygon_classifications_treeWidget.setHorizontalScrollBarPolicy(Qt.ScrollBarAsNeeded)
        self.polygon_classifications_treeWidget.setDragEnabled(False)
        self.polygon_classifications_treeWidget.setDragDropMode(QAbstractItemView.NoDragDrop)
        self.polygon_classifications_treeWidget.setAlternatingRowColors(True)
        self.polygon_classifications_treeWidget.setTextElideMode(Qt.ElideLeft)

        self.refresh_classification_colors_treewidget()

        okButton = QPushButton("&OK")
        cancelButton = QPushButton("Cancel")

        buttonLayout = QHBoxLayout()
        buttonLayout.addStretch()
        buttonLayout.addWidget(okButton)
        buttonLayout.addWidget(cancelButton)

        layout = QGridLayout()

        layout.addWidget(self.polygon_classifications_treeWidget, 0, 0, 1, 3)
        layout.addLayout(buttonLayout, 1, 0, 1, 3)

        self.setLayout(layout)

        self.connect(okButton, SIGNAL("clicked()"),
                     self, SLOT("accept()"))
        self.connect(cancelButton, SIGNAL("clicked()"),
                     self, SLOT("reject()"))

        self.setWindowTitle("Polygon intersection colors")

    def refresh_classification_colors_treewidget(self):

        if len(PolygonIntersectionRepresentationDialog.colors) < len(self.polygon_classifications):
            dupl_factor = 1 + int(len(self.polygon_classifications) / len(PolygonIntersectionRepresentationDialog.colors))
            curr_colors = dupl_factor * PolygonIntersectionRepresentationDialog.colors
        else:
            curr_colors = PolygonIntersectionRepresentationDialog.colors

        self.polygon_classifications_treeWidget.clear()

        for classification_id, color in zip(self.polygon_classifications, curr_colors):
            tree_item = QTreeWidgetItem(self.polygon_classifications_treeWidget)
            tree_item.setText(0, unicode(classification_id))

            color_QgsColorButtonV2 = QgsColorButtonV2()
            color_QgsColorButtonV2.setColor(QColor(color))
            self.polygon_classifications_treeWidget.setItemWidget(tree_item, 1, color_QgsColorButtonV2)


def get_profile_plot_params(dialog):

    profile_params = {}

    # get profile plot parameters

    try:
        profile_params['plot_min_elevation_user'] = float(dialog.plot_min_value_QLineedit.text())
    except:
        profile_params['plot_min_elevation_user'] = None

    try:
        profile_params['plot_max_elevation_user'] = float(dialog.plot_max_value_QLineedit.text())
    except:
        profile_params['plot_max_elevation_user'] = None

    try:
        profile_params['vertical_exaggeration'] = float(dialog.DEM_exageration_ratio_Qlineedit.text())
        assert profile_params['vertical_exaggeration'] > 0
    except:
        profile_params['vertical_exaggeration'] = 1

    profile_params['filled_height'] = dialog.plotProfile_height_filled_checkbox.isChecked()
    profile_params['filled_slope'] = dialog.plotProfile_slope_filled_checkbox.isChecked()
    profile_params['plot_height_choice'] = dialog.plotProfile_height_checkbox.isChecked()
    profile_params['plot_slope_choice'] = dialog.plotProfile_slope_checkbox.isChecked()
    profile_params['plot_slope_absolute'] = dialog.plotProfile_slope_absolute_qradiobutton.isChecked()
    profile_params['plot_slope_directional'] = dialog.plotProfile_slope_directional_qradiobutton.isChecked()
    profile_params['invert_xaxis'] = dialog.plotProfile_invert_xaxis_checkbox.isChecked()

    return profile_params


class PlotTopoProfileDialog(QDialog):

    def __init__(self, plugin_name, profile_length, natural_elev_min, natural_elev_max, parent=None):

        super(PlotTopoProfileDialog, self).__init__(parent)

        self.plugin_name = plugin_name

        # pre-process elevation values

        # suggested plot elevation range
        z_padding = 0.5
        delta_z = natural_elev_max - natural_elev_min
        if delta_z < 0.0:
            warn(self,
                 self.plugin_name,
                 "Error: min elevation larger then max elevation")
            return
        elif delta_z == 0.0:
            plot_z_min = floor(natural_elev_min) - 10
            plot_z_max = ceil(natural_elev_max) + 10
        else:
            plot_z_min = floor(natural_elev_min - delta_z * z_padding)
            plot_z_max = ceil(natural_elev_max + delta_z * z_padding)
        delta_plot_z = plot_z_max - plot_z_min

        # suggested exxageration value
        w_to_h_rat = float(profile_length) / float(delta_plot_z)
        sugg_ve = 0.2*w_to_h_rat

        layout = QVBoxLayout()

        # profile options

        plotProfile_Layout = QVBoxLayout()

        plotsettings_groupbox = QGroupBox("Plot settings")

        plotsetting_layout = QGridLayout()

        plotsetting_layout.addWidget(QLabel(self.tr("Vertical exaggeration")), 0, 0, 1, 1)
        self.DEM_exageration_ratio_Qlineedit = QLineEdit()
        self.DEM_exageration_ratio_Qlineedit.setText("%f" % sugg_ve)
        plotsetting_layout.addWidget(self.DEM_exageration_ratio_Qlineedit, 0, 1, 1, 1)

        plotsetting_layout.addWidget(QLabel(self.tr("Plot z min value")), 1, 0, 1, 1)
        self.plot_min_value_QLineedit = QLineEdit()
        self.plot_min_value_QLineedit.setText("%f" % plot_z_min)
        plotsetting_layout.addWidget(self.plot_min_value_QLineedit, 1, 1, 1, 1)

        plotsetting_layout.addWidget(QLabel(self.tr("Plot z max value")), 1, 2, 1, 1)
        self.plot_max_value_QLineedit = QLineEdit()
        self.plot_max_value_QLineedit.setText("%f" % plot_z_max)
        plotsetting_layout.addWidget(self.plot_max_value_QLineedit, 1, 3, 1, 1)

        plotsettings_groupbox.setLayout(plotsetting_layout)

        plotProfile_Layout.addWidget(plotsettings_groupbox)

        plotVariables_groupbox = QGroupBox("Plot variables")

        plotvariables_layout = QGridLayout()

        self.plotProfile_height_filled_checkbox = QCheckBox(self.tr("filled"))
        plotvariables_layout.addWidget(self.plotProfile_height_filled_checkbox, 0, 0, 1, 1)

        self.plotProfile_height_checkbox = QCheckBox(self.tr("Height"))
        self.plotProfile_height_checkbox.setChecked(True)
        plotvariables_layout.addWidget(self.plotProfile_height_checkbox, 0, 1, 1, 1)

        self.plotProfile_slope_filled_checkbox = QCheckBox(self.tr("filled"))
        plotvariables_layout.addWidget(self.plotProfile_slope_filled_checkbox, 1, 0, 1, 1)

        self.plotProfile_slope_checkbox = QCheckBox(self.tr("Slope (degrees)"))
        plotvariables_layout.addWidget(self.plotProfile_slope_checkbox, 1, 1, 1, 1)

        self.plotProfile_slope_absolute_qradiobutton = QRadioButton(self.tr("absolute"))
        self.plotProfile_slope_absolute_qradiobutton.setChecked(True);
        plotvariables_layout.addWidget(self.plotProfile_slope_absolute_qradiobutton, 1, 2, 1, 1)

        self.plotProfile_slope_directional_qradiobutton = QRadioButton(self.tr("directional"))
        plotvariables_layout.addWidget(self.plotProfile_slope_directional_qradiobutton, 1, 3, 1, 1)

        plotvariables_layout.addWidget(QLabel("Note: to  calculate correctly the slope, the project must have a CRS set or the DEM(s) must not be in lon-lat"), 2, 1, 2, 3)
        plotVariables_groupbox.setLayout(plotvariables_layout)

        plotProfile_Layout.addWidget(plotVariables_groupbox)

        advanced_parameters_groupbox = QGroupBox("Advanced parameters")

        advparams_layout = QHBoxLayout()

        self.plotProfile_invert_xaxis_checkbox = QCheckBox(self.tr("Flip x-axis direction"))
        advparams_layout.addWidget(self.plotProfile_invert_xaxis_checkbox)

        advanced_parameters_groupbox.setLayout(advparams_layout)

        plotProfile_Layout.addWidget(advanced_parameters_groupbox)

        layout.addLayout(plotProfile_Layout)

        # ok/cancel section

        okButton = QPushButton("&OK")
        cancelButton = QPushButton("Cancel")

        buttonLayout = QHBoxLayout()
        buttonLayout.addStretch()
        buttonLayout.addWidget(okButton)
        buttonLayout.addWidget(cancelButton)

        layout.addLayout(buttonLayout)

        self.setLayout(layout)

        self.connect(okButton, SIGNAL("clicked()"),
                     self, SLOT("accept()"))
        self.connect(cancelButton, SIGNAL("clicked()"),
                     self, SLOT("reject()"))

        self.setWindowTitle("Topographic plot parameters")


class FigureExportDialog(QDialog):

    def __init__(self, plugin_name, parent=None):

        super(FigureExportDialog, self).__init__(parent)

        self.plugin_name = plugin_name

        layout = QVBoxLayout()

        # main parameters gropbox

        main_params_groupBox = QGroupBox("Main graphic parameters")

        main_params_layout = QGridLayout()

        main_params_layout.addWidget(QLabel(self.tr("Figure width (inches)")), 0, 0, 1, 1)
        self.figure_width_inches_QLineEdit = QLineEdit("10")
        main_params_layout.addWidget(self.figure_width_inches_QLineEdit, 0, 1, 1, 1)

        main_params_layout.addWidget(QLabel(self.tr("Resolution (dpi)")), 0, 2, 1, 1)
        self.figure_resolution_dpi_QLineEdit = QLineEdit("200")
        main_params_layout.addWidget(self.figure_resolution_dpi_QLineEdit, 0, 3, 1, 1)

        main_params_layout.addWidget(QLabel(self.tr("Font size (pts)")), 0, 4, 1, 1)
        self.figure_fontsize_pts_QLineEdit = QLineEdit("12")
        main_params_layout.addWidget(self.figure_fontsize_pts_QLineEdit, 0, 5, 1, 1)

        main_params_groupBox.setLayout(main_params_layout)

        layout.addWidget(main_params_groupBox)

        # additional parameters groupbox

        add_params_groupBox = QGroupBox(self.tr("Subplot configuration tools parameters"))

        add_params_layout = QGridLayout()

        add_params_layout.addWidget(QLabel("Top space"), 0, 2, 1, 1)
        self.top_space_value_QDoubleSpinBox = QDoubleSpinBox()
        self.top_space_value_QDoubleSpinBox.setRange(0.0, 1.0)
        self.top_space_value_QDoubleSpinBox.setDecimals(2)
        self.top_space_value_QDoubleSpinBox.setSingleStep(0.01)
        self.top_space_value_QDoubleSpinBox.setValue(0.96)
        add_params_layout.addWidget(self.top_space_value_QDoubleSpinBox, 0, 3, 1, 1)

        add_params_layout.addWidget(QLabel("Left space"), 1, 0, 1, 1)
        self.left_space_value_QDoubleSpinBox = QDoubleSpinBox()
        self.left_space_value_QDoubleSpinBox.setRange(0.0, 1.0)
        self.left_space_value_QDoubleSpinBox.setDecimals(2)
        self.left_space_value_QDoubleSpinBox.setSingleStep(0.01)
        self.left_space_value_QDoubleSpinBox.setValue(0.1)
        add_params_layout.addWidget(self.left_space_value_QDoubleSpinBox, 1, 1, 1, 1)

        add_params_layout.addWidget(QLabel("Right space"), 1, 4, 1, 1)
        self.right_space_value_QDoubleSpinBox = QDoubleSpinBox()
        self.right_space_value_QDoubleSpinBox.setRange(0.0, 1.0)
        self.right_space_value_QDoubleSpinBox.setDecimals(2)
        self.right_space_value_QDoubleSpinBox.setSingleStep(0.01)
        self.right_space_value_QDoubleSpinBox.setValue(0.96)
        add_params_layout.addWidget(self.right_space_value_QDoubleSpinBox, 1, 5, 1, 1)

        add_params_layout.addWidget(QLabel("Bottom space"), 2, 2, 1, 1)
        self.bottom_space_value_QDoubleSpinBox = QDoubleSpinBox()
        self.bottom_space_value_QDoubleSpinBox.setRange(0.0, 1.0)
        self.bottom_space_value_QDoubleSpinBox.setDecimals(2)
        self.bottom_space_value_QDoubleSpinBox.setSingleStep(0.01)
        self.bottom_space_value_QDoubleSpinBox.setValue(0.06)
        add_params_layout.addWidget(self.bottom_space_value_QDoubleSpinBox, 2, 3, 1, 1)

        add_params_layout.addWidget(QLabel("Blank width space between subplots"), 3, 0, 1, 2)
        self.blank_width_space_value_QDoubleSpinBox = QDoubleSpinBox()
        self.blank_width_space_value_QDoubleSpinBox.setRange(0.0, 1.0)
        self.blank_width_space_value_QDoubleSpinBox.setDecimals(2)
        self.blank_width_space_value_QDoubleSpinBox.setSingleStep(0.01)
        self.blank_width_space_value_QDoubleSpinBox.setValue(0.1)
        add_params_layout.addWidget(self.blank_width_space_value_QDoubleSpinBox, 3, 2, 1, 1)

        add_params_layout.addWidget(QLabel("Blank height space between subplots"), 3, 3, 1, 2)
        self.blank_height_space_value_QDoubleSpinBox = QDoubleSpinBox()
        self.blank_height_space_value_QDoubleSpinBox.setRange(0.0, 1.0)
        self.blank_height_space_value_QDoubleSpinBox.setDecimals(2)
        self.blank_height_space_value_QDoubleSpinBox.setSingleStep(0.01)
        self.blank_height_space_value_QDoubleSpinBox.setValue(0.1)
        add_params_layout.addWidget(self.blank_height_space_value_QDoubleSpinBox, 3, 5, 1, 1)

        add_params_layout.setRowMinimumHeight(3, 50)

        add_params_groupBox.setLayout(add_params_layout)

        layout.addWidget(add_params_groupBox)

        # graphic parameters import and export

        graphic_params_io_groupBox = QGroupBox("Graphic parameters save/load")

        graphic_params_io_layout = QHBoxLayout()

        self.graphic_params_save_QPushButton = QPushButton("Save")
        self.graphic_params_save_QPushButton.clicked.connect(self.output_graphic_params_save)
        graphic_params_io_layout.addWidget(self.graphic_params_save_QPushButton)

        self.graphic_params_load_QPushButton = QPushButton("Load")
        self.graphic_params_load_QPushButton.clicked.connect(self.output_graphic_params_load)
        graphic_params_io_layout.addWidget(self.graphic_params_load_QPushButton)

        graphic_params_io_groupBox.setLayout(graphic_params_io_layout)

        layout.addWidget(graphic_params_io_groupBox)

        # output file parameters

        output_file_groupBox = QGroupBox(self.tr("Output file - available formats: tif, pdf, svg"))

        output_file_layout = QGridLayout()

        self.figure_outpath_QLineEdit = QLineEdit()
        output_file_layout.addWidget(self.figure_outpath_QLineEdit, 3, 0, 1, 1)

        self.figure_outpath_QPushButton = QPushButton(self.tr("Choose"))
        self.figure_outpath_QPushButton.clicked.connect(self.define_figure_outpath)
        output_file_layout.addWidget(self.figure_outpath_QPushButton, 3, 1, 1, 1)

        output_file_groupBox.setLayout(output_file_layout)

        layout.addWidget(output_file_groupBox)

        # execution buttons

        decide_QWiget = QWidget()

        buttonLayout = QHBoxLayout()
        buttonLayout.addStretch()

        okButton = QPushButton("&OK")
        cancelButton = QPushButton("Cancel")

        buttonLayout.addWidget(okButton)
        buttonLayout.addWidget(cancelButton)

        decide_QWiget.setLayout(buttonLayout)

        layout.addWidget(decide_QWiget)

        self.setLayout(layout)

        self.connect(okButton, SIGNAL("clicked()"),
                     self, SLOT("accept()"))
        self.connect(cancelButton, SIGNAL("clicked()"),
                     self, SLOT("reject()"))

        self.setWindowTitle("Export figure")

    def output_graphic_params_save(self):

        output_file_path = new_file_path(self, "Define output configuration file", "*.txt", "txt")

        if not output_file_path:
            return

        out_configuration_string = """figure width = %f
resolution (dpi) = %d
font size (pts) = %f
top space = %f
left space = %f        
right space = %f        
bottom space = %f  
blank width space = %f
blank height space = %f""" % (float(self.figure_width_inches_QLineEdit.text()),
                              int(self.figure_resolution_dpi_QLineEdit.text()),
                              float(self.figure_fontsize_pts_QLineEdit.text()),
                              float(self.top_space_value_QDoubleSpinBox.value()),
                              float(self.left_space_value_QDoubleSpinBox.value()),
                              float(self.right_space_value_QDoubleSpinBox.value()),
                              float(self.bottom_space_value_QDoubleSpinBox.value()),
                              float(self.blank_width_space_value_QDoubleSpinBox.value()),
                              float(self.blank_height_space_value_QDoubleSpinBox.value()))

        with open(output_file_path, "w") as ofile:
            ofile.write(out_configuration_string)

        info(self,
             self.plugin_name,
             "Graphic parameters saved")

    def output_graphic_params_load(self):

        input_file_path = old_file_path(self, "Choose input configuration file", "*.txt", "txt")

        if not input_file_path:
            return

        with open(input_file_path, "r") as ifile:
            config_lines = ifile.readlines()

        try:
            figure_width_inches = float(config_lines[0].split("=")[1])
            figure_resolution_dpi = int(config_lines[1].split("=")[1])
            figure_fontsize_pts = float(config_lines[2].split("=")[1])
            top_space_value = float(config_lines[3].split("=")[1])
            left_space_value = float(config_lines[4].split("=")[1])
            right_space_value = float(config_lines[5].split("=")[1])
            bottom_space_value = float(config_lines[6].split("=")[1])
            blank_width_space = float(config_lines[7].split("=")[1])
            blank_height_space = float(config_lines[8].split("=")[1])
        except:
            warn(self,
                 self.plugin_name,
                 "Error in configuration file")
            return

        self.figure_width_inches_QLineEdit.setText(str(figure_width_inches))
        self.figure_resolution_dpi_QLineEdit.setText(str(figure_resolution_dpi))
        self.figure_fontsize_pts_QLineEdit.setText(str(figure_fontsize_pts))
        self.top_space_value_QDoubleSpinBox.setValue(top_space_value)
        self.left_space_value_QDoubleSpinBox.setValue(left_space_value)
        self.right_space_value_QDoubleSpinBox.setValue(right_space_value)
        self.bottom_space_value_QDoubleSpinBox.setValue(bottom_space_value)
        self.blank_width_space_value_QDoubleSpinBox.setValue(blank_width_space)
        self.blank_height_space_value_QDoubleSpinBox.setValue(blank_height_space)

    def define_figure_outpath(self):

        outfile_path = new_file_path(self, "Create", "", "Images (*.svg *.pdf *.tif)")

        self.figure_outpath_QLineEdit.setText(outfile_path)


class TopographicProfileExportDialog(QDialog):

    def __init__(self, plugin_name, selected_dem_params, parent=None):

        super(TopographicProfileExportDialog, self).__init__(parent)

        self.plugin_name = plugin_name

        layout = QVBoxLayout()

        ##
        # Profile source

        source_groupBox = QGroupBox(self.tr("Profile sources"))

        source_layout = QGridLayout()

        self.src_allselecteddems_QRadioButton = QRadioButton(self.tr("All selected DEMs"))
        source_layout.addWidget(self.src_allselecteddems_QRadioButton, 1, 0, 1, 2)
        self.src_allselecteddems_QRadioButton.setChecked(True)

        self.src_singledem_QRadioButton = QRadioButton(self.tr("Single DEM"))
        source_layout.addWidget(self.src_singledem_QRadioButton, 2, 0, 1, 1)

        self.src_singledemlist_QComboBox = QComboBox()
        selected_dem_layers = [dem_param.layer for dem_param in selected_dem_params]
        for qgsRasterLayer in selected_dem_layers:
            self.src_singledemlist_QComboBox.addItem(qgsRasterLayer.name())
        source_layout.addWidget(self.src_singledemlist_QComboBox, 2, 1, 1, 1)

        self.src_singlegpx_QRadioButton = QRadioButton(self.tr("GPX file"))
        source_layout.addWidget(self.src_singlegpx_QRadioButton, 3, 0, 1, 1)

        source_groupBox.setLayout(source_layout)

        layout.addWidget(source_groupBox)

        ##
        # Output type

        output_type_groupBox = QGroupBox(self.tr("Output format"))

        output_type_layout = QGridLayout()

        self.outtype_shapefile_point_QRadioButton = QRadioButton(self.tr("shapefile - point"))
        output_type_layout.addWidget(self.outtype_shapefile_point_QRadioButton, 0, 0, 1, 1)
        self.outtype_shapefile_point_QRadioButton.setChecked(True)

        self.outtype_shapefile_line_QRadioButton = QRadioButton(self.tr("shapefile - line"))
        output_type_layout.addWidget(self.outtype_shapefile_line_QRadioButton, 1, 0, 1, 1)

        self.outtype_csv_QRadioButton = QRadioButton(self.tr("csv"))
        output_type_layout.addWidget(self.outtype_csv_QRadioButton, 2, 0, 1, 1)

        output_type_groupBox.setLayout(output_type_layout)

        layout.addWidget(output_type_groupBox)

        ##
        # Output name/path

        output_path_groupBox = QGroupBox(self.tr("Output file"))

        output_path_layout = QGridLayout()

        self.outpath_QLineEdit = QLineEdit()
        output_path_layout.addWidget(self.outpath_QLineEdit, 0, 0, 1, 1)

        self.outpath_QPushButton = QPushButton("....")
        self.outpath_QPushButton.clicked.connect(self.define_outpath)
        output_path_layout.addWidget(self.outpath_QPushButton, 0, 1, 1, 1)

        self.load_output_checkBox = QCheckBox("load output shapefile in project")
        self.load_output_checkBox.setChecked(True)
        output_path_layout.addWidget(self.load_output_checkBox, 1, 0, 1, 2)

        output_path_groupBox.setLayout(output_path_layout)

        layout.addWidget(output_path_groupBox)

        decide_QWiget = QWidget()

        buttonLayout = QHBoxLayout()
        buttonLayout.addStretch()

        okButton = QPushButton("&OK")
        cancelButton = QPushButton("Cancel")

        buttonLayout.addWidget(okButton)
        buttonLayout.addWidget(cancelButton)

        decide_QWiget.setLayout(buttonLayout)

        layout.addWidget(decide_QWiget)

        self.setLayout(layout)

        self.connect(okButton, SIGNAL("clicked()"),
                     self, SLOT("accept()"))
        self.connect(cancelButton, SIGNAL("clicked()"),
                     self, SLOT("reject()"))

        self.setWindowTitle("Export topographic profile")

    def define_outpath(self):

        if self.outtype_shapefile_line_QRadioButton.isChecked() or self.outtype_shapefile_point_QRadioButton.isChecked():
            outfile_path = new_file_path(self, "Save file", "", "Shapefiles (*.shp)")
        elif self.outtype_csv_QRadioButton.isChecked():
            outfile_path = new_file_path(self, "Save file", "", "Csv (*.csv)")
        else:
            warn(self,
                 self.plugin_name,
                 self.tr("Output type definiton error"))
            return

        self.outpath_QLineEdit.setText(outfile_path)


class PointDataExportDialog(QDialog):

    def __init__(self, plugin_name, parent=None):

        super(PointDataExportDialog, self).__init__(parent)

        self.plugin_name = plugin_name

        layout = QVBoxLayout()

        ##
        # Output type

        output_type_groupBox = QGroupBox(self.tr("Output format"))

        output_type_layout = QGridLayout()

        self.outtype_shapefile_point_QRadioButton = QRadioButton(self.tr("shapefile - point"))
        output_type_layout.addWidget(self.outtype_shapefile_point_QRadioButton, 0, 0, 1, 1)
        self.outtype_shapefile_point_QRadioButton.setChecked(True)

        self.outtype_csv_QRadioButton = QRadioButton(self.tr("csv"))
        output_type_layout.addWidget(self.outtype_csv_QRadioButton, 1, 0, 1, 1)

        output_type_groupBox.setLayout(output_type_layout)

        layout.addWidget(output_type_groupBox)

        ##
        # Output name/path

        output_path_groupBox = QGroupBox(self.tr("Output path"))

        output_path_layout = QGridLayout()

        self.outpath_QLineEdit = QLineEdit()
        output_path_layout.addWidget(self.outpath_QLineEdit, 0, 0, 1, 1)

        self.outpath_QPushButton = QPushButton(self.tr("Choose"))
        self.outpath_QPushButton.clicked.connect(self.define_outpath)
        output_path_layout.addWidget(self.outpath_QPushButton, 0, 1, 1, 1)

        self.load_output_checkBox = QCheckBox("load output shapefile in project")
        self.load_output_checkBox.setChecked(True)
        output_path_layout.addWidget(self.load_output_checkBox, 1, 0, 1, 2)

        output_path_groupBox.setLayout(output_path_layout)

        layout.addWidget(output_path_groupBox)

        decide_QWiget = QWidget()

        buttonLayout = QHBoxLayout()
        buttonLayout.addStretch()

        okButton = QPushButton("&OK")
        cancelButton = QPushButton("Cancel")

        buttonLayout.addWidget(okButton)
        buttonLayout.addWidget(cancelButton)

        decide_QWiget.setLayout(buttonLayout)

        layout.addWidget(decide_QWiget)

        self.setLayout(layout)

        self.connect(okButton, SIGNAL("clicked()"),
                     self, SLOT("accept()"))
        self.connect(cancelButton, SIGNAL("clicked()"),
                     self, SLOT("reject()"))

        self.setWindowTitle("Export geological attitudes")

    def define_outpath(self):

        if self.outtype_shapefile_point_QRadioButton.isChecked():
            outfile_path = new_file_path(self, "Path", "*.shp", "Shapefile")
        elif self.outtype_csv_QRadioButton.isChecked():
            outfile_path = new_file_path(self, "Path", "*.csv", "Csv")
        else:
            warn(self,
                 self.plugin_name,
                 self.tr("Output type definiton error"))
            return

        self.outpath_QLineEdit.setText(outfile_path)


class LineDataExportDialog(QDialog):

    def __init__(self, plugin_name, parent=None):

        super(LineDataExportDialog, self).__init__(parent)

        self.plugin_name = plugin_name

        layout = QVBoxLayout()

        ##
        # Output type

        output_type_groupBox = QGroupBox(self.tr("Output format"))

        output_type_layout = QGridLayout()

        self.outtype_shapefile_line_QRadioButton = QRadioButton(self.tr("shapefile - line"))
        output_type_layout.addWidget(self.outtype_shapefile_line_QRadioButton, 0, 0, 1, 1)
        self.outtype_shapefile_line_QRadioButton.setChecked(True)

        self.outtype_csv_QRadioButton = QRadioButton(self.tr("csv"))
        output_type_layout.addWidget(self.outtype_csv_QRadioButton, 0, 1, 1, 1)

        output_type_groupBox.setLayout(output_type_layout)

        layout.addWidget(output_type_groupBox)

        ##
        # Output name/path

        output_path_groupBox = QGroupBox(self.tr("Output file"))

        output_path_layout = QGridLayout()

        self.outpath_QLineEdit = QLineEdit()
        output_path_layout.addWidget(self.outpath_QLineEdit, 0, 0, 1, 1)

        self.outpath_QPushButton = QPushButton("....")
        self.outpath_QPushButton.clicked.connect(self.define_outpath)
        output_path_layout.addWidget(self.outpath_QPushButton, 0, 1, 1, 1)

        self.load_output_checkBox = QCheckBox("load output in project")
        self.load_output_checkBox.setChecked(True)
        output_path_layout.addWidget(self.load_output_checkBox, 1, 0, 1, 2)

        output_path_groupBox.setLayout(output_path_layout)

        layout.addWidget(output_path_groupBox)

        decide_QWiget = QWidget()

        buttonLayout = QHBoxLayout()
        buttonLayout.addStretch()

        okButton = QPushButton("&OK")
        cancelButton = QPushButton("Cancel")

        buttonLayout.addWidget(okButton)
        buttonLayout.addWidget(cancelButton)

        decide_QWiget.setLayout(buttonLayout)

        layout.addWidget(decide_QWiget)

        self.setLayout(layout)

        self.connect(okButton, SIGNAL("clicked()"),
                     self, SLOT("accept()"))
        self.connect(cancelButton, SIGNAL("clicked()"),
                     self, SLOT("reject()"))

        self.setWindowTitle("Export")

    def define_outpath(self):

        if self.outtype_shapefile_line_QRadioButton.isChecked():
            outfile_path = new_file_path(self, "Save file", "", "Shapefiles (*.shp)")
        elif self.outtype_csv_QRadioButton.isChecked():
            outfile_path = new_file_path(self, "Save file", "", "Csv (*.csv)")
        else:
            warn(self,
                 self.plugin_name,
                 self.tr("Output type definiton error"))
            return

        self.outpath_QLineEdit.setText(outfile_path)


class StatisticsDialog(QDialog):

    def __init__(self, plugin_name, topo_profiles, parent=None):

        super(StatisticsDialog, self).__init__(parent)

        self.plugin_name = plugin_name

        profiles_stats = zip(topo_profiles.names,
                             zip(topo_profiles.statistics_elev,
                                 topo_profiles.statistics_dirslopes,
                                 topo_profiles.statistics_slopes))

        layout = QVBoxLayout()

        self.text_widget = QTextEdit()
        self.text_widget.setReadOnly(True)

        stat_report = "General statistics\n"
        stat_report += "\nProfile length: %f\n" % topo_profiles.profile_length
        stat_report += "\nTopographic elevations\n"
        stat_report += " - min: {}\n".format(topo_profiles.natural_elev_range[0])
        stat_report += " - max: {}\n\n".format(topo_profiles.natural_elev_range[1])
        stat_report += self.report_stats(profiles_stats)

        self.text_widget.setPlainText(stat_report)

        layout.addWidget(self.text_widget)

        self.setLayout(layout)

        self.setWindowTitle("Statistics")

    def report_stats(self, profiles_stats):

        def type_report(values):

            type_report = 'min: %s\n' % (values['min'])
            type_report += 'max: %s\n' % (values['max'])
            type_report += 'mean: %s\n' % (values['mean'])
            type_report += 'variance: %s\n' % (values['var'])
            type_report += 'standard deviation: %s\n\n' % (values['std'])

            return type_report

        report = 'Dataset statistics\n'
        types = ['elevations', 'directional slopes', 'absolute slopes']
        for name, stats in profiles_stats:
            report += '\ndataset name\n%s\n\n' % name
            for type, stat_val in zip(types, stats):
                report += '%s\n\n' % type
                report += type_report(stat_val)

        return report
