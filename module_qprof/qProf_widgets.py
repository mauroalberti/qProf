
from qgis.PyQt import uic

from ..qygsf.geometries.grids.statistics import *
from ..qygsf.io.profiles import *
from ..qygsf.utils.string_utils.utils_string import *
from ..qygsf.utils.qgis_utils.canvas import *
from ..qygsf.utils.qt_utils.filesystem import *
from ..qygsf.utils.qgis_utils.lines import *
from ..qygsf.utils.qgis_utils.messages import *
from ..qygsf.utils.qgis_utils.rasters import *

from .qProf_plotting import *
from .qProf_export import *

max_parallel_profiles_number = 20


class ActionWidget(QWidget):

    def __init__(self,
                 current_directory,
                 plugin_name,
                 canvas
                 ):

        super(ActionWidget, self).__init__()

        self.plugin_name = plugin_name
        self.canvas = canvas

        self.current_directory = current_directory
        uic.loadUi(f"{self.current_directory}/ui/choices_treewidget.ui", self)

        self.settings = QSettings("alberese", self.plugin_name)
        self.settings_gpxdir_key = "gpx/last_used_dir"

        self.actions_qtreewidget = self.actionsTreeWidget

        self.profile_operations = {
            "Load from line layer": self.define_profile_lines_from_line_layer,
            "digitize trace": self.digitize_rubberband_line,
            "clear trace": self.clear_rubberband_line,
            "save trace": self.save_rubberband_line,
            "Define in text window": self.define_track_source_from_text_window,
            "Read from GPX file": self.define_track_source_from_gpx_file,
            "From DEMs/grids": self.elevations_from_dems,
            "From GPX file": self.elevations_from_gpx,
            "create single profile": self.plot_grids_profile,
            "define parallel profiles parameters": self.define_parallel_profiles,
            "plot profiles": self.plot_parallel_profiles,
        }

        self.actions_qtreewidget.itemDoubleClicked.connect(self.activate_action_window)

        self.init_parameters()

        QgsProject.instance().crsChanged.connect(self.init_parameters)
        QgsProject.instance().crsChanged.connect(self.warn_of_reset)

    def warn_of_reset(self):

        warn(
            self.plugin_name,
            "Project CRS has been changed. Please redefine parameters"
        )

    def init_parameters(self):

        self.parallel_profiles_params = None

        self.profile_track_source_type = TrackSource.UNDEFINED
        self.invert_line_profile = False
        self.line_from_points_list = None
        self.input_gpx_file_path = None
        self.gpx_choice = GPXElevationUsage.NOT_USED
        self.gpx_track_name = None

        self.profile_name = None
        self.profile_line_list = None  # list of Lines, in the project CRS, undensified
        self.input_geoprofiles_set = GeoProfilesSet_()  # main instance for the geoprofiles
        self.profile_windows = []  # used to maintain alive the plots, i.e. to avoid the C++ objects being destroyed
        self.selected_dems = None
        self.selected_dem_parameters = None
        self.gpx_name = ''

    def activate_action_window(self):

        current_item_text = self.actions_qtreewidget.currentItem().text(0)

        operation = self.profile_operations.get(current_item_text)

        if operation is not None:
            operation()

    def define_profile_lines_from_line_layer(self
                                             ):
        """
        Should define:
         - source type -> self.profile_track_source = TrackSource.LINE_LAYER
         - list of undensified, inverted-in-case, CRS-projected lines

        """

        self.clear_rubberband_line()

        current_line_layers = loaded_line_layers()

        if len(current_line_layers) == 0:
            warn(
                self.plugin_name,
                "No available line layers"
            )
            return

        dialog = SourceLineLayerDialog(
            self.plugin_name,
            current_line_layers
        )

        if dialog.exec_():
            line_qgsvectorlayer, invert_profile, order_field_ndx = self.line_layer_params(dialog)
        else:
            warn(
                self.plugin_name,
                "No defined line source"
            )
            return

        line_order_fld_ndx = int(order_field_ndx) - 1 if order_field_ndx else None

        success, result = try_load_line_layer(
            line_layer=line_qgsvectorlayer,
            project_crs=projectCrs(),
            line_order_fld_ndx=line_order_fld_ndx,
            invert_direction=invert_profile
        )

        if not success:
            msg = result
            warn(
                self.plugin_name,
                f"Line layer not read: {msg}"
            )
            return

        self.profile_name = line_qgsvectorlayer.sourceName()
        self.profile_line_list = result
        self.profile_track_source_type = TrackSource.LINE_LAYER

        ok(
            self.plugin_name,
            "Line layer read"
        )

    def try_get_point_list(self,
                           dialog
                           ) -> Tuple[bool, Union[str, Line2D]]:

        try:

            raw_point_string = dialog.point_list_qtextedit.toPlainText()
            raw_point_list = raw_point_string.split("\n")
            raw_point_list = [clean_string(str(unicode_txt)) for unicode_txt in raw_point_list]
            data_list = [rp for rp in raw_point_list if rp != ""]

            point_list = [to_float(xy_pair.split(",")) for xy_pair in data_list]
            line_2d = xytuple_list_to_line2d(point_list)

            return True, line_2d

        except Exception as e:

            return False, str(e)

    def define_track_source_from_text_window(self):

        self.clear_rubberband_line()

        self.init_topo_labels()

        dialog = LoadPointListDialog(self.plugin_name)

        if dialog.exec_():

            success, result = self.try_get_point_list(dialog)
            if not success:
                msg = result
                warn(
                    self.plugin_name,
                    msg
                )
                return
            line2d = result

        else:

            warn(
                self.plugin_name,
                "No defined text_window source"
            )
            return

        try:

            npts = line2d.num_pts()
            if npts < 2:
                warn(
                    self.plugin_name,
                    "Defined text-window source with less than two points"
                )
                return

        except:

            warn(
                self.plugin_name,
                "No defined text-window source"
            )
            return

        self.profile_name = "Text input"
        self.profile_line_list = [line2d]
        self.profile_track_source_type = TrackSource.POINT_LIST

        ok(
            self.plugin_name,
            "Point list read"
        )

    def define_track_source_from_gpx_file(self):

        self.clear_rubberband_line()

        dialog = GpxInputDialog(
            self.plugin_name,
            self.settings,
            self.settings_gpxdir_key
        )

        if dialog.exec_():

            self.input_gpx_file_path = str(dialog.input_gpx_file_path.text())
            self.invert_line_profile = dialog.invert_track_direction.isChecked()

            success, results = try_extract_track_from_gpxfile(
                source_gpx_path=self.input_gpx_file_path,
                invert_profile=self.invert_line_profile
            )

            if not success:
                msg = results
                warn(
                    self.plugin_name,
                    msg
                )
                return

            name, line4d = results

            self.profile_name = os.path.basename(self.input_gpx_file_path)
            self.profile_line_list = [line4d]
            self.gpx_track_name = name
            self.profile_track_source_type = TrackSource.GPX_FILE

            ok(
                self.plugin_name,
                "GPX file read"
            )

        else:

            return

    def get_dem_parameters(self,
                           dem: QgsRasterLayer):

        success, results = try_raster_qgis_params(dem)

        if not success:
            msg = f"Error with {dem.name()} as source"
            warn(
                self.plugin_name,
                msg
            )
            return

        return QGisRasterParameters(*results)

    def get_selected_dems_params(self,
                                 dialog):

        selected_dems = []
        for dem_qgis_ndx in range(dialog.listDEMs_treeWidget.topLevelItemCount()):
            curr_DEM_item = dialog.listDEMs_treeWidget.topLevelItem(dem_qgis_ndx)
            if curr_DEM_item.checkState(0) == 2:
                selected_dems.append(dialog.singleband_raster_layers_in_project[dem_qgis_ndx])

        return selected_dems

    def elevations_from_dems(self):

        self.selected_dems = None
        self.selected_dem_parameters = []

        current_raster_layers = loaded_monoband_raster_layers()
        if len(current_raster_layers) == 0:
            warn(
                self.plugin_name,
                "No loaded DEM"
            )
            return

        dialog = SourceDEMsDialog(
            self.plugin_name,
            current_raster_layers
        )

        if dialog.exec_():
            selected_dems = self.get_selected_dems_params(dialog)
        else:
            warn(
                self.plugin_name,
                "No chosen DEM"
            )
            return

        if len(selected_dems) == 0:
            warn(
                self.plugin_name,
                "No selected DEM",
            )
            return
        else:
            self.selected_dems = selected_dems

        # get geodata

        self.selected_dem_parameters = [self.get_dem_parameters(dem) for dem in selected_dems]

        ok(
            self.plugin_name,
            "DEMs read" if len(selected_dems) >= 1 else "DEM read"
        )

    def elevations_from_gpx(self):

        dialog = GPXElevationDialog(
            self.plugin_name
        )

        if dialog.exec_():
            self.gpx_choice = dialog.get_gpx_choice()

    def line_layer_params(
            self,
            dialog):

        line_layer = dialog.line_shape
        invert_profile = dialog.qcbxInvertProfile.isChecked()
        order_field_ndx = dialog.Trace2D_order_field_comboBox.currentIndex()

        return line_layer, invert_profile, order_field_ndx

    def init_topo_labels(self):
        """
        Initialize topographic label and order parameters.

        :return:
        """

        self.profiles_labels = None
        self.profiles_order = None

    def check_pre_profile(self):

        for geoprofile in self.input_geoprofiles_set.geoprofiles:
            if not geoprofile._named_lines.statistics_calculated:
                warn(
                    self.plugin_name,
                    "Profile statistics not yet calculated"
                )
                return False

        return True

    def calculate_profile_statistics(self,
        geoprofiles
    ):

        for geoprofile in geoprofiles:

            for name, line3d in geoprofile._named_lines:

                statistics_elev = [get_statistics(p) for p in line3d.z_array()]
                statistics_dirslopes = [get_statistics(p) for p in line3d.dir_slopes()]
                statistics_slopes = [get_statistics(p) for p in np.absolute(line3d.dir_slopes())]

                profile_length = line3d.incremental_length_2d()[-1] - line3d.incremental_length_2d()[0]
                natural_elev_range = (
                    np.nanmin(np.array([ds_stats["min"] for ds_stats in statistics_elev])),
                    np.nanmax(np.array([ds_stats["max"] for ds_stats in statistics_elev])))

                statistics_calculated = True

        """
        dialog = StatisticsDialog(
            self.plugin_name,
            self.input_geoprofiles
        )

        dialog.exec_()
        """

    """
    def try_prepare_single_topo_profiles(self
        ) -> Tuple[bool, str]:

        self.input_geoprofiles = GeoProfilesSet()  # reset any previous created profiles

        self.demline_source = "dem_source"
        topo_source_type = self.demline_source
        if topo_source_type == self.demline_source:

            try:

                selected_dems = self.selected_dems
                selected_dem_parameters = self.selected_dem_parameters

            except Exception as e:

                return False, f"Input DEMs definition not correct: {e}"

            try:

                print(f"DEBUG: I am at the start")

                # get DEMs resolutions in project CRS and choose the min value

                dem_resolutions_prj_crs_list = []
                for dem, dem_params in zip(self.selected_dems, self.selected_dem_parameters):
                    dem_resolutions_prj_crs_list.append(
                        get_dem_resolution_in_prj_crs(
                            dem,
                            dem_params,
                            projectCrs())
                    )

                min_dem_resolution = min(dem_resolutions_prj_crs_list)

                if min_dem_resolution > 1:
                    sample_distance = round(min_dem_resolution)
                else:
                    sample_distance = min_dem_resolution

            except Exception as e:

                return False, f"Sample distance value not correct: {e}"

            if self.qcbxDigitizeLineSource.isChecked():
                if self.digitized_profile_line2dt is None or \
                   self.digitized_profile_line2dt.num_pts < 2:
                    wrn(self,
                         self.plugin_name,
                         "No digitized line available")
                    return
                else:
                    source_profile_lines = [self.digitized_profile_line2dt]
            else:

            try:

                source_profile_lines = self.profile_lines

            except:

                return False, "DEM-line profile source not correctly created [1]"

            if source_profile_lines is None:

                return False, "DEM-line profile source not correctly created [2]"

        elif topo_source_type == self.gpxfile_source:
            stop_rubberband()
            try:
                source_gpx_path = str(self.qlneInputGPXFile.text())
                if source_gpx_path == '':
                    wrn(self,
                         self.plugin_name,
                         "Source GPX file is not set")
                    return
            except Exception as e:
                wrn(self,
                     self.plugin_name,
                     "Source GPX file not correctly set: {}".format(e))
                return

        else:
            wrn(self,
                 self.plugin_name,
                 "Debug: uncorrect type source for topo sources def")
            return

        # calculates profiles

        if topo_source_type == self.demline_source:  # sources are DEM(s) and line

            # check total number of points in line(s) to create
            estimated_total_num_pts = 0
            for profile_line in source_profile_lines:

                profile_length = profile_line.length_2d
                profile_num_pts = profile_length / sample_distance
                estimated_total_num_pts += profile_num_pts

            estimated_total_num_pts = int(ceil(estimated_total_num_pts))

            if estimated_total_num_pts > pt_num_threshold:

                return False, f"There are {estimated_total_num_pts} estimated points (limit is {pt_num_threshold}) in profile(s) to create.\nTry increasing sample distance value"

            for profile_line in source_profile_lines:

                try:

                    topo_profiles = topoprofiles_from_dems(
                        source_profile_line=profile_line,
                        sample_distance=sample_distance,
                        selected_dems=selected_dems,
                        selected_dem_parameters=selected_dem_parameters
                    )

                except Exception as e:

                    return False, f"Error with data source read: {e}"

                if topo_profiles is None:

                    return False, "Debug: profile not created"

                geoprofile = GeoProfile()
                geoprofile.source_data_type = topo_source_type
                geoprofile.original_line = profile_line
                geoprofile.sample_distance = sample_distance
                geoprofile.set_topo_profiles(topo_profiles)

                self.input_geoprofiles.append(geoprofile)

        elif topo_source_type == self.gpxfile_source:  # source is GPX file

            try:
                topo_profiles = topoprofiles_from_gpxfile(source_gpx_path,
                                                          invert_profile,
                                                          self.gpxfile_source)
            except Exception as e:
                wrn(self,
                     self.plugin_name,
                     "Error with profile calculation from GPX file: {}".format(e))
                return

            if topo_profiles is None:
                wrn(self,
                     self.plugin_name,
                     "Debug: profile not created")
                return

            geoprofile = GeoProfile()
            geoprofile.source_data_type = topo_source_type
            geoprofile.original_line = source_profile_lines
            geoprofile.sample_distance = sample_distance
            geoprofile.set_topo_profiles(topo_profiles)

            self.input_geoprofiles.append(geoprofile)

        else:  # source error
            error(self,
                  self.plugin_name,
                 "Debug: profile calculation not defined")
            return

        return True, ""
    """

    def get_profile_plot_params(self,
        dialog,
        surface_names: List[str]
    ) -> Dict:

        profile_params = {}

        # get profile plot parameters

        try:
            profile_params['plot_min_elevation_user'] = float(dialog.qledtPlotMinValue.text())
        except:
            profile_params['plot_min_elevation_user'] = None

        try:
            profile_params['plot_max_elevation_user'] = float(dialog.qledtPlotMaxValue.text())
        except:
            profile_params['plot_max_elevation_user'] = None

        profile_params['set_vertical_exaggeration'] = dialog.qcbxSetVerticalExaggeration.isChecked()
        try:
            profile_params['vertical_exaggeration'] = float(dialog.qledtDemExagerationRatio.text())
            assert profile_params['vertical_exaggeration'] > 0
        except:
            profile_params['vertical_exaggeration'] = 1

        profile_params['filled_height'] = dialog.qcbxPlotFilledHeight.isChecked()
        profile_params['filled_slope'] = dialog.qcbxPlotFilledSlope.isChecked()
        profile_params['plot_height_choice'] = dialog.qcbxPlotProfileHeight.isChecked()
        profile_params['plot_slope_choice'] = dialog.qcbxPlotProfileSlope.isChecked()
        profile_params['plot_slope_absolute'] = dialog.qrbtPlotAbsoluteSlope.isChecked()
        profile_params['plot_slope_directional'] = dialog.qrbtPlotDirectionalSlope.isChecked()
        profile_params['invert_xaxis'] = dialog.qcbxInvertXAxisProfile.isChecked()

        if hasattr(dialog, 'visible_elevation_layers') and dialog.visible_elevation_layers is not None:
            profile_params['visible_elev_lyrs'] = dialog.visible_elevation_layers
        else:
            profile_params['visible_elev_lyrs'] = surface_names

        if hasattr(dialog, 'elevation_layer_colors') and dialog.elevation_layer_colors is not None:
            profile_params['elev_lyr_colors'] = dialog.elevation_layer_colors
        else:
            profile_params['elev_lyr_colors'] = [QColor('red')] * len(surface_names)

        return profile_params

    def plot_grids_profile(self):

        if self.profile_track_source_type == TrackSource.UNDEFINED:
            warn(
                self.plugin_name,
                "No profile track source defined"
            )
            return

        success, result = try_prepare_grids_profile(
            profile_line=self.profile_line_list[0],
            track_source=self.profile_track_source_type,
            gpx_elevation_usage=self.gpx_choice,
            selected_grids=self.selected_dems,
            selected_grids_parameters=self.selected_dem_parameters,
            gpx_track_name=self.gpx_track_name
        )

        if not success:
            msg = result
            error(
                self.plugin_name,
                msg
            )
            return

        grids_profile = result

        print(f"DEBUG: grids_profile: {grids_profile}")

        profiles_min_elevs = []
        profiles_max_elevs = []
        profiles_lengths = []

        for _, line3d in grids_profile._named_lines:
            profiles_min_elevs.append(line3d.z_min())
            profiles_max_elevs.append(line3d.z_max())
            profiles_lengths.append(line3d.length_2d())

        surface_names = [name for name, _ in grids_profile._named_lines]

        if self.input_geoprofiles_set.plot_params is None:
            surface_colors = None
        else:
            surface_colors = self.input_geoprofiles_set.plot_params.get('elev_lyr_colors')

        # pre-process input data to account for multi.profiles

        profile_length = np.nanmax(profiles_lengths)
        natural_elev_min = np.nanmin(profiles_min_elevs)
        natural_elev_max = np.nanmax(profiles_max_elevs)

        print(f"Profiles max length: {profile_length}")
        print(f"Profiles min elevation: {natural_elev_min}")
        print(f"Profiles max elevation: {natural_elev_max}")

        if np.isnan(profile_length) or profile_length == 0.0:
            error(
                self.plugin_name,
                f"Max profile length is {profile_length}.\nCheck profile trace."
            )
            return

        if np.isnan(natural_elev_min) or np.isnan(natural_elev_max):
            error(
                self.plugin_name,
                f"Max elevation in profile(s) is {natural_elev_max} and min is {natural_elev_min}.\nCheck profile trace location vs. DEM(s). [#1]"
            )
            return

        if natural_elev_max <= natural_elev_min:
            error(
                self.plugin_name,
                "Error: min elevation larger then max elevation"
            )
            return

        dialog = PlotTopoProfileDialog(self.plugin_name,
                                       profiles_lengths,
                                       profiles_min_elevs,
                                       profiles_max_elevs,
                                       surface_names,
                                       surface_colors)

        if dialog.exec_():
            plot_params = self.get_profile_plot_params(
                dialog,
                surface_names
            )
        else:
            return

        profiles_created = True

        # plot profiles

        profile_window = plot_gridsprofile(
            named_grids_profile=grids_profile,
            plot_params=plot_params
        )

        self.profile_windows.append(profile_window)

    def define_parallel_profiles(self):

        dialog = ParallelProfilesParametersDialog(
            self.plugin_name
        )

        if dialog.exec_():

            profiles_spacing = dialog.spacing_wdgt.value()
            num_left_profiles = dialog.num_left_profiles_wdgt.value()
            num_right_profiles = dialog.num_right_profiles_wdgt.value()

        else:

            warn(
                self.plugin_name,
                "No parameters defined for parallel profiles "
            )
            return

        if profiles_spacing <= 0.0:
            warn(
                self.plugin_name,
                "Spacing between parallel profiles cannot be zero"
            )
            return
        elif num_left_profiles + num_right_profiles == 0:
            warn(
                self.plugin_name,
                "Total number of parallel profiles cannot be zero"
            )
            return
        elif num_left_profiles + num_right_profiles >= max_parallel_profiles_number:
            warn(
                self.plugin_name,
                f"Total number of parallel profiles ({num_left_profiles + num_right_profiles}) cannot be greater than {max_parallel_profiles_number}"
            )
            return
        else:
            ok(
                self.plugin_name,
                "Parallel profiles parameters defined"
            )

            self.parallel_profiles_params = {
                "spacing": profiles_spacing,
                "num_left_profiles": num_left_profiles,
                "num_right_profiles": num_right_profiles
            }

            print(self.parallel_profiles_params)

    def plot_parallel_profiles(self):

        if self.profile_track_source_type == TrackSource.UNDEFINED:
            warn(
                self.plugin_name,
                "No profile track source defined"
            )
            return

        if self.profile_track_source_type == TrackSource.GPX_FILE:
            warn(
                self.plugin_name,
                "Parallel profiles are not implemented for GPX-derived profile source"
            )
            return

        success, result = try_prepare_grids_profile(
            profile_line=self.profile_line_list[0],
            track_source=self.profile_track_source_type,
            gpx_elevation_usage=self.gpx_choice,
            selected_grids=self.selected_dems,
            selected_grids_parameters=self.selected_dem_parameters,
            gpx_track_name=self.gpx_track_name
        )

        if not success:
            msg = result
            error(
                self.plugin_name,
                msg
            )
            return

        grids_profile = result

        profiles_min_elevs = []
        profiles_max_elevs = []
        profiles_lengths = []

        for _, line3d in grids_profile._named_lines:
            profiles_min_elevs.append(line3d.z_min())
            profiles_max_elevs.append(line3d.z_max())
            profiles_lengths.append(line3d.length_2d())

        surface_names = [name for name, _ in grids_profile._named_lines]

        if self.input_geoprofiles_set.plot_params is None:
            surface_colors = None
        else:
            surface_colors = self.input_geoprofiles_set.plot_params.get('elev_lyr_colors')

        # pre-process input data to account for multi.profiles

        profile_length = np.nanmax(profiles_lengths)
        natural_elev_min = np.nanmin(profiles_min_elevs)
        natural_elev_max = np.nanmax(profiles_max_elevs)

        print(f"DEBUG: Profiles max length: {profile_length}")
        print(f"DEBUG: Profiles min elevation: {natural_elev_min}")
        print(f"DEBUG: Profiles max elevation: {natural_elev_max}")

        if np.isnan(profile_length) or profile_length == 0.0:
            error(
                self.plugin_name,
                f"Max profile length is {profile_length}.\nCheck profile trace."
            )
            return

        if np.isnan(natural_elev_min) or np.isnan(natural_elev_max):
            error(
                self.plugin_name,
                f"Max elevation in profile(s) is {natural_elev_max} and min is {natural_elev_min}.\nCheck profile trace location vs. DEM(s). [#2]"
            )
            return

        if natural_elev_max <= natural_elev_min:
            error(
                self.plugin_name,
                "Error: min elevation larger then max elevation"
            )
            return

        dialog = PlotTopoProfileDialog(self.plugin_name,
                                       profiles_lengths,
                                       profiles_min_elevs,
                                       profiles_max_elevs,
                                       surface_names,
                                       surface_colors)

        if dialog.exec_():
            plot_params = self.get_profile_plot_params(
                dialog,
                surface_names
            )
        else:
            return

        profiles_created = True

        # plot profiles

        plot_addit_params = dict()
        """
        plot_addit_params["add_trendplunge_label"] = self.plot_prj_add_trendplunge_label.isChecked()
        plot_addit_params["add_ptid_label"] = self.plot_prj_add_pt_id_label.isChecked()
        """
        plot_addit_params["add_trendplunge_label"] = False
        plot_addit_params["add_ptid_label"] = False

        plot_addit_params["polygon_class_colors"] = None  # self.polygon_classification_colors
        plot_addit_params["plane_attitudes_colors"] = None  # self.plane_attitudes_colors

        """
        num_subplots = plot_addit_params['plot_height_choice'] + plot_addit_params['plot_slope_choice']
        grid_spec = gridspec.GridSpec(num_subplots, 1)
        """

        geoprofile = GeoProfile()
        if grids_profile.num_lines() == 1:
            sz_profile = grids_profile.to_sz_arrays(0)
        else:
            raise Exception("Currently doesn't manage the plot of multiple grids in a geological profile")

        print(f"DEBUG: sz_profile min z: {sz_profile.y_min()}")
        print(f"DEBUG: sz_profile max z: {sz_profile.y_max()}")
        geoprofile.topo_profile = sz_profile
        profile_window = plot_gridsprofile(
            named_grids_profile=grids_profile,
            plot_params=plot_params
        )

        self.profile_windows.append(profile_window)

    def digitize_rubberband_line(self):

        self.clear_rubberband_line()

        self.previous_maptool = self.canvas.mapTool()  # Save the standard map tool for restoring it at the end

        info(
            self.plugin_name,
            "Now you can digitize the trace on the map.\nLeft click: add point\nRight click: end adding point"
        )

        self.rubberband = QgsRubberBand(self.canvas)
        self.rubberband.setWidth(2)
        self.rubberband.setColor(QColor(Qt.red))

        self.digitize_maptool = MapDigitizeTool(self.canvas)
        self.canvas.setMapTool(self.digitize_maptool)

        self.digitize_maptool.moved.connect(self.canvas_refresh_profile_line)
        self.digitize_maptool.leftClicked.connect(self.profile_add_point)
        self.digitize_maptool.rightClicked.connect(self.canvas_end_profile_line)

    def canvas_refresh_profile_line(self, position):

        """
        if len(self.profile_canvas_points) == 0:
            return
        """

        x, y = xy_from_canvas(self.canvas, position)
        #print(f"DEBUG: canvas_refresh_profile_line -> {x}, {y}")

        self.refresh_rubberband(self.profile_canvas_points + [[x, y]])

    def profile_add_point(self, position):

        x, y = xy_from_canvas(self.canvas, position)
        #print(f"DEBUG: profile_add_point -> {x} {y}")

        self.profile_canvas_points.append([x, y])

    def canvas_end_profile_line(self):

        self.refresh_rubberband(self.profile_canvas_points)

        self.line_from_digitation = None

        if len(self.profile_canvas_points) <= 1:
            warn(
                self.plugin_name,
                "At least two non-coincident points are required"
            )
            return

        raw_line = Line2D(
            [Point2D(x, y) for x, y in self.profile_canvas_points]).remove_coincident_points()

        if raw_line.num_pts() <= 1:
            warn(
                self.plugin_name,
                "Just one non-coincident point"
            )
            return

        self.profile_canvas_points = []
        self.restore_previous_map_tool()

        self.line_from_digitation = raw_line
        self.profile_name = "Digitized line"
        self.profile_line_list = [raw_line]
        self.profile_track_source_type = TrackSource.DIGITATION

    def restore_previous_map_tool(self):

        self.canvas.unsetMapTool(self.digitize_maptool)
        self.canvas.setMapTool(self.previous_maptool)

    def refresh_rubberband(self,
                           xy_list
                           ):

        self.rubberband.reset(QgsWkbTypes.LineGeometry)
        for x, y in xy_list:
            self.rubberband.addPoint(QgsPointXY(x, y))

    def clear_rubberband_line(self):

        self.profile_track_source_type = TrackSource.UNDEFINED

        self.profile_canvas_points = []
        self.line_from_digitation = None

        try:
            self.rubberband.reset()
        except:
            pass

    def save_rubberband_line(self):

        def output_profile_line(
                output_format,
                output_filepath,
                pts2dt,
                proj_sr
        ):

            points = [[n, pt2dt.x, pt2dt.y] for n, pt2dt in enumerate(pts2dt)]
            if output_format == "csv":
                success, msg = write_generic_csv(
                    output_filepath,
                    ['id', 'x', 'y'],
                    points
                )
                if not success:
                    warn(
                        self.plugin_name,
                        msg
                    )
            elif output_format == "shapefile - line":
                success, msg = write_rubberband_profile_lnshp(
                    output_filepath,
                    ['id'],
                    points,
                    proj_sr)
                if not success:
                    warn(
                        self.plugin_name,
                        msg
                    )
            else:
                error(
                    self.plugin_name,
                    "Debug: error in export format"
                )
                return

            if success:
                info(
                    self.plugin_name,
                    "Line saved"
                )

        def get_format_type():

            if dialog.outtype_shapefile_line_QRadioButton.isChecked():
                return "shapefile - line"
            elif dialog.outtype_csv_QRadioButton.isChecked():
                return "csv"
            else:
                return ""

        if self.line_from_digitation is None:

            warn(
                self.plugin_name,
                "No available line to save [1]"
            )
            return

        elif self.line_from_digitation.num_pts() < 2:

            warn(
                self.plugin_name,
                "No available line to save [2]"
            )
            return

        dialog = LineDataExportDialog(self.plugin_name)
        if dialog.exec_():
            output_format = get_format_type()
            if output_format == "":
                warn(
                    self.plugin_name,
                    "Error in output format"
                )
                return
            output_filepath = dialog.outpath_QLineEdit.text()
            if len(output_filepath) == 0:
                warn(
                    self.plugin_name,
                    "Error in output path"
                )
                return
            add_to_project = dialog.load_output_checkBox.isChecked()
        else:
            warn(
                self.plugin_name,
                "No export defined"
            )
            return

        # get project CRS information
        project_crs_osr = proj4str()

        output_profile_line(
            output_format,
            output_filepath,
            self.line_from_digitation.pts,
            project_crs_osr)

        # add theme to QGis project
        if output_format == "shapefile - line" and add_to_project:
            try:

                digitized_line_layer = QgsVectorLayer(output_filepath,
                                                      QFileInfo(output_filepath).baseName(),
                                                      "ogr")
                QgsProject.instance().addMapLayer(digitized_line_layer)

            except:

                error(
                    self.plugin_name,
                    "Unable to load layer in project"
                )
                return

    def closeEvent(self, evnt):

        self.clear_rubberband_line()

        """
        try:
            
            self.digitize_maptool.moved.disconnect(self.canvas_refresh_profile_line)
            self.digitize_maptool.leftClicked.disconnect(self.profile_add_point)
            self.digitize_maptool.rightClicked.disconnect(self.canvas_end_profile_line)
            self.restore_previous_map_tool()

        except:

            pass
        """


class SourceDEMsDialog(QDialog):

    def __init__(self,
                 plugin_name,
                 raster_layers,
                 parent=None
                 ):

        super(SourceDEMsDialog, self).__init__(parent)

        self.plugin_name = plugin_name

        self.singleband_raster_layers_in_project = raster_layers

        self.listDEMs_treeWidget = QTreeWidget()
        self.listDEMs_treeWidget.setColumnCount(2)
        self.listDEMs_treeWidget.headerItem().setText(0, "Select")
        self.listDEMs_treeWidget.headerItem().setText(1, "Name")
        self.listDEMs_treeWidget.setHorizontalScrollBarPolicy(Qt.ScrollBarAsNeeded)
        self.listDEMs_treeWidget.setDragEnabled(False)
        self.listDEMs_treeWidget.setDragDropMode(QAbstractItemView.NoDragDrop)
        self.listDEMs_treeWidget.setAlternatingRowColors(True)
        self.listDEMs_treeWidget.setSelectionMode(QAbstractItemView.SingleSelection)
        self.listDEMs_treeWidget.setTextElideMode(Qt.ElideLeft)

        self.populate_raster_layer_treewidget()

        self.listDEMs_treeWidget.resizeColumnToContents(0)
        self.listDEMs_treeWidget.resizeColumnToContents(1)

        okButton = QPushButton("&OK")
        cancelButton = QPushButton("Cancel")

        buttonLayout = QHBoxLayout()
        buttonLayout.addStretch()
        buttonLayout.addWidget(okButton)
        buttonLayout.addWidget(cancelButton)

        layout = QGridLayout()

        layout.addWidget(
            self.listDEMs_treeWidget,
            0, 0, 1, 3)
        layout.addLayout(
            buttonLayout,
            1, 0, 1, 3)

        self.setLayout(layout)

        okButton.clicked.connect(self.accept)
        cancelButton.clicked.connect(self.reject)

        self.setWindowTitle("Define source DEMs")

    def populate_raster_layer_treewidget(self):

        self.listDEMs_treeWidget.clear()

        for raster_layer in self.singleband_raster_layers_in_project:
            tree_item = QTreeWidgetItem(self.listDEMs_treeWidget)
            tree_item.setText(1, raster_layer.name())
            tree_item.setFlags(tree_item.flags() | Qt.ItemIsUserCheckable)
            tree_item.setCheckState(0, 0)


class GPXElevationDialog(QDialog):

    def __init__(self,
                 plugin_name,
                 parent=None
                 ):

        super(GPXElevationDialog, self).__init__(parent)

        self.plugin_name = plugin_name

        gpx_choices_layout = QVBoxLayout()

        self.doNotUseGPXElevations_qradiobutton = QRadioButton(
            "Do not use GPX elevations")
        self.useGPXElevationsWithDEMs_qradiobutton = QRadioButton(
            "Use GPX elevations together with DEM elevations")
        self.useOnlyGPXElevations_qradiobutton = QRadioButton(
            "Use only GPX elevations")
        self.doNotUseGPXElevations_qradiobutton.setChecked(True)

        gpx_choices_layout.addWidget(self.doNotUseGPXElevations_qradiobutton)
        gpx_choices_layout.addWidget(self.useGPXElevationsWithDEMs_qradiobutton)
        gpx_choices_layout.addWidget(self.useOnlyGPXElevations_qradiobutton)

        verticalSpacer = QSpacerItem(20, 30, QSizePolicy.Minimum, QSizePolicy.Expanding)
        gpx_choices_layout.addItem(verticalSpacer)

        okButton = QPushButton("&OK")
        cancelButton = QPushButton("Cancel")

        buttonLayout = QHBoxLayout()
        buttonLayout.addStretch()
        buttonLayout.addWidget(okButton)
        buttonLayout.addWidget(cancelButton)

        gpx_choices_layout.addLayout(
            buttonLayout)

        self.setLayout(gpx_choices_layout)

        okButton.clicked.connect(self.accept)
        cancelButton.clicked.connect(self.reject)

        self.setWindowTitle("Set GPX elevation role")

    def get_gpx_choice(self):

        if self.doNotUseGPXElevations_qradiobutton.isChecked():
            return GPXElevationUsage.NOT_USED
        elif self.useGPXElevationsWithDEMs_qradiobutton.isChecked():
            return GPXElevationUsage.USE_WITH_DEMS
        elif self.useOnlyGPXElevations_qradiobutton.isChecked():
            return GPXElevationUsage.ONLY_TO_USE
        else:
            raise Exception("Debug: not correctly set radio buttons for GPX elevation usage.")


class SourceLineLayerDialog(QDialog):

    def __init__(self,
                 plugin_name,
                 current_line_layers,
                 parent=None
                 ):

        super(SourceLineLayerDialog, self).__init__(parent)

        self.plugin_name = plugin_name
        self.current_line_layers = current_line_layers

        layout = QGridLayout()

        layout.addWidget(
            QLabel(self.tr("Input line layer:")),
            0, 0, 1, 1)

        self.LineLayers_comboBox = QComboBox()
        layout.addWidget(
            self.LineLayers_comboBox,
            0, 1, 1, 3)
        self.refresh_input_profile_layer_combobox()

        '''
        self.qrbtLineIsMultiProfile = QCheckBox(self.tr("Layer with multiple profiles ->"))
        layout.addWidget(
            self.qrbtLineIsMultiProfile,
            1, 0, 1, 2)

        layout.addWidget(
            QLabel(self.tr("category field:")),
            1, 2, 1, 1)
        self.Trace2D_label_field_comboBox = QComboBox()
        layout.addWidget(
            self.Trace2D_label_field_comboBox,
            1, 3, 1, 1)

        self.refresh_label_field_combobox()
        self.LineLayers_comboBox.currentIndexChanged[int].connect(self.refresh_label_field_combobox)
        '''

        self.qcbxInvertProfile = QCheckBox("Invert orientation")
        layout.addWidget(
            self.qcbxInvertProfile,
            1, 0, 1, 1)

        layout.addWidget(
            QLabel(self.tr("Line order field:")),
            2, 0, 1, 1)

        self.Trace2D_order_field_comboBox = QComboBox()
        layout.addWidget(
            self.Trace2D_order_field_comboBox,
            2, 1, 1, 3)

        self.refresh_order_field_combobox()

        self.LineLayers_comboBox.currentIndexChanged[int].connect(self.refresh_order_field_combobox)

        okButton = QPushButton("&OK")
        cancelButton = QPushButton("Cancel")

        buttonLayout = QHBoxLayout()
        buttonLayout.addStretch()
        buttonLayout.addWidget(okButton)
        buttonLayout.addWidget(cancelButton)

        layout.addLayout(
            buttonLayout,
            3, 0, 1, 3)

        self.setLayout(layout)

        okButton.clicked.connect(self.accept)
        cancelButton.clicked.connect(self.reject)

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

    def refresh_label_field_combobox(self):

        self.Trace2D_label_field_comboBox.clear()
        self.Trace2D_label_field_comboBox.addItem('--optional--')

        shape_qgis_ndx = self.LineLayers_comboBox.currentIndex()
        self.line_shape = self.current_line_layers[shape_qgis_ndx]

        line_layer_field_list = self.line_shape.dataProvider().fields().toList()
        for field in line_layer_field_list:
            self.Trace2D_label_field_comboBox.addItem(field.name())


class LoadPointListDialog(QDialog):

    def __init__(self, plugin_name, parent=None):

        super(LoadPointListDialog, self).__init__(parent)

        self.plugin_name = plugin_name
        layout = QGridLayout()

        layout.addWidget(
            QLabel(self.tr("Point list, with at least two points.")),
            0, 0, 1, 1)
        layout.addWidget(
            QLabel(self.tr("Each point is defined by a comma-separated, x-y coordinate pair (same CRS as project), one for each row")), 1, 0,
            1, 1)
        layout.addWidget(
            QLabel(self.tr("Example:\n549242.7, 242942.2\n578370.3, 322634.5")),
            2, 0, 1, 1)

        self.point_list_qtextedit = QTextEdit()
        layout.addWidget(
            self.point_list_qtextedit,
            3, 0, 1, 1)

        okButton = QPushButton("&OK")
        cancelButton = QPushButton("Cancel")

        buttonLayout = QHBoxLayout()
        buttonLayout.addStretch()
        buttonLayout.addWidget(okButton)
        buttonLayout.addWidget(cancelButton)

        layout.addLayout(
            buttonLayout,
            4, 0, 1, 3)

        self.setLayout(layout)

        okButton.clicked.connect(self.accept)
        cancelButton.clicked.connect(self.reject)

        self.setWindowTitle("Point list")


class ElevationLineStyleDialog(QDialog):

    def __init__(self,
                 plugin_name,
                 layer_names,
                 layer_colors,
                 parent=None
                 ):

        super(ElevationLineStyleDialog, self).__init__(parent)

        self.plugin_name = plugin_name

        self.qtwdElevationLayers = QTreeWidget()
        self.qtwdElevationLayers.setColumnCount(3)
        self.qtwdElevationLayers.headerItem().setText(0, "View")
        self.qtwdElevationLayers.headerItem().setText(1, "Name")
        self.qtwdElevationLayers.headerItem().setText(2, "Color")
        self.qtwdElevationLayers.setHorizontalScrollBarPolicy(Qt.ScrollBarAsNeeded)
        self.qtwdElevationLayers.setDragEnabled(False)
        self.qtwdElevationLayers.setDragDropMode(QAbstractItemView.NoDragDrop)
        self.qtwdElevationLayers.setAlternatingRowColors(True)
        self.qtwdElevationLayers.setSelectionMode(QAbstractItemView.SingleSelection)
        self.qtwdElevationLayers.setTextElideMode(Qt.ElideLeft)

        self.populate_elevation_layer_treewidget(layer_names, layer_colors)

        self.qtwdElevationLayers.resizeColumnToContents(0)
        self.qtwdElevationLayers.resizeColumnToContents(1)
        self.qtwdElevationLayers.resizeColumnToContents(2)

        okButton = QPushButton("&OK")
        cancelButton = QPushButton("Cancel")

        buttonLayout = QHBoxLayout()
        buttonLayout.addStretch()
        buttonLayout.addWidget(okButton)
        buttonLayout.addWidget(cancelButton)

        layout = QGridLayout()

        layout.addWidget(
            self.qtwdElevationLayers,
            0, 0, 1, 3)
        layout.addLayout(
            buttonLayout,
            1, 0, 1, 3)

        self.setLayout(layout)

        okButton.clicked.connect(self.accept)
        cancelButton.clicked.connect(self.reject)

        self.setWindowTitle("Define elevation line style")

    def populate_elevation_layer_treewidget(self,
                                            layer_names,
                                            layer_colors
                                            ):

        self.qtwdElevationLayers.clear()

        if layer_colors is None:
            num_available_colors = 0
        else:
            num_available_colors = len(layer_colors)

        for ndx, layer_name in enumerate(layer_names):
            tree_item = QTreeWidgetItem(self.qtwdElevationLayers)
            tree_item.setText(1, layer_name)
            color_button = QgsColorButton()
            if ndx < num_available_colors:
                color_button.setColor(layer_colors[ndx])
            else:
                color_button.setColor(QColor('red'))
            self.qtwdElevationLayers.setItemWidget(tree_item, 2, color_button)
            tree_item.setFlags(tree_item.flags() | Qt.ItemIsUserCheckable)
            tree_item.setCheckState(0, 2)


class PolygonIntersectionRepresentationDialog(QDialog):

    colors = [
        "darkseagreen",
        "darkgoldenrod",
        "darkviolet",
        "hotpink",
        "powderblue",
        "yellowgreen",
        "palevioletred",
        "seagreen",
        "darkturquoise",
        "beige",
        "darkkhaki",
        "red",
        "yellow",
        "magenta",
        "blue",
        "cyan",
        "chartreuse"
    ]

    def __init__(self,
                 plugin_name,
                 polygon_classification_set,
                 parent=None
                 ):

        super(PolygonIntersectionRepresentationDialog, self).__init__(parent)

        self.plugin_name = plugin_name
        self.polygon_classifications = list(polygon_classification_set)

        self.polygon_classifications_treeWidget = QTreeWidget()
        self.polygon_classifications_treeWidget.setColumnCount(2)
        self.polygon_classifications_treeWidget.headerItem().setText(0, "Name")
        self.polygon_classifications_treeWidget.headerItem().setText(1, "Color")
        self.polygon_classifications_treeWidget.setHorizontalScrollBarPolicy(Qt.ScrollBarAsNeeded)
        self.polygon_classifications_treeWidget.setDragEnabled(False)
        self.polygon_classifications_treeWidget.setDragDropMode(QAbstractItemView.NoDragDrop)
        self.polygon_classifications_treeWidget.setAlternatingRowColors(True)
        self.polygon_classifications_treeWidget.setTextElideMode(Qt.ElideLeft)

        self.update_classification_colors_treewidget()

        self.polygon_classifications_treeWidget.resizeColumnToContents(0)
        self.polygon_classifications_treeWidget.resizeColumnToContents(1)

        okButton = QPushButton("&OK")
        cancelButton = QPushButton("Cancel")

        buttonLayout = QHBoxLayout()
        buttonLayout.addStretch()
        buttonLayout.addWidget(okButton)
        buttonLayout.addWidget(cancelButton)

        layout = QGridLayout()

        layout.addWidget(
            self.polygon_classifications_treeWidget,
            0, 0, 1, 3)
        layout.addLayout(
            buttonLayout,
            1, 0, 1, 3)

        self.setLayout(layout)

        okButton.clicked.connect(self.accept)
        cancelButton.clicked.connect(self.reject)

        self.setWindowTitle("Polygon intersection colors")

    def update_classification_colors_treewidget(self):

        if len(PolygonIntersectionRepresentationDialog.colors) < len(self.polygon_classifications):
            dupl_factor = 1 + int(len(self.polygon_classifications) / len(PolygonIntersectionRepresentationDialog.colors))
            curr_colors = dupl_factor * PolygonIntersectionRepresentationDialog.colors
        else:
            curr_colors = PolygonIntersectionRepresentationDialog.colors

        self.polygon_classifications_treeWidget.clear()

        for classification_id, color in zip(self.polygon_classifications, curr_colors):
            tree_item = QTreeWidgetItem(self.polygon_classifications_treeWidget)
            tree_item.setText(0, str(classification_id))

            color_QgsColorButton = QgsColorButton()
            color_QgsColorButton.setColor(QColor(color))
            self.polygon_classifications_treeWidget.setItemWidget(tree_item, 1, color_QgsColorButton)


class PlotTopoProfileDialog(QDialog):

    def __init__(self,
                 plugin_name,
                 profile_length_set,
                 natural_elev_min_set,
                 natural_elev_max_set,
                 elevation_layer_names,
                 elevation_layer_colors,
                 parent=None
                 ):

        super(PlotTopoProfileDialog, self).__init__(parent)

        self.plugin_name = plugin_name
        self.elevation_layer_names = elevation_layer_names
        self.elevation_layer_colors = elevation_layer_colors

        # pre-process input data to account for multi.profiles

        profile_length = np.nanmax(profile_length_set)
        natural_elev_min = np.nanmin(natural_elev_min_set)
        natural_elev_max = np.nanmax(natural_elev_max_set)

        '''
        if np.isnan(profile_length) or profile_length == 0.0:
            error(
                self,
                self.plugin_name,
                f"Max profile length is {profile_length}.\nCheck profile trace."
            )
            return

        if np.isnan(natural_elev_min) or np.isnan(natural_elev_max):
            error(
                self,
                self.plugin_name,
                f"Max elevation in profile(s) is {natural_elev_max} and min is {natural_elev_min}.\nCheck profile trace."
            )
            return
        '''

        # pre-process elevation values

        # suggested plot elevation range

        z_padding = 0.5
        delta_z = natural_elev_max - natural_elev_min
        if delta_z < 0.0:
            warn(
                self.plugin_name,
                "Error: min elevation larger then max elevation"
            )
            return
        elif delta_z == 0.0:
            plot_z_min = floor(natural_elev_min) - 10
            plot_z_max = ceil(natural_elev_max) + 10
        else:
            plot_z_min = floor(natural_elev_min - delta_z * z_padding)
            plot_z_max = ceil(natural_elev_max + delta_z * z_padding)
        delta_plot_z = plot_z_max - plot_z_min

        # suggested exaggeration value

        w_to_h_rat = float(profile_length) / float(delta_plot_z)
        sugg_ve = 0.2*w_to_h_rat

        # Prepare the dialog

        layout = QVBoxLayout()

        # Axes

        qlytProfilePlot = QVBoxLayout()

        qgbxPlotSettings = QGroupBox("X axis")

        XAxisSettings_qgridlayout = QGridLayout()

        self.qcbxInvertXAxisProfile = QCheckBox(self.tr("Flip x-axis direction"))
        XAxisSettings_qgridlayout.addWidget(
            self.qcbxInvertXAxisProfile,
            0, 0, 1, 1)

        # Y variables

        qgbxYVariables = QGroupBox("Y axis")

        YAxis_qgridlayout = QGridLayout()

        self.qcbxPlotProfileHeight = QCheckBox(self.tr("Elevation"))
        self.qcbxPlotProfileHeight.setChecked(True)
        YAxis_qgridlayout.addWidget(
            self.qcbxPlotProfileHeight,
            0, 0, 1, 1)

        self.qcbxSetVerticalExaggeration = QCheckBox("Fixed vertical exaggeration")
        self.qcbxSetVerticalExaggeration.setChecked(True)
        YAxis_qgridlayout.addWidget(
            self.qcbxSetVerticalExaggeration,
            0, 1, 1, 1)

        self.qledtDemExagerationRatio = QLineEdit()
        self.qledtDemExagerationRatio.setText("%f" % sugg_ve)
        YAxis_qgridlayout.addWidget(
            self.qledtDemExagerationRatio,
            0, 2, 1, 1)

        YAxis_qgridlayout.addWidget(
            QLabel(self.tr("Plot z max value")),
            1, 1, 1, 1)

        self.qledtPlotMaxValue = QLineEdit()
        self.qledtPlotMaxValue.setText("%f" % plot_z_max)
        YAxis_qgridlayout.addWidget(
            self.qledtPlotMaxValue,
            1, 2, 1, 1)

        YAxis_qgridlayout.addWidget(
            QLabel(self.tr("Plot z min value")),
            2, 1, 1, 1)

        self.qledtPlotMinValue = QLineEdit()
        self.qledtPlotMinValue.setText("%f" % plot_z_min)
        YAxis_qgridlayout.addWidget(
            self.qledtPlotMinValue,
            2, 2, 1, 1)

        qgbxPlotSettings.setLayout(XAxisSettings_qgridlayout)

        qlytProfilePlot.addWidget(qgbxPlotSettings)

        self.qcbxPlotProfileSlope = QCheckBox(self.tr("Slope (degrees)"))
        YAxis_qgridlayout.addWidget(
            self.qcbxPlotProfileSlope,
            3, 0, 1, 1)

        self.qrbtPlotAbsoluteSlope = QRadioButton(self.tr("absolute"))
        self.qrbtPlotAbsoluteSlope.setChecked(True)
        YAxis_qgridlayout.addWidget(
            self.qrbtPlotAbsoluteSlope,
            3, 1, 1, 1)

        self.qrbtPlotDirectionalSlope = QRadioButton(self.tr("directional"))
        YAxis_qgridlayout.addWidget(
            self.qrbtPlotDirectionalSlope,
            3, 2, 1, 1)

        YAxis_qgridlayout.addWidget(
            QLabel("Note: to  calculate correctly the slope, the project must have\na planar CRS set or the DEM(s) must not be in lon-lat"),
            4, 1, 1, 3)

        qgbxYVariables.setLayout(YAxis_qgridlayout)

        qlytProfilePlot.addWidget(qgbxYVariables)

        # Style parameters

        qgbxStyleParameters = QGroupBox("Plot style")

        qlyStyleParameters = QGridLayout()

        self.qcbxPlotFilledHeight = QCheckBox(self.tr("Filled height"))
        qlyStyleParameters.addWidget(
            self.qcbxPlotFilledHeight,
            0, 0, 1, 1)

        self.qcbxPlotFilledSlope = QCheckBox(self.tr("Filled slope"))
        qlyStyleParameters.addWidget(
            self.qcbxPlotFilledSlope,
            0, 1, 1, 1)

        self.qpbtDefineTopoColors = QPushButton(self.tr("Elevation line visibility and colors"))
        self.qpbtDefineTopoColors.clicked.connect(self.define_profile_colors)
        qlyStyleParameters.addWidget(
            self.qpbtDefineTopoColors,
            1, 0, 1, 3)

        qgbxStyleParameters.setLayout(qlyStyleParameters)

        qlytProfilePlot.addWidget(qgbxStyleParameters)

        layout.addLayout(qlytProfilePlot)

        # ok/cancel section

        okButton = QPushButton("&OK")
        cancelButton = QPushButton("Cancel")

        buttonLayout = QHBoxLayout()
        buttonLayout.addStretch()
        buttonLayout.addWidget(okButton)
        buttonLayout.addWidget(cancelButton)

        layout.addLayout(buttonLayout)

        self.setLayout(layout)

        okButton.clicked.connect(self.accept)
        cancelButton.clicked.connect(self.reject)

        self.setWindowTitle("Topographic plot parameters")

    def define_profile_colors(self):

        def layer_styles(dialog):

            layer_visibilities = []
            layer_colors = []

            for layer_ndx in range(len(self.elevation_layer_names)):
                curr_item = dialog.qtwdElevationLayers.topLevelItem(layer_ndx)
                if curr_item.checkState(0) == 2:
                    layer_visibilities.append(True)
                else:
                    layer_visibilities.append(False)
                layer_colors.append(dialog.qtwdElevationLayers.itemWidget(curr_item, 2).color())

            return layer_visibilities, layer_colors

        if len(self.elevation_layer_names) == 0:
            warn(
                self.plugin_name,
                "No loaded elevation layer"
            )
            return

        dialog = ElevationLineStyleDialog(
            self.plugin_name,
            self.elevation_layer_names,
            self.elevation_layer_colors
        )

        if dialog.exec_():
            visible_elevation_layers, layer_colors = layer_styles(dialog)
        else:
            return

        if len(visible_elevation_layers) == 0:
            warn(
                self.plugin_name,
                "No visible layer"
            )
            return
        else:
            self.visible_elevation_layers = visible_elevation_layers
            self.elevation_layer_colors = layer_colors


class FigureExportDialog(QDialog):

    def __init__(self,
                 plugin_name,
                 parent=None
                 ):

        super(FigureExportDialog, self).__init__(parent)

        self.plugin_name = plugin_name

        layout = QVBoxLayout()

        # main parameters gropbox

        main_params_groupBox = QGroupBox("Main graphic parameters")

        main_params_layout = QGridLayout()

        main_params_layout.addWidget(
            QLabel(self.tr("Figure width (inches)")),
            0, 0, 1, 1)
        self.figure_width_inches_QLineEdit = QLineEdit("10")
        main_params_layout.addWidget(
            self.figure_width_inches_QLineEdit,
            0, 1, 1, 1)

        main_params_layout.addWidget(
            QLabel(self.tr("Resolution (dpi)")),
            0, 2, 1, 1)
        self.figure_resolution_dpi_QLineEdit = QLineEdit("200")
        main_params_layout.addWidget(
            self.figure_resolution_dpi_QLineEdit,
            0, 3, 1, 1)

        main_params_layout.addWidget(
            QLabel(self.tr("Font size (pts)")),
            0, 4, 1, 1)
        self.figure_fontsize_pts_QLineEdit = QLineEdit("12")
        main_params_layout.addWidget(
            self.figure_fontsize_pts_QLineEdit,
            0, 5, 1, 1)

        main_params_groupBox.setLayout(main_params_layout)

        layout.addWidget(main_params_groupBox)

        # additional parameters groupbox

        add_params_groupBox = QGroupBox(self.tr("Subplot configuration tools parameters"))

        add_params_layout = QGridLayout()

        add_params_layout.addWidget(
            QLabel("Top space"),
            0, 2, 1, 1)
        self.top_space_value_QDoubleSpinBox = QDoubleSpinBox()
        self.top_space_value_QDoubleSpinBox.setRange(0.0, 1.0)
        self.top_space_value_QDoubleSpinBox.setDecimals(2)
        self.top_space_value_QDoubleSpinBox.setSingleStep(0.01)
        self.top_space_value_QDoubleSpinBox.setValue(0.96)
        add_params_layout.addWidget(
            self.top_space_value_QDoubleSpinBox,
            0, 3, 1, 1)

        add_params_layout.addWidget(
            QLabel("Left space"),
            1, 0, 1, 1)
        self.left_space_value_QDoubleSpinBox = QDoubleSpinBox()
        self.left_space_value_QDoubleSpinBox.setRange(0.0, 1.0)
        self.left_space_value_QDoubleSpinBox.setDecimals(2)
        self.left_space_value_QDoubleSpinBox.setSingleStep(0.01)
        self.left_space_value_QDoubleSpinBox.setValue(0.1)
        add_params_layout.addWidget(
            self.left_space_value_QDoubleSpinBox,
            1, 1, 1, 1)

        add_params_layout.addWidget(
            QLabel("Right space"),
            1, 4, 1, 1)
        self.right_space_value_QDoubleSpinBox = QDoubleSpinBox()
        self.right_space_value_QDoubleSpinBox.setRange(0.0, 1.0)
        self.right_space_value_QDoubleSpinBox.setDecimals(2)
        self.right_space_value_QDoubleSpinBox.setSingleStep(0.01)
        self.right_space_value_QDoubleSpinBox.setValue(0.96)
        add_params_layout.addWidget(
            self.right_space_value_QDoubleSpinBox,
            1, 5, 1, 1)

        add_params_layout.addWidget(
            QLabel("Bottom space"),
            2, 2, 1, 1)
        self.bottom_space_value_QDoubleSpinBox = QDoubleSpinBox()
        self.bottom_space_value_QDoubleSpinBox.setRange(0.0, 1.0)
        self.bottom_space_value_QDoubleSpinBox.setDecimals(2)
        self.bottom_space_value_QDoubleSpinBox.setSingleStep(0.01)
        self.bottom_space_value_QDoubleSpinBox.setValue(0.06)
        add_params_layout.addWidget(
            self.bottom_space_value_QDoubleSpinBox,
            2, 3, 1, 1)

        add_params_layout.addWidget(
            QLabel("Blank width space between subplots"),
            3, 0, 1, 2)
        self.blank_width_space_value_QDoubleSpinBox = QDoubleSpinBox()
        self.blank_width_space_value_QDoubleSpinBox.setRange(0.0, 1.0)
        self.blank_width_space_value_QDoubleSpinBox.setDecimals(2)
        self.blank_width_space_value_QDoubleSpinBox.setSingleStep(0.01)
        self.blank_width_space_value_QDoubleSpinBox.setValue(0.1)
        add_params_layout.addWidget(
            self.blank_width_space_value_QDoubleSpinBox,
            3, 2, 1, 1)

        add_params_layout.addWidget(
            QLabel("Blank height space between subplots"),
            3, 3, 1, 2)
        self.blank_height_space_value_QDoubleSpinBox = QDoubleSpinBox()
        self.blank_height_space_value_QDoubleSpinBox.setRange(0.0, 1.0)
        self.blank_height_space_value_QDoubleSpinBox.setDecimals(2)
        self.blank_height_space_value_QDoubleSpinBox.setSingleStep(0.01)
        self.blank_height_space_value_QDoubleSpinBox.setValue(0.1)
        add_params_layout.addWidget(
            self.blank_height_space_value_QDoubleSpinBox,
            3, 5, 1, 1)

        add_params_layout.setRowMinimumHeight(3, 50)

        add_params_groupBox.setLayout(add_params_layout)

        layout.addWidget(add_params_groupBox)

        # graphic parameters import and export

        graphic_params_io_groupBox = QGroupBox("Graphic parameters save/load")

        graphic_params_io_layout = QHBoxLayout()

        self.graphic_params_save_QPushButton = QPushButton("Save")
        self.graphic_params_save_QPushButton.clicked.connect(self.output_graphic_params_save)
        graphic_params_io_layout.addWidget(
            self.graphic_params_save_QPushButton)

        self.graphic_params_load_QPushButton = QPushButton("Load")
        self.graphic_params_load_QPushButton.clicked.connect(self.output_graphic_params_load)
        graphic_params_io_layout.addWidget(
            self.graphic_params_load_QPushButton)

        graphic_params_io_groupBox.setLayout(graphic_params_io_layout)

        layout.addWidget(
            graphic_params_io_groupBox)

        # output file parameters

        output_file_groupBox = QGroupBox(self.tr("Output file - available formats: tif, pdf, svg"))

        output_file_layout = QGridLayout()

        self.figure_outpath_QLineEdit = QLineEdit()
        output_file_layout.addWidget(
            self.figure_outpath_QLineEdit,
            3, 0, 1, 1)

        self.figure_outpath_QPushButton = QPushButton(self.tr("Choose"))
        self.figure_outpath_QPushButton.clicked.connect(self.define_figure_outpath)
        output_file_layout.addWidget(
            self.figure_outpath_QPushButton,
            3, 1, 1, 1)

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

        okButton.clicked.connect(self.accept)
        cancelButton.clicked.connect(self.reject)

        self.setWindowTitle("Export figure")

    def output_graphic_params_save(self):

        output_file_path = new_file_path(
            self,
            "Define output configuration file",
            "*.txt",
            "txt"
        )

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

        ok(
            self.plugin_name,
            "Graphic parameters saved"
        )

    def output_graphic_params_load(self):

        input_file_path = old_file_path(
            self,
            "Choose input configuration file",
            "*.txt",
            "txt"
        )

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

            error(
                self.plugin_name,
                "Error in configuration file"
            )
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

        outfile_path = new_file_path(
            self,
            "Create",
            "",
            "Images (*.svg *.pdf *.tif)"
        )

        self.figure_outpath_QLineEdit.setText(outfile_path)


class TopographicProfileExportDialog(QDialog):

    def __init__(self,
                 plugin_name,
                 selected_dem_params,
                 parent=None
                 ):

        super(TopographicProfileExportDialog, self).__init__(parent)

        self.plugin_name = plugin_name

        layout = QVBoxLayout()

        # Profile source

        source_groupBox = QGroupBox(self.tr("Profile sources"))

        source_layout = QGridLayout()

        self.src_allselecteddems_QRadioButton = QRadioButton(self.tr("All selected DEMs"))
        source_layout.addWidget(self.src_allselecteddems_QRadioButton, 1, 0, 1, 2)
        self.src_allselecteddems_QRadioButton.setChecked(True)

        self.src_singledem_QRadioButton = QRadioButton(self.tr("Single DEM"))
        source_layout.addWidget(
            self.src_singledem_QRadioButton,
            2, 0, 1, 1)

        self.src_singledemlist_QComboBox = QComboBox()
        selected_dem_layers = [dem_param.layer for dem_param in selected_dem_params]
        for qgsRasterLayer in selected_dem_layers:
            self.src_singledemlist_QComboBox.addItem(qgsRasterLayer.name())
        source_layout.addWidget(
            self.src_singledemlist_QComboBox,
            2, 1, 1, 1)

        self.src_singlegpx_QRadioButton = QRadioButton(self.tr("GPX file"))
        source_layout.addWidget(
            self.src_singlegpx_QRadioButton,
            3, 0, 1, 1)

        source_groupBox.setLayout(source_layout)

        layout.addWidget(source_groupBox)

        # Output type

        output_type_groupBox = QGroupBox(self.tr("Output format"))

        output_type_layout = QGridLayout()

        self.outtype_shapefile_point_QRadioButton = QRadioButton(self.tr("shapefile - point"))
        output_type_layout.addWidget(self.outtype_shapefile_point_QRadioButton, 0, 0, 1, 1)
        self.outtype_shapefile_point_QRadioButton.setChecked(True)

        self.outtype_shapefile_line_QRadioButton = QRadioButton(self.tr("shapefile - line"))
        output_type_layout.addWidget(
            self.outtype_shapefile_line_QRadioButton,
            1, 0, 1, 1)

        self.outtype_csv_QRadioButton = QRadioButton(self.tr("csv"))
        output_type_layout.addWidget(
            self.outtype_csv_QRadioButton,
            2, 0, 1, 1)

        output_type_groupBox.setLayout(output_type_layout)

        layout.addWidget(output_type_groupBox)

        # Output name/path

        output_path_groupBox = QGroupBox(self.tr("Output file"))

        output_path_layout = QGridLayout()

        self.outpath_QLineEdit = QLineEdit()
        output_path_layout.addWidget(
            self.outpath_QLineEdit,
            0, 0, 1, 1)

        self.outpath_QPushButton = QPushButton("....")
        self.outpath_QPushButton.clicked.connect(self.define_outpath)
        output_path_layout.addWidget(
            self.outpath_QPushButton,
            0, 1, 1, 1)

        self.load_output_checkBox = QCheckBox("load output shapefile in project")
        self.load_output_checkBox.setChecked(True)
        output_path_layout.addWidget(
            self.load_output_checkBox,
            1, 0, 1, 2)

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

        okButton.clicked.connect(self.accept)
        cancelButton.clicked.connect(self.reject)

        self.setWindowTitle("Export topographic profile")

    def define_outpath(self):

        if self.outtype_shapefile_line_QRadioButton.isChecked() or self.outtype_shapefile_point_QRadioButton.isChecked():
            outfile_path = new_file_path(
                self,
                "Save file",
                "",
                "Shapefiles (*.shp)"
            )
        elif self.outtype_csv_QRadioButton.isChecked():
            outfile_path = new_file_path(
                self,
                "Save file",
                "",
                "Csv (*.csv)"
            )
        else:
            error(
                self.plugin_name,
                self.tr("Output type definition error")
            )
            return

        self.outpath_QLineEdit.setText(outfile_path)


class PointDataExportDialog(QDialog):

    def __init__(self,
                 plugin_name,
                 parent=None
                 ):

        super(PointDataExportDialog, self).__init__(parent)

        self.plugin_name = plugin_name

        layout = QVBoxLayout()

        # Output type

        output_type_groupBox = QGroupBox(self.tr("Output format"))

        output_type_layout = QGridLayout()

        self.outtype_shapefile_point_QRadioButton = QRadioButton(self.tr("shapefile - point"))
        output_type_layout.addWidget(
            self.outtype_shapefile_point_QRadioButton,
            0, 0, 1, 1)
        self.outtype_shapefile_point_QRadioButton.setChecked(True)

        self.outtype_csv_QRadioButton = QRadioButton(self.tr("csv"))
        output_type_layout.addWidget(
            self.outtype_csv_QRadioButton,
            1, 0, 1, 1)

        output_type_groupBox.setLayout(output_type_layout)

        layout.addWidget(output_type_groupBox)

        # Output name/path

        output_path_groupBox = QGroupBox(self.tr("Output path"))

        output_path_layout = QGridLayout()

        self.outpath_QLineEdit = QLineEdit()
        output_path_layout.addWidget(
            self.outpath_QLineEdit,
            0, 0, 1, 1)

        self.outpath_QPushButton = QPushButton(self.tr("Choose"))
        self.outpath_QPushButton.clicked.connect(self.define_outpath)
        output_path_layout.addWidget(
            self.outpath_QPushButton,
            0, 1, 1, 1)

        self.load_output_checkBox = QCheckBox("load output shapefile in project")
        self.load_output_checkBox.setChecked(True)
        output_path_layout.addWidget(
            self.load_output_checkBox,
            1, 0, 1, 2)

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

        okButton.clicked.connect(self.accept)
        cancelButton.clicked.connect(self.reject)

        self.setWindowTitle("Export geological attitudes")

    def define_outpath(self):

        if self.outtype_shapefile_point_QRadioButton.isChecked():
            outfile_path = new_file_path(
                self,
                "Path",
                "*.shp",
                "Shapefile"
            )
        elif self.outtype_csv_QRadioButton.isChecked():
            outfile_path = new_file_path(
                self,
                "Path",
                "*.csv",
                "Csv"
            )
        else:
            error(
                self.plugin_name,
                self.tr("Output type definition error")
            )
            return

        self.outpath_QLineEdit.setText(outfile_path)


class LineDataExportDialog(QDialog):

    def __init__(self,
                 plugin_name,
                 parent=None
                 ):

        super(LineDataExportDialog, self).__init__(parent)

        self.plugin_name = plugin_name

        layout = QVBoxLayout()

        # Output type

        output_type_groupBox = QGroupBox(self.tr("Output format"))

        output_type_layout = QGridLayout()

        self.outtype_shapefile_line_QRadioButton = QRadioButton(self.tr("shapefile - line"))
        output_type_layout.addWidget(self.outtype_shapefile_line_QRadioButton, 0, 0, 1, 1)
        self.outtype_shapefile_line_QRadioButton.setChecked(True)

        self.outtype_csv_QRadioButton = QRadioButton(self.tr("csv"))
        output_type_layout.addWidget(
            self.outtype_csv_QRadioButton,
            0, 1, 1, 1)

        output_type_groupBox.setLayout(output_type_layout)

        layout.addWidget(output_type_groupBox)

        # Output name/path

        output_path_groupBox = QGroupBox(self.tr("Output file"))

        output_path_layout = QGridLayout()

        self.outpath_QLineEdit = QLineEdit()
        output_path_layout.addWidget(
            self.outpath_QLineEdit,
            0, 0, 1, 1)

        self.outpath_QPushButton = QPushButton("....")
        self.outpath_QPushButton.clicked.connect(self.define_outpath)
        output_path_layout.addWidget(
            self.outpath_QPushButton,
            0, 1, 1, 1)

        self.load_output_checkBox = QCheckBox("load output in project")
        self.load_output_checkBox.setChecked(True)
        output_path_layout.addWidget(
            self.load_output_checkBox,
            1, 0, 1, 2)

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

        okButton.clicked.connect(self.accept)
        cancelButton.clicked.connect(self.reject)

        self.setWindowTitle("Export")

    def define_outpath(self):

        if self.outtype_shapefile_line_QRadioButton.isChecked():
            outfile_path = new_file_path(
                self,
                "Save file",
                "",
                "Shapefiles (*.shp)"
            )
        elif self.outtype_csv_QRadioButton.isChecked():
            outfile_path = new_file_path(
                self,
                "Save file",
                "",
                "Csv (*.csv)"
            )
        else:
            error(
                self.plugin_name,
                self.tr("Output type definition error")
            )
            return

        self.outpath_QLineEdit.setText(outfile_path)


class StatisticsDialog(QDialog):

    def __init__(
            self,
            plugin_name,
            geoprofile_set,
            parent=None
    ):

        super(StatisticsDialog, self).__init__(parent)

        self.plugin_name = plugin_name

        layout = QVBoxLayout()

        self.text_widget = QTextEdit()
        self.text_widget.setReadOnly(True)

        num_profiles = geoprofile_set.geoprofiles_num
        stat_report = f"\nGeneral statistics for {num_profiles} profiles\n"

        for ndx in range(num_profiles):

            profile_elevations = geoprofile_set.geoprofile(ndx)._named_lines

            profiles_stats = list(
                zip(
                    profile_elevations.surface_names,
                    list(
                        zip(
                            profile_elevations.statistics_elev,
                            profile_elevations.statistics_dirslopes,
                            profile_elevations.statistics_slopes
                        )
                    )
                )
            )

            stat_report += f"\nStatistics for profile # {ndx+1}"
            stat_report += f"\n\tLength: {profile_elevations.x_length}"
            stat_report += "\n\tTopographic elevations"
            stat_report += f"\n\t - min: {profile_elevations.natural_elev_range[0]}"
            stat_report += f"\n\t - max: {profile_elevations.natural_elev_range[1]}"
            stat_report += self.report_stats(profiles_stats)

        for ndx in range(num_profiles):

            topo_profiles = geoprofile_set.geoprofile(ndx)._named_lines
            resampled_line_xs = topo_profiles.x_array
            resampled_line_ys = topo_profiles.y_array

            if resampled_line_xs is not None:

                stat_report += f"\nSampling points ({len(resampled_line_xs)}) for profile # {ndx + 1}"

                for ln_ndx, (x, y) in enumerate(zip(resampled_line_xs, resampled_line_ys)):
                   stat_report += f"\n{ln_ndx+1}, {x}, {y}"

        self.text_widget.insertPlainText(stat_report)

        layout.addWidget(self.text_widget)

        self.setLayout(layout)

        self.setWindowTitle("Statistics")

    def report_stats(self, profiles_stats):

        def type_report(values):

            type_report = f"min: {values['min']}\n"
            type_report += f"max: {values['max']}\n"
            type_report += f"mean: {values['mean']}\n"
            type_report += f"variance: {values['var']}\n"
            type_report += f"standard deviation: {values['std']}\n\n"

            return type_report

        report = 'Dataset statistics\n'
        types = [
            'elevations',
            'directional slopes',
            'absolute slopes'
        ]

        for name, stats in profiles_stats:
            report += f"\ndataset name\n{name}\n\n"
            for tp, stat_val in zip(types, stats):
                report += f"{tp}\n\n"
                report += type_report(stat_val)

        return report


class DigitizeLineDialog(QDialog):

    def __init__(
            self,
            plugin_name,
            canvas,
            parent=None
    ):

        super(DigitizeLineDialog, self).__init__(parent)

        self.plugin_name = plugin_name
        self.canvas = canvas
        self.previous_maptool = self.canvas.mapTool()  # Save the standard map tool for restoring it at the end

        layout = QVBoxLayout()

        self.qpbtDigitizeLine = QPushButton(self.tr("Digitize trace"))
        self.qpbtDigitizeLine.setToolTip(
            "Digitize the trace on the map.\n"
            "Left click: add point\n"
            "Right click: end adding point"
        )
        self.qpbtDigitizeLine.clicked.connect(self.digitize_line)

        layout.addWidget(self.qpbtDigitizeLine)

        self.qpbtClearLine = QPushButton(self.tr("Clear"))
        self.qpbtClearLine.clicked.connect(self.clear_rubberband)
        layout.addWidget(self.qpbtClearLine)

        self.qpbtClearLine = QPushButton(self.tr("Save"))
        self.qpbtClearLine.clicked.connect(self.save_rubberband)
        layout.addWidget(self.qpbtClearLine)

        self.setLayout(layout)

        self.setWindowTitle("Digitize line")

    '''
    def connect_digitize_maptool(self):

        self.digitize_maptool.moved.connect(self.canvas_refresh_profile_line)
        self.digitize_maptool.leftClicked.connect(self.profile_add_point)
        self.digitize_maptool.rightClicked.connect(self.canvas_end_profile_line)
    '''

    def digitize_line(self):

        self.clear_rubberband()

        info(
            self.plugin_name,
            "Now you can digitize a line on the map.\nLeft click: add point\nRight click: end adding point"
        )

        self.rubberband = QgsRubberBand(self.canvas)
        self.rubberband.setWidth(2)
        self.rubberband.setColor(QColor(Qt.red))

        self.digitize_maptool = MapDigitizeTool(self.canvas)
        self.canvas.setMapTool(self.digitize_maptool)

        self.digitize_maptool.moved.connect(self.canvas_refresh_profile_line)
        self.digitize_maptool.leftClicked.connect(self.profile_add_point)
        self.digitize_maptool.rightClicked.connect(self.canvas_end_profile_line)

        #print(f"DEBUG: exiting digitize_line")

    def canvas_refresh_profile_line(self, position):

        """
        if len(self.profile_canvas_points) == 0:
            return
        """

        x, y = xy_from_canvas(self.canvas, position)
        #print(f"DEBUG: canvas_refresh_profile_line -> {x}, {y}")

        self.refresh_rubberband(self.profile_canvas_points + [[x, y]])

    def profile_add_point(self, position):

        x, y = xy_from_canvas(self.canvas, position)
        #print(f"DEBUG: profile_add_point -> {x} {y}")

        self.profile_canvas_points.append([x, y])

    def canvas_end_profile_line(self):

        self.refresh_rubberband(self.profile_canvas_points)

        self.digitized_profile_line2dt = None
        if len(self.profile_canvas_points) > 1:
            raw_line = Line2D(
                [Point2D(x, y) for x, y in self.profile_canvas_points]).remove_coincident_points()
            if raw_line.num_pts() > 1:
                self.digitized_profile_line2dt = raw_line

        self.profile_canvas_points = []

        self.restore_previous_map_tool()

    def restore_previous_map_tool(self):

        self.canvas.unsetMapTool(self.digitize_maptool)
        self.canvas.setMapTool(self.previous_maptool)

    def refresh_rubberband(self,
                           xy_list
                           ):

        self.rubberband.reset(QgsWkbTypes.LineGeometry)
        for x, y in xy_list:
            self.rubberband.addPoint(QgsPointXY(x, y))

    def clear_rubberband(self):

        self.profile_canvas_points = []
        self.digitized_profile_line2dt = None

        try:

            self.rubberband.reset()

        except:

            pass

    def save_rubberband(self):

        def output_profile_line(
                output_format,
                output_filepath,
                pts2dt,
                proj_sr
        ):

            points = [[n, pt2dt.x, pt2dt.y] for n, pt2dt in enumerate(pts2dt)]
            if output_format == "csv":
                success, msg = write_generic_csv(output_filepath,
                                                 ['id', 'x', 'y'],
                                                 points)
                if not success:
                    warn(
                        self.plugin_name,
                        msg
                    )
            elif output_format == "shapefile - line":
                success, msg = write_rubberband_profile_lnshp(
                    output_filepath,
                    ['id'],
                    points,
                    proj_sr)
                if not success:
                    warn(
                        self.plugin_name,
                        msg
                    )
            else:
                error(
                    self.plugin_name,
                    "Debug: error in export format"
                )
                return

            if success:
                ok(
                    self.plugin_name,
                    "Line saved"
                )

        def get_format_type():

            if dialog.outtype_shapefile_line_QRadioButton.isChecked():
                return "shapefile - line"
            elif dialog.outtype_csv_QRadioButton.isChecked():
                return "csv"
            else:
                return ""

        if self.digitized_profile_line2dt is None:

            warn(
                self.plugin_name,
                "No available line to save [1]"
            )
            return

        elif self.digitized_profile_line2dt.num_pts() < 2:

            warn(
                self.plugin_name,
                "No available line to save [2]"
            )
            return

        dialog = LineDataExportDialog(self.plugin_name)
        if dialog.exec_():
            output_format = get_format_type()
            if output_format == "":
                warn(
                    self.plugin_name,
                    "Error in output format"
                )
                return
            output_filepath = dialog.outpath_QLineEdit.text()
            if len(output_filepath) == 0:
                warn(
                    self.plugin_name,
                    "Error in output path"
                )
                return
            add_to_project = dialog.load_output_checkBox.isChecked()
        else:
            warn(
                self.plugin_name,
                "No export defined"
            )
            return

        # get project CRS information
        project_crs_osr = proj4str()

        output_profile_line(
            output_format,
            output_filepath,
            self.digitized_profile_line2dt.pts,
            project_crs_osr)

        # add theme to QGis project
        if output_format == "shapefile - line" and add_to_project:
            try:

                digitized_line_layer = QgsVectorLayer(output_filepath,
                                                      QFileInfo(output_filepath).baseName(),
                                                      "ogr")
                QgsProject.instance().addMapLayer(digitized_line_layer)

            except:

                error(
                    self.plugin_name,
                    "Unable to load layer in project"
                )
                return


class GpxInputDialog(QDialog):

    def __init__(self,
                 plugin_name,
                 settings,
                 settings_gpxdir_key,
                 parent=None
                 ):

        super(GpxInputDialog, self).__init__(parent)

        self.plugin_name = plugin_name
        self.settings = settings
        self.settings_gpxdir_key = settings_gpxdir_key

        layout = QGridLayout()

        layout.addWidget(
            QLabel(self.tr("Choose input file:")),
            0, 0, 1, 1)

        self.input_gpx_file_path = QLineEdit()
        self.input_gpx_file_path.setPlaceholderText("my_track.gpx")
        layout.addWidget(
            self.input_gpx_file_path,
            0, 1, 1, 1)

        self.define_input_gpx_file = QPushButton("...")
        self.define_input_gpx_file.clicked.connect(self.select_input_gpx_file)
        layout.addWidget(
            self.define_input_gpx_file,
            0, 2, 1, 1)

        self.invert_track_direction = QCheckBox(
            "Invert direction"
        )

        layout.addWidget(
            self.invert_track_direction,
            1, 0, 1, 2)

        okButton = QPushButton("&OK")
        cancelButton = QPushButton("Cancel")
        okButton.clicked.connect(self.accept)
        cancelButton.clicked.connect(self.reject)

        layout.addWidget(okButton,
                         2, 1, 1, 1)
        layout.addWidget(cancelButton,
                         2, 2, 1, 1)

        self.setLayout(layout)

        self.setWindowTitle("Input GPX")

    def select_input_gpx_file(self):

        gpx_last_used_dir = self.settings.value(self.settings_gpxdir_key,
                                                "")
        file_name, __ = QFileDialog.getOpenFileName(
            self,
            self.tr("Open GPX file"),
            gpx_last_used_dir,
            "GPX (*.gpx *.GPX)"
        )
        if not file_name:
            return
        else:
            update_directory_key(self.settings,
                                 self.settings_gpxdir_key,
                                 file_name)
            self.input_gpx_file_path.setText(file_name)


class ParallelProfilesParametersDialog(QDialog):

    def __init__(
            self,
            plugin_name,
            parent=None
    ):

        super(ParallelProfilesParametersDialog, self).__init__(parent)

        self.plugin_name = plugin_name

        layout = QGridLayout()

        layout.addWidget(QLabel("Profiles spacing"), 0, 0, 1, 1)

        self.spacing_wdgt = QDoubleSpinBox()
        self.spacing_wdgt.setValue(100.0)
        self.spacing_wdgt.setMinimum(0.0)
        self.spacing_wdgt.setMaximum(10000.0)
        self.spacing_wdgt.setSingleStep(50.0)

        layout.addWidget(self.spacing_wdgt, 0, 1, 1, 1)

        layout.addWidget(QLabel("Num. left profiles"), 1, 0, 1, 1)

        self.num_left_profiles_wdgt = QSpinBox()
        self.num_left_profiles_wdgt.setValue(5)
        self.num_left_profiles_wdgt.setMinimum(0)

        layout.addWidget(self.num_left_profiles_wdgt, 1, 1, 1, 1)

        layout.addWidget(QLabel("Num. right profiles"), 2, 0, 1, 1)

        self.num_right_profiles_wdgt = QSpinBox()
        self.num_right_profiles_wdgt.setValue(5)
        self.num_right_profiles_wdgt.setMinimum(0)

        layout.addWidget(self.num_right_profiles_wdgt, 2, 1, 1, 1)

        okButton = QPushButton("&OK")
        cancelButton = QPushButton("Cancel")
        okButton.clicked.connect(self.accept)
        cancelButton.clicked.connect(self.reject)

        verticalSpacer = QSpacerItem(20, 30, QSizePolicy.Minimum, QSizePolicy.Expanding)
        layout.addItem(verticalSpacer)

        layout.addWidget(okButton,
                         4, 0, 1, 1)
        layout.addWidget(cancelButton,
                         4, 1, 1, 1)

        self.setLayout(layout)

        self.setWindowTitle("Define parallel profiles parameters")
        self.setMinimumSize(300, 150)







