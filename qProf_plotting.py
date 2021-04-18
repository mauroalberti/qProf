
import matplotlib as mpl

from .qygsf.profiles.geoprofiles import *
from .qygsf.utils.mpl_utils.mpl_widget import *
from .qygsf.utils.qt_utils.colors import *

from .config.settings import *


def plot_structural_attitude(
        plot_addit_params,
        axes,
        section_length,
        vertical_exaggeration,
        structural_attitude_list,
        color
):

    # TODO:  manage case for possible nan z values
    projected_z = [structural_attitude.pt_3d.z for structural_attitude in structural_attitude_list if
                   0.0 <= structural_attitude.sign_hor_dist <= section_length]

    # TODO:  manage case for possible nan z values
    projected_s = [structural_attitude.sign_hor_dist for structural_attitude in structural_attitude_list if
                   0.0 <= structural_attitude.sign_hor_dist <= section_length]

    projected_ids = [structural_attitude.id for structural_attitude in structural_attitude_list if
                     0.0 <= structural_attitude.sign_hor_dist <= section_length]

    axes.plot(projected_s, projected_z, 'o', color=color)

    # plot segments representing structural data
    for structural_attitude in structural_attitude_list:
        if 0.0 <= structural_attitude.sign_hor_dist <= section_length:
            structural_segment_s, structural_segment_z = define_plot_structural_segment(structural_attitude,
                                                                                        section_length,
                                                                                        vertical_exaggeration)

            axes.plot(structural_segment_s, structural_segment_z, '-', color=color)

    if plot_addit_params["add_trendplunge_label"] or plot_addit_params["add_ptid_label"]:

        src_dip_dirs = [structural_attitude.src_geol_plane.dd for structural_attitude in
                        structural_attitude_list if 0.0 <= structural_attitude.sign_hor_dist <= section_length]
        src_dip_angs = [structural_attitude.src_geol_plane.da for structural_attitude in
                        structural_attitude_list if 0.0 <= structural_attitude.sign_hor_dist <= section_length]

        for rec_id, src_dip_dir, src_dip_ang, s, z in zip(projected_ids, src_dip_dirs, src_dip_angs, projected_s,
                                                          projected_z):

            if plot_addit_params["add_trendplunge_label"] and plot_addit_params["add_ptid_label"]:
                label = "%s-%03d/%02d" % (rec_id, src_dip_dir, src_dip_ang)
            elif plot_addit_params["add_ptid_label"]:
                label = "%s" % rec_id
            elif plot_addit_params["add_trendplunge_label"]:
                label = "%03d/%02d" % (src_dip_dir, src_dip_ang)

            axes.annotate(label, (s + 15, z + 15))


def plot_projected_line_set(
        axes,
        curve_set,
        labels
):

    colors = colors_addit * (int(len(curve_set) / len(colors_addit)) + 1)
    for multiline_2d, label, color in zip(curve_set, labels, colors):
        for line_2d in multiline_2d.lines:
            plot_line(axes, line_2d.x_array, line_2d.y_array, color, name=label)


def plot_profile_lines_intersection_points(
        axes,
        profile_lines_intersection_points
):

    for s, pt3d, intersection_id, color in profile_lines_intersection_points:
        axes.plot(s, pt3d.z, 'o', color=color)
        if str(intersection_id).upper() != "NULL" or str(intersection_id) != '':
            axes.annotate(str(intersection_id), (s + 25, pt3d.z + 25))


def plot_profile_polygon_intersection_line(
        plot_addit_params,
        axes,
        intersection_line_value
):

    classification, line3d, s_list = intersection_line_value
    z_list = [pt3d.z for pt3d in line3d.pts]

    if plot_addit_params["polygon_class_colors"] is None:
        color = "red"
    else:
        color = plot_addit_params["polygon_class_colors"][str(classification)]

    plot_line(axes, s_list, z_list, color, linewidth=3.0, name=classification)


def create_axes(
    profile_window,
    plot_x_range,
    plot_y_range,
    grid_spec,
    ndx_subplot
) -> 'mpl.axes._axes.Axes.AxesSubplot':

    x_min, x_max = plot_x_range
    y_min, y_max = plot_y_range
    axes = profile_window.canvas.fig.add_subplot(grid_spec[ndx_subplot])
    #print(x_min, x_max)
    axes.set_xlim(x_min, x_max)
    axes.set_ylim(y_min, y_max)

    axes.grid(True, color='lightgrey', linestyle='-.', linewidth=0.25, alpha=50)
    #print(f"DOCUMENTATION: axes: {type(axes)}")
    return axes


def plot_topo_profile_lines(
    geoprofile,
    plot_params,
    profile_window,
    grid_spec,
    ndx_subplot,
    topo_type,
    plot_s_range,
    plot_z_range,
    filled_choice
):

    topo_profiles = [line3d for _, line3d in geoprofile.named_topoprofiles]
    topoline_colors = plot_params['elev_lyr_colors']
    topoline_visibilities = plot_params['visible_elev_lyrs']

    axes = create_axes(
        profile_window,
        plot_s_range,
        plot_z_range,
        grid_spec,
        ndx_subplot
    )

    if plot_params['invert_xaxis']:
        axes.invert_xaxis()

    s_values_list = [line3d.incremental_length_2d() for line3d in topo_profiles]

    if topo_type == 'elevation':
        z_values_list = [line3d.z_array() for line3d in topo_profiles]
        plot_z_filled = plot_z_range[0]
    else:
        if plot_params['plot_slope_absolute']:
            z_values_list = [line3d.absolute_slopes() for line3d in topo_profiles]
        else:
            z_values_list = [line3d.dir_slopes() for line3d in topo_profiles]
        plot_z_filled = 0.0

    for s, z, topoline_color, topoline_visibility in zip(
            s_values_list,
            z_values_list,
            topoline_colors,
            topoline_visibilities
    ):

        if topoline_visibility:

            if filled_choice:
                plot_filled_line(
                    axes,
                    s,
                    z,
                    plot_z_filled,
                    qcolor2rgbmpl(topoline_color))

            plot_line(
                axes,
                s,
                z,
                qcolor2rgbmpl(topoline_color))

    return axes


def plot_geoprofile(
    geoprofile: GeoProfile_,
    plot_params,
    plot_addit_params: Dict,
    slope_padding=0.2
):

    # extract/define plot parameters

    #plot_params = input_geoprofiles_set.plot_params

    set_vertical_exaggeration = plot_params["set_vertical_exaggeration"]
    vertical_exaggeration = plot_params['vertical_exaggeration']

    plot_height_choice = plot_params['plot_height_choice']
    plot_slope_choice = plot_params['plot_slope_choice']

    if plot_height_choice:
        # defines plot min and max values
        plot_z_min = plot_params['plot_min_elevation_user']
        plot_z_max = plot_params['plot_max_elevation_user']

    # populate the plot

    profile_window = MplWidget('Profile')

    num_subplots = plot_height_choice + plot_slope_choice
    grid_spec = mpl.gridspec.GridSpec(num_subplots, 1)

    ndx_subplot = -1

    plot_s_min, plot_s_max = 0, geoprofile.s_max()
    #print(f"A: plot_s_min, plot_s_max")

    # if slopes to be calculated and plotted
    if plot_slope_choice:
        plot_slope_max = 90.0
        # defines slope value lists and the min and max values
        if plot_params['plot_slope_absolute']:
            #slopes = geoprofile.topo_profiles.absolute_slopes
            plot_slope_min = 0.0
        else:
            #slopes = geoprofile.topo_profiles.profile_dirslopes
            plot_slope_min = -90.0

        """
        profiles_slope_min = np.nanmin(np.array(list(map(np.nanmin, slopes))))
        profiles_slope_max = np.nanmax(np.array(list(map(np.nanmax, slopes))))

        delta_slope = profiles_slope_max - profiles_slope_min
        plot_slope_min = profiles_slope_min - delta_slope * slope_padding
        plot_slope_max = profiles_slope_max + delta_slope * slope_padding
        """

    # plot topographic profile elevations

    if plot_height_choice:
        ndx_subplot += 1
        axes_elevation = plot_topo_profile_lines(
            geoprofile,
            plot_params,
            profile_window,
            grid_spec,
            ndx_subplot,
            'elevation',
            (plot_s_min, plot_s_max),
            (plot_z_min, plot_z_max),
            plot_params['filled_height'])
        if set_vertical_exaggeration:
            axes_elevation.set_aspect(vertical_exaggeration)
        axes_elevation.set_anchor('W')  # align left

    # plot topographic profile slopes

    if plot_slope_choice:
        ndx_subplot += 1
        axes_slopes = plot_topo_profile_lines(
            geoprofile,
            plot_params,
            profile_window,
            grid_spec,
            ndx_subplot,
            'slope',
            (plot_s_min, plot_s_max),
            (plot_slope_min, plot_slope_max),
            plot_params['filled_slope'])
        axes_slopes.set_anchor('W')  # align left

    # plot geological outcrop intersections

    if len(geoprofile.polygons_intersections) > 0:
        for line_intersection_value in geoprofile.polygons_intersections:
            plot_profile_polygon_intersection_line(plot_addit_params,
                                                   axes_elevation,
                                                   line_intersection_value)

    # plot geological attitudes intersections

    if len(geoprofile.profile_attitudes) > 0:
        for plane_attitude_set, color in zip(geoprofile.profile_attitudes, plot_addit_params["plane_attitudes_colors"]):
            plot_structural_attitude(plot_addit_params,
                                     axes_elevation,
                                     plot_s_max,
                                     vertical_exaggeration,
                                     plane_attitude_set,
                                     color)

    # plot geological traces projections

    if len(geoprofile.projected_lines) > 0:
        for curve_set, labels in zip(geoprofile.projected_lines, geoprofile.projected_lines_ids):
            plot_projected_line_set(axes_elevation,
                                    curve_set,
                                    labels)

    # plot line-profile intersections

    if len(geoprofile.line_intersections) > 0:
        plot_profile_lines_intersection_points(axes_elevation,
                                               geoprofile.line_intersections)

    profile_window.canvas.fig.tight_layout()
    profile_window.canvas.draw()

    return profile_window


def plot_geoprofiles(
        input_geoprofiles_set: GeoProfilesSet_,
        plot_addit_params,
        slope_padding=0.2
):

    #print(f"DOC: geoprofiles in plot_geoprofiles: {type(input_geoprofiles_set)}")

    # extract/define plot parameters

    plot_params = input_geoprofiles_set.plot_params

    set_vertical_exaggeration = plot_params["set_vertical_exaggeration"]
    vertical_exaggeration = plot_params['vertical_exaggeration']

    plot_height_choice = plot_params['plot_height_choice']
    plot_slope_choice = plot_params['plot_slope_choice']

    if plot_height_choice:
        # defines plot min and max values
        plot_z_min = plot_params['plot_min_elevation_user']
        plot_z_max = plot_params['plot_max_elevation_user']

    # populate the plot

    profile_window = MplWidget('Profile')

    #print(f"DEBUG: plot_height_choice: {plot_height_choice}")
    #print(f"DEBUG: plot_slope_choice: {plot_slope_choice}")
    #print(f"DEBUG: geoprofiles.geoprofiles_num: {input_geoprofiles_set.geoprofiles_num}")

    num_subplots = (plot_height_choice + plot_slope_choice) * input_geoprofiles_set.geoprofiles_num
    grid_spec = mpl.gridspec.GridSpec(num_subplots, 1)

    ndx_subplot = -1
    for ndx in range(input_geoprofiles_set.geoprofiles_num):

        geoprofile = input_geoprofiles_set.geoprofile(ndx)
        plot_s_min, plot_s_max = 0, geoprofile.named_topoprofiles.profile_length()
        #print(f"B: plot_s_min, plot_s_max")

        # if slopes to be calculated and plotted
        if plot_slope_choice:
            plot_slope_max = 90.0
            # defines slope value lists and the min and max values
            if plot_params['plot_slope_absolute']:
                #slopes = geoprofile.topo_profiles.absolute_slopes
                plot_slope_min = 0.0
            else:
                #slopes = geoprofile.topo_profiles.profile_dirslopes
                plot_slope_min = -90.0

            """
            profiles_slope_min = np.nanmin(np.array(list(map(np.nanmin, slopes))))
            profiles_slope_max = np.nanmax(np.array(list(map(np.nanmax, slopes))))

            delta_slope = profiles_slope_max - profiles_slope_min
            plot_slope_min = profiles_slope_min - delta_slope * slope_padding
            plot_slope_max = profiles_slope_max + delta_slope * slope_padding
            """

        # plot topographic profile elevations

        if plot_height_choice:
            ndx_subplot += 1
            axes_elevation = plot_topo_profile_lines(
                geoprofile,
                plot_params,
                profile_window,
                grid_spec,
                ndx_subplot,
                'elevation',
                (plot_s_min, plot_s_max),
                (plot_z_min, plot_z_max),
                plot_params['filled_height'])
            if set_vertical_exaggeration:
                axes_elevation.set_aspect(vertical_exaggeration)
            axes_elevation.set_anchor('W')  # align left

        # plot topographic profile slopes

        if plot_slope_choice:
            ndx_subplot += 1
            axes_slopes = plot_topo_profile_lines(
                geoprofile,
                plot_params,
                profile_window,
                grid_spec,
                ndx_subplot,
                'slope',
                (plot_s_min, plot_s_max),
                (plot_slope_min, plot_slope_max),
                plot_params['filled_slope'])
            axes_slopes.set_anchor('W')  # align left

        # plot geological outcrop intersections

        if len(geoprofile.polygons_intersections) > 0:
            for line_intersection_value in geoprofile.polygons_intersections:
                plot_profile_polygon_intersection_line(plot_addit_params,
                                                       axes_elevation,
                                                       line_intersection_value)

        # plot geological attitudes intersections

        if len(geoprofile.profile_attitudes) > 0:
            for plane_attitude_set, color in zip(geoprofile.profile_attitudes, plot_addit_params["plane_attitudes_colors"]):
                plot_structural_attitude(plot_addit_params,
                                         axes_elevation,
                                         plot_s_max,
                                         vertical_exaggeration,
                                         plane_attitude_set,
                                         color)

        # plot geological traces projections

        if len(geoprofile.projected_lines) > 0:
            for curve_set, labels in zip(geoprofile.projected_lines, geoprofile.projected_lines_ids):
                plot_projected_line_set(axes_elevation,
                                        curve_set,
                                        labels)

        # plot line-profile intersections

        if len(geoprofile.line_intersections) > 0:
            plot_profile_lines_intersection_points(axes_elevation,
                                                   geoprofile.line_intersections)

    profile_window.canvas.fig.tight_layout()
    profile_window.canvas.draw()

    return profile_window







