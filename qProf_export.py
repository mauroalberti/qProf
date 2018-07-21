# -*- coding: utf-8 -*-

from builtins import zip
from builtins import str
from builtins import map
from builtins import range
import os

import ogr


def write_rubberband_profile_lnshp(fileName, header_list, points, sr):

    shape_driver_name = "ESRI Shapefile"
    shape_driver = ogr.GetDriverByName(shape_driver_name)
    if shape_driver is None:
        return False, "%s driver is not available" % shape_driver_name

    try:
        datasource = shape_driver.CreateDataSource(str(fileName))
    except TypeError:
        datasource = shape_driver.CreateDataSource(str(fileName))

    if datasource is None:
        return False, "Creation of %s shapefile failed" % os.path.split(fileName)[1]

    layer = datasource.CreateLayer('profile', sr, geom_type=ogr.wkbLineString)
    if layer is None:
        return False, "Output layer creation failed"

    # creates required fields
    layer.CreateField(ogr.FieldDefn(header_list[0], ogr.OFTInteger))
    featureDefn = layer.GetLayerDefn()

    # loops through output records

    for ndx in range(len(points)-1):

        _, x0, y0 = points[ndx]
        _, x1, y1 = points[ndx+1]

        ln_feature = ogr.Feature(featureDefn)
        segment_2d = ogr.CreateGeometryFromWkt('LINESTRING(%f %f, %f %f)' % (x0, y0, x1, y1))
        ln_feature.SetGeometry(segment_2d)

        ln_feature.SetField(header_list[0], 1)
        layer.CreateFeature(ln_feature)

        ln_feature.Destroy()

    datasource.Destroy()

    return True, "done"


def write_generic_csv(output_filepath, header_list, parsed_results, sep=","):

    try:
        with open(str(output_filepath), 'w') as f:
            f.write(sep.join(header_list) + '\n')
            for rec in parsed_results:
                out_rec_string = ''
                for val in rec:
                    out_rec_string += str(val) + sep
                f.write(out_rec_string[:-1] + '\n')
        return True, "done"
    except Exception as e:
        return False, e.message


def write_topography_singledem_csv(fileName, header_list, multiprofile_dem_data, current_dem_ndx, sep=","):

    try:
        with open(str(fileName), 'w') as f:
            f.write(sep.join(header_list) + '\n')
            for prof_ndx, profile_data in enumerate(multiprofile_dem_data):
                for rec in profile_data:
                    rec_id, x, y, cum2ddist = rec[:4]
                    z = rec[3 + current_dem_ndx * 3 + 1]
                    cum3ddist = rec[3 + current_dem_ndx * 3 + 2]
                    slope = rec[3 + current_dem_ndx * 3 + 3]
                    outdata_list = list(map(str, [prof_ndx+1, rec_id, x, y, cum2ddist, z, cum3ddist, slope]))
                    f.write(sep.join(outdata_list) + "\n")
        return True, "done"
    except Exception as e:
        return False, e.message


def write_topography_singledem_ptshp(out_file_path, header_list, multiprofile_dem_data, current_dem_ndx, sr):

    shape_driver_name = "ESRI Shapefile"
    shape_driver = ogr.GetDriverByName(shape_driver_name)
    if shape_driver is None:
        return False, "%s driver is not available" % shape_driver_name

    try:
        datasource = shape_driver.CreateDataSource(str(out_file_path))
    except TypeError:
        datasource = shape_driver.CreateDataSource(str(out_file_path))

    if datasource is None:
        return False, "Creation of %s shapefile failed" % os.path.split(out_file_path)[1]

    layer = datasource.CreateLayer('profile', sr, geom_type=ogr.wkbPoint)
    if layer is None:
        return False, "Output layer creation failed"

        # creates required fields
    layer.CreateField(ogr.FieldDefn(header_list[0], ogr.OFTInteger))
    layer.CreateField(ogr.FieldDefn(header_list[1], ogr.OFTInteger))
    layer.CreateField(ogr.FieldDefn(header_list[2], ogr.OFTReal))
    layer.CreateField(ogr.FieldDefn(header_list[3], ogr.OFTReal))
    layer.CreateField(ogr.FieldDefn(header_list[4], ogr.OFTReal))
    layer.CreateField(ogr.FieldDefn(header_list[5], ogr.OFTReal))
    layer.CreateField(ogr.FieldDefn(header_list[6], ogr.OFTReal))
    layer.CreateField(ogr.FieldDefn(header_list[7], ogr.OFTReal))

    featureDefn = layer.GetLayerDefn()

    # loops through geoprofiles and output records

    for prof_ndx, profile_data in enumerate(multiprofile_dem_data):
        for rec in profile_data:
            rec_id, x, y, cumdist2D = rec[0], rec[1], rec[2], rec[3]
            z = rec[3 + current_dem_ndx * 3 + 1]
            cumdist3D = rec[3 + current_dem_ndx * 3 + 2]
            slopedegr = rec[3 + current_dem_ndx * 3 + 3]

            if z == "":
                continue

            pt_feature = ogr.Feature(featureDefn)

            pt = ogr.Geometry(ogr.wkbPoint25D)
            pt.SetPoint(0, x, y, z)
            pt_feature.SetGeometry(pt)

            pt_feature.SetField(header_list[0], prof_ndx)
            pt_feature.SetField(header_list[1], rec_id)
            pt_feature.SetField(header_list[2], x)
            pt_feature.SetField(header_list[3], y)
            pt_feature.SetField(header_list[4], cumdist2D)
            pt_feature.SetField(header_list[5], z)
            if cumdist3D != '':
                pt_feature.SetField(header_list[6], cumdist3D)
            if slopedegr != '':
                pt_feature.SetField(header_list[7], slopedegr)

            layer.CreateFeature(pt_feature)

            pt_feature.Destroy()

    datasource.Destroy()

    return True, "done"


def write_topography_singledem_lnshp(fileName, header_list, multiprofile_dem_data, current_dem_ndx, sr):

    shape_driver_name = "ESRI Shapefile"
    shape_driver = ogr.GetDriverByName(shape_driver_name)
    if shape_driver is None:
        return False, "%s driver is not available" % shape_driver_name

    try:
        datasource = shape_driver.CreateDataSource(str(fileName))
    except TypeError:
        datasource = shape_driver.CreateDataSource(str(fileName))

    if datasource is None:
        return False, "Creation of %s shapefile failed" % os.path.split(fileName)[1]

    layer = datasource.CreateLayer('profile', sr, geom_type=ogr.wkbLineString25D)
    if layer is None:
        return False, "Output layer creation failed"

        # creates required fields
    layer.CreateField(ogr.FieldDefn(header_list[0], ogr.OFTInteger))
    layer.CreateField(ogr.FieldDefn(header_list[1], ogr.OFTInteger))
    layer.CreateField(ogr.FieldDefn(header_list[4], ogr.OFTReal))
    layer.CreateField(ogr.FieldDefn(header_list[6], ogr.OFTReal))
    layer.CreateField(ogr.FieldDefn(header_list[7], ogr.OFTReal))

    featureDefn = layer.GetLayerDefn()

    # loops through output records

    for prof_ndx, profile_data in enumerate(multiprofile_dem_data):

        for ndx in range(len(profile_data) - 1):

            rec_a = profile_data[ndx]
            rec_b = profile_data[ndx + 1]

            rec_id = rec_a[0]
            x0, y0, z0 = rec_a[1], rec_a[2], rec_a[3 + current_dem_ndx * 3 + 1]
            x1, y1, z1 = rec_b[1], rec_b[2], rec_b[3 + current_dem_ndx * 3 + 1]
            cum3ddist = rec_b[3 + current_dem_ndx * 3 + 2]
            slope_degr = rec_b[3 + current_dem_ndx * 3 + 3]

            if z0 == '' or z1 == '':
                continue

            ln_feature = ogr.Feature(featureDefn)
            segment_3d = ogr.CreateGeometryFromWkt('LINESTRING(%f %f %f, %f %f %f)' % (x0, y0, z0, x1, y1, z1))
            ln_feature.SetGeometry(segment_3d)

            ln_feature.SetField(header_list[0], prof_ndx)
            ln_feature.SetField(header_list[1], rec_id)
            ln_feature.SetField(header_list[4], rec_b[3])

            if cum3ddist != '':
                ln_feature.SetField(header_list[6], cum3ddist)

            if slope_degr != '':
                ln_feature.SetField(header_list[7], slope_degr)

            layer.CreateFeature(ln_feature)

            ln_feature.Destroy()

    datasource.Destroy()

    return True, "done"


def write_topography_multidems_csv(fileName, multi_dems_headers, multiprofile_dem_data, sep=","):

    try:
        with open(str(fileName), 'w') as f:
            f.write(sep.join(multi_dems_headers) + '\n')
            for prof_ndx, profile_data in enumerate(multiprofile_dem_data):
                for rec in profile_data:
                    out_rec_string = sep.join(map(str, [prof_ndx+1]+rec))
                    f.write(out_rec_string + '\n')
        return True, "done"
    except Exception as e:
        return False, e.message


def write_topography_multidems_ptshp(fileName, multidems_headers, dem_names, multiprofile_dem_data, sr):

    shape_driver_name = "ESRI Shapefile"
    shape_driver = ogr.GetDriverByName(shape_driver_name)
    if shape_driver is None:
        return False, "%s driver is not available" % shape_driver_name

    try:
        datasource = shape_driver.CreateDataSource(str(fileName))
    except TypeError:
        datasource = shape_driver.CreateDataSource(str(fileName))

    if datasource is None:
        return False, "Creation of %s shapefile failed" % os.path.split(fileName)[1]

    layer = datasource.CreateLayer('profile', sr, geom_type=ogr.wkbPoint)
    if layer is None:
        return False, "Output layer creation failed"

    # creates required fields
    layer.CreateField(ogr.FieldDefn(multidems_headers[0], ogr.OFTInteger))
    layer.CreateField(ogr.FieldDefn(multidems_headers[1], ogr.OFTInteger))
    layer.CreateField(ogr.FieldDefn(multidems_headers[2], ogr.OFTReal))
    layer.CreateField(ogr.FieldDefn(multidems_headers[3], ogr.OFTReal))
    layer.CreateField(ogr.FieldDefn(multidems_headers[4], ogr.OFTReal))

    for dem_ndx in range(len(dem_names)):
        layer.CreateField(ogr.FieldDefn(multidems_headers[4 + dem_ndx * 3 + 1], ogr.OFTReal))
        layer.CreateField(ogr.FieldDefn(multidems_headers[4 + dem_ndx * 3 + 2], ogr.OFTReal))
        layer.CreateField(ogr.FieldDefn(multidems_headers[4 + dem_ndx * 3 + 3], ogr.OFTReal))

    featureDefn = layer.GetLayerDefn()

    field_names = []
    for i in range(featureDefn.GetFieldCount()):
        field_names.append(featureDefn.GetFieldDefn(i).GetName())

    assert len(multidems_headers) == len(field_names)

    # loops through output records

    for prof_ndx, profile_data in enumerate(multiprofile_dem_data):

        for rec in profile_data:

            pt_feature = ogr.Feature(featureDefn)

            pt = ogr.Geometry(ogr.wkbPoint)
            pt.SetPoint(0, rec[1], rec[2])
            pt_feature.SetGeometry(pt)

            pt_feature.SetField(field_names[0], prof_ndx)
            pt_feature.SetField(field_names[1], rec[0])
            pt_feature.SetField(field_names[2], rec[1])
            pt_feature.SetField(field_names[3], rec[2])
            pt_feature.SetField(field_names[4], rec[3])
            for dem_ndx in range(len(dem_names)):
                dem_height = rec[3 + dem_ndx * 3 + 1]
                if dem_height != '':
                    pt_feature.SetField(field_names[4 + dem_ndx * 3 + 1], dem_height)
                cum3ddist = rec[3 + dem_ndx * 3 + 2]
                if cum3ddist != '':
                    pt_feature.SetField(field_names[4 + dem_ndx * 3 + 2], cum3ddist)
                slope = rec[3 + dem_ndx * 3 + 3]
                if slope != '':
                    pt_feature.SetField(field_names[4 + dem_ndx * 3 + 3], slope)

            layer.CreateFeature(pt_feature)

            pt_feature.Destroy()

    datasource.Destroy()

    return True, "Done"


def write_topography_multidems_lnshp(fileName, header_list, dem_names, multiprofile_dem_data, sr):

    shape_driver_name = "ESRI Shapefile"
    shape_driver = ogr.GetDriverByName(shape_driver_name)
    if shape_driver is None:
        return False, "%s driver is not available" % shape_driver_name

    try:
        datasource = shape_driver.CreateDataSource(str(fileName))
    except TypeError:
        datasource = shape_driver.CreateDataSource(str(fileName))

    if datasource is None:
        return False, "Creation of %s shapefile failed" % os.path.split(fileName)[1]

    layer = datasource.CreateLayer('profile', sr, geom_type=ogr.wkbLineString)
    if layer is None:
        return False, "Output layer creation failed"

    # creates required fields

    layer.CreateField(ogr.FieldDefn(header_list[0], ogr.OFTInteger))  # prof ndx
    layer.CreateField(ogr.FieldDefn(header_list[1], ogr.OFTInteger))  # rec ndx
    layer.CreateField(ogr.FieldDefn(header_list[4], ogr.OFTReal))  # cum dist 2d
    for dem_ndx in range(len(dem_names)):
        layer.CreateField(ogr.FieldDefn(header_list[4 + dem_ndx * 3 + 1], ogr.OFTReal))
        layer.CreateField(ogr.FieldDefn(header_list[4 + dem_ndx * 3 + 2], ogr.OFTReal))
        layer.CreateField(ogr.FieldDefn(header_list[4 + dem_ndx * 3 + 3], ogr.OFTReal))

    featureDefn = layer.GetLayerDefn()

    field_names = []
    for i in range(featureDefn.GetFieldCount()):
        field_names.append(featureDefn.GetFieldDefn(i).GetName())

    # loops through output records
    for prof_ndx, profile_data in enumerate(multiprofile_dem_data):
        for ndx in range(len(profile_data) - 1):

            rec_a = profile_data[ndx]
            rec_b = profile_data[ndx + 1]

            rec_a_x, rec_a_y = rec_a[1], rec_a[2]
            rec_b_x, rec_b_y = rec_b[1], rec_b[2]

            ln_feature = ogr.Feature(featureDefn)

            segment_2d = ogr.CreateGeometryFromWkt('LINESTRING(%f %f, %f %f)' % (rec_a_x, rec_a_y, rec_b_x, rec_b_y))
            ln_feature.SetGeometry(segment_2d)

            ln_feature.SetField(field_names[0], prof_ndx)
            ln_feature.SetField(field_names[1], rec_a[0])
            ln_feature.SetField(field_names[4], rec_b[3])
            for dem_ndx, dem_name in enumerate(dem_names):
                dem_height = rec_b[3 + dem_ndx * 3 + 1]
                if dem_height != '':
                    ln_feature.SetField(field_names[2 + dem_ndx * 3 + 1], dem_height)
                cum3ddist = rec_b[3 + dem_ndx * 3 + 2]
                if cum3ddist != '':
                    ln_feature.SetField(field_names[2 + dem_ndx * 3 + 2], cum3ddist)
                slope = rec_b[3 + dem_ndx * 3 + 3]
                if slope != '':
                    ln_feature.SetField(field_names[2 + dem_ndx * 3 + 3], slope)

            layer.CreateFeature(ln_feature)
            ln_feature.Destroy()

    datasource.Destroy()

    return True, "Done"


def write_topography_gpx_ptshp(output_filepath, header_list, gpx_parsed_results, sr):

    shape_driver_name = "ESRI Shapefile"
    shape_driver = ogr.GetDriverByName(shape_driver_name)
    if shape_driver is None:
        return False, "%s driver is not available" % shape_driver_name

    try:
        datasource = shape_driver.CreateDataSource(str(output_filepath))
    except TypeError:
        datasource = shape_driver.CreateDataSource(str(output_filepath))

    if datasource is None:
        return False, "Creation of %s shapefile failed" % os.path.split(output_filepath)[1]

    layer = datasource.CreateLayer('profile', sr, geom_type=ogr.wkbPoint)
    if layer is None:
        return False, "Point layer creation failed"

        # creates required fields
    layer.CreateField(ogr.FieldDefn(header_list[0], ogr.OFTInteger))
    layer.CreateField(ogr.FieldDefn(header_list[1], ogr.OFTReal))
    layer.CreateField(ogr.FieldDefn(header_list[2], ogr.OFTReal))
    time_field = ogr.FieldDefn(header_list[3], ogr.OFTString)
    time_field.SetWidth(20)
    layer.CreateField(time_field)
    layer.CreateField(ogr.FieldDefn(header_list[4], ogr.OFTReal))
    layer.CreateField(ogr.FieldDefn(header_list[5], ogr.OFTReal))
    layer.CreateField(ogr.FieldDefn(header_list[6], ogr.OFTReal))
    layer.CreateField(ogr.FieldDefn(header_list[7], ogr.OFTReal))

    featureDefn = layer.GetLayerDefn()

    # loops through output records
    for rec in gpx_parsed_results:

        pt_feature = ogr.Feature(featureDefn)

        pt = ogr.Geometry(ogr.wkbPoint)
        pt.SetPoint(0, rec[2], rec[1])
        pt_feature.SetGeometry(pt)

        pt_feature.SetField(header_list[0], rec[0])
        pt_feature.SetField(header_list[1], rec[1])
        pt_feature.SetField(header_list[2], rec[2])

        pt_feature.SetField(header_list[3], str(rec[3]))
        if rec[4] != '':
            pt_feature.SetField(header_list[4], str(rec[4]))
        pt_feature.SetField(header_list[5], rec[5])
        if rec[6] != '':
            pt_feature.SetField(header_list[6], rec[6])
        if rec[7] != '':
            pt_feature.SetField(header_list[7], rec[7])

        layer.CreateFeature(pt_feature)

        pt_feature.Destroy()

    datasource.Destroy()

    return True, "done"


def write_topography_gpx_lnshp(output_filepath, header_list, gpx_parsed_results, sr):

    shape_driver_name = "ESRI Shapefile"
    shape_driver = ogr.GetDriverByName(shape_driver_name)
    if shape_driver is None:
        return False, "%s driver is not available" % shape_driver_name

    try:
        datasource = shape_driver.CreateDataSource(str(output_filepath))
    except TypeError:
        datasource = shape_driver.CreateDataSource(str(output_filepath))

    if datasource is None:
        return False, "Creation of %s shapefile failed" % os.path.split(output_filepath)[1]

    layer = datasource.CreateLayer('profile', sr, geom_type=ogr.wkbLineString25D)
    if layer is None:
        return False, "Output layer creation failed"

        # creates required fields
    layer.CreateField(ogr.FieldDefn(header_list[0], ogr.OFTInteger))
    time_beg_field = ogr.FieldDefn('time_beg', ogr.OFTString)
    time_beg_field.SetWidth(20)
    layer.CreateField(time_beg_field)
    time_end_field = ogr.FieldDefn('time_end', ogr.OFTString)
    time_end_field.SetWidth(20)
    layer.CreateField(time_end_field)
    layer.CreateField(ogr.FieldDefn(header_list[5], ogr.OFTReal))
    layer.CreateField(ogr.FieldDefn(header_list[6], ogr.OFTReal))
    layer.CreateField(ogr.FieldDefn(header_list[7], ogr.OFTReal))

    featureDefn = layer.GetLayerDefn()

    # loops through output records
    for ndx in range(len(gpx_parsed_results) - 1):

        rec_a = gpx_parsed_results[ndx]
        rec_b = gpx_parsed_results[ndx + 1]

        lon0, lat0, z0 = rec_a[2], rec_a[1], rec_a[4]
        lon1, lat1, z1 = rec_b[2], rec_b[1], rec_b[4]

        if z0 == '' or z1 == '':
            continue

        ln_feature = ogr.Feature(featureDefn)
        segment_3d = ogr.CreateGeometryFromWkt('LINESTRING(%f %f %f, %f %f %f)' % (lon0, lat0, z0, lon1, lat1, z1))
        ln_feature.SetGeometry(segment_3d)

        ln_feature.SetField(header_list[0], rec_a[0])
        ln_feature.SetField('time_beg', str(rec_a[3]))
        ln_feature.SetField('time_end', str(rec_b[3]))
        ln_feature.SetField(header_list[5], rec_b[5])
        if rec_b[6] != '':
            ln_feature.SetField(header_list[6], rec_b[6])
        if rec_b[7] != '':
            ln_feature.SetField(header_list[7], rec_b[7])

        layer.CreateFeature(ln_feature)

        ln_feature.Destroy()

    datasource.Destroy()

    return True, "done"


def write_geological_attitudes_ptshp(fileName, parsed_crosssect_results, sr):

    shape_driver_name = "ESRI Shapefile"
    shape_driver = ogr.GetDriverByName(shape_driver_name)
    if shape_driver is None:
        return False, "%s driver is not available" % shape_driver_name

    try:
        datasource = shape_driver.CreateDataSource(str(fileName))
    except TypeError:
        datasource = shape_driver.CreateDataSource(str(fileName))

    if datasource is None:
        return False, "Creation of %s shapefile failed" % os.path.split(fileName)[1]

    layer = datasource.CreateLayer('profile', sr, geom_type=ogr.wkbPoint25D)
    if layer is None:
        return False, "Output layer creation failed"

    # creates required fields
    layer.CreateField(ogr.FieldDefn('id', ogr.OFTString))
    layer.CreateField(ogr.FieldDefn('or_pt_x', ogr.OFTReal))
    layer.CreateField(ogr.FieldDefn('or_pt_y', ogr.OFTReal))
    layer.CreateField(ogr.FieldDefn('or_pt_z', ogr.OFTReal))
    layer.CreateField(ogr.FieldDefn('prj_pt_x', ogr.OFTReal))
    layer.CreateField(ogr.FieldDefn('prj_pt_y', ogr.OFTReal))
    layer.CreateField(ogr.FieldDefn('prj_pt_z', ogr.OFTReal))
    layer.CreateField(ogr.FieldDefn('s', ogr.OFTReal))
    layer.CreateField(ogr.FieldDefn('or_dpdir', ogr.OFTReal))
    layer.CreateField(ogr.FieldDefn('or_dpang', ogr.OFTReal))
    layer.CreateField(ogr.FieldDefn('tr_dpang', ogr.OFTReal))
    layer.CreateField(ogr.FieldDefn('tr_dpdir', ogr.OFTString))

    featureDefn = layer.GetLayerDefn()

    # loops through output records
    for rec in parsed_crosssect_results:
        pt_id, or_pt_x, or_pt_y, or_pt_z, pr_pt_x, pr_pt_y, pr_pt_z, s, or_dipdir, or_dipangle, tr_dipangle, tr_dipdir = rec

        pt_feature = ogr.Feature(featureDefn)

        pt = ogr.Geometry(ogr.wkbPoint25D)
        pt.SetPoint(0, pr_pt_x, pr_pt_y, pr_pt_z)
        pt_feature.SetGeometry(pt)

        pt_feature.SetField('id', str(pt_id))
        pt_feature.SetField('or_pt_x', or_pt_x)
        pt_feature.SetField('or_pt_y', or_pt_y)
        pt_feature.SetField('or_pt_z', or_pt_z)
        pt_feature.SetField('prj_pt_x', pr_pt_x)
        pt_feature.SetField('prj_pt_y', pr_pt_y)
        pt_feature.SetField('prj_pt_z', pr_pt_z)
        pt_feature.SetField('s', s)
        pt_feature.SetField('or_dpdir', or_dipdir)
        pt_feature.SetField('or_dpang', or_dipangle)
        pt_feature.SetField('tr_dpang', tr_dipangle)
        pt_feature.SetField('tr_dpdir', str(tr_dipdir))

        layer.CreateFeature(pt_feature)

        pt_feature.Destroy()

    datasource.Destroy()

    return True, "done"


def write_intersection_line_csv(output_filepath, header_list, parsed_results, sep=","):

    try:
        with open(str(output_filepath), 'w') as f:
            f.write(sep.join(header_list) + '\n')
            for classification, line3d, s_list in parsed_results:
                for pt, s in zip(line3d.pts, s_list):
                    out_values = [classification, s, pt.x, pt.y, pt.z]
                    out_val_strings = [str(val) for val in out_values]
                    f.write(sep.join(out_val_strings) + '\n')
        return True, "done"
    except Exception as e:
        return False, e.message


def write_intersection_line_ptshp(fileName, header_list, intersline_results, sr):

    shape_driver_name = "ESRI Shapefile"
    shape_driver = ogr.GetDriverByName(shape_driver_name)
    if shape_driver is None:
        return False, "%s driver is not available" % shape_driver_name

    try:
        datasource = shape_driver.CreateDataSource(str(fileName))
    except TypeError:
        datasource = shape_driver.CreateDataSource(str(fileName))

    if datasource is None:
        return False, "Creation of %s shapefile failed" % os.path.split(fileName)[1]

    layer = datasource.CreateLayer('profile', sr, geom_type=ogr.wkbPoint25D)
    if layer is None:
        return False, "Output layer creation failed"

    # creates required fields
    layer.CreateField(ogr.FieldDefn(header_list[0], ogr.OFTString))
    layer.CreateField(ogr.FieldDefn(header_list[1], ogr.OFTReal))
    layer.CreateField(ogr.FieldDefn(header_list[2], ogr.OFTReal))
    layer.CreateField(ogr.FieldDefn(header_list[3], ogr.OFTReal))
    layer.CreateField(ogr.FieldDefn(header_list[4], ogr.OFTReal))

    featureDefn = layer.GetLayerDefn()

    # loops through output records
    for rec_id, s, x, y, z in intersline_results:
        pt_feature = ogr.Feature(featureDefn)

        pt = ogr.Geometry(ogr.wkbPoint25D)
        pt.SetPoint(0, x, y, z)
        pt_feature.SetGeometry(pt)

        pt_feature.SetField(header_list[0], str(rec_id))
        pt_feature.SetField(header_list[1], s)
        pt_feature.SetField(header_list[2], x)
        pt_feature.SetField(header_list[3], y)
        pt_feature.SetField(header_list[4], z)

        layer.CreateFeature(pt_feature)

        pt_feature.Destroy()

    datasource.Destroy()

    return True, "done"


def write_intersection_polygon_lnshp(fileName, header_list, intersline_results, sr):

    shape_driver_name = "ESRI Shapefile"
    shape_driver = ogr.GetDriverByName(shape_driver_name)
    if shape_driver is None:
        return False, "%s driver is not available" % shape_driver_name

    try:
        datasource = shape_driver.CreateDataSource(str(fileName))
    except TypeError:
        datasource = shape_driver.CreateDataSource(str(fileName))

    if datasource is None:
        return False, "Creation of %s shapefile failed" % os.path.split(fileName)[1]

    layer = datasource.CreateLayer('profile', sr, geom_type=ogr.wkbLineString25D)
    if layer is None:
        return False, "Output layer creation failed"

    # creates required fields
    layer.CreateField(ogr.FieldDefn(header_list[0], ogr.OFTString))
    layer.CreateField(ogr.FieldDefn(header_list[1], ogr.OFTReal))

    featureDefn = layer.GetLayerDefn()

    # loops through output records

    for classification, line3d, s_list in intersline_results:

        assert len(line3d.pts) == len(s_list)

        # loops through output records

        for ndx in range(len(line3d.pts) - 1):
            rec_a = line3d.pts[ndx]
            rec_b = line3d.pts[ndx + 1]

            x0, y0, z0 = rec_a.x, rec_a.y, rec_a.z
            x1, y1, z1 = rec_b.x, rec_b.y, rec_b.z
            s = s_list[ndx + 1]

            ln_feature = ogr.Feature(featureDefn)
            segment_3d = ogr.CreateGeometryFromWkt('LINESTRING(%f %f %f, %f %f %f)' % (x0, y0, z0, x1, y1, z1))
            ln_feature.SetGeometry(segment_3d)

            ln_feature.SetField(header_list[0], str(classification))
            ln_feature.SetField(header_list[1], s)

            layer.CreateFeature(ln_feature)

            ln_feature.Destroy()

    datasource.Destroy()

    return True, "Output layer created"

