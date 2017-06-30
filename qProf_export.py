# -*- coding: utf-8 -*-

import os

import ogr


def write_generic_csv(output_filepath, header_list, parsed_results, sep=","):

    with open(unicode(output_filepath), 'w') as f:
        f.write(sep.join(header_list) + '\n')
        for rec in parsed_results:
            out_rec_string = ''
            for val in rec:
                out_rec_string += str(val) + sep
            f.write(out_rec_string[:-1] + '\n')


def write_topography_singleDEM_csv(fileName, header_list, export_data, current_dem_ndx, sep=","):

    with open(unicode(fileName), 'w') as f:
        f.write(sep.join(header_list) + '\n')
        for rec in export_data:
            rec_id, x, y, cum2ddist = rec[0], rec[1], rec[2], rec[3]
            z = rec[3 + current_dem_ndx * 3 + 1]
            cum3ddist = rec[3 + current_dem_ndx * 3 + 2]
            slope = rec[3 + current_dem_ndx * 3 + 3]

            outdata_list = [str(val) for val in [rec_id, x, y, cum2ddist, z, cum3ddist, slope]]
            f.write(sep.join(outdata_list) + "\n")


def write_topography_allDEMs_csv(fileName, header_list, export_data, sep=","):

    with open(unicode(fileName), 'w') as f:
        f.write(sep.join(header_list) + '\n')
        for rec in export_data:
            out_rec_string = ''
            for val in rec:
                out_rec_string += str(val) + sep
            f.write(out_rec_string[:-1] + '\n')


def write_intersection_polygon_lnshp(fileName, header_list, intersline_results, sr):

    shape_driver_name = "ESRI Shapefile"
    shape_driver = ogr.GetDriverByName(shape_driver_name)
    if shape_driver is None:
        return False, "%s driver is not available" % shape_driver_name

    try:
        shp_datasource = shape_driver.CreateDataSource(unicode(fileName))
    except TypeError:
        shp_datasource = shape_driver.CreateDataSource(str(fileName))

    if shp_datasource is None:
        return False, "Creation of %s shapefile failed" % os.path.split(fileName)[1]

    lnshp_layer = shp_datasource.CreateLayer('profile', sr, geom_type=ogr.wkbLineString25D)
    if lnshp_layer is None:
        return False, "Output layer creation failed"

    # creates required fields
    lnshp_layer.CreateField(ogr.FieldDefn(header_list[0], ogr.OFTString))
    lnshp_layer.CreateField(ogr.FieldDefn(header_list[1], ogr.OFTReal))

    lnshp_featureDefn = lnshp_layer.GetLayerDefn()

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

            ln_feature = ogr.Feature(lnshp_featureDefn)
            segment_3d = ogr.CreateGeometryFromWkt('LINESTRING(%f %f %f, %f %f %f)' % (x0, y0, z0, x1, y1, z1))
            ln_feature.SetGeometry(segment_3d)

            ln_feature.SetField(header_list[0], str(classification))
            ln_feature.SetField(header_list[1], s)

            lnshp_layer.CreateFeature(ln_feature)

            ln_feature.Destroy()

    shp_datasource.Destroy()

    return True, "Output layer created"


def write_line_csv(output_filepath, header_list, parsed_results, sep=","):

    with open(unicode(output_filepath), 'w') as f:
        f.write(sep.join(header_list) + '\n')
        for classification, line3d, s_list in parsed_results:
            for pt, s in zip(line3d.pts, s_list):
                out_values = [classification, s, pt.x, pt.y, pt.z]
                out_val_strings = [str(val) for val in out_values]
                f.write(sep.join(out_val_strings) + '\n')


def write_topography_allDEMs_ptshp(fileName, header_list, dem_names, export_data, sr):

    shape_driver_name = "ESRI Shapefile"
    shape_driver = ogr.GetDriverByName(shape_driver_name)
    if shape_driver is None:
        return False, "%s driver is not available" % shape_driver_name

    try:
        shp_datasource = shape_driver.CreateDataSource(unicode(fileName))
    except TypeError:
        shp_datasource = shape_driver.CreateDataSource(str(fileName))

    if shp_datasource is None:
        return False, "Creation of %s shapefile failed" % os.path.split(fileName)[1]

    ptshp_layer = shp_datasource.CreateLayer('profile', sr, geom_type=ogr.wkbPoint)
    if ptshp_layer is None:
        return False, "Output layer creation failed"

        # creates required fields
    ptshp_layer.CreateField(ogr.FieldDefn(header_list[0], ogr.OFTInteger))
    ptshp_layer.CreateField(ogr.FieldDefn(header_list[1], ogr.OFTReal))
    ptshp_layer.CreateField(ogr.FieldDefn(header_list[2], ogr.OFTReal))
    ptshp_layer.CreateField(ogr.FieldDefn(header_list[3], ogr.OFTReal))

    for dem_ndx in range(len(dem_names)):
        ptshp_layer.CreateField(ogr.FieldDefn(header_list[3 + dem_ndx * 3 + 1], ogr.OFTReal))
        ptshp_layer.CreateField(ogr.FieldDefn(header_list[3 + dem_ndx * 3 + 2], ogr.OFTReal))
        ptshp_layer.CreateField(ogr.FieldDefn(header_list[3 + dem_ndx * 3 + 3], ogr.OFTReal))

    ptshp_featureDefn = ptshp_layer.GetLayerDefn()

    field_names = []
    for i in range(ptshp_featureDefn.GetFieldCount()):
        field_names.append(ptshp_featureDefn.GetFieldDefn(i).GetName())

    assert len(header_list) == len(field_names)

    # loops through output records

    for rec in export_data:

        pt_feature = ogr.Feature(ptshp_featureDefn)

        pt = ogr.Geometry(ogr.wkbPoint)
        pt.SetPoint_2D(0, rec[1], rec[2])
        pt_feature.SetGeometry(pt)

        pt_feature.SetField(field_names[0], rec[0])
        pt_feature.SetField(field_names[1], rec[1])
        pt_feature.SetField(field_names[2], rec[2])
        pt_feature.SetField(field_names[3], rec[3])
        for dem_ndx in range(len(dem_names)):
            dem_height = rec[3 + dem_ndx * 3 + 1]
            if dem_height != '':
                pt_feature.SetField(field_names[3 + dem_ndx * 3 + 1], dem_height)
            cum3ddist = rec[3 + dem_ndx * 3 + 2]
            if cum3ddist != '':
                pt_feature.SetField(field_names[3 + dem_ndx * 3 + 2], cum3ddist)
            slope = rec[3 + dem_ndx * 3 + 3]
            if slope != '':
                pt_feature.SetField(field_names[3 + dem_ndx * 3 + 3], slope)

        ptshp_layer.CreateFeature(pt_feature)

        pt_feature.Destroy()
    shp_datasource.Destroy()

    return True, "Done"


def write_topography_allDEMs_lnshp(fileName, header_list, dem_names, export_data, sr):

    shape_driver_name = "ESRI Shapefile"
    shape_driver = ogr.GetDriverByName(shape_driver_name)
    if shape_driver is None:
        return False, "%s driver is not available" % shape_driver_name

    try:
        shp_datasource = shape_driver.CreateDataSource(unicode(fileName))
    except TypeError:
        shp_datasource = shape_driver.CreateDataSource(str(fileName))

    if shp_datasource is None:
        return False, "Creation of %s shapefile failed" % os.path.split(fileName)[1]

    lnshp_layer = shp_datasource.CreateLayer('profile', sr, geom_type=ogr.wkbLineString)
    if lnshp_layer is None:
        return False, "Output layer creation failed"

        # creates required fields
    lnshp_layer.CreateField(ogr.FieldDefn(header_list[0], ogr.OFTInteger))
    lnshp_layer.CreateField(ogr.FieldDefn(header_list[3], ogr.OFTReal))
    for dem_ndx in range(len(dem_names)):
        lnshp_layer.CreateField(ogr.FieldDefn(header_list[3 + dem_ndx * 3 + 1], ogr.OFTReal))
        lnshp_layer.CreateField(ogr.FieldDefn(header_list[3 + dem_ndx * 3 + 2], ogr.OFTReal))
        lnshp_layer.CreateField(ogr.FieldDefn(header_list[3 + dem_ndx * 3 + 3], ogr.OFTReal))

    lnshp_featureDefn = lnshp_layer.GetLayerDefn()

    field_names = []
    for i in range(lnshp_featureDefn.GetFieldCount()):
        field_names.append(lnshp_featureDefn.GetFieldDefn(i).GetName())

    # loops through output records
    for ndx in range(len(export_data) - 1):

        rec_a = export_data[ndx]
        rec_b = export_data[ndx + 1]

        rec_a_x, rec_a_y = rec_a[1], rec_a[2]
        rec_b_x, rec_b_y = rec_b[1], rec_b[2]

        ln_feature = ogr.Feature(lnshp_featureDefn)

        segment_2d = ogr.CreateGeometryFromWkt('LINESTRING(%f %f, %f %f)' % (rec_a_x, rec_a_y, rec_b_x, rec_b_y))
        ln_feature.SetGeometry(segment_2d)

        ln_feature.SetField(field_names[0], rec_a[0])
        ln_feature.SetField(field_names[1], rec_b[3])
        for dem_ndx, dem_name in enumerate(dem_names):
            dem_height = rec_b[3 + dem_ndx * 3 + 1]
            if dem_height != '':
                ln_feature.SetField(field_names[1 + dem_ndx * 3 + 1], dem_height)
            cum3ddist = rec_b[3 + dem_ndx * 3 + 2]
            if cum3ddist != '':
                ln_feature.SetField(field_names[1 + dem_ndx * 3 + 2], cum3ddist)
            slope = rec_b[3 + dem_ndx * 3 + 3]
            if slope != '':
                ln_feature.SetField(field_names[1 + dem_ndx * 3 + 3], slope)

        lnshp_layer.CreateFeature(ln_feature)
        ln_feature.Destroy()
    shp_datasource.Destroy()

    return True, "Done"
