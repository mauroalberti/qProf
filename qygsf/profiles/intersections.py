from .geoprofiles import *
from ..geometries.shapes.space3d import *
from ..geometries.shapes.space4d import *


def calculate_distance_with_sign(
        projected_point,
        section_init_pt,
        section_vector
):

    assert projected_point.z != np.nan
    assert projected_point.z is not None

    projected_vector = Segment4D(section_init_pt, projected_point).vector()
    cos_alpha = section_vector.cos_angle(projected_vector)

    return projected_vector.len_3d * cos_alpha


def get_intersection_slope(intersection_versor_3d, section_vector):

    slope_radians = abs(radians(intersection_versor_3d.slope))
    scalar_product_for_downward_sense = section_vector.sp(intersection_versor_3d.downward)
    if scalar_product_for_downward_sense > 0.0:
        intersection_downward_sense = "right"
    elif scalar_product_for_downward_sense == 0.0:
        intersection_downward_sense = "vertical"
    else:
        intersection_downward_sense = "left"

    return slope_radians, intersection_downward_sense


def calculate_intersection_versor(section_cartes_plane, structural_cartes_plane):

    return section_cartes_plane.inters_versor(structural_cartes_plane)


def calculate_nearest_intersection(intersection_versor_3d, section_cartes_plane, structural_cartes_plane,
                                   structural_pt):

    dummy_inters_point = section_cartes_plane.inters_point(structural_cartes_plane)
    dummy_structural_vector = Segment4D(dummy_inters_point, structural_pt).vector()
    dummy_distance = dummy_structural_vector.sp(intersection_versor_3d)
    offset_vector = intersection_versor_3d.scale(dummy_distance)

    return Point4D(dummy_inters_point.x + offset_vector.x,
                   dummy_inters_point.y + offset_vector.y,
                   dummy_inters_point.z + offset_vector.z)


def calculate_axis_intersection(map_axis, section_cartes_plane, structural_pt):

    axis_versor = map_axis.as_gvect().versor()
    l, m, n = axis_versor.x, axis_versor.y, axis_versor.z
    axis_param_line = ParamLine3D(structural_pt, l, m, n)
    return axis_param_line.intersect_cartes_plane(section_cartes_plane)


def map_measure_to_section(structural_rec, section_data, map_axis=None):

    # extract source data
    structural_pt, structural_plane, structural_pt_id = structural_rec
    section_init_pt, section_cartes_plane, section_vector = section_data['init_pt'], section_data['cartes_plane'], \
                                                            section_data['vector']

    # transform geological plane attitude into Cartesian plane      
    structural_cartes_plane = structural_plane.plane(structural_pt)

    ## intersection versor
    intersection_versor_3d = calculate_intersection_versor(section_cartes_plane, structural_cartes_plane)

    # calculate slope of geological plane onto section plane
    slope_radians, intersection_downward_sense = get_intersection_slope(intersection_versor_3d, section_vector)

    # intersection point
    if map_axis is None:
        intersection_point_3d = calculate_nearest_intersection(intersection_versor_3d, section_cartes_plane,
                                                               structural_cartes_plane, structural_pt)
    else:
        intersection_point_3d = calculate_axis_intersection(map_axis, section_cartes_plane, structural_pt)

    # horizontal spat_distance between projected structural point and profile start
    signed_distance_from_section_start = calculate_distance_with_sign(intersection_point_3d, section_init_pt,
                                                                      section_vector)

    # solution for current structural point
    return PlaneAttitude(structural_pt_id,
                         structural_pt,
                         structural_plane,
                         intersection_point_3d,
                         slope_radians,
                         intersection_downward_sense,
                         signed_distance_from_section_start)


def map_struct_pts_on_section(structural_data, section_data, mapping_method):
    """
    defines:
        - 2D x-y location in section
        - plane-plane segment intersection
    """

    if mapping_method['method'] == 'nearest':
        return [map_measure_to_section(structural_rec, section_data) for structural_rec in structural_data]

    if mapping_method['method'] == 'common axis':
        map_axis = Axis(mapping_method['trend'], mapping_method['plunge'])
        return [map_measure_to_section(structural_rec, section_data, map_axis) for structural_rec in structural_data]

    if mapping_method['method'] == 'individual axes':
        assert len(mapping_method['individual_axes_values']) == len(structural_data)
        result = []
        for structural_rec, (trend, plunge) in zip(structural_data, mapping_method['individual_axes_values']):
            try:
                map_axis = Axis(trend, plunge)
                result.append(map_measure_to_section(structural_rec, section_data, map_axis))
            except:
                continue
        return result


class IntersectionParameters(object):
    """
    IntersectionParameters class.
    Manages the metadata for spdata results (DEM source filename, source point, plane attitude.
    
    """

    def __init__(self, sourcename, srcPt, srcPlaneAttitude):
        self._sourcename = sourcename
        self._srcPt = srcPt
        self._srcPlaneAttitude = srcPlaneAttitude


class Intersections(object):

    def __init__(self):

        self.parameters = None

        self.xcoords_x = []
        self.xcoords_y = []
        self.ycoords_x = []
        self.ycoords_y = []

        self.links = None
        self.networks = {}

    def get_intersections(self):
        """
        Initialize a structured array of the possible and found links for each intersection.
        It will store a list of the possible connections for each intersection,
        together with the found connections.       
        """

        # data type for structured array storing intersection parameters
        dt = np.dtype([('id', np.uint16),
                       ('i', np.uint16),
                       ('j', np.uint16),
                       ('pi_dir', np.str_, 1),
                       ('conn_from', np.uint16),
                       ('conn_to', np.uint16),
                       ('start', np.bool_)
                       ])

        # number of valid intersections
        num_intersections = len(list(self.xcoords_x[np.logical_not(np.isnan(self.xcoords_x))])) + \
                            len(list(self.ycoords_y[np.logical_not(np.isnan(self.ycoords_y))]))

        # creation and initialization of structured array of valid intersections in the x-direction
        links = np.zeros(num_intersections, dtype=dt)

        # filling array with values

        curr_ndx = 0
        for i in range(self.xcoords_x.shape[0]):
            for j in range(self.xcoords_x.shape[1]):
                if not isnan(self.xcoords_x[i, j]):
                    links[curr_ndx] = (curr_ndx + 1, i, j, 'x', 0, 0, False)
                    curr_ndx += 1

        for i in range(self.ycoords_y.shape[0]):
            for j in range(self.ycoords_y.shape[1]):
                if not isnan(self.ycoords_y[i, j]):
                    links[curr_ndx] = (curr_ndx + 1, i, j, 'y', 0, 0, False)
                    curr_ndx += 1

        return links

    def set_neighbours(self):

        # shape of input arrays (equal shapes)
        num_rows, num_cols = self.xcoords_x.shape

        # dictionary storing intersection links
        neighbours = {}

        # search and connect intersection points   
        for curr_ndx in range(self.links.shape[0]):

            # get current point location (i, j) and direction type (pi_dir)
            curr_id = self.links[curr_ndx]['id']
            curr_i = self.links[curr_ndx]['i']
            curr_j = self.links[curr_ndx]['j']
            curr_dir = self.links[curr_ndx]['pi_dir']

            # check possible connected spdata
            near_intersections = []

            if curr_dir == 'x':

                if curr_i < num_rows - 1 and curr_j < num_cols - 1:

                    try:  # -- A
                        id_link = self.links[(self.links['i'] == curr_i + 1) & \
                                             (self.links['j'] == curr_j + 1) & \
                                             (self.links['pi_dir'] == 'y')]['id']
                        if len(list(id_link)) == 1:
                            near_intersections.append(id_link[0])
                    except:
                        pass
                    try:  # -- B
                        id_link = self.links[(self.links['i'] == curr_i + 1) & \
                                             (self.links['j'] == curr_j) & \
                                             (self.links['pi_dir'] == 'x')]['id']
                        if len(list(id_link)) == 1:
                            near_intersections.append(id_link[0])
                    except:
                        pass
                    try:  # -- C
                        id_link = self.links[(self.links['i'] == curr_i + 1) & \
                                             (self.links['j'] == curr_j) & \
                                             (self.links['pi_dir'] == 'y')]['id']
                        if len(list(id_link)) == 1:
                            near_intersections.append(id_link[0])
                    except:
                        pass

                if curr_i > 0 and curr_j < num_cols - 1:

                    try:  # -- E
                        id_link = self.links[(self.links['i'] == curr_i) & \
                                             (self.links['j'] == curr_j) & \
                                             (self.links['pi_dir'] == 'y')]['id']
                        if len(list(id_link)) == 1:
                            near_intersections.append(id_link[0])
                    except:
                        pass
                    try:  # -- F
                        id_link = self.links[(self.links['i'] == curr_i - 1) & \
                                             (self.links['j'] == curr_j) & \
                                             (self.links['pi_dir'] == 'x')]['id']
                        if len(list(id_link)) == 1:
                            near_intersections.append(id_link[0])
                    except:
                        pass
                    try:  # -- G
                        id_link = self.links[(self.links['i'] == curr_i) & \
                                             (self.links['j'] == curr_j + 1) & \
                                             (self.links['pi_dir'] == 'y')]['id']
                        if len(list(id_link)) == 1:
                            near_intersections.append(id_link[0])
                    except:
                        pass

            if curr_dir == 'y':

                if curr_i > 0 and curr_j < num_cols - 1:

                    try:  # -- D
                        id_link = self.links[(self.links['i'] == curr_i) & \
                                             (self.links['j'] == curr_j) & \
                                             (self.links['pi_dir'] == 'x')]['id']
                        if len(list(id_link)) == 1:
                            near_intersections.append(id_link[0])
                    except:
                        pass
                    try:  # -- F
                        id_link = self.links[(self.links['i'] == curr_i - 1) & \
                                             (self.links['j'] == curr_j) & \
                                             (self.links['pi_dir'] == 'x')]['id']
                        if len(list(id_link)) == 1:
                            near_intersections.append(id_link[0])
                    except:
                        pass
                    try:  # -- G
                        id_link = self.links[(self.links['i'] == curr_i) & \
                                             (self.links['j'] == curr_j + 1) & \
                                             (self.links['pi_dir'] == 'y')]['id']
                        if len(list(id_link)) == 1:
                            near_intersections.append(id_link[0])
                    except:
                        pass

                if curr_i > 0 and curr_j > 0:

                    try:  # -- H
                        id_link = self.links[(self.links['i'] == curr_i) & \
                                             (self.links['j'] == curr_j - 1) & \
                                             (self.links['pi_dir'] == 'x')]['id']
                        if len(list(id_link)) == 1:
                            near_intersections.append(id_link[0])
                    except:
                        pass
                    try:  # -- I
                        id_link = self.links[(self.links['i'] == curr_i) & \
                                             (self.links['j'] == curr_j - 1) & \
                                             (self.links['pi_dir'] == 'y')]['id']
                        if len(list(id_link)) == 1:
                            near_intersections.append(id_link[0])
                    except:
                        pass
                    try:  # -- L
                        id_link = self.links[(self.links['i'] == curr_i - 1) & \
                                             (self.links['j'] == curr_j - 1) & \
                                             (self.links['pi_dir'] == 'x')]['id']
                        if len(list(id_link)) == 1:
                            near_intersections.append(id_link[0])
                    except:
                        pass

            neighbours[curr_id] = near_intersections

        return neighbours

    def follow_path(self, start_id):
        """
        Creates a path of connected intersections from a given start intersection.
        
        """
        from_id = start_id

        while self.links[from_id - 1]['conn_to'] == 0:

            conns = self.neighbours[from_id]
            num_conn = len(conns)
            if num_conn == 0:
                raise Exception('no connected intersection')
            elif num_conn == 1:
                if self.links[conns[0] - 1]['conn_from'] == 0 and self.links[conns[0] - 1]['conn_to'] != from_id:
                    to_id = conns[0]
                else:
                    raise Exception('no free connection')
            elif num_conn == 2:
                if self.links[conns[0] - 1]['conn_from'] == 0 and self.links[conns[0] - 1]['conn_to'] != from_id:
                    to_id = conns[0]
                elif self.links[conns[1] - 1]['conn_from'] == 0 and self.links[conns[1] - 1]['conn_to'] != from_id:
                    to_id = conns[1]
                else:
                    raise Exception('no free connection')
            else:
                raise Exception('multiple connection')

            # set connection
            self.links[to_id - 1]['conn_from'] = from_id
            self.links[from_id - 1]['conn_to'] = to_id

            # prepare for next step
            from_id = to_id

    def path_closed(self, start_id):

        from_id = start_id

        while self.links[from_id - 1]['conn_to'] != 0:

            to_id = self.links[from_id - 1]['conn_to']

            if to_id == start_id: return True

            from_id = to_id

        return False

    def invert_path(self, start_id):

        self.links[start_id - 1]['start'] = False

        curr_id = start_id

        while curr_id != 0:

            prev_from_id = self.links[curr_id - 1]['conn_from']
            prev_to_id = self.links[curr_id - 1]['conn_to']

            self.links[curr_id - 1]['conn_from'] = prev_to_id
            self.links[curr_id - 1]['conn_to'] = prev_from_id

            if self.links[curr_id - 1]['conn_from'] == 0:
                self.links[curr_id - 1]['start'] = True

            curr_id = prev_to_id

        return

    def patch_path(self, start_id):

        if self.path_closed(start_id): return

        from_id = start_id

        conns = self.neighbours[from_id]
        try:
            conns.remove(self.links[from_id - 1]['conn_to'])
        except:
            pass

        num_conn = len(conns)

        if num_conn != 1: return

        new_toid = self.links[conns[0] - 1]

        if self.links[new_toid]['conn_to'] > 0 \
                and self.links[new_toid]['conn_to'] != from_id \
                and self.links[new_toid]['conn_from'] == 0:

            if self.path_closed(new_toid): return
            self.invert_path(from_id)
            self.links[from_id - 1]['conn_to'] = new_toid
            self.links[new_toid - 1]['conn_from'] = from_id
            self.links[new_toid - 1]['start'] = False

    def define_paths(self):

        # simple networks starting from border
        for ndx in range(self.links.shape[0]):

            if len(self.neighbours[ndx + 1]) != 1 or \
                            self.links[ndx]['conn_from'] > 0 or \
                            self.links[ndx]['conn_to'] > 0:
                continue

            try:
                self.follow_path(ndx + 1)
            except:
                continue

        # inner, simple networks

        for ndx in range(self.links.shape[0]):

            if len(self.neighbours[ndx + 1]) != 2 or \
                            self.links[ndx]['conn_to'] > 0 or \
                            self.links[ndx]['start'] == True:
                continue

            try:
                self.links[ndx]['start'] = True
                self.follow_path(ndx + 1)
            except:
                continue

        # inner, simple networks, connection of FROM

        for ndx in range(self.links.shape[0]):

            if len(self.neighbours[ndx + 1]) == 2 and \
                            self.links[ndx]['conn_from'] == 0:
                try:
                    self.patch_path(ndx + 1)
                except:
                    continue

    def define_networks(self):
        """
        Creates list of connected intersections,
        to output as line shapefile
        
        
        """

        pid = 0
        networks = {}

        # open, simple networks
        for ndx in range(self.links.shape[0]):

            if len(self.neighbours[ndx + 1]) != 1: continue

            network_list = []

            to_ndx = ndx + 1

            while to_ndx != 0:
                network_list.append(to_ndx)

                to_ndx = self.links[to_ndx - 1]['conn_to']

            if len(network_list) > 1:
                pid += 1

                networks[pid] = network_list

        # closed, simple networks
        for ndx in range(self.links.shape[0]):

            if len(self.neighbours[ndx + 1]) != 2 or \
                            self.links[ndx]['start'] == False:
                continue

            start_id = ndx + 1

            network_list = []

            to_ndx = ndx + 1

            while to_ndx != 0:

                network_list.append(to_ndx)

                to_ndx = self.links[to_ndx - 1]['conn_to']

                if to_ndx == start_id:
                    network_list.append(to_ndx)
                    break

            if len(network_list) > 1:
                pid += 1

                networks[pid] = network_list

        return networks
