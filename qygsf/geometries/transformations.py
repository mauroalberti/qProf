
from affine import Affine

from ..mathematics.scalars import *
from ..mathematics.deformations import *
from ..orientations.orientations import Axis


def gdal_to_affine(
    geotransform: Tuple[numbers.Real, numbers.Real, numbers.Real, numbers.Real, numbers.Real, numbers.Real]
) -> Affine:
    """
    Create an affine transformation
    from a GDAL geotransform tuple.

    """

    return Affine.from_gdal(*geotransform)


def forward_transformation(
    trans: Affine,
    row: numbers.Real,
    col: numbers.Real
) -> Tuple[numbers.Real, numbers.Real]:
    """
    Calculate the x, y coordinates given an affine transformation
    and the row, col values.

    """

    return trans * (col, row)


def simple_shear_horiz_matrix(phi_angle_degr, alpha_angle_degr):

    phi_angle_rad = radians(phi_angle_degr)
    alpha_angle_rad = radians(alpha_angle_degr)

    gamma = tan(phi_angle_rad)
    sin_a = sin(alpha_angle_rad)
    cos_a = cos(alpha_angle_rad)

    return np.array([(1.0 - gamma * sin_a * cos_a, gamma * cos_a * cos_a, 0.0),
                     (-gamma * sin_a * sin_a, 1.0 + gamma * sin_a * cos_a, 0.0),
                     (0.0, 0.0, 1.0)])


def simple_shear_vert_matrix(
        phi_angle_degr,
        alpha_angle_degr
):

    phi_angle_rad = radians(phi_angle_degr)
    alpha_angle_rad = radians(alpha_angle_degr)

    gamma = tan(phi_angle_rad)
    sin_a = sin(alpha_angle_rad)
    cos_a = cos(alpha_angle_rad)

    return np.array([(1.0, 0.0, gamma * cos_a),
                     (0.0, 1.0, gamma * sin_a),
                     (0.0, 0.0, 1.0)])


def deformation_matrices(deform_params):

    def_matrices = []

    for deform_param in deform_params:
        if deform_param['type'] == 'displacement':
            displ_x = deform_param['parameters']['delta_x']
            displ_y = deform_param['parameters']['delta_y']
            displ_z = deform_param['parameters']['delta_z']
            deformation = {'increment': 'additive',
                           'matrix': np.array([displ_x, displ_y, displ_z])}
        elif deform_param['type'] == 'rotation':
            rot_matr = rotation_matrix(deform_param['parameters']['rotation axis trend'],
                                       deform_param['parameters']['rotation axis plunge'],
                                       deform_param['parameters']['rotation angle'])
            deformation = {'increment': 'multiplicative',
                           'matrix': rot_matr,
                           'shift_pt': np.array([deform_param['parameters']['center x'],
                                                 deform_param['parameters']['center y'],
                                                 deform_param['parameters']['center z']])}
        elif deform_param['type'] == 'scaling':
            scal_matr = matrScaling(deform_param['parameters']['x factor'],
                                       deform_param['parameters']['y factor'],
                                       deform_param['parameters']['z factor'])
            deformation = {'increment': 'multiplicative',
                           'matrix': scal_matr,
                           'shift_pt': np.array([deform_param['parameters']['center x'],
                                                 deform_param['parameters']['center y'],
                                                 deform_param['parameters']['center z']])}
        elif deform_param['type'] == 'simple shear - horizontal':
            simple_shear_horiz_matr = simple_shear_horiz_matrix(deform_param['parameters']['psi angle (degr.)'],
                                                                deform_param['parameters']['alpha angle (degr.)'])
            deformation = {'increment': 'multiplicative',
                           'matrix': simple_shear_horiz_matr,
                           'shift_pt': np.array([deform_param['parameters']['center x'],
                                                 deform_param['parameters']['center y'],
                                                 deform_param['parameters']['center z']])}
        elif deform_param['type'] == 'simple shear - vertical':
            simple_shear_vert_matr = simple_shear_vert_matrix(deform_param['parameters']['psi angle (degr.)'],
                                                              deform_param['parameters']['alpha angle (degr.)'])
            deformation = {'increment': 'multiplicative',
                           'matrix': simple_shear_vert_matr,
                           'shift_pt': np.array([deform_param['parameters']['center x'],
                                                 deform_param['parameters']['center y'],
                                                 deform_param['parameters']['center z']])}
        else:
            continue

        def_matrices.append(deformation)

    return def_matrices


def rotation_matrix(
        rot_axis_trend,
        rot_axis_plunge,
        rot_angle
):

    phi = radians(rot_angle)

    rotation_versor = Axis(rot_axis_trend, rot_axis_plunge).as_direction().as_versor()

    l = rotation_versor.x
    m = rotation_versor.y
    n = rotation_versor.z

    cos_phi = cos(phi)
    sin_phi = sin(phi)

    a11 = cos_phi + ((l * l) * (1 - cos_phi))
    a12 = ((l * m) * (1 - cos_phi)) - (n * sin_phi)
    a13 = ((l * n) * (1 - cos_phi)) + (m * sin_phi)

    a21 = ((l * m) * (1 - cos_phi)) + (n * sin_phi)
    a22 = cos_phi + ((m * m) * (1 - cos_phi))
    a23 = ((m * n) * (1 - cos_phi)) - (l * sin_phi)

    a31 = ((l * n) * (1 - cos_phi)) - (m * sin_phi)
    a32 = ((m * n) * (1 - cos_phi)) + (l * sin_phi)
    a33 = cos_phi + ((n * n) * (1 - cos_phi))

    return np.array([(a11, a12, a13),
                     (a21, a22, a23),
                     (a31, a32, a33)])

