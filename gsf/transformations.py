# -*- coding: utf-8 -*-

from __future__ import division
from math import radians, sin, cos, tan

import numpy as np

from .array_utils import almost_zero
from .geometry import GAxis


class RefFrame(object):

    def __init__(self, versor_x, versor_y, versor_z):

        assert almost_zero(versor_x.scalar_product(versor_y))
        assert almost_zero(versor_x.scalar_product(versor_z))
        assert almost_zero(versor_y.scalar_product(versor_z))

        self.axes = [versor_x, versor_y, versor_z]

    def rotation_matrix(self, rotated_frame):

        for frame_axis in self.axes:
            assert almost_zero(frame_axis.lenght_3d() - 1.0)
        for frame_axis in rotated_frame.axes:
            assert almost_zero(frame_axis.lenght_3d() - 1.0)

        rot_matrix = np.zeros((3, 3)) * np.nan

        for i, rot_frame_versor in enumerate(rotated_frame.axes):
            for j, init_frame_versor in enumerate(self.axes):
                rot_matrix[i, j] = init_frame_versor.scalar_product(rot_frame_versor)

        return rot_matrix


def rotation_matrix(rot_axis_trend, rot_axis_plunge, rot_angle):

    phi = radians(rot_angle)

    rotation_versor = GAxis(rot_axis_trend, rot_axis_plunge).versor_3d()

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


def scaling_matrix(scale_factor_x, scale_factor_y, scale_factor_z):

    return np.array([(scale_factor_x, 0.0, 0.0),
                     (0.0, scale_factor_y, 0.0),
                     (0.0, 0.0, scale_factor_z)])


def simple_shear_horiz_matrix(phi_angle_degr, alpha_angle_degr):

    phi_angle_rad = radians(phi_angle_degr)
    alpha_angle_rad = radians(alpha_angle_degr)

    gamma = tan(phi_angle_rad)
    sin_a = sin(alpha_angle_rad)
    cos_a = cos(alpha_angle_rad)

    return np.array([(1.0 - gamma * sin_a * cos_a, gamma * cos_a * cos_a, 0.0),
                     (-gamma * sin_a * sin_a, 1.0 + gamma * sin_a * cos_a, 0.0),
                     (0.0, 0.0, 1.0)])


def simple_shear_vert_matrix(phi_angle_degr, alpha_angle_degr):

    phi_angle_rad = radians(phi_angle_degr)
    alpha_angle_rad = radians(alpha_angle_degr)

    gamma = tan(phi_angle_rad)
    sin_a = sin(alpha_angle_rad)
    cos_a = cos(alpha_angle_rad)

    return np.array([(1.0, 0.0, gamma * cos_a),
                     (0.0, 1.0, gamma * sin_a),
                     (0.0, 0.0, 1.0)])


def deformation_matrices(deform_params):

    deformation_matrices = []

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
            scal_matr = scaling_matrix(deform_param['parameters']['x factor'],
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

        deformation_matrices.append(deformation)

    return deformation_matrices

