
from math import radians, tan, sin, cos

import numpy as np


def matrScaling(
        scale_factor_x,
        scale_factor_y,
        scale_factor_z
):
    """
    
    :param scale_factor_x: 
    :param scale_factor_y: 
    :param scale_factor_z: 
    :return:

    Examples:
    """

    return np.array([(scale_factor_x, 0.0, 0.0),
                     (0.0, scale_factor_y, 0.0),
                     (0.0, 0.0, scale_factor_z)])


def matrHorizSimpleShear(
        phi_angle_degr,
        alpha_angle_degr
):
    """
    
    :param phi_angle_degr: 
    :param alpha_angle_degr: 
    :return:

    Examples:
    """

    phi_angle_rad = radians(phi_angle_degr)
    alpha_angle_rad = radians(alpha_angle_degr)

    gamma = tan(phi_angle_rad)
    sin_a = sin(alpha_angle_rad)
    cos_a = cos(alpha_angle_rad)

    return np.array([(1.0 - gamma * sin_a * cos_a, gamma * cos_a * cos_a, 0.0),
                     (-gamma * sin_a * sin_a, 1.0 + gamma * sin_a * cos_a, 0.0),
                     (0.0, 0.0, 1.0)])


def matrVertSimpleShear(
        phi_angle_degr,
        alpha_angle_degr
):
    """
    
    :param phi_angle_degr: 
    :param alpha_angle_degr: 
    :return:

    Examples:
    """

    phi_angle_rad = radians(phi_angle_degr)
    alpha_angle_rad = radians(alpha_angle_degr)

    gamma = tan(phi_angle_rad)
    sin_a = sin(alpha_angle_rad)
    cos_a = cos(alpha_angle_rad)

    return np.array([(1.0, 0.0, gamma * cos_a),
                     (0.0, 1.0, gamma * sin_a),
                     (0.0, 0.0, 1.0)])


if __name__ == "__main__":

    import doctest
    doctest.testmod()
