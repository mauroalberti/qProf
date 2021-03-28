
from pygsf.orientations.orientations import *
from pygsf.geology.ptbaxes import *


"""

Dataset from: 
Kagan Y.Y., 1991
3-D rotation of double-couple earthquake sources
 
   values in tuples are (trend, plunge)
   
"""

# Original values

k91_fs_T_axis_vals, k91_fs_P_axis_vals = (120, 24), (232, 41)
k91_ss_T_axis_vals, k91_ss_P_axis_vals = (295, 55), (51, 17)

# Kagan (1991) provided rotation solutions

# provided values for each solution are azimuth, colatitude angle, rotation angle (all in degrees)

k91_rot_sol_1 = dict(az=24.8, colat_b=101.2, rot_ang=102.8)
k91_rot_sol_2 = dict(az=257.5, colat_b=79.7, rot_ang=104.3)
k91_rot_sol_3 = dict(az=144.8, colat_b=105.2, rot_ang=124.1)
k91_rot_sol_4 = dict(az=96.8, colat_b=16.7, rot_ang=165.9)

k91_rot_sols = [
    k91_rot_sol_1,
    k91_rot_sol_2,
    k91_rot_sol_3,
    k91_rot_sol_4]


def sols2rotaxis(rot_sol: dict) -> RotationAxis:

    az = rot_sol["az"]
    colat = rot_sol["colat_b"]
    rot_ang = rot_sol["rot_ang"]

    pl_from_N = plng2colatBottom(colat)

    gv = Direct(az, pl_from_N)

    return RotationAxis.fromDirect(gv, rot_ang)


# gsf-scope derived parameters

k91_fs_T_gaxis, k91_fs_P_gaxis = Axis.fromAzPl(*k91_fs_T_axis_vals), Axis.fromAzPl(*k91_fs_P_axis_vals)
k91_ss_T_gaxis, k91_ss_P_gaxis = Axis.fromAzPl(*k91_ss_T_axis_vals), Axis.fromAzPl(*k91_ss_P_axis_vals)

k91_fs_PTBaxes = PTBAxes(
    p_axis=k91_fs_P_gaxis,
    t_axis=k91_fs_T_gaxis)

k91_ss_PTBaxes = PTBAxes(
    p_axis=k91_ss_P_gaxis,
    t_axis=k91_ss_T_gaxis)

k91_fs_quater = k91_fs_PTBaxes.toQuatern()
k91_ss_quater = k91_ss_PTBaxes.toQuatern()

