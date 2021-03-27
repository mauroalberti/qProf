
import numbers

from pygsf.geometries.shapes.space3d import *


class ProfileAttitude:
    """
    Represent a geological attitude projected onto a vertical profile.
    """

    def __init__(
        self,
        rec_id: int,
        s: numbers.Real,
        z: numbers.Real,
        slope_degr: numbers.Real,
        down_sense: str,
        dist: numbers.Real,
        src_dip_dir: numbers.Real,
        src_dip_ang: numbers.Real
    ):
        """
        :param rec_id: the identifier of the observation.
        :type rec_id: int.
        :param s: the signed horizontal distance along the profile.
        :type s: numbers.Real (note: may exceed a profile range, both before, negative values, or after its end.
        :param z: the height of the attitude in the profile.
        :type z: numbers.Real.
        :param slope_degr: the slope of the attitude in the profile. Unit is degrees.
        :type slope_degr: numbers.Real.
        :param down_sense: downward sense, to the right or to the profile left.
        :type down_sense: str.
        :param dist: the distance between the attitude point and the point projection on the profile.
        :type: dist: numbers.Real
        :param src_dip_dir: source dip direction
        :type src_dip_dir: numbers.Real
        :param src_dip_ang: source dip angle
        :type src_dip_ang: numbers.Real
        """

        self.id = rec_id
        self.s = s
        self.z = z
        self.slope_degr = slope_degr
        self.down_sense = down_sense
        self.dist = dist
        self.src_dip_dir = src_dip_dir
        self.src_dip_ang = src_dip_ang

    def __repr__(self) -> str:
        """
        Creates the representation of a ProfileAttitude instance.

        :return: the representation of a ProfileAttitude instance.
        :rtype: str.
        """

        return f"ProfileAttitude(id={self.id}, s={self.s}, z={self.z}, " \
               f"slope_degr={self.slope_degr}, down_sense={self.down_sense}, " \
               f"dist={self.dist}, src_dip_dir={self.src_dip_dir}, src_dip_ang={self.src_dip_ang})"

    def create_segment_for_plot(
            self,
            profile_length: numbers.Real,
            vertical_exaggeration: numbers.Real = 1,
            segment_scale_factor: numbers.Real = 70.0):
        """

        :param profile_length:
        :param vertical_exaggeration:
        :param segment_scale_factor: the scale factor controlling the attitude segment length in the plot.
        :return:
        """

        ve = float(vertical_exaggeration)

        z0 = self.z

        h_dist = self.s
        slope_rad = radians(self.slope_degr)
        intersection_downward_sense = self.down_sense
        length = profile_length / segment_scale_factor

        s_slope = sin(float(slope_rad))
        c_slope = cos(float(slope_rad))

        if c_slope == 0.0:
            height_corr = length / ve
            structural_segment_s = [h_dist, h_dist]
            structural_segment_z = [z0 + height_corr, z0 - height_corr]
        else:
            t_slope = s_slope / c_slope
            width = length * c_slope

            length_exag = width * sqrt(1 + ve*ve * t_slope*t_slope)

            corr_width = width * length / length_exag
            corr_height = corr_width * t_slope

            structural_segment_s = [h_dist - corr_width, h_dist + corr_width]
            structural_segment_z = [z0 + corr_height, z0 - corr_height]

            if intersection_downward_sense == "left":
                structural_segment_z = [z0 - corr_height, z0 + corr_height]

        return structural_segment_s, structural_segment_z


