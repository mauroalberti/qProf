
from .faults import *

from ..orientations.orientations import RotationAxis, sortRotations


class PTBAxes(object):
    """
    Represent the triad of P, T and B kinematic axes.
    It can also calculate the M plane.
    """

    def __init__(self, p_axis: Axis, t_axis: Axis):
        """
        Create a new PTBAxes instances, given the two
        P and T axes (provided as Axis instances).
        T and P axes are recalculated to be strictly orthogonal,
        based on fixed T orientation.

        Example:
          >>> PTBAxes(p_axis=Axis(0, 0), t_axis=Axis(90, 0))
          PTBAxes(P Axis(az: 0.00°, pl: -0.00°), T Axis(az: 90.00°, pl: 0.00°))
          >>> PTBAxes(p_axis=Axis(0, 0), t_axis=Axis(80, 0))
          Traceback (most recent call last):
          ...
          Exception: P and T axes must be sub-orthogonal
        """

        if not (isinstance(p_axis, Axis) and isinstance(t_axis, Axis)):
            raise Exception("P and T axes must be of type Axis")

        if not p_axis.is_sub_orthogonal(t_axis):
            raise Exception("P and T axes must be sub-orthogonal")

        b_vect = t_axis.normal_versor(p_axis)

        self._t_versor = t_axis.as_versor()
        self._p_versor = b_vect.cross_product(self._t_versor).versor()

    @classmethod
    def fromVects(cls,
        p_vector: Vect3D,
        t_vector: Vect3D
    ) -> 'PTBAxes':
        """
        Class method to create a PTBAxes instance from T and P axis vectors.
        Vectors are not required to be normalized but are required to be
        sub-orthogonal.

        :param t_vector: the asVect representing the T axis (Vect instance).
        :param p_vector: the asVect representing the P axis (Vect instance).
        :return: a PTBAxes instance.

        Example:
          >>> PTBAxes.fromVects(p_vector=Vect3D(0,1,0), t_vector=Vect3D(1,0,0))
          PTBAxes(P Axis(az: 0.00°, pl: -0.00°), T Axis(az: 90.00°, pl: -0.00°))
          >>> PTBAxes.fromVects(p_vector=Vect3D(1,0,0), t_vector=Vect3D(0,0,-1))
          PTBAxes(P Axis(az: 90.00°, pl: -0.00°), T Axis(az: 0.00°, pl: 90.00°))
          >>> PTBAxes.fromVects(p_vector=Vect3D(-1,1,0), t_vector=Vect3D(1,1,0))
          PTBAxes(P Axis(az: 315.00°, pl: -0.00°), T Axis(az: 45.00°, pl: -0.00°))
          >>> PTBAxes.fromVects(p_vector=Vect3D(0.5, 1, 0), t_vector=Vect3D(1, 1, 0))
          Traceback (most recent call last):
          ...
          Exception: P and T axes must be sub-orthogonal
        """

        if not isinstance(t_vector, Vect3D):
            raise Exception("T vector must be of type Vect")

        if t_vector.is_close_to_zero:
            raise Exception("T vector cannot be zero-valued")

        if not isinstance(p_vector, Vect3D):
            raise Exception("P vector must be of type Vect")

        if p_vector.is_close_to_zero:
            raise Exception("P vector cannot be zero-valued")

        return cls(
            p_axis=Axis.fromVect(p_vector),
            t_axis=Axis.fromVect(t_vector))

    @classmethod
    def fromFaultSlick(cls, fault: Fault, slick_ndx: numbers.Integral=0) -> 'PTBAxes':
        """
        Class method to calculate P-T axes from a FaultSlick instance.
        Return P and T axes and a third Boolean variable,
        indicating whether the P-T derivation is from a slickenline with a known movement sense (True)
        or with unknown/uncertain movement sense (False).

        Example:
          >>> PTBAxes.fromFaultSlick(Fault(90, 45, slickenlines=[Slick(90, 45)]))
          PTBAxes(P Axis(az: 0.00°, pl: -90.00°), T Axis(az: 90.00°, pl: -0.00°))
        """

        if not isinstance(fault, Fault):
            raise Exception("First argument for fault input must have type Fault")

        slick = fault.slick(slick_ndx)

        if slick.hasUnknownSense:
            raise Exception("Slickenline must have knonw movement sense")

        s_versor = slick.geom.as_versor()
        f_versor = fault.plane.normDirectFrwrd().as_versor()

        p_versor = (f_versor - s_versor).versor()
        t_versor = (f_versor + s_versor).versor()

        return cls.fromVects(
            t_vector=t_versor,
            p_vector=p_versor)

    @classmethod
    def fromQuatern(cls, quaternion: Quaternion) -> 'PTBAxes':
        """
        Creates a PTBAxes instance from a given quaternion.
        Formula extracted from eq. 10 in:
        Kagan, Y.Y, 1991. 3-D rotation of double-couple earthquake sources.

        :param quaternion: a Quaternion instance.
        :return:a PTBAxes instance.
        """

        if not isinstance(quaternion, Quaternion):
            raise Exception("Input argument must be of Quaternion type")

        q0, q1, q2, q3 = quaternion.normalize().components()

        q0q0 = q0*q0
        q0q1 = q0*q1
        q0q2 = q0*q2
        q0q3 = q0*q3

        q1q1 = q1*q1
        q1q2 = q1*q2
        q1q3 = q1*q3

        q2q2 = q2*q2
        q2q3 = q2*q3

        q3q3 = q3*q3

        t1 = q0q0 + q1q1 - q2q2 - q3q3
        t2 = 2*(q1q2 + q0q3)
        t3 = 2*(q1q3 - q0q2)

        p1 = 2*(q1q2 - q0q3)
        p2 = q0q0 - q1q1 + q2q2 - q3q3
        p3 = 2*(q2q3 + q0q1)

        t_vector = Vect3D(t1, t2, t3)
        p_vector = Vect3D(p1, p2, p3)

        return PTBAxes.fromVects(t_vector=t_vector, p_vector=p_vector)

    def __repr__(self):

        return "PTBAxes(P {}, T {})".format(
            Axis.fromVect(self._p_versor),
            Axis.fromVect(self._t_versor))

    @property
    def PVersor(self) -> Vect3D:
        """
        Return the P versor component of the PTBAxes instance.

        :return: P versor instance

        Example:
          >>> PTBAxes(p_axis=Axis(0, 0), t_axis=Axis(90, 0)).PVersor
          Vect3D(-0.0000, 1.0000, 0.0000)
        """

        return self._p_versor

    @property
    def TVersor(self) -> Vect3D:
        """
        Return the T versor component of the PTBAxes instance.

        :return: T versor instance

        Example:
          >>> PTBAxes(p_axis=Axis(0, 0), t_axis=Axis(90, 0)).TVersor
          Vect3D(1.0000, 0.0000, -0.0000)
        """

        return self._t_versor

    @property
    def BVersor(self) -> Vect3D:
        """
        Return the B versor component of the PTBAxes instance.

        :return: B versor instance

        Example:
          >>> PTBAxes(p_axis=Axis(0, 0), t_axis=Axis(90, 0)).BVersor
          Vect3D(0.0000, 0.0000, 1.0000)
        """

        return self.TVersor.cross_product(self.PVersor)

    @property
    def PAxis(self) -> Axis:
        """
        Return the P axis.

        Example:
          >>> PTBAxes(p_axis=Axis(0, 0), t_axis=Axis(90, 0)).PAxis
          Axis(az: 0.00°, pl: -0.00°)
          >>> PTBAxes(p_axis=Axis(0, 90), t_axis=Axis(90, 0)).PAxis
          Axis(az: 0.00°, pl: 90.00°)
        """

        return Axis.fromVect(self.PVersor)

    @property
    def TAxis(self) -> Axis:
        """
        Return the T axis.

        Example:
          >>> PTBAxes(p_axis=Axis(0, 0), t_axis=Axis(90, 0)).TAxis
          Axis(az: 90.00°, pl: 0.00°)
          >>> PTBAxes(p_axis=Axis(0, -90), t_axis=Axis(90, 0)).TAxis
          Axis(az: 90.00°, pl: 0.00°)
        """

        return Axis.fromVect(self.TVersor)

    @property
    def BAxis(self) -> Axis:
        """
        Calculate the B axis.

        Example:
          >>> PTBAxes(p_axis=Axis(0, 0), t_axis=Axis(90, 0)).BAxis
          Axis(az: 0.00°, pl: -90.00°)
          >>> PTBAxes(p_axis=Axis(0, 90), t_axis=Axis(0, 0)).BAxis
          Axis(az: 270.00°, pl: -0.00°)
        """

        return Axis.fromVect(self.BVersor)

    @property
    def MPlane(self) -> Plane:
        """
        Calculate M plane.

        Example:
          >>> PTBAxes(p_axis=Axis(0, 90), t_axis=Axis(90, 0)).MPlane.is_sub_parallel(Plane(0.0, 90.0))
          True
          >>> PTBAxes(p_axis=Axis(45, 45), t_axis=Axis(225, 45)).MPlane.is_sub_parallel(Plane(315.00, 90.00))
          True
        """

        return self.PAxis.common_plane(self.TAxis)

    def almostEqual(self, another: 'PTBAxes', tolerance_angle: numbers.Real=VECTOR_ANGLE_THRESHOLD) -> bool:
        """
        Checks for equivalence between two PTBAXes instances
        intersect a given tolerance angle (default is VECTOR_ANGLE_THRESHOLD)

        :param another: a PTBAXes instance.
        :param tolerance_angle: the tolerance angle for the equality check (numbers.Real)
        :return: Boolean.

        Examples:
          >>> fm1 = PTBAxes(p_axis=Axis(0, 0), t_axis=Axis(90, 0))
          >>> fm2 = PTBAxes(p_axis=Axis(0, 0.5), t_axis=Axis(90, 0))
          >>> fm1.almostEqual(fm2)
          True
          >>> fm3 = PTBAxes(p_axis=Axis(180.5, 0), t_axis=Axis(90.5, 0))
          >>> fm1.almostEqual(fm3)
          True
          >>> fm3.almostEqual(fm2)
          True
          >>> fm4 = PTBAxes(p_axis=Axis(181.5, 0), t_axis=Axis(91.5, 0))
          >>> fm1.almostEqual(fm4)
          False
        """

        if not isinstance(another, PTBAxes):
            raise Exception("Argument must be of PTBAxes type")

        if not self.PAxis.is_sub_parallel(another.PAxis, tolerance_angle):
            return False

        if not self.TAxis.is_sub_parallel(another.TAxis, tolerance_angle):
            return False

        return True

    def toMatrix(self) -> np.ndarray:
        """
        Creates a rotation matrix from the PTB as_vect xyz.
        Formula as in:
        - eq. 3 in Kagan, Y. Y., 2007. Simplified algorithms for calculating double-couple rotation.
        - eq. 31 in Kagan, Y. Y., 2008. On geometric complexity of earthquake focal zone and fault system.

        :return: a 3x3 numpy arrays fo floats.

        Example:
          >>> arraysAreClose(PTBAxes(p_axis=Axis(0, 0), t_axis=Axis(90, 0)).toMatrix(), np.identity(3))
          True
        """

        t = self.TVersor
        p = self.PVersor
        b = self.BVersor

        return np.array([
            [t.x, p.x, b.x],
            [t.y, p.y, b.y],
            [t.z, p.z, b.z]
        ])

    def toQuatern(self) -> Quaternion:
        """
        Transforms the focal mechanism into a quaternion.

        :return: a Quaternion instance.

        Example:
          >>> PTBAxes(p_axis=Axis(232, 41), t_axis=Axis(120, 24)).toQuatern()
          Quaternion(-0.41567, 0.85017, -0.31120, -0.08706)
          >>> PTBAxes(p_axis=Axis(51, 17), t_axis=Axis(295, 55)).toQuatern()
          Quaternion(0.38380, 0.30459, 0.80853, -0.32588)
        """

        return Quaternion.fromRotMatr(self.toMatrix())


def focmech_rotate(
    fm: PTBAxes,
    ra: RotationAxis
) -> PTBAxes:
    """
    Rotate a fochal mechanism (a PTBAxes instance) to a new orientation
    via a rotation axis.

    :param fm: the focal mechanism to rotVectByAxis
    :param ra: the rotation axis
    :return: the rotated focal mechanism
    """

    qfm = fm.toQuatern()
    qra = ra.toRotQuater()

    qrot = qra * qfm

    return PTBAxes.fromQuatern(qrot)


def focmechs_invert_rotations(fm1: PTBAxes, fm2: PTBAxes) -> List[RotationAxis]:
    """
    Calculate the rotations between two focal mechanisms, sensu Kagan.
    See Kagan Y.Y. papers for theoretical basis.
    Practical implementation derive from Alberti, 2010:
    Analysis of kinematic correlations in faults and focal mechanisms with GIS and Fortran programs.

    :param fm1: a PTBAxes instance
    :param fm2: another PTBAxes instance
    :return:a list of 4 rotation axes, sorted by increasing rotation angle
    """

    # processing of equal focal mechanisms

    t_axes_angle = fm1.TAxis.angle_as_degrees(fm2.TAxis)
    p_axes_angle = fm1.PAxis.angle_as_degrees(fm2.PAxis)

    if t_axes_angle < 0.5 and p_axes_angle < 0.5:
        return []

    # transformation of XYZ axes cartesian xyz (fm1,2) into quaternions q1,2

    focmec1_matrix = fm1.toMatrix()
    focmec2_matrix = fm2.toMatrix()

    fm1_quaternion = Quaternion.fromRotMatr(focmec1_matrix)
    fm2_quaternion = Quaternion.fromRotMatr(focmec2_matrix)

    # calculation of quaternion inverse q1,2[-1]

    fm1_inversequatern = fm1_quaternion.inverse
    fm2_inversequatern = fm2_quaternion.inverse

    # calculation of rotation quaternion : q' = q2*q1[-1]

    base_rot_quater = fm2_quaternion * fm1_inversequatern

    # calculation of secondary rotation pure quaternions: a(i,j,k) = q2*(i,j,k)*q2[-1]

    axes_quats = [
        Quaternion.i(),
        Quaternion.j(),
        Quaternion.k()]

    suppl_prod2quat = map(lambda ax_quat: fm2_quaternion * (ax_quat * fm2_inversequatern), axes_quats)

    # calculation of the other 3 rotation quaternions: q'(i,j,k) = a(i,j,k)*q'

    rotations_quaternions = list(map(lambda quat: quat * base_rot_quater, suppl_prod2quat))
    rotations_quaternions.append(base_rot_quater)

    rotation_axes = list(map(lambda quat: RotationAxis.fromQuater(quat).toMinRotAxis(), rotations_quaternions))

    return sortRotations(rotation_axes)


if __name__ == "__main__":

    import doctest
    doctest.testmod()
