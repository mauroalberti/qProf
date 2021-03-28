from math import radians, cos, sin

from numpy import amin

from qygsf import CPlane3D, Vect
from qygsf.mathematics.defaults import MIN_ANGLE_DEGR_DISORIENTATION


class Plane:
    """
    Geological plane.
    Defined by dip direction and dip angle (both in degrees):
     - dip direction: [0.0, 360.0[ clockwise, from 0 (North);
     - dip angle: [0, 90.0]: downward-pointing.
    """

    def __init__(self, srcAzimuth, srcDipAngle, isRHRStrike=False):
        """
        Geological plane constructor.

        @param  srcAzimuth:  Azimuth of the plane (RHR strike or dip direction).
        @type  srcAzimuth:  number or string convertible to float.
        @param  srcDipAngle:  Dip angle of the plane (0-90°).
        @type  srcDipAngle:  number or string convertible to float.

        @return:  GPlane.

        Example:
          >>> Plane(0, 90)
          GPlane(000.00, +90.00)
          >>> Plane(0, 90, isRHRStrike=True)
          GPlane(090.00, +90.00)
          >>> Plane(0, 90, True)
          GPlane(090.00, +90.00)
        """

        def rhrstrk2dd(rhr_strk):
            """Converts RHR strike value to dip direction value.

            Example:
                >>> rhrstrk2dd(285.5)
                15.5
            """

            return (rhr_strk + 90.0) % 360.0

        if isRHRStrike:
            self._dipdir = rhrstrk2dd(srcAzimuth)
        else:
            self._dipdir = srcAzimuth % 360.0
        self._dipangle = float(srcDipAngle)

    @property
    def dd(self):
        """
        Return the dip direction of the geological plane.

        Example:
          >>> Plane(34.2, 89.7).dd
          34.2
        """

        return self._dipdir

    @property
    def da(self):
        """
        Return the dip angle of the geological plane.

        Example:
          >>> Plane(183, 77).da
          77.0

        """

        return self._dipangle

    @property
    def dda(self):
        """
        Return a tuple storing the dip direction and dip angle values of a geological plane.

        Example:
          >>> gp = Plane(89.4, 17.2)
          >>> gp.dda
          (89.4, 17.2)
        """

        return self.dd, self.da

    @property
    def strike_rhr(self):
        """
        Return the strike according to the right-hand-rule.

        Examples:
          >>> Plane(90, 45).strike_rhr
          0.0
          >>> Plane(45, 89).strike_rhr
          315.0
          >>> Plane(275, 38).strike_rhr
          185.0
          >>> Plane(0, 38).strike_rhr
          270.0
        """

        return (self.dd - 90.0) % 360.0

    @property
    def strike_lhr(self):
        """
        Return the strike according to the left-hand-rule.

        Examples:
          >>> Plane(90, 45).strike_lhr
          180.0
          >>> Plane(45, 89).strike_lhr
          135.0
          >>> Plane(275, 38).strike_lhr
          5.0
          >>> Plane(0, 38).strike_lhr
          90.0
        """

        return (self.dd + 90.0) % 360.0

    def __repr__(self):

        return "GPlane({:06.2f}, {:+06.2f})".format(*self.dda)

    @property
    def normal(self):
        """
        Return the geological vector normal to the geological plane,
        pointing in the same direction as the geological plane.

        Example:
            >>> Plane(90, 55).normal
            GVect(090.00, -35.00)
        """

        trend = self.dd % 360.0
        plunge = self.da - 90.0

        return Direct(trend, plunge)

    def plane(self, point):
        """
        Given a GPlane instance and a provided Point instance,
        calculate the corresponding Plane instance.

        Example:
          >>> Plane(0, 0).plane(Point4D(0, 0, 0))
          Plane(0.0000, 0.0000, 1.0000, -0.0000)
          >>> Plane(90, 45).plane(Point4D(0, 0, 0))
          Plane(0.7071, 0.0000, 0.7071, -0.0000)
          >>> Plane(0, 90).plane(Point4D(0, 0, 0))
          Plane(0.0000, 1.0000, -0.0000, -0.0000)
        """

        normal_versor = self.normal.versor()
        a, b, c = normal_versor.x, normal_versor.y, normal_versor.z
        d = - (a * point.x + b * point.y + c * point.z)
        return CPlane3D(a, b, c, d)

    def angle(self, another):
        """
        Calculate angle (in degrees) between two geoplanes.

        >>> p1 = Plane(100.0, 50.0)
        >>> p1.angle(p1)
        0.0

        >>> p2 = Plane(300.0, 10.0)
        >>> p3 = Plane(300.0, 90.0)
        >>> p2.angle(p3)
        80.0
        """

        return self.normal.as_axis().angle(another.normal.as_axis())

    def rake_to_gv(self, rake):
        """
        Calculate GVect given a GPlane instance and a rake value.
        The rake is defined according to the Aki and Richards, 1980 conventions:
        rake = 0° -> left-lateral
        rake = 90° -> reverse
        rake = +/- 180° -> right-lateral
        rake = -90° -> normal

        Examples:
          >>> Plane(180, 45).rake_to_gv(0.0)
          GVect(090.00, +00.00)
          >>> Plane(180, 45).rake_to_gv(90.0)
          GVect(000.00, -45.00)
          >>> Plane(180, 45).rake_to_gv(-90.0)
          GVect(180.00, +45.00)
          >>> Plane(180, 45).rake_to_gv(180.0)
          GVect(270.00, -00.00)
          >>> Plane(180, 45).rake_to_gv(-180.0)
          GVect(270.00, +00.00)
        """

        rk = radians(rake)
        strk = radians(self.strike_rhr)
        dip = radians(self.da)

        x = cos(rk)*sin(strk)-sin(rk)*cos(dip)*cos(strk)
        y = cos(rk)*cos(strk)+sin(rk)*cos(dip)*sin(strk)
        z = sin(rk) * sin(dip)

        return Vect(x, y, z).gvect


class Direct:
    """
    Geological vector.
    Defined by trend and plunge (both in degrees):
     - trend: [0.0, 360.0[ clockwise, from 0 (North):
     - plunge: [-90.0, 90.0].
    """

    def __init__(self, srcTrend, srcPlunge):
        """
        Geological vector constructor.
        srcTrend: Trend range: [0.0, 360.0[ clockwise, from 0 (North)
        srcPlunge: Plunge: [-90.0, 90.0],
        negative value: upward pointing axis, positive values: downward axis;

        Example:
          >>> a = Direct(120, -27)
          >>> b = Direct(54, -320)
          Traceback (most recent call last):
          ...
          AssertionError: plunge must be between -90° and +90° (comprised)
        """

        assert -90.0 <= float(srcPlunge) <= 90.0, "plunge must be between -90° and +90° (comprised)"
        self._trend = float(srcTrend) % 360.0
        self._plunge = float(srcPlunge)

    @property
    def tr(self):
        """
        Return trend of the geological direction.
        Range is [0, 360[

        Example:
          >>> Direct(420, -17).tr
          60.0
          >>> Direct(-20, 49).tr
          340.0
        """

        return self._trend

    @property
    def pl(self):
        """
        Return plunge of the geological direction.
        Range is [-90, 90]

        Example:
          >>> Direct(420, -17).pl
          -17.0
        """

        return self._plunge

    @property
    def tp(self):
        """
        Return trend and plunge of the geological direction.

        Example:
          >>> Direct(-90, -45).tp
          (270.0, -45.0)
        """

        return self.tr, self.pl

    def __repr__(self):

        return "GVect({:06.2f}, {:+06.2f})".format(*self.tp)

    def copy(self):
        """
        Return a copy of the GVect instance.

        Example:
          >>> Direct(10, 20).copy()
          GVect(010.00, +20.00)
        """

        return Direct(*self.tp)

    def opposite(self):
        """
        Return the opposite GVect.

        Example:
          >>> Direct(0, 30).opposite()
          GVect(180.00, -30.00)
          >>> Direct(315, 10).opposite()
          GVect(135.00, -10.00)
        """

        trend = (self.tr + 180.) % 360.
        plunge = -self.pl

        return Direct(trend, plunge)

    def versor(self):
        """
        Return the Vect corresponding to the geological vector.

        Examples:
          >>> Direct(0, 90).versor()
          Vect(0.0000, 0.0000, -1.0000)
          >>> Direct(0, -90).versor()
          Vect(0.0000, 0.0000, 1.0000)
          >>> Direct(90, 90).versor()
          Vect(0.0000, 0.0000, -1.0000)

        """

        north_coord = cos(radians(self.pl)) * cos(radians(self.tr))
        east_coord = cos(radians(self.pl)) * sin(radians(self.tr))
        down_coord = sin(radians(self.pl))

        return Vect(east_coord, north_coord, -down_coord)

    @property
    def is_upward(self):
        """
        Check whether the instance is pointing upward or horizontal.

        Examples:
          >>> Direct(10, 15).is_upward
          False
          >>> Direct(257.4, 0.0).is_upward
          False
          >>> Direct(90, -45).is_upward
          True
        """

        if self.versor().z > 0.0:
            return True
        else:
            return False

    @property
    def is_downward(self):
        """
        Check whether the instance is pointing downward or horizontal.

        Examples:
          >>> Direct(10, 15).is_downward
          True
          >>> Direct(257.4, 0.0).is_downward
          False
          >>> Direct(90, -45).is_downward
          False
        """

        if self.versor().z < 0.0:
            return True
        else:
            return False

    @property
    def upward(self):
        """
        Return upward-point geological vector.

        Examples:
          >>> Direct(90, -45).upward
          GVect(090.00, -45.00)
          >>> Direct(180, 45).upward
          GVect(000.00, -45.00)
          >>> Direct(0, 0).upward
          GVect(180.00, -00.00)
          >>> Direct(0, 90).upward
          GVect(180.00, -90.00)
        """

        if self.is_upward:
            return self.copy()
        else:
            return self.opposite()

    @property
    def downward(self):
        """
        Return downward-pointing geological vector.

        Examples:
          >>> Direct(90, -45).downward
          GVect(270.00, +45.00)
          >>> Direct(180, 45).downward
          GVect(180.00, +45.00)
          >>> Direct(0, 0).downward
          GVect(180.00, -00.00)
          >>> Direct(0, 90).downward
          GVect(000.00, +90.00)
        """

        if self.is_downward:
            return self.copy()
        else:
            return self.opposite()

    def angle(self, another):
        """
        Calculate angle (in degrees) between the two
        GVect instances.

        Examples:
          >>> Direct(0, 90).angle(Direct(90, 0)) # doctest: +NUMBER
          90.0000000
          >>> Direct(0, 0).angle(Direct(270, 0)) # doctest: +NUMBER
          90.0000000
          >>> Direct(0, 0).angle(Direct(0, 0)) # doctest: +NUMBER
          0.0000000
          >>> Direct(0, 0).angle(Direct(180, 0)) # doctest: +NUMBER
          180.0000000
        """

        return self.versor().angle(another.versor())

    @property
    def normal_gplane(self):
        """
        Return the geological plane that is normal to the geological vector.

        Examples:
          >>> Direct(0, 45).normal_gplane
          GPlane(180.00, +45.00)
          >>> Direct(0, -45).normal_gplane
          GPlane(000.00, +45.00)
          >>> Direct(0, 90).normal_gplane
          GPlane(180.00, +00.00)
        """

        down_axis = self.downward
        dipdir = (down_axis.tr + 180.0) % 360.0
        dipangle = 90.0 - down_axis.pl

        return Plane(dipdir, dipangle)

    def common_plane(self, another):
        """
        Calculate GPlane instance defined by the two GVect instances.

        Examples:
          >>> Direct(0, 0).common_plane(Direct(90, 0))
          GPlane(180.00, +00.00)
          >>> Direct(0, 0).common_plane(Direct(90, 90))
          GPlane(090.00, +90.00)
          >>> Direct(45, 0).common_plane(Direct(135, 45))
          GPlane(135.00, +45.00)
          >>> Direct(315, 45).common_plane(Direct(135, 45))
          GPlane(225.00, +90.00)
        """

        normal = self.versor().vp(another.versor())
        return normal.gvect.normal_gplane

    def as_axis(self):
        """
        Create GAxis instance with the same attitude as the self instance.

        Example:
          >>> Direct(220, 32).as_axis()
          GAxis(220.00, +32.00)
        """

        return Axis(*self.tp)

    def vp(self, another):
        """
        Calculate the GVect instance that is normal to the two provided sources.
        Angle between sources must be larger than MIN_ANGLE_DEGR_DISORIENTATION,
        otherwise a SubparallelLineationException will be raised.

        Example:
          >>> Direct(0, 0).vp(Direct(4, 0))
          Traceback (most recent call last):
          ...
          SubparallelLineationException: Sources must not be sub- or anti-parallel
          >>> Direct(0, 0).vp(Direct(179, 0))
          Traceback (most recent call last):
          ...
          SubparallelLineationException: Sources must not be sub- or anti-parallel
          >>> Direct(0, 0).vp(Direct(5.1, 0))
          GVect(000.00, +90.00)
          >>> Direct(90, 45).vp(Direct(90, 0))
          GVect(180.00, +00.00)
        """

        if not MIN_ANGLE_DEGR_DISORIENTATION <= self.angle(another) <= 180. - MIN_ANGLE_DEGR_DISORIENTATION:
            raise Exception("Sources must not be sub- or anti-parallel")

        return self.versor().vp(another.versor()).gvect


class Axis(Direct):

    def __init__(self, srcTrend, srcPlunge):

        super(Axis, self).__init__(srcTrend, srcPlunge)

    def __repr__(self):

        return "GAxis({:06.2f}, {:+06.2f})".format(*self.tp)

    def as_gvect(self):
        """
        Create GVect instance with the same attitude as the self instance.

        Example:
          >>> Axis(220, 32).as_gvect()
          GVect(220.00, +32.00)
        """

        return Direct(*self.tp)

    def versor(self):
        """
        Create a versor parallel to the geological axis.

        :return:
        """

        return self.as_gvect().versor()

    def angle(self, another):
        """
        Calculate angle (in degrees) between the two
        GAxis instances.

        Examples:
          >>> Axis(0, 90).angle(Axis(90, 0)) # doctest: +NUMBER
          90.0000000
          >>> Axis(0, 0).angle(Axis(270, 0)) # doctest: +NUMBER
          90.0000000
          >>> Axis(0, 0).angle(Axis(0, 0)) # doctest: +NUMBER
          0.0000000
          >>> Axis(0, 0).angle(Axis(180, 0)) # doctest: +NUMBER
          0.0000000
        """

        angle_vers = self.versor().angle(another.versor())
        return min(angle_vers, 180. - angle_vers)

    @property
    def normal_gplane(self):
        """
        Calculate the geological plane that is normal to the given GAxis instance.

        Example:
          >>> Axis(0, 90).normal_gplane
          GPlane(180.00, +00.00)
          >>> Axis(0, 0).normal_gplane
          GPlane(000.00, +90.00)
          >>> Axis(45, 45).normal_gplane
          GPlane(225.00, +45.00)
        """

        return self.as_gvect().normal_gplane

    @property
    def is_upward(self):
        """
        Check whether the instance is pointing upward or horizontal.

        Examples:
          >>> Axis(10, 15).is_upward
          False
          >>> Axis(257.4, 0.0).is_upward
          False
          >>> Axis(90, -45).is_upward
          True
        """

        return self.as_gvect().is_upward

    @property
    def is_downward(self):
        """
        Check whether the instance is pointing downward or horizontal.

        Examples:
          >>> Axis(10, 15).is_downward
          True
          >>> Axis(257.4, 0.0).is_downward
          False
          >>> Axis(90, -45).is_downward
          False
        """

        return self.as_gvect().is_downward

    @property
    def upward(self):
        """
        Return upward-point geological axis.

        Examples:
          >>> Axis(90, -45).upward
          GAxis(090.00, -45.00)
          >>> Axis(180, 45).upward
          GAxis(000.00, -45.00)
          >>> Axis(0, 0).upward
          GAxis(180.00, -00.00)
          >>> Axis(0, 90).upward
          GAxis(180.00, -90.00)
        """

        return self.as_gvect().upward.as_axis()

    @property
    def downward(self):
        """
        Return downward-pointing geological vector.

        Examples:
          >>> Axis(90, -45).downward
          GAxis(270.00, +45.00)
          >>> Axis(180, 45).downward
          GAxis(180.00, +45.00)
          >>> Axis(0, 0).downward
          GAxis(180.00, -00.00)
          >>> Axis(0, 90).downward
          GAxis(000.00, +90.00)
        """

        return self.as_gvect().downward.as_axis()

    def common_plane(self, another):
        """
        Calculate GPlane instance defined by the two GAxis instances.

        Examples:
          >>> Axis(0, 0).common_plane(Axis(90, 0))
          GPlane(180.00, +00.00)
          >>> Axis(0, 0).common_plane(Axis(90, 90))
          GPlane(090.00, +90.00)
          >>> Axis(45, 0).common_plane(Axis(135, 45))
          GPlane(135.00, +45.00)
          >>> Axis(315, 45).common_plane(Axis(135, 45))
          GPlane(225.00, +90.00)
        """

        return self.as_gvect().common_plane(another.as_gvect())

    def vp(self, another):
        """
        Calculate the GAxis instance that is perpendicular to the two provided.
        The two source GAxis must not be subparallel (threshold is MIN_ANGLE_DEGR_DISORIENTATION),
        otherwise a SubparallelLineationException will be raised.

        Example:
          >>> Axis(0, 0).vp(Axis(4, 0))
          Traceback (most recent call last):
          ...
          SubparallelLineationException: Sources must not be sub- or anti-parallel
          >>> Axis(0, 0).vp(Axis(180, 0))
          Traceback (most recent call last):
          ...
          SubparallelLineationException: Sources must not be sub- or anti-parallel
          >>> Axis(90, 0).vp(Axis(180, 0))
          GAxis(000.00, +90.00)
          >>> Axis(90, 45).vp(Axis(180, 0))
          GAxis(270.00, +45.00)
          >>> Axis(270, 45).vp(Axis(180, 90))
          GAxis(180.00, -00.00)
        """

        return self.as_gvect().vp(another.as_gvect()).as_axis()