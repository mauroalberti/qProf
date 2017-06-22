# -*- coding: utf-8 -*-


from .geometry import *
from .errors import SlickelineTypeException, SlickelineSenseException
from .math_utils import isclose

class Slickenline(object):
    """
    Slickeline.
    It can be defined through a GVect instance, in which case
    it has a movement sense, or via a GAxis, when the movement sense 
    is unknown or not sure.
    When the movement sense is known, the GVect instance indicates
    the displacement of the block that is:
    for a horizontal or a dipping, non vertical fault: the upper block
    for a vertical fault: the block individuated by the (formal) dip direction. 
    """

    def __init__(self, mov_lin):
        """"
        Slickenline constructor.
        The 'mov_lin' argument is a GVect or a GAxis instance. 
        Depending on that, the movement sense will be known
        when a GVect provided, and unknown/uncertain when a GAxis provided.

        Example:
          >>> Slickenline(GVect(90, 10))
          Slickenline(090.00, +10.00, True)
        """

        assert isinstance(mov_lin, (GVect, GAxis)), "Movement is not of the correct type"
        self._mov_lin = mov_lin

    def has_known_sense(self):
        """
        Check whether the slickenline has known movement sense.

        Example:
          >>> Slickenline(GVect(90, 45)).has_known_sense()
          True
          >>> Slickenline(GAxis(90, 45)).has_known_sense()
          False
        """

        if isinstance(self._mov_lin, GAxis):
            return False
        elif isinstance(self._mov_lin, GVect):
            return True
        else:
            raise SlickelineTypeException("Error with provided slickeline type")

    def has_unknown_sense(self):
        """
        Check whether the slickenline has unknown/uncertain movement sense.

        Example:
          >>> Slickenline(GAxis(90, 45)).has_unknown_sense()
          True
          >>> Slickenline(GVect(90, 45)).has_unknown_sense()
          False
        """

        return not self.has_known_sense()

    def set_known_sense(self):
        """
        Set (formal) movement sense to Slickline instance without known/certain
        movemment sense. Raise SlickelineSenseException when already known.

        Example:
          >>> Slickenline(GVect(180, -30)).set_known_sense()
          Traceback (most recent call last):
          ...
          SlickelineSenseException: Slickenline must have unknown movement sense
          >>> Slickenline(GAxis(180, -30)).set_known_sense() 
          Slickenline(180.00, -30.00, True)
        """

        if self.has_known_sense():
            raise SlickelineSenseException("Slickenline must have unknown movement sense")

        return Slickenline(self.lin.as_vect())

    def set_unknown_sense(self):
        """
        Set to unknown/uncertain the movement sense for the current Slickline instance. 
        Raise SlickelineSenseException when already unknown.

        Example:
          >>> Slickenline(GAxis(180, -30)).set_unknown_sense()
          Traceback (most recent call last):
          ...
          SlickelineSenseException: Slickenline must have known movement sense
          >>> Slickenline(GVect(180, -30)).set_unknown_sense() 
          Slickenline(180.00, -30.00, False)
        """

        if self.has_unknown_sense():
            raise SlickelineSenseException("Slickenline must have known movement sense")

        return Slickenline(self.lin.as_axis())

    @property
    def lin(self):
        """
        Return the slickenline orientation value,
        as a GVect (known movement sense) or a GAxis instance (unknown movement sense).

        Example:
          >>> Slickenline(GVect(90, 45)).lin
          GVect(090.00, +45.00)
          >>> Slickenline(GAxis(90, 45)).lin
          GAxis(090.00, +45.00)
        """

        return self._mov_lin

    @property
    def vals(self):
        """
        The slickenline parameters.
        """

        known_mov = self.has_known_sense()

        return self._mov_lin.tr, self._mov_lin.pl, known_mov

    def __repr__(self):

        return "Slickenline({:06.2f}, {:+06.2f}, {})".format(*self.vals)

    def invert(self):
        """
        Invert the slickenline sense, when known, otherwise raise SlickelineSenseException.

        Example:
         >>> Slickenline(GAxis(30, 45)).invert()
         Traceback (most recent call last):
         ...
         SlickelineSenseException: Slickenline must have know movement sense
         >>> Slickenline(GVect(30, 45)).invert()
         Slickenline(210.00, -45.00, True)
        """

        if not self.has_known_sense():
            raise SlickelineSenseException("Slickenline must have know movement sense")

        return Slickenline(self.lin.opposite())


class FaultSlick(object):
    """
    Represent a couple of geological observations made up by a fault plane, 
    represented by a GPlane instance, and a slickenline observation, 
    represented by a Slickenline instance.
    """

    def __init__(self, fault_plane, slickenline):
        """
        Create an instance of a FaultSlick.

        Example:
          >>> FaultSlick(GPlane(90, 45), Slickenline(GAxis(90, 45)))
          FaultSlick(GPlane(090.00, +45.00), Slickenline(090.00, +45.00, False))
        """

        assert isinstance(fault_plane, GPlane), "Provided fault plane must be a GPlane instance"
        assert isinstance(slickenline, Slickenline), "Provided slickenline must be a Slickenline instance"

        assert isclose(fault_plane.normal.angle(slickenline.lin), 90.), "Slickenline is not within fault plane"
        self._fltpln = fault_plane
        self._slick = slickenline

    @property
    def fp(self):
        """
        Return fault plane, as a GPlane instance.

        Example:
          >>> FaultSlick(GPlane(90, 45), Slickenline(GAxis(90, 45))).fp
          GPlane(090.00, +45.00)
        """

        return self._fltpln

    @property
    def sl(self):
        """
        Return the slickenline associated with the fault. 

        Example:
          >>> FaultSlick(GPlane(90, 45), Slickenline(GAxis(90, 45))).sl
          Slickenline(090.00, +45.00, False)
        """

        return self._slick

    @property
    def known_sense(self):
        """
        Check if Slickenline in FautlSlick instance has known
        movement sense.

        Example: 
          >>> FaultSlick(GPlane(90, 45), Slickenline(GAxis(90, 45))).known_sense
          False
          >>> FaultSlick(GPlane(90, 45), Slickenline(GVect(90, 45))).known_sense
          True
        """

        if self.sl.has_known_sense():
            return True
        else:
            return False

    def set_known_sense(self):
        """
        Create FaultSlick instance with known movement sense
         from another instance, raising SlickelineSenseException if already known.

        Example: 
          >>> FaultSlick(GPlane(90, 45), Slickenline(GVect(90, 45))).set_known_sense()
          Traceback (most recent call last):
          ...
          SlickelineSenseException: Fault slickenline sense must be unknown
          >>> FaultSlick(GPlane(0, 45), Slickenline(GAxis(0, 45))).set_known_sense()
          FaultSlick(GPlane(000.00, +45.00), Slickenline(000.00, +45.00, True))
        """

        if self.known_sense:
            raise SlickelineSenseException("Fault slickenline sense must be unknown")

        return FaultSlick(self.fp, self.sl.set_known_sense())

    def set_unknown_sense(self):
        """
        Create FaultSlick instance with unknown/uncertain movement sense.
        Raise SlickelineSenseException if source instance has already unknown movement sense.

        Example: 
          >>> FaultSlick(GPlane(90, 45), Slickenline(GAxis(90, 45))).set_unknown_sense()
          Traceback (most recent call last):
          ...
          SlickelineSenseException: Fault slickenline sense must be known
          >>> FaultSlick(GPlane(0, 45), Slickenline(GVect(0, 45))).set_unknown_sense()
          FaultSlick(GPlane(000.00, +45.00), Slickenline(000.00, +45.00, False))
        """

        if not self.known_sense:
            raise SlickelineSenseException("Fault slickenline sense must be known")

        return FaultSlick(self.fp, self.sl.set_unknown_sense())

    def __repr__(self):

        return "FaultSlick({}, {})".format(self.fp, self.sl)

    def opposite_mov(self):
        """
        Create FaultSlick instance with opposite movement, when the source instance
        has defined movement sense, otherwise raise SlickelineSenseException.

        Example:
          >>> FaultSlick(GPlane(90, 45), Slickenline(GAxis(90, 45))).opposite_mov()
          Traceback (most recent call last):
          ...
          SlickelineSenseException: Fault slickenline must have known movement sense
          >>> FaultSlick(GPlane(90, 45), Slickenline(GVect(90, 45))).opposite_mov()
          FaultSlick(GPlane(090.00, +45.00), Slickenline(270.00, -45.00, True))
        """

        if not self.known_sense:
            raise SlickelineSenseException("Fault slickenline must have known movement sense")

        return FaultSlick(self.fp, self.sl.invert())

    def PTaxes(self):
        """
        Calculate P-T axes. 
        Return P axis, T axis and a third variable, boolean, 
        indicating if the P-T derivation is from a slickenline
        with a known movement sense (True) or with
        unknown/uncertain movement sense (False).
        
        Example:
          >>> FaultSlick(GPlane(90, 45), Slickenline(GVect(90, 45))).PTaxes()
          PTBAxes(P: GAxis(000.00, -90.00), T: GAxis(090.00, +00.00), True)
        """

        s_versor = self.sl.lin.versor_full
        f_versor = self.fp.normal.versor_full
        T_axis = (f_versor + s_versor).gaxis
        P_axis = (f_versor - s_versor).gaxis
        known = self.known_sense

        return PTBAxes(P_axis, T_axis, known)


class PTBAxes(object):
    """
    Represent the triad of P, T and B kinematic axes.
    It can also calculate the M plane.
    """

    def __init__(self, p_axis, t_axis, known=True):
        """
        Create a new PTBAxes instances, given the two
        P and T axes.
        
        Example:
          >>> PTBAxes(GVect(0, 0), GVect(90, 0))
          Traceback (most recent call last):
          ...
          AssertionError: P axis must be an instance of GAxis
          >>> PTBAxes(GAxis(0, 0), GAxis(80, 0))
          Traceback (most recent call last):
          ...
          AssertionError: P and T axes must be perpendicular
          >>> PTBAxes(GAxis(0, 0), GAxis(90, 0))
          PTBAxes(P: GAxis(000.00, +00.00), T: GAxis(090.00, +00.00), True)
        """

        assert isinstance(p_axis, GAxis), "P axis must be an instance of GAxis"
        assert isinstance(t_axis, GAxis), "T axis must be an instance of GAxis"
        assert isclose(p_axis.angle(t_axis), 90.), "P and T axes must be perpendicular"

        self._pax = p_axis
        self._tax = t_axis
        self._known = known

    @property
    def Paxis(self):
        """
        Return the P axis.
        
        Example:
          >>> PTBAxes(GAxis(0, 0), GAxis(90, 0)).Paxis
          GAxis(000.00, +00.00)
        """

        return self._pax

    @property
    def Taxis(self):
        """
        Return the T axis.

        Example:
          >>> PTBAxes(GAxis(0, 0), GAxis(90, 0)).Taxis
          GAxis(090.00, +00.00)
        """

        return self._tax

    @property
    def known(self):
        """
        Indicate if the movement sense is known.

        Example:
          >>> PTBAxes(GAxis(0, 0), GAxis(90, 0)).known
          True
        """

        return self._known

    def __repr__(self):

        return "PTBAxes(P: {}, T: {}, {})".format(self.Paxis, self.Taxis, self.known)

    @property
    def Baxis(self):
        """
        Calculate the B axis.
        
        Example:
          >>> PTBAxes(GAxis(0, 0), GAxis(90, 0)).Baxis
          GAxis(000.00, +90.00)
        """

        return self.Paxis.vp(self.Taxis).as_axis()

    @property
    def Mplane(self):
        """
        Calculate M plane.
        
        Example:
          >>> PTBAxes(GAxis(0, 90), GAxis(90, 0)).Mplane
          GPlane(000.00, +90.00)
          >>> PTBAxes(GAxis(45, 45), GAxis(225, 45)).Mplane
          GPlane(315.00, +90.00)
        """

        return self.Paxis.common_plane(self.Taxis)

if __name__ == "__main__":

    import doctest
    import numtest  # external module, used in doctest float checks
    doctest.testmod()
