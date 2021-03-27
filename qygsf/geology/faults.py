from qygsf import PTBAxes
from qygsf.geometries.shapes.geometry import *
from qygsf.errors import *
from qygsf.mathematics.math_utils import *


class Slick(object):
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
          >>> Slick(GVect(90, 10))
          Slickenline(090.00, +10.00, True)
        """

        assert isinstance(mov_lin, (GVect, GAxis)), "Movement is not of the correct type"
        self._mov_lin = mov_lin

    def hasKnownSense(self):
        """
        Check whether the slickenline has known movement sense.

        Example:
          >>> Slick(GVect(90, 45)).hasKnownSense()
          True
          >>> Slick(GAxis(90, 45)).hasKnownSense()
          False
        """

        if isinstance(self._mov_lin, GAxis):
            return False
        elif isinstance(self._mov_lin, GVect):
            return True
        else:
            raise Exception("Error with provided slickeline type")

    def hasUnknownSense(self):
        """
        Check whether the slickenline has unknown/uncertain movement sense.

        Example:
          >>> Slick(GAxis(90, 45)).hasUnknownSense()
          True
          >>> Slick(GVect(90, 45)).hasUnknownSense()
          False
        """

        return not self.hasKnownSense()

    def setKnownSense(self):
        """
        Set (formal) movement sense to Slickline instance without known/certain
        movemment sense. Raise SlickelineSenseException when already known.

        Example:
          >>> Slick(GVect(180, -30)).setKnownSense()
          Traceback (most recent call last):
          ...
          SlickelineSenseException: Slickenline must have unknown movement sense
          >>> Slick(GAxis(180, -30)).setKnownSense() 
          Slickenline(180.00, -30.00, True)
        """

        if self.hasKnownSense():
            raise Exception("Slickenline must have unknown movement sense")

        return Slick(self.lin.as_gvect())

    def setUnknownSense(self):
        """
        Set to unknown/uncertain the movement sense for the current Slickline instance. 
        Raise SlickelineSenseException when already unknown.

        Example:
          >>> Slick(GAxis(180, -30)).setUnknownSense()
          Traceback (most recent call last):
          ...
          SlickelineSenseException: Slickenline must have known movement sense
          >>> Slick(GVect(180, -30)).setUnknownSense() 
          Slickenline(180.00, -30.00, False)
        """

        if self.hasUnknownSense():
            raise Exception("Slickenline must have known movement sense")

        return Slick(self.lin.as_axis())

    @property
    def lin(self):
        """
        Return the slickenline orientation value,
        as a GVect (known movement sense) or a GAxis instance (unknown movement sense).

        Example:
          >>> Slick(GVect(90, 45)).lin
          GVect(090.00, +45.00)
          >>> Slick(GAxis(90, 45)).lin
          GAxis(090.00, +45.00)
        """

        return self._mov_lin

    @property
    def vals(self):
        """
        The slickenline parameters.
        """

        known_mov = self.hasKnownSense()

        return self._mov_lin.tr, self._mov_lin.pl, known_mov

    def __repr__(self):

        return "Slickenline({:06.2f}, {:+06.2f}, {})".format(*self.vals)

    def invert(self):
        """
        Invert the slickenline sense, when known, otherwise raise SlickelineSenseException.

        Example:
         >>> Slick(GAxis(30, 45)).invert()
         Traceback (most recent call last):
         ...
         SlickelineSenseException: Slickenline must have know movement sense
         >>> Slick(GVect(30, 45)).invert()
         Slickenline(210.00, -45.00, True)
        """

        if not self.hasKnownSense():
            raise Exception("Slickenline must have know movement sense")

        return Slick(self.lin.opposite())


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
          >>> FaultSlick(GPlane(90, 45), Slick(GAxis(90, 45)))
          FaultSlick(GPlane(090.00, +45.00), Slickenline(090.00, +45.00, False))
        """

        assert isinstance(fault_plane, GPlane), "Provided fault plane must be a GPlane instance"
        assert isinstance(slickenline, Slick), "Provided slickenline must be a Slickenline instance"

        assert isclose(fault_plane.normal.angle(slickenline.lin), 90.), "Slickenline is not within fault plane"
        self._fltpln = fault_plane
        self._slick = slickenline

    @property
    def fp(self):
        """
        Return fault plane, as a GPlane instance.

        Example:
          >>> FaultSlick(GPlane(90, 45), Slick(GAxis(90, 45))).fp
          GPlane(090.00, +45.00)
        """

        return self._fltpln

    @property
    def sl(self):
        """
        Return the slickenline associated with the fault. 

        Example:
          >>> FaultSlick(GPlane(90, 45), Slick(GAxis(90, 45))).sl
          Slickenline(090.00, +45.00, False)
        """

        return self._slick

    @property
    def known_sense(self):
        """
        Check if Slickenline in FautlSlick instance has known
        movement sense.

        Example: 
          >>> FaultSlick(GPlane(90, 45), Slick(GAxis(90, 45))).known_sense
          False
          >>> FaultSlick(GPlane(90, 45), Slick(GVect(90, 45))).known_sense
          True
        """

        if self.sl.hasKnownSense():
            return True
        else:
            return False

    def set_known_sense(self):
        """
        Create FaultSlick instance with known movement sense
         from another instance, raising SlickelineSenseException if already known.

        Example: 
          >>> FaultSlick(GPlane(90, 45), Slick(GVect(90, 45))).set_known_sense()
          Traceback (most recent call last):
          ...
          SlickelineSenseException: Fault slickenline sense must be unknown
          >>> FaultSlick(GPlane(0, 45), Slick(GAxis(0, 45))).set_known_sense()
          FaultSlick(GPlane(000.00, +45.00), Slickenline(000.00, +45.00, True))
        """

        if self.known_sense:
            raise Exception("Fault slickenline sense must be unknown")

        return FaultSlick(self.fp, self.sl.setKnownSense())

    def set_unknown_sense(self):
        """
        Create FaultSlick instance with unknown/uncertain movement sense.
        Raise SlickelineSenseException if source instance has already unknown movement sense.

        Example: 
          >>> FaultSlick(GPlane(90, 45), Slick(GAxis(90, 45))).set_unknown_sense()
          Traceback (most recent call last):
          ...
          SlickelineSenseException: Fault slickenline sense must be known
          >>> FaultSlick(GPlane(0, 45), Slick(GVect(0, 45))).set_unknown_sense()
          FaultSlick(GPlane(000.00, +45.00), Slickenline(000.00, +45.00, False))
        """

        if not self.known_sense:
            raise Exception("Fault slickenline sense must be known")

        return FaultSlick(self.fp, self.sl.setUnknownSense())

    def __repr__(self):

        return "FaultSlick({}, {})".format(self.fp, self.sl)

    def opposite_mov(self):
        """
        Create FaultSlick instance with opposite movement, when the source instance
        has defined movement sense, otherwise raise SlickelineSenseException.

        Example:
          >>> FaultSlick(GPlane(90, 45), Slick(GAxis(90, 45))).opposite_mov()
          Traceback (most recent call last):
          ...
          SlickelineSenseException: Fault slickenline must have known movement sense
          >>> FaultSlick(GPlane(90, 45), Slick(GVect(90, 45))).opposite_mov()
          FaultSlick(GPlane(090.00, +45.00), Slickenline(270.00, -45.00, True))
        """

        if not self.known_sense:
            raise Exception("Fault slickenline must have known movement sense")

        return FaultSlick(self.fp, self.sl.invert())

    def PTaxes(self):
        """
        Calculate P-T axes. 
        Return P axis, T axis and a third variable, boolean, 
        indicating if the P-T derivation is from a slickenline
        with a known movement sense (True) or with
        unknown/uncertain movement sense (False).
        
        Example:
          >>> FaultSlick(GPlane(90, 45), Slick(GVect(90, 45))).PTaxes()
          PTBAxes(P: GAxis(000.00, -90.00), T: GAxis(090.00, +00.00), True)
        """

        s_versor = self.sl.lin.versor()
        f_versor = self.fp.normal.versor()
        T_axis = (f_versor + s_versor).gaxis
        P_axis = (f_versor - s_versor).gaxis
        known = self.known_sense

        return PTBAxes(P_axis, T_axis, known)


if __name__ == "__main__":

    import doctest

    doctest.testmod()
