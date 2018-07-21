

from builtins import object
class RasterParametersException(Exception):
    """
    Exception for raster parameters.
    """
    pass


class VectorInputException(Exception):
    """
    Exception for vector input parameters.
    """
    pass


class FunInputException(Exception):
    """
    Exception for function input errors.
    """
    pass


class OutputException(Exception):
    """
    Exception for output errors.
    """
    pass


class ConnectionException(object):
    pass


class AnaliticSurfaceIOException(Exception):
    pass


class AnaliticSurfaceCalcException(Exception):
    pass


class GPXIOException(Exception):
    pass


class VectorIOException(Exception):

    pass


class OGRIOException(Exception):
    """
    Exception for raster parameters.
    """
    pass
