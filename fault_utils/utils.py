
from .errors import RakeInputException


def rake_to_apsg_movsense(rake):

    if rake > 0.0 and rake < 180.0:  # reverse faults according to Aki & Richards, 1980 convention
        return 1
    elif rake < 0.0 and rake > -180.0:  # normal faults according to Aki & Richards, 1980 convention
        return -1
    elif abs(rake) == 0.0 or abs(rake) == 180.0:
        raise RakeInputException("Currently transcurrent data (rake = +/-180 or = 0.0) are not handled in plots")
    else:
        raise RakeInputException("Input rake value not acceptable")

def movsense_to_apsg_movsense(str_val):

    if str_val == "R":
        return 1
    elif str_val == "N":
        return -1
    else:
        raise RakeInputException("Input rake value not acceptable")

