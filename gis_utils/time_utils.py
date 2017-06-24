
import calendar


def standard_gpstime_to_seconds(time_str):
    """
    Naive and experimental implementation of a "normal" date-time string from GPX
    into a float value that should represent seconds since January 1, 1970.
    Use just as a first , it-could-be-correct, tool.

    Example:
      >>> standard_gpstime_to_seconds("1970-01-01T00:00:00Z")
      0
      >>> standard_gpstime_to_seconds("1970-01-01T01:0:00Z")
      3600
    """

    date, hhmmss = time_str.split("T")
    if hhmmss.endswith("Z"):
        hhmmss = hhmmss[:-1]

    year, month, day = map(int, date.split("-"))
    hour, minutes, seconds = hhmmss.split(":")
    hour, minutes, seconds = int(hour), int(minutes), float(seconds)

    # modified from:
    # https://stackoverflow.com/questions/7852855/how-to-convert-a-python-datetime-object-to-seconds

    t = year, month, day, hour, minutes, seconds
    secs = calendar.timegm(t)

    return secs


if __name__ == "__main__":

    import doctest
    import numtest  # external module, used in doctest float checks
    doctest.testmod()