
def qcolor2rgbmpl(
        qcolor
):

    red = qcolor.red() / 255.0
    green = qcolor.green() / 255.0
    blue = qcolor.blue() / 255.0
    return red, green, blue

