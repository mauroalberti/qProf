
import matplotlib.pyplot as plt

import mplstereonet as ms

from .type_checks import *
from ...utils.dicts import *


default_gvect_marker_upward_symbol = "x"
default_gvect_marker_downward_symbol = "o"
default_gvect_marker_color = "blue"

default_gaxis_marker_upward_symbol = "+"
default_gaxis_marker_downward_symbol = "s"
default_gaxis_marker_color = "orange"

default_gplane_downward_linestyle = "solid"
default_gplane_upward_linestyle = "dashed"
default_gplane_line_color = "blue"


def plot(data, force=''):
    """
    Plot geological data with matplotlib and mplstereonet
    """

    def params_gvect(gvect, kwargs, force_emisphere):

        if force_emisphere == 'lower':
            default_marker = default_gvect_marker_downward_symbol
            if gvect.is_upward:
                plot_orien = gvect.opposite()
            else:
                plot_orien = gvect
        elif force_emisphere == 'upper':
            default_marker = default_gvect_marker_upward_symbol
            if gvect.is_downward:
                plot_orien = gvect.opposite()
            else:
                plot_orien = gvect
        elif not force_emisphere:
            plot_orien = gvect
            default_marker = default_gvect_marker_downward_symbol if not plot_orien.is_upward else default_gvect_marker_upward_symbol
        else:
            raise PlotException("Invalid force emisphere parameter")

        if plot_orien.is_upward:  # apparently mplstereonet does not handle negative plunges
            plot_orien = plot_orien.mirrorHoriz()

        bearing, plunge = plot_orien.d
        symbol = kwargs.get("m", default_marker)
        color = kwargs.get("c", default_gvect_marker_color)

        return plunge, bearing, symbol, color

    def params_gaxis(gaxis, kwargs, force_emisphere):

        if (not force_emisphere) or (force_emisphere == 'lower'):
            default_marker = default_gaxis_marker_downward_symbol
            if gaxis.is_upward:
                plot_gaxis = gaxis.opposite()
            else:
                plot_gaxis = gaxis
        elif force_emisphere == 'upper':
            default_marker = default_gaxis_marker_upward_symbol
            if gaxis.is_downward:
                plot_gaxis = gaxis.opposite()
            else:
                plot_gaxis = gaxis
        else:
            raise PlotException("Invalid force emisphere parameter")

        if plot_gaxis.is_upward:  # apparently mplstereonet does not handle negative plunges
            plot_gaxis = plot_gaxis.mirrorHoriz()

        bearing, plunge = plot_gaxis.d
        symbol = kwargs.get("m", default_marker)
        color = kwargs.get("c", default_gaxis_marker_color)

        return plunge, bearing, symbol, color

    def params_gplane(gplane, kwargs, force_emisphere):

        if (not force_emisphere) or (force_emisphere == 'lower'):
            default_line_style = default_gplane_downward_linestyle
            plot_gplane = gplane
        elif force_emisphere == 'upper':
            default_line_style = default_gplane_upward_linestyle
            plot_gplane = gplane.mirrorVertPPlane()
        else:
            raise PlotException("Invalid force emisphere parameter")

        strike, dip = plot_gplane.srda

        line_style = kwargs.get("m", default_line_style)
        color = kwargs.get("c", default_gplane_line_color)

        return strike, dip, line_style, color

    if not isinstance(data, list):
        data = [data]

    if force not in ('', 'upper', 'lower'):
        raise PlotException("Force parameter not isValid")

    fig, ax = ms.subplots()

    for rec in data:

        if isinstance(rec, tuple):
            if isinstance(rec[-1], str):
                params = rec[-1]
                objs = rec[:-1]
            else:
                params = None
                objs = rec
        else:
            objs = [rec]
            params = None

        if params:
            kwargs = string2dict(params)
        else:
            kwargs = dict()

        for obj in objs:
            if isDirect(obj):
                plunge, bearing, symbol, color = params_gvect(obj, kwargs, force_emisphere=force)
                ax.line(
                    plunge,
                    bearing,
                    marker=symbol,
                    markerfacecolor=color,
                    markeredgecolor=color)
            elif isAxis(obj):
                plunge, bearing, symbol, color = params_gaxis(obj, kwargs, force_emisphere=force)
                ax.line(
                    plunge,
                    bearing,
                    marker=symbol,
                    markerfacecolor=color,
                    markeredgecolor=color)
            elif isPlane(obj):
                strike, dip, linestyle, color = params_gplane(obj, kwargs, force_emisphere=force)
                ax.plane(
                    strike,
                    dip,
                    linestyle=linestyle,
                    color=color)

    ax.grid()

    return fig, ax


def splot(data, force=''):
    """
    Plot geological data with matplotlib and mplstereonet
    """

    def params_gvect(gvect, kwargs, force_emisphere):

        if force_emisphere == 'lower':
            default_marker = default_gvect_marker_downward_symbol
            if gvect.is_upward:
                plot_orien = gvect.opposite()
            else:
                plot_orien = gvect
        elif force_emisphere == 'upper':
            default_marker = default_gvect_marker_upward_symbol
            if gvect.is_downward:
                plot_orien = gvect.opposite()
            else:
                plot_orien = gvect
        elif not force_emisphere:
            plot_orien = gvect
            default_marker = default_gvect_marker_downward_symbol if not plot_orien.is_upward else default_gvect_marker_upward_symbol
        else:
            raise PlotException("Invalid force emisphere parameter")

        if plot_orien.is_upward:  # apparently mplstereonet does not handle negative plunges
            plot_orien = plot_orien.mirrorHoriz()

        bearing, plunge = plot_orien.d
        symbol = kwargs.get("m", default_marker)
        color = kwargs.get("c", default_gvect_marker_color)

        return plunge, bearing, symbol, color

    def params_gaxis(gaxis, kwargs, force_emisphere):

        if (not force_emisphere) or (force_emisphere == 'lower'):
            default_marker = default_gaxis_marker_downward_symbol
            if gaxis.is_upward:
                plot_gaxis = gaxis.opposite()
            else:
                plot_gaxis = gaxis
        elif force_emisphere == 'upper':
            default_marker = default_gaxis_marker_upward_symbol
            if gaxis.is_downward:
                plot_gaxis = gaxis.opposite()
            else:
                plot_gaxis = gaxis
        else:
            raise PlotException("Invalid force emisphere parameter")

        if plot_gaxis.is_upward:  # apparently mplstereonet does not handle negative plunges
            plot_gaxis = plot_gaxis.mirrorHoriz()

        bearing, plunge = plot_gaxis.d
        symbol = kwargs.get("m", default_marker)
        color = kwargs.get("c", default_gaxis_marker_color)

        return plunge, bearing, symbol, color

    def params_gplane(gplane, kwargs, force_emisphere):

        if (not force_emisphere) or (force_emisphere == 'lower'):
            default_line_style = default_gplane_downward_linestyle
            plot_gplane = gplane
        elif force_emisphere == 'upper':
            default_line_style = default_gplane_upward_linestyle
            plot_gplane = gplane.mirrorVertPPlane()
        else:
            raise PlotException("Invalid force emisphere parameter")

        strike, dip = plot_gplane.srda

        line_style = kwargs.get("m", default_line_style)
        color = kwargs.get("c", default_gplane_line_color)

        return strike, dip, line_style, color

    if not isinstance(data, list):
        data = [data]

    if force not in ('', 'upper', 'lower'):
        raise PlotException("Force parameter not isValid")

    fig, ax = ms.subplots()

    for rec in data:

        if isinstance(rec, tuple):
            if isinstance(rec[-1], str):
                params = rec[-1]
                objs = rec[:-1]
            else:
                params = None
                objs = rec
        else:
            objs = [rec]
            params = None

        if params:
            kwargs = string2dict(params)
        else:
            kwargs = dict()

        for obj in objs:
            if isDirect(obj):
                plunge, bearing, symbol, color = params_gvect(
                    obj,
                    kwargs,
                    force_emisphere=force)
                ax.line(
                    plunge,
                    bearing,
                    marker=symbol,
                    markerfacecolor=color,
                    markeredgecolor=color)
            elif isAxis(obj):
                plunge, bearing, symbol, color = params_gaxis(
                    obj,
                    kwargs,
                    force_emisphere=force)
                ax.line(
                    plunge,
                    bearing,
                    marker=symbol,
                    markerfacecolor=color,
                    markeredgecolor=color)
            elif isPlane(obj):
                strike, dip, linestyle, color = params_gplane(
                    obj,
                    kwargs,
                    force_emisphere=force)
                ax.plane(
                    strike,
                    dip,
                    linestyle=linestyle,
                    color=color)

    ax.grid()
    plt.show()


def splotDirect(gvect):

    if gvect.is_upward:
        plunge, bearing = gvect.mirrorHoriz().pt
        symbol = "x"
    else:
        plunge, bearing = gvect.pt
        symbol = "o"

    fig, ax = ms.subplots()

    ax.line(plunge, bearing, marker=symbol)
    ax.grid()
    plt.show()


class PlotException(Exception):

    pass

