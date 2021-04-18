
from qgis.core import Qgis
from qgis.utils import iface

success_dur = 1  # sec
info_dur = 2  # secs
warn_dur = 4  # secs
err_dur = 6  # secs


def ok(
        header: str,
        msg: str
):
    iface.messageBar().pushMessage(
        header,
        msg,
        level=Qgis.Success,
        duration=success_dur
    )


def info(
        header: str,
        msg: str
):
    iface.messageBar().pushMessage(
        header,
        msg,
        level=Qgis.Info,
        duration=info_dur
    )


def warn(
    header: str,
    msg: str
):

    iface.messageBar().pushMessage(
        header,
        msg,
        level=Qgis.Warning,
        duration=warn_dur
    )


def error(
        header: str,
        msg: str
):
    iface.messageBar().pushMessage(
        header,
        msg,
        level=Qgis.Critical,
        duration=err_dur
    )