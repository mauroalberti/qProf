
from __future__ import division

from PyQt4.QtCore import QSettings, QFileInfo
from PyQt4.QtGui import QFileDialog



def lastUsedDir(plugin_lastDir_QSettStr):
    """
    modified from module RASTERCALC by Barry Rowlingson

    Args:
        plugin_lastDir_QSettStr: string, e.g. "/qProf/lastDir"

    Returns:
        last plugin directory
    """

    settings = QSettings()
    return settings.value(plugin_lastDir_QSettStr, "", type=str)



def setLastUsedDir(plugin_lastDir_QSettStr, lastDir):
    """
    modified from module RASTERCALC by Barry Rowlingson
    """

    path = QFileInfo(lastDir).absolutePath()
    settings = QSettings()
    settings.setValue(plugin_lastDir_QSettStr, str(path))


def new_file_path(parent, show_msg, path, filter_text):

    output_filename = QFileDialog.getSaveFileName(parent,
                                                  show_msg,
                                                  path,
                                                  filter_text)
    if not output_filename:
        return ''
    else:
        return output_filename


def old_file_path(parent, show_msg, filter_extension, filter_text):

    input_filename = QFileDialog.getOpenFileName(parent,
                                                 parent.tr(show_msg),
                                                 filter_extension,
                                                 filter_text)
    if not input_filename:
        return ''
    else:
        return input_filename
