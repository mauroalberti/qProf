
from builtins import str

from qgis.PyQt.QtCore import QFileInfo
from qgis.PyQt.QtWidgets import QFileDialog



def update_directory_key(settings, settings_dir_key, fileName):
    """
    modified from module RASTERCALC by Barry Rowlingson
    """

    path = QFileInfo(fileName).absolutePath()
    settings.setValue(settings_dir_key,
                      str(path))


def new_file_path(parent, show_msg, path, filter_text):

    output_filename, __ = QFileDialog.getSaveFileName(
        parent,
        show_msg,
        path,
        filter_text
    )
    if not output_filename:
        return ''
    else:
        return output_filename


def old_file_path(parent, show_msg, filter_extension, filter_text):

    input_filename, __ = QFileDialog.getOpenFileName(parent,
                                                 parent.tr(show_msg),
                                                 filter_extension,
                                                 filter_text)
    if not input_filename:
        return ''
    else:
        return input_filename
