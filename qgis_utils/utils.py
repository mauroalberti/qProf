from typing import Tuple

from qgis.core import QgsProject
from qgis.PyQt.QtGui import QIcon
from qgis.PyQt.QtWidgets import QAction


def check_project_planar_crs() -> Tuple[bool, str]:

    project = QgsProject.instance()

    if project.count() == 0:
        msg = "Is a project open or has layers loaded?"
        return False, msg

    curr_proj_crs = project.crs()

    if not curr_proj_crs.isValid():
        msg = "Current project crs is not valid.\nPlease apply a valid crs to the current project."
        return False, msg

    if curr_proj_crs.isGeographic():
        msg = "Current project crs is geographic.\nPlease apply a planar crs to the current project."
        return False, msg

    msg = "Project open and with planar crs defined."
    return True, msg


def create_action(
        icon_path,
        text,
        callback,
        enabled_flag=True,
        status_tip=None,
        whats_this=None,
        parent=None,
        object_name=None):

    """
    # adapted from RedLayers by E. Ferreguti
    Create an action.

    :param icon_path: Path to the icon for this action. Can be a resource
        path (e.g. ':/plugins/foo/bar.png') or a normal file system path.
    :type icon_path: str

    :param text: Text that should be shown in menu items for this action.
    :type text: str

    :param callback: Function to be called when the action is triggered.
    :type callback: function

    :param enabled_flag: A flag indicating if the action should be enabled
        by default. Defaults to True.
    :type enabled_flag: bool

    :param status_tip: Optional text to show in a popup when mouse pointer
        hovers over the action.
    :type status_tip: str

    :param parent: Parent widget for the new action. Defaults None.
    :type parent: QWidget

    :param whats_this: Optional text to show in the status bar when the
        mouse pointer hovers over the action.

    :param object_name: Optional name to identify objects during customization
    :type object_name: str

    :returns: The action that was created.
    :rtype: QAction
    """

    icon = QIcon(icon_path)

    action = QAction(icon, text, parent)
    if callback:
        action.triggered.connect(callback)
    action.setEnabled(enabled_flag)

    if status_tip:
        action.setStatusTip(status_tip)

    if whats_this:
        action.setWhatsThis(whats_this)

    if object_name:
        action.setObjectName(object_name)

    return action

