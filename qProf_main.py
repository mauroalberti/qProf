
from .qygsf.utils.qt_utils.actions import *
from .qygsf.utils.qgis_utils.project import check_project_planar_crs

from .config.params import *

from .qProf_widgets import *


class qProf_main(object):

    def __init__(self, interface):

        self.plugin_name = plugin_name
        self.interface = interface
        self.main_window = self.interface.mainWindow()
        self.canvas = self.interface.mapCanvas()
        self.current_directory = os.path.dirname(__file__)
        self.actions = []

    def initGui(self):

        self.qactOpenMainWin = create_action(
            icon_path=f":/plugins/{self.plugin_name}/icons/qProf_main.svg",
            text=self.plugin_name,
            callback=self.open_main_widget,
            whats_this="Topographic and geological profiles",
            parent=self.interface.mainWindow())
        self.interface.addPluginToMenu(self.plugin_name,
                                       self.qactOpenMainWin)
        self.actions.append(self.qactOpenMainWin)

    def unload(self):

        self.interface.removePluginMenu(self.plugin_name, self.qactOpenMainWin)

    def open_main_widget(self):

        correct, msg = check_project_planar_crs()

        if not correct:

            wrn(
                None,
                self.plugin_name,
                msg
            )
            return

        self.qProf_QWidget = ActionWidget(
            self.current_directory,
            self.plugin_name,
            self.canvas
        )

        # show dialog
        self.qProf_QWidget.show()

