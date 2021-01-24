
from . import resources

from .qgis_utils.utils import *
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
            f":/plugins/{self.plugin_name}/icons/qProf_main.svg",
            self.plugin_name,
            self.open_main_widget,
            whats_this="Topographic and geological profiles",
            parent=self.interface.mainWindow())
        self.interface.addPluginToMenu(self.plugin_name,
                                       self.qactOpenMainWin)
        self.actions.append(self.qactOpenMainWin)

    def unload(self):

        self.interface.removePluginMenu(self.plugin_name, self.qactOpenMainWin)

    def open_main_widget(self):

        project = QgsProject.instance()

        if project.count() == 0:
            warn(
                parent=None,
                header=plugin_name,
                msg="No project/layer available.\nPlease open project and add layers.")
            return

        self.qProf_QWidget = ActionWidget(
            self.current_directory,
            self.plugin_name,
            self.canvas
        )

        # show dialog
        self.qProf_QWidget.show()

        '''
        qprof_DockWidget = QDockWidget(self.plugin_name,
                                       self.interface.mainWindow())
        qprof_DockWidget.setAttribute(Qt.WA_DeleteOnClose)
        qprof_DockWidget.setAllowedAreas(Qt.LeftDockWidgetArea | Qt.RightDockWidgetArea)
        self.qProf_QWidget = qprof_QWidget(self.plugin_name,
                                           self.canvas)
        qprof_DockWidget.setWidget(self.qProf_QWidget)
        qprof_DockWidget.destroyed.connect(self.qProf_QWidget.closeEvent)
        self.interface.addDockWidget(Qt.RightDockWidgetArea, qprof_DockWidget)
        '''

