from PyQt4.QtCore import *
from PyQt4.QtGui import *

import resources

from qProf_QWidget import qprof_QWidget

_plugin_name_ = "qProf"


class qProf_main(object):

    def __init__(self, interface):

        self.plugin_name = _plugin_name_
        self.interface = interface
        self.main_window = self.interface.mainWindow()
        self.canvas = self.interface.mapCanvas()


    def initGui(self):

        self.qprof_QAction = QAction(QIcon(":/plugins/{}/icon.png".format(self.plugin_name)),
                                     self.plugin_name,
                                     self.interface.mainWindow())
        self.qprof_QAction.setWhatsThis("Topographic and geological profiles")
        self.qprof_QAction.triggered.connect(self.open_qprof)
        self.interface.addPluginToMenu(self.plugin_name,
                                       self.qprof_QAction)

    def unload(self):

        self.interface.removePluginMenu(self.plugin_name, self.qprof_QAction)

    def open_qprof(self):

        qprof_DockWidget = QDockWidget(self.plugin_name,
                                       self.interface.mainWindow())
        qprof_DockWidget.setAttribute(Qt.WA_DeleteOnClose)
        qprof_DockWidget.setAllowedAreas(Qt.LeftDockWidgetArea | Qt.RightDockWidgetArea)
        self.qProf_QWidget = qprof_QWidget(self.plugin_name,
                                           self.canvas)
        qprof_DockWidget.setWidget(self.qProf_QWidget)
        qprof_DockWidget.destroyed.connect(self.qProf_QWidget.closeEvent)
        self.interface.addDockWidget(Qt.RightDockWidgetArea, qprof_DockWidget)
