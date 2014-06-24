
from PyQt4.QtCore import *
from PyQt4.QtGui import *

from qgis.core import *

import resources

from qProf_QWidget import qprof_QWidget



class qProf_gui( object ):

    def __init__(self, interface):

        self.interface = interface
        self.main_window = self.interface.mainWindow()
        self.canvas = self.interface.mapCanvas()        


    def initGui(self):
        
        self.qprof_QAction = QAction(QIcon(":/plugins/qProf/icon.png"), "qProf", self.interface.mainWindow())
        self.qprof_QAction.setWhatsThis( "Topographic and geological profiles" ) 
        self.qprof_QAction.triggered.connect( self.open_qprof )
        self.interface.addPluginToMenu("QProf", self.qprof_QAction)


    def unload(self):

        self.interface.removePluginMenu( "qProf", self.qprof_QAction )


    def open_qprof(self):

        qprof_DockWidget = QDockWidget( 'qProf', self.interface.mainWindow() )        
        qprof_DockWidget.setAttribute(Qt.WA_DeleteOnClose)
        qprof_DockWidget.setAllowedAreas( Qt.LeftDockWidgetArea | Qt.RightDockWidgetArea )        
        self.qProf_QWidget = qprof_QWidget( self.canvas )        
        qprof_DockWidget.setWidget( self.qProf_QWidget ) 
        qprof_DockWidget.destroyed.connect( self.qProf_QWidget.closeEvent )       
        self.interface.addDockWidget( Qt.RightDockWidgetArea, qprof_DockWidget )
        

