
# Import the PyQt and QGIS libraries
from PyQt4.QtCore import *
from PyQt4.QtGui import *

from qgis.core import *

# Initialize Qt resources from file resources.py
import resources

# Import the code for the dialog
from qProf_dialog import qProfDialog


class qProf_gui(object):

    def __init__(self, iface):
        # Save reference to the QGIS interface
        self.iface = iface

    def initGui(self):
        # Create action that will start plugin configuration
        self.action = QAction(QIcon(":/plugins/qProf/icon.png"), "qProf", self.iface.mainWindow())
                   
        # connect the action to the run method
        QObject.connect(self.action, SIGNAL("triggered()"), self.run)

        # Add toolbar button and menu item
        self.iface.addToolBarIcon(self.action)
        self.iface.addPluginToMenu("QProf", self.action)


    def unload(self):
        # Remove the plugin menu item and icon
        self.iface.removePluginMenu("qProf",self.action)
        self.iface.removeToolBarIcon(self.action)


    # run method that performs all the real work
    def run(self):

        # create and show the dialog
        dlg = qProfDialog()        
        dlg.show()
        dlg.exec_()
        

