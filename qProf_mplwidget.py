"""
        
        This file contains modified code from Tosi book - Matplotlib for Python Developers
        
"""

# Python Qt4 bindings for GUI objects
from PyQt4.QtCore import *
from PyQt4.QtGui import *

from matplotlib import rcParams


# import the Qt4Agg FigureCanvas object, that binds Figure to
# Qt4Agg backend. It also inherits from QWidget
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas

# import the NavigationToolbar Qt4Agg widget
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QTAgg as NavigationToolbar

# Matplotlib Figure object
from matplotlib.figure import Figure

       

class MplCanvas(FigureCanvas):
    """
    Class to represent the FigureCanvas widget.
    """
    
    def __init__( self ):
 
       
        self.set_rcParams()    
         
        self.fig = Figure()
        #self.ax = self.fig.add_subplot(111)
        FigureCanvas.__init__(self, self.fig)
            

    def set_rcParams(self):
    
        rcParams["font.size"] = 9.0
        rcParams["xtick.direction"] = 'out'        
        rcParams["ytick.direction"] = 'out'  
        
        rcParams["figure.subplot.left"] = 0.1  
        rcParams["figure.subplot.right"] = 0.96  
        rcParams["figure.subplot.bottom"] = 0.06  
        rcParams["figure.subplot.top"] = 0.96  
        rcParams["figure.subplot.wspace"] = 0.1  
        rcParams["figure.subplot.hspace"] = 0.1  
        
        rcParams["figure.facecolor"] = 'white' 

       
        
class MplWidget( QWidget ):
    """
    Widget defined in Qt Designer.
    """
    
    def __init__(self, window_title = 'Profile', parent = None):
        
        # initialization of Qt MainWindow widget
        QWidget.__init__(self)
        self.setWindowTitle ( window_title )

        # set the canvas and the navigation toolbar
        self.canvas = MplCanvas( )
        self.ntb = NavigationToolbar(self.canvas, self) 

        inputWidget = QWidget() 
        inputLayout = QHBoxLayout() 
        inputLayout.addWidget( QLabel( self.tr("Set profile colors") ) )        
        inputWidget.setLayout(inputLayout) 
        
        # manage the navigation toolbar
        self.window_tabs = QTabWidget()
        self.window_tabs.setSizePolicy(QSizePolicy.Minimum, QSizePolicy.Fixed)        
        self.window_tabs.addTab( self.ntb, "Map" )
        
        #TO BE ADDED IN SUSEQUENT RELEASE
        #self.window_tabs.addTab( inputWidget, "Rendering" )         

        # create a vertical box layout
        self.vbl = QVBoxLayout()       
        
        # add widgets to the vertical box
        self.vbl.addWidget(self.window_tabs)
        self.vbl.addWidget(self.canvas)
                                                             
        # set the layout to the vertical box
        self.setLayout(self.vbl)
        
        self.show()
        
