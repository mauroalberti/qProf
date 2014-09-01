"""
        
        This file contains modified code from Tosi book - Matplotlib for Python Developers
        
"""

import numpy as np

from PyQt4.QtCore import *
from PyQt4.QtGui import *

from matplotlib import rcParams
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QTAgg as NavigationToolbar
from matplotlib.figure import Figure

from .utils import valid_intervals

       

class MplCanvas(FigureCanvas):
    """
    Class to represent the FigureCanvas widget.
    """
    
    def __init__( self ): 
       
        self.set_rcParams()    
         
        self.fig = Figure()
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
        self.window_tabs.addTab( self.ntb, "View" )
        
        #TO BE ADDED IN SUBSEQUENT RELEASE
        #self.window_tabs.addTab( inputWidget, "Rendering" )         

        # create a vertical box layout
        self.vbl = QVBoxLayout()       
        
        # add widgets to the vertical box
        self.vbl.addWidget(self.window_tabs)
        self.vbl.addWidget(self.canvas)
                                                             
        # set the layout to the vertical box
        self.setLayout(self.vbl)
        
        self.show()
       
            
def plot_line( axes, x_list, y_list, linecolor, name='' ):
    
    if name != '':                       
        axes.plot( x_list, y_list,'-', color=linecolor, label = unicode(name) )
    else:
        axes.plot( x_list, y_list,'-', color=linecolor )
        
            
def plot_filled_line( axes, x_list, y_list, plot_y_min, facecolor, alpha = 0.1 ):

    y_values_array = np.array( y_list )
    x_values_array = np.array( x_list )
    for val_int in valid_intervals( y_values_array ):               
        axes.fill_between( x_values_array[ val_int['start'] : val_int['end']+1 ], 
                          plot_y_min, 
                          y_values_array[ val_int['start'] : val_int['end']+1 ], 
                          facecolor = facecolor, 
                          alpha = alpha )
    
         
