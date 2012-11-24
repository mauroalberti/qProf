"""
        
        This file contains modified code from Tosi book - Matplotlib for Python Developers
        
"""


import os
import webbrowser

# Python Qt4 bindings for GUI objects
from PyQt4 import QtCore, QtGui

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
    
    def __init__( self, map_extent_x, map_extent_y ):
        
        self.set_rcParams()        
        self.setup_Figure(map_extent_x, map_extent_y)        
        self.setup_FigureCanvas()
            

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


    def setup_Figure(self, map_extent_x, map_extent_y):        
        
        # setup Matplotlib Figure and Axis
        
        self.fig = Figure()
        
        self.ax = self.fig.add_subplot(111)
        
        #self.ax.set_aspect(1.)

        x_min_init,x_max_init = map_extent_x 
        y_min_init,y_max_init = map_extent_y
                
        self.ax.set_xlim(x_min_init,x_max_init)
        self.ax.set_ylim(y_min_init,y_max_init) 
        
        
    def setup_FigureCanvas(self):        
        
        FigureCanvas.__init__(self, self.fig)
        
        """
        # we define the widget as expandable
        FigureCanvas.setSizePolicy(self,
                                   QtGui.QSizePolicy.Expanding,
                                   QtGui.QSizePolicy.Expanding)
        
        # notify the system of updated policy
        FigureCanvas.updateGeometry(self) 
        """
        
        
class MplWidget(QtGui.QWidget):
    """
    Widget defined in Qt Designer.
    """
    
    def __init__(self, map_extent_x, map_extent_y, parent = None):
        
        # initialization of Qt MainWindow widget
        QtGui.QWidget.__init__(self, parent)

        # set the canvas to the Matplotlib widget
        self.canvas = MplCanvas( map_extent_x, map_extent_y )
                
        program_folder = os.path.join(os.path.dirname(__file__), "ims" )
        
        # manage the navigation toolbar
        self.ntb = NavigationToolbar(self.canvas, self)        
        #self.ntb.removeAction(self.ntb.buttons[0])
        self.ntb.clear()

        a = self.ntb.addAction(self.ntb._icon(os.path.join(program_folder, "world.png")), 'Home', self.zoom2fullextent)
        a.setToolTip('Reset original view')
        a = self.ntb.addAction(self.ntb._icon(os.path.join(program_folder, "arrow_left.png")), 'Back', self.ntb.back)
        a.setToolTip('Back to previous view')
        a = self.ntb.addAction(self.ntb._icon(os.path.join(program_folder, "arrow_right.png")), 'Forward', self.ntb.forward)
        a.setToolTip('Forward to next view')

        a = self.ntb.addAction(self.ntb._icon(os.path.join(program_folder, "arrow_out.png")), 'Pan', self.ntb.pan)
        a.setToolTip('Pan axes with left mouse, zoom with right')
        a = self.ntb.addAction(self.ntb._icon(os.path.join(program_folder, "zoom.png")), 'Zoom', self.ntb.zoom)
        a.setToolTip('Zoom to rectangle')


                                          
        # create a vertical box layout
        self.vbl = QtGui.QVBoxLayout()
        
        # add widgets to the vertical box
        self.vbl.addWidget(self.ntb)
        self.vbl.addWidget(self.canvas)
                                                             
        # set the layout to the vertical box
        self.setLayout(self.vbl)
        
           
                
    def zoom2fullextent(self):
        """
        Emit the signal for updating the map view to the extent of the DEM, in alternative of
        the shapefile, or at the standard extent.
        
        """
        
        self.canvas.emit( QtCore.SIGNAL("zoom_to_full_view") )
 
 

        
               
        
