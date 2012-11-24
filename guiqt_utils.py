
import os

from PyQt4.QtCore import *
from PyQt4.QtGui import *


def push_button_with_icon( icon_path, button_size, callback_function ):

    image_pixmap =  QPixmap( icon_path )
    image_icon =  QIcon( image_pixmap )
   
    push_button = QPushButton() 
    push_button.setIcon( image_icon )
    push_button.setIconSize( image_pixmap.rect().size())
    push_button.setFixedSize ( button_size[0], button_size[1] )        
    QObject.connect(push_button, 
                 SIGNAL("clicked()"),
                 callback_function)
    
    return push_button         