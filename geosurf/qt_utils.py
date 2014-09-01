
from PyQt4.QtCore import QSettings, QFileInfo
from PyQt4.QtGui import QFileDialog


# from module RASTERCALC by Barry Rowlingson    
def lastUsedDir():
    settings = QSettings()
    return settings.value( "/qProf/lastDir", "", type=str )


# from module RASTERCALC by Barry Rowlingson
def setLastUsedDir(lastDir):
    path = QFileInfo( lastDir ).absolutePath()
    settings = QSettings()
    settings.setValue( "/qProf/lastDir", str(path) )
    
 
def define_save_file_name( parent, show_msg, filter_extension, filter_text ):
        
    output_filename = QFileDialog.getSaveFileName(parent, 
                                                  show_msg, 
                                                  filter_extension, 
                                                  filter_text )        
    if not output_filename: 
        return ''
    else:
        return output_filename 
    
    
def define_existing_file_name( parent, show_msg, filter_extension, filter_text ):
        
    input_filename = QFileDialog.getOpenFileName( parent, 
                                                  show_msg, 
                                                  filter_extension, 
                                                  filter_text )        
    if not input_filename: 
        return ''
    else:
        return input_filename   
    
        