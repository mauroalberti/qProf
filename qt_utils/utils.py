
from PyQt4.QtCore import QSettings, QFileInfo


# from module RASTERCALC by Barry Rowlingson    
def lastUsedDir():
    settings = QSettings()
    return settings.value( "/qProf/lastDir", "", type=str )


# from module RASTERCALC by Barry Rowlingson
def setLastUsedDir(lastDir):
    path = QFileInfo( lastDir ).absolutePath()
    settings = QSettings()
    settings.setValue( "/qProf/lastDir", str(path) )
    
    