
from PyQt4.QtCore import *
from PyQt4.QtGui import *


def info(parent, header, msg):
    
    QMessageBox.information(parent, header, msg)


def warn(parent, header, msg):

    QMessageBox.warning(parent, header, msg)


def error(parent, header, msg):

    QMessageBox.error(parent, header, msg)
    
    
def update_ComboBox(combobox, init_choice, names):

    combobox.clear()

    if len(names) == 0:
        return

    if init_choice:
        combobox.addItem(init_choice)

    combobox.addItems(names)
