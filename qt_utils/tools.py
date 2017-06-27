
from PyQt4.QtCore import *
from PyQt4.QtGui import *


def update_ComboBox(combobox, init_choice, names):

    combobox.clear()

    if len(names) == 0:
        return

    if init_choice:
        combobox.addItem(init_choice)

    combobox.addItems(names)
