
from qgis.PyQt.QtWidgets import *


def info(parent, header, msg):

    QMessageBox.information(
        parent,
        header,
        str(msg)
    )


def warn(parent, header, msg):

    QMessageBox.warning(
        parent,
        header,
        str(msg)
    )


def error(
        parent,
        header,
        msg
):

    QMessageBox.critical(
        parent,
        header,
        str(msg)
    )

    
def update_ComboBox(combobox, init_choice, names):

    combobox.clear()

    if len(names) == 0:
        return

    if init_choice:
        combobox.addItem(init_choice)

    combobox.addItems(names)
