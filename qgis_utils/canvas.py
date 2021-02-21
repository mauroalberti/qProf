
from qgis.gui import *

from qgis.PyQt.QtCore import *
from qgis.PyQt.QtGui import *


class MapDigitizeTool(QgsMapTool):

    moved = pyqtSignal(dict)
    leftClicked = pyqtSignal(dict)
    rightClicked = pyqtSignal(dict)

    def __init__(self, canvas):

        QgsMapTool.__init__(self, canvas)
        self.canvas = canvas
        self.cursor = QCursor(Qt.CrossCursor)

    def canvasMoveEvent(self, event):

        print(f"DEBUG: {event.pos().x()}, {event.pos().y()}")

        self.moved.emit({'x': event.pos().x(), 'y': event.pos().y()})

    def canvasReleaseEvent(self, event):

        if event.button() == Qt.RightButton:
            self.rightClicked.emit({'x': event.pos().x(), 'y': event.pos().y()})
        elif event.button() == Qt.LeftButton:
            self.leftClicked.emit({'x': event.pos().x(), 'y': event.pos().y()})
        else:
            return

    def canvasDoubleClickEvent(self, event):

        self.doubleClicked.emit({'x': event.pos().x(), 'y': event.pos().y()})

    def activate(self):

        QgsMapTool.activate(self)
        self.canvas.setCursor(self.cursor)

    def deactivate(self):

        QgsMapTool.deactivate(self)

    def isZoomTool(self):

        return False

    def setCursor(self, cursor):

        self.cursor = QCursor(cursor)

"""
Modified from: profiletool, script: tools/ptmaptool.py

#-----------------------------------------------------------
# 
# Profile
# Copyright (C) 2008  Borys Jurgiel
# Copyright (C) 2012  Patrice Verchere
#-----------------------------------------------------------
# 
# licensed under the terms of GNU GPL 2
# 
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License along
# with this program; if not, print to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#
#---------------------------------------------------------------------
"""


class PointMapToolEmitPoint(QgsMapToolEmitPoint):

    def __init__(self,
                 canvas,
                 button):

        super(PointMapToolEmitPoint, self).__init__(canvas)
        self.canvas = canvas
        self.cursor = QCursor(Qt.CrossCursor)
        self.button = button

    def setCursor(self, cursor):

        self.cursor = QCursor(cursor)


def xy_from_canvas(
        canvas,
        position
):

    mapPos = canvas.getCoordinateTransform().toMapCoordinates(position["x"], position["y"])

    return mapPos.x(), mapPos.y()


