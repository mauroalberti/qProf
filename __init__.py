"""
/***************************************************************************
 qProf - plugin for Quantum GIS

 geoProfiles
                              -------------------
        begin                : 2012-11-06
        version              : 0.1.3 for QuantumGIS - 2013-04-02
        copyright            : (C) 2012-2013 M. Alberti (implementation) and M. Zaniero (concept)
        email                : alberti.m65@gmail.com
        
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 3 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
"""


from qProf_gui import qProf_gui


def name():
    return "qProf"

def description():
    return "Compute topographic profiles from DEM and paths"

def version():
    return "0.1.3"

def icon():
    return "icon.png"

def qgisMinimumVersion():
    return "1.8"

def classFactory(iface):    
    # create qgSurf_gui class   
    return qProf_gui(iface)



