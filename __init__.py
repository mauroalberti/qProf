"""
/***************************************************************************
 qProf - plugin for Quantum GIS

 geoProfiles
                              -------------------
        begin                : 2012-11-06
        version              : 0.2.1 for QuantumGIS 2.0 - 22/10/2013
        copyright            : (C) 2012-2013 M. Alberti (implementation and concept) and M. Zaniero (concept)
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


def classFactory(iface):    
    # create qgSurf_gui class   
    return qProf_gui(iface)



