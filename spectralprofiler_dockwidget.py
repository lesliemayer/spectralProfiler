# -*- coding: utf-8 -*-
"""
/***************************************************************************
 SpectralProfilerDockWidget
                                 A QGIS plugin
 SP2
                             -------------------
        begin                : 2016-08-02
        git sha              : $Format:%H$
        copyright            : (C) 2016 by J
        email                : j
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
"""

from collections import defaultdict
import os
import logging

from PyQt4 import QtCore, QtGui, uic
from PyQt4.QtCore import pyqtSignal

from plot_dialog import PlotDialog

FORM_CLASS, _ = uic.loadUiType(os.path.join(
    os.path.dirname(__file__), 'spectralprofiler_dockwidget_base.ui'))


class SpectralProfilerDockWidget(QtGui.QDockWidget, FORM_CLASS):

    # wavelengths = [  512.6,   518.4,   524.7,   530.4,   536.5,   542.8,   548.7,
    #      554.5,   560.5,   566.7,   572.6,   578.5,   584.5,   590.6,
    #      596.7,   602.5,   608.6,   614.6,   620.5,   626.7,   632.7,
    #      638.6,   644.6,   650.6,   656.6,   662.6,   668.8,   674.7,
    #      680.6,   686.7,   692.6,   698.6,   704.7,   710.8,   716.7,
    #      722.7,   728.7,   734.7,   740.7,   746.8,   752.8,   758.7,
    #      764.8,   770.7,   776.7,   782.7,   788.8,   794.7,   800.7,
    #      806.8,   812.7,   818.7,   824.8,   830.8,   836.8,   842.8,
    #      848.8,   854.6,   860.7,   866.7,   872.7,   878.7,   884.6,
    #      890.7,   896.6,   902.7,   908.7,   914.6,   920.6,   926.6,
    #      932.6,   938.6,   944.6,   950.6,   955.4,   963.5,   971.4,
    #      979.7,   987.6,   993.7,  1013.1,  1019.5,  1027.7,  1035.5,
    #     1043.6,  1051.7,  1059.7,  1067.8,  1075.8,  1083.6,  1091.8,
    #     1099.7,  1107.7,  1115.9,  1123.8,  1131.8,  1139.7,  1147.8,
    #     1155.7,  1163.8,  1171.8,  1179.8,  1187.8,  1195.8,  1203.9,
    #     1211.9,  1219.8,  1227.9,  1235.9,  1244. ,  1252. ,  1259.8,
    #     1267.8,  1275.9,  1284.2,  1292. ,  1299.8,  1307.8,  1315.9,
    #     1323.8,  1331.8,  1339.8,  1347.8,  1355.8,  1363.8,  1371.8,
    #     1379.8,  1387.8,  1395.9,  1403.8,  1411.8,  1419.8,  1427.9,
    #     1435.7,  1443.8,  1451.9,  1459.8,  1467.8,  1475.8,  1483.9,
    #     1491.8,  1499.8,  1507.8,  1515.7,  1523.8,  1531.7,  1539.7,
    #     1547.7,  1555.5,  1563.7,  1571.7,  1579.6,  1587.7,  1595.7,
    #     1603.7,  1611.7,  1620.1,  1628.1,  1636.1,  1644.2,  1717.6,
    #     1725.6,  1733.7,  1742. ,  1749.7,  1757.7,  1766.3,  1773.6,
    #     1782.2,  1789.8,  1797.6,  1805.8,  1813.7,  1822. ,  1830. ,
    #     1837.6,  1845.6,  1853.7,  1861.8,  1870.1,  1877.3,  1885.7,
    #     1893.7,  1901.5,  1910. ,  1918. ,  1925.3,  1934.3,  1948.8,
    #     1957.6,  1965.9,  1973.3,  1981.3,  1989.4,  1997.7,  2005.8,
    #     2013. ,  2021.5,  2029.3,  2037.4,  2045.8,  2053.3,  2061.3,
    #     2069.4,  2077. ,  2085.5,  2093. ,  2101.9,  2109.2,  2117. ,
    #     2125.4,  2132.9,  2141.5,  2149. ,  2156.8,  2165.2,  2172.8,
    #     2181. ,  2189.4,  2196.8,  2204.7,  2213. ,  2221.2,  2228.7,
    #     2236.8,  2245. ,  2252.5,  2260.7,  2269.2,  2276.6,  2284.7,
    #     2292.7,  2300.4,  2308.9,  2316.4,  2324. ,  2332.6,  2340.6,
    #     2348.3,  2356.2,  2364.6,  2372.2,  2380.2,  2388.5,  2396.2,
    #     2404.2,  2412.2,  2420.2,  2428. ,  2436.3,  2444.3,  2451.9,
    #     2460.1,  2467.9,  2476. ,  2484.1,  2492.6,  2500.1,  2508.1,
    #     2516.1,  2524.1,  2532.1,  2540. ,  2548. ,  2556. ,  2564. ,
    #     2572. ,  2579.9,  2587.9]

    closingPlugin = pyqtSignal()

    def __init__(self, parent=None):
        """Constructor."""
        super(SpectralProfilerDockWidget, self).__init__(parent)
        # Set up the user interface from Designer.
        # After setupUI you can access any designer object by doing
        # self.<objectname>, and you can use autoconnect slots - see
        # http://qt-project.org/doc/qt-4.8/designer-using-a-ui-file.html
        # #widgets-and-dialogs-with-auto-connect
        self.setupUi(self)

        self.initgui()
        self.window_key = 0
        self.plot_windows = {}

    def initgui(self):


        # #Spectral Smoothing

        # plugin won't load if this is commented out
        self.plot_selected.clicked.connect(self.plot)

        self.radioButton_2.setChecked(True)  # set Mare to true


    # This is called when "Plot Selected" button is hit
    def plot(self):
        print("spectralprofiler_dockwidget : plot")

        # pcorrect = None
        # for i in range(self.correction_vlb.count()):
        #     widget = self.correction_vlb.itemAt(i).widget()
        #     if isinstance(widget, QtGui.QRadioButton):
        #         if widget.isChecked():
        #             pcorrect = widget.text()

        # get the value of the radio buttons
        # pcorrect = None
        # for i in range(self.groupBox.count()):
        #     widget = self.groupBox.itemAt(i).widget()
        #     if isinstance(widget, QtGui.QRadioButton):
        #        if widget.isChecked():
        #           pcorrect = widget.text()




        # Create the plot dialog
        dialog = PlotDialog()

        # not using this window anymore - lrm
        # Set the title of the window to a number
        # dialog.setWindowTitle('{}'.format(self.window_key))
        # self.plot_windows[self.window_key] = dialog
        # self.window_key += 1

        # Get the selected spectra and add them to the dialog
        v_layer = self.parent.v_layer
        spectra = self.parent.spectra

        # lrm ----------------------------------------------
        # selected = v_layer.selectedFeatures()
        # for i in selected:
        #     #attrs = i.attributeMap()
        #     fields = i.fields()
        #     #print "fields = " + str(fields)
        #     # for (k, attr) in attrs.iteritems():
        #     #      print "%d: %s" % (k, attr.toString())
        # --------------------------------------------------

        # Get filename and selected sp data :
        selected = v_layer.selectedFeatures()

        # set up a default dictionary (defaultdict is from collections library)
        #  will look like :
        # {u'SP_2B2_01_00896_N233_E3127.spc': [75, 76, 77, 78, 79, 80]})
        d = defaultdict(list)
        for s in selected:
            fname = s.attribute('FILENAME')
            id = s.attribute('OBSERVATION_ID')
            d[fname].append(id)



        # d is filename and observation numbers
        # u'SP_2C_02_03860_S136_E3557.spc': [20, 21, 22, 23, 24, 25, 26, 27]
        print("spectralprofiler_dockwidget : plot : d = {}".format(d))

        # get the angles info for selected observations - lrm
        emission_angle = []
        for s in selected:
            emission_angle.append(s.attribute('EMISSION_ANGLE'))
        logging.debug("spectralprofiler_dockwidget : emission_angle = %s", emission_angle)
        
        incidence_angle = []
        for s in selected:
            incidence_angle.append(s.attribute('INCIDENCE_ANGLE'))
        logging.debug("spectralprofiler_dockwidget : incidence_angle = %s", incidence_angle)
        
        phase_angle = []
        for s in selected:
            phase_angle.append(s.attribute('PHASE_ANGLE'))
        logging.debug("spectralprofiler_dockwidget : phase_angle = %s", phase_angle)


        selected_spectra = {}  # u'SP_2C_02_03860_S136_E3557.spc' - a pandas class ***

        for k, obs in d.iteritems():
            selected_spectra[k] = spectra[k].spectra.iloc[obs]
        print("spectralprofiler_dockwidget : plot : selected_spectra = {}".format(selected_spectra))
        # lrm : need this or crashes
        dialog.set_spectra(selected_spectra)


        # lrm :
        #dialog.plot()
        dialog.plot(emission_angle, incidence_angle, phase_angle)

        # CALL sp_extract_plugin here????????
        # spectraPlot = SP_EXTRACT(spectra, emission_angle, incidence_angle, phase_angle)
        # spectraPlot.make_plots()