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


from PyQt4 import QtCore, QtGui, uic
from PyQt4.QtCore import pyqtSignal

from plot_dialog import PlotDialog

FORM_CLASS, _ = uic.loadUiType(os.path.join(
    os.path.dirname(__file__), 'spectralprofiler_dockwidget_base.ui'))


class SpectralProfilerDockWidget(QtGui.QDockWidget, FORM_CLASS):

    wavelengths = [  512.6,   518.4,   524.7,   530.4,   536.5,   542.8,   548.7,
         554.5,   560.5,   566.7,   572.6,   578.5,   584.5,   590.6,
         596.7,   602.5,   608.6,   614.6,   620.5,   626.7,   632.7,
         638.6,   644.6,   650.6,   656.6,   662.6,   668.8,   674.7,
         680.6,   686.7,   692.6,   698.6,   704.7,   710.8,   716.7,
         722.7,   728.7,   734.7,   740.7,   746.8,   752.8,   758.7,
         764.8,   770.7,   776.7,   782.7,   788.8,   794.7,   800.7,
         806.8,   812.7,   818.7,   824.8,   830.8,   836.8,   842.8,
         848.8,   854.6,   860.7,   866.7,   872.7,   878.7,   884.6,
         890.7,   896.6,   902.7,   908.7,   914.6,   920.6,   926.6,
         932.6,   938.6,   944.6,   950.6,   955.4,   963.5,   971.4,
         979.7,   987.6,   993.7,  1013.1,  1019.5,  1027.7,  1035.5,
        1043.6,  1051.7,  1059.7,  1067.8,  1075.8,  1083.6,  1091.8,
        1099.7,  1107.7,  1115.9,  1123.8,  1131.8,  1139.7,  1147.8,
        1155.7,  1163.8,  1171.8,  1179.8,  1187.8,  1195.8,  1203.9,
        1211.9,  1219.8,  1227.9,  1235.9,  1244. ,  1252. ,  1259.8,
        1267.8,  1275.9,  1284.2,  1292. ,  1299.8,  1307.8,  1315.9,
        1323.8,  1331.8,  1339.8,  1347.8,  1355.8,  1363.8,  1371.8,
        1379.8,  1387.8,  1395.9,  1403.8,  1411.8,  1419.8,  1427.9,
        1435.7,  1443.8,  1451.9,  1459.8,  1467.8,  1475.8,  1483.9,
        1491.8,  1499.8,  1507.8,  1515.7,  1523.8,  1531.7,  1539.7,
        1547.7,  1555.5,  1563.7,  1571.7,  1579.6,  1587.7,  1595.7,
        1603.7,  1611.7,  1620.1,  1628.1,  1636.1,  1644.2,  1717.6,
        1725.6,  1733.7,  1742. ,  1749.7,  1757.7,  1766.3,  1773.6,
        1782.2,  1789.8,  1797.6,  1805.8,  1813.7,  1822. ,  1830. ,
        1837.6,  1845.6,  1853.7,  1861.8,  1870.1,  1877.3,  1885.7,
        1893.7,  1901.5,  1910. ,  1918. ,  1925.3,  1934.3,  1948.8,
        1957.6,  1965.9,  1973.3,  1981.3,  1989.4,  1997.7,  2005.8,
        2013. ,  2021.5,  2029.3,  2037.4,  2045.8,  2053.3,  2061.3,
        2069.4,  2077. ,  2085.5,  2093. ,  2101.9,  2109.2,  2117. ,
        2125.4,  2132.9,  2141.5,  2149. ,  2156.8,  2165.2,  2172.8,
        2181. ,  2189.4,  2196.8,  2204.7,  2213. ,  2221.2,  2228.7,
        2236.8,  2245. ,  2252.5,  2260.7,  2269.2,  2276.6,  2284.7,
        2292.7,  2300.4,  2308.9,  2316.4,  2324. ,  2332.6,  2340.6,
        2348.3,  2356.2,  2364.6,  2372.2,  2380.2,  2388.5,  2396.2,
        2404.2,  2412.2,  2420.2,  2428. ,  2436.3,  2444.3,  2451.9,
        2460.1,  2467.9,  2476. ,  2484.1,  2492.6,  2500.1,  2508.1,
        2516.1,  2524.1,  2532.1,  2540. ,  2548. ,  2556. ,  2564. ,
        2572. ,  2579.9,  2587.9]

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

        # lrm
        # # Continuum Correction
        # self.continuum_bottom.valueChanged.connect(lambda i: self.continuum_bottom_label.setText('{}'.format(
        #     self.wavelengths[i])))
        # self.continuum_top.valueChanged.connect(lambda i: self.continuum_top_label.setText('{}'.format(
        #     self.wavelengths[i])))

        # lrm
        # # Add + icon to button to add a wavelength
        # self.add_endpoint.setIcon(QtGui.QIcon(os.path.join(os.path.dirname(os.path.abspath(__file__)),
        #                                                    'icons/add.png')))
        # # Add - icon to button for removing a wavelength
        # self.remove_endpoint.setIcon(QtGui.QIcon(os.path.join(os.path.dirname(os.path.abspath(__file__)),
        #                                                       'icons/remove.png')))

        # lrm
        # self.add_endpoint.clicked.connect(self.add_continuum_endpoint)
        # self.remove_endpoint.clicked.connect(self.remove_continuum_endpoint)

        # Clipping
        self.clipping_lower.valueChanged.connect(lambda i: self.clipping_lower_label.setText('{}'.format(
            self.wavelengths[i])))
        self.clipping_upper.valueChanged.connect(lambda i: self.clipping_upper_label.setText('{}'.format(
            self.wavelengths[i])))

        # #Spectral Smoothing

        self.plot_selected.clicked.connect(self.plot)


    def closeEvent(self, event):
        self.closingPlugin.emit()
        event.accept()

    def add_continuum_endpoint(self):
        widget = QtGui.QWidget()
        hbox = QtGui.QHBoxLayout()
        hbox.setContentsMargins(0,0,0,0)
        widget.setLayout(hbox)

        slider = QtGui.QSlider(QtCore.Qt.Horizontal)
        slider.setRange(0, 268)
        label = QtGui.QLabel('{}'.format(self.wavelengths[0]))
        slider.valueChanged.connect(lambda i: label.setText('{}'.format(self.wavelengths[i])))

        hbox.addWidget(slider)
        hbox.addWidget(label)

        self.cc_slider_vlb.addWidget(widget)

    def remove_continuum_endpoint(self):
        row = self.cc_slider_vlb.count()
        widget = self.cc_slider_vlb.takeAt(row - 1)

        if widget is not None:
            widget.widget().deleteLater()

    # This is called when "Plot Selected" button is hit
    def plot(self):
        print("spectralprofiler_dockwidget : plot")

        # Get the continuum endpoints
        # continuum_endpoints = []
        # for i in range(self.cc_slider_vlb.count()):
        #     row_widget = self.cc_slider_vlb.itemAt(i).widget()
        #     for element in row_widget.children():
        #         if isinstance(element, QtGui.QLabel):
        #             continuum_endpoints.append(float(element.text()))
        # lrm
        continuum_endpoints = [500., 1700.]  # This is for the correction

        # lrm
        #correction_method = self.correction_method.currentText()
        #print("type of correction_method = {}".format(type(correction_method)))
        # smoothing_method = self.smoothing_method.currentText()
        # lrm
        correction_method = unicode("Linear")
        #smoothing_method = unicode("Gaussian")
        smoothing_method = unicode("None")
        print("spectralprofiler_dockwidget : plot : correction_method = {}".format(correction_method))
        print("spectralprofiler_dockwidget : plot : smoothing_method = {}".format(smoothing_method))


        # lrm
        #smoothing_window_size = int(self.smoothing_window_size.value())
        smoothing_window_size = 7


        #offset = float(self.offset.value())
        #print("spectralprofiler_dockwidget : plot : offset = {}".format(offset))
        offset = float(0.)  # lrm
        print("spectralprofiler_dockwidget : plot : offset = {}".format(offset))

        # The plotting x-axis range
        # clipping_lower = float(self.clipping_lower_label.text())
        # clipping_upper = float(self.clipping_upper_label.text())
        clipping_lower = float(500.)
        clipping_upper = float(1700.)
        print("spectralprofiler_dockwidget : plot : clipping_lower, clipping_upper = {} {}".format(clipping_lower, clipping_upper))

        # pcorrect = None
        # for i in range(self.correction_vlb.count()):
        #     widget = self.correction_vlb.itemAt(i).widget()
        #     if isinstance(widget, QtGui.QRadioButton):
        #         if widget.isChecked():
        #             pcorrect = widget.text()
        pcorrect = unicode("Highlands")
        print("spectralprofiler_dockwidget : plot : pcorrect = {}".format(pcorrect))

        # Create the plot dialog
        dialog = PlotDialog()

        dialog.wavelengths = self.wavelengths
        dialog.setWindowTitle('{}'.format(self.window_key))
        self.plot_windows[self.window_key] = dialog
        self.window_key += 1

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
        d = defaultdict(list)
        for s in selected:
            fname = s.attribute('FILENAME')
            id = s.attribute('OBSERVATION_ID')
            d[fname].append(id)

        # d is filename and observation numbers
        # u'SP_2C_02_03860_S136_E3557.spc': [20, 21, 22, 23, 24, 25, 26, 27]
        print("spectralprofiler_dockwidget : plot : d = {}".format(d))


        selected_spectra = {}  # u'SP_2C_02_03860_S136_E3557.spc' - a pandas class ***

        for k, obs in d.iteritems():
            selected_spectra[k] = spectra[k].spectra.iloc[obs]
        print("spectralprofiler_dockwidget : plot : selected_spectra = {}".format(selected_spectra))
        dialog.set_spectra(selected_spectra)

        # THIS IS WHERE WE GET THE TKINTER ERROR :
        dialog.plot(continuum_endpoints, correction_method, smoothing_method, smoothing_window_size, offset,
                    clipping_lower, clipping_upper, pcorrect)