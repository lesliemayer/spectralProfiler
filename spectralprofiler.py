# -*- coding: utf-8 -*-
"""
/***************************************************************************
 SpectralProfiler
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
from PyQt4.QtCore import QSettings, QTranslator, qVersion, QCoreApplication, Qt, QVariant
from PyQt4.QtGui import QAction, QIcon

from PyQt4 import QtCore, QtGui

# Initialize Qt resources from file resources.py
import resources

# Import dependencies
import io_spectral_profiler as sp
import qgis
from qgis.core import *

# Import the code for the DockWidget
from spectralprofiler_dockwidget import SpectralProfilerDockWidget
import os.path

# for writing to the log file
import logging

# lrm : set up the logging file
logging.basicConfig(level=logging.DEBUG,
                    filename=r'C:\Users\lrmayer\Documents\Mayer\QGIS_Plugin\sp_extract_plugin.log', filemode='w')

FIELDS = [QgsField('SPACECRAFT_CLOCK_COUNT', QVariant.Double),
          QgsField('VIS_FOCAL_PLANE_TEMPERATURE', QVariant.Double),
          QgsField('NIR1_FOCAL_PLANE_TEMPERATURE', QVariant.Double),
          QgsField('NIR2_FOCAL_PLANE_TEMPERATURE', QVariant.Double),
          QgsField('SPECTROMETER_TEMPERATURE_1', QVariant.Double),
          QgsField('SPECTROMETER_TEMPERATURE_2', QVariant.Double),
          QgsField('SPECTROMETER_TEMPERATURE_3', QVariant.Double),
          QgsField('SPECTROMETER_TEMPERATURE_4', QVariant.Double),
          QgsField('HALOGEN_BULB_RADIANCE', QVariant.Double),
          QgsField('HALOGEN_BULB_VOLTAGE1', QVariant.Double),
          QgsField('HALOGEN_BULB_VOLTAGE2', QVariant.Double),
          QgsField('HALOGEN_BULB_TEMPERATURE1', QVariant.Double),
          QgsField('HALOGEN_BULB_TEMPERATURE2', QVariant.Double),
          QgsField( 'SPACECRAFT_ALTITUDE', QVariant.Double),
          QgsField('SPACECRAFT_GROUND_SPEED', QVariant.Double),
          QgsField('SUB_SPACECRAFT_LATITUDE', QVariant.Double),
          QgsField('SUB_SPACECRAFT_LONGITUDE', QVariant.Double),
          QgsField('CENTER_LATITUDE', QVariant.Double),
          QgsField('CENTER_LONGITUDE', QVariant.Double),
          QgsField('EMISSION_ANGLE', QVariant.Double),
          QgsField('SPACECRAFT_AZIMUTH', QVariant.Double),
          QgsField('INCIDENCE_ANGLE', QVariant.Double),
          QgsField('SOLAR_AZIMUTH_ANGLE', QVariant.Double),
          QgsField('PHASE_ANGLE', QVariant.Double),
          QgsField('SP_TEMPERATURE', QVariant.Double),
          QgsField('SP_PELTIER_HOT_TEMPERATURE', QVariant.Double),
          QgsField('SP_N2_RADIATOR_TEMPERATURE', QVariant.Double),
          QgsField('SP_CAL_VIS_TEMPERATURE', QVariant.Double),
          QgsField('SP_CAL_NIR_TEMPERATURE', QVariant.Double),
          QgsField('DPU_TEMPERATURE', QVariant.Double),
          QgsField('SP_POWER_P5V', QVariant.Double),
              QgsField('SP_POWER_M15V', QVariant.Double),
                       QgsField('SP_POWER_P15V', QVariant.Double),
                                QgsField('CALIBRATION', QVariant.Double),
          QgsField('SP_PELTIER', QVariant.Double),
              QgsField('TC_MI_STATUS', QVariant.Double),
                       QgsField('CLOCK_COUNT_ERR_FLAG', QVariant.Double),
          QgsField('SPATIAL_RESOLUTION_FLAG', QVariant.Double),
              QgsField('GEOMETRIC_INFO_RECAL_FLAG', QVariant.Double),
          QgsField('SUPPORT_IMAGE_LINE_POSITION', QVariant.Double),
              QgsField('SUPPORT_IMAGE_COLUMN_POSITION', QVariant.Double),
          QgsField('THUMBNAIL_LINE_POSITION', QVariant.Double),
              QgsField('THUMBNAIL_COLUMN_POSITION', QVariant.Double),
          QgsField('FILENAME', QVariant.String),
          QgsField('OBSERVATION_ID', QVariant.String)]

FIELD_LOOKUP = {i: FIELDS[i] for i in range(45)}


class SpectralProfiler:
    """QGIS Plugin Implementation."""

    def __init__(self, iface):
        """Constructor.

        :param iface: An interface instance that will be passed to this class
            which provides the hook by which you can manipulate the QGIS
            application at run time.
        :type iface: QgsInterface
        """
        # Save reference to the QGIS interface
        self.iface = iface

        # initialize plugin directory
        self.plugin_dir = os.path.dirname(__file__)

        # initialize locale
        locale = QSettings().value('locale/userLocale')[0:2]
        locale_path = os.path.join(
            self.plugin_dir,
            'i18n',
            'SpectralProfiler_{}.qm'.format(locale))

        if os.path.exists(locale_path):
            self.translator = QTranslator()
            self.translator.load(locale_path)

            if qVersion() > '4.3.3':
                QCoreApplication.installTranslator(self.translator)

        # Declare instance attributes
        self.actions = []
        self.menu = self.tr(u'&Spectral Profiler')
        # TODO: We are going to let the user set this up in a future iteration
        self.toolbar = self.iface.addToolBar(u'SpectralProfiler')
        self.toolbar.setObjectName(u'SpectralProfiler')


        #print "** INITIALIZING SpectralProfiler"

        self.pluginIsActive = False
        self.dockwidget = None
        self.v_layer = None
        self.spectra = {}
        self.plots = {}
        self.canvas = qgis.utils.iface.mapCanvas()

        self.run()
    # noinspection PyMethodMayBeStatic
    def tr(self, message):
        """Get the translation for a string using Qt translation API.

        We implement this ourselves since we do not inherit QObject.

        :param message: String for translation.
        :type message: str, QString

        :returns: Translated version of message.
        :rtype: QString
        """
        # noinspection PyTypeChecker,PyArgumentList,PyCallByClass
        return QCoreApplication.translate('SpectralProfiler', message)


    def add_action(
        self,
        icon_path,
        text,
        callback,
        enabled_flag=True,
        add_to_menu=True,
        add_to_toolbar=True,
        status_tip=None,
        whats_this=None,
        parent=None):
        """Add a toolbar icon to the toolbar.

        :param icon_path: Path to the icon for this action. Can be a resource
            path (e.g. ':/plugins/foo/bar.png') or a normal file system path.
        :type icon_path: str

        :param text: Text that should be shown in menu items for this action.
        :type text: str

        :param callback: Function to be called when the action is triggered.
        :type callback: function

        :param enabled_flag: A flag indicating if the action should be enabled
            by default. Defaults to True.
        :type enabled_flag: bool

        :param add_to_menu: Flag indicating whether the action should also
            be added to the menu. Defaults to True.
        :type add_to_menu: bool

        :param add_to_toolbar: Flag indicating whether the action should also
            be added to the toolbar. Defaults to True.
        :type add_to_toolbar: bool

        :param status_tip: Optional text to show in a popup when mouse pointer
            hovers over the action.
        :type status_tip: str

        :param parent: Parent widget for the new action. Defaults None.
        :type parent: QWidget

        :param whats_this: Optional text to show in the status bar when the
            mouse pointer hovers over the action.

        :returns: The action that was created. Note that the action is also
            added to self.actions list.
        :rtype: QAction
        """

        base = os.path.dirname(os.path.abspath(__file__))
        icon_path = os.path.join(base, icon_path)
        icon = QIcon(icon_path)

        action = QAction(icon, text, parent)
        action.triggered.connect(callback)
        action.setEnabled(enabled_flag)

        if status_tip is not None:
            action.setStatusTip(status_tip)

        if whats_this is not None:
            action.setWhatsThis(whats_this)

        if add_to_toolbar:
            self.toolbar.addAction(action)

        if add_to_menu:
            self.iface.addPluginToMenu(
                self.menu,
                action)

        self.actions.append(action)

        return action

    def initGui(self):
        """Create the menu entries and toolbar icons inside the QGIS GUI."""

        self.add_action("icons/load.png",
                        "Load .spc File",
                        self.openspc)

        self.add_action("icons/plot.png",
                        "Open the plotting dock.",
                        self.reopen)



    #--------------------------------------------------------------------------

    def onClosePlugin(self):
        """Cleanup necessary items here when plugin dockwidget is closed"""

        #print "** CLOSING SpectralProfiler"

        # disconnects
        self.dockwidget.closingPlugin.disconnect(self.onClosePlugin)

        # remove this statement if dockwidget is to remain
        # for reuse if plugin is reopened
        # Commented next statement since it causes QGIS crashe
        # when closing the docked window:
        # self.dockwidget = None
        for layer in self.iface.legendInterface().layers():
            if layer.name() == 'sp_observations':
                try:
                    QgsMapLayerRegistry.instance().removeMapLayer(layer.id())
                except: pass

        self.pluginIsActive = False
        self.dockwidget = None

    def unload(self):
        """Removes the plugin menu item and icon from QGIS GUI."""

        #print "** UNLOAD SpectralProfiler"

        for action in self.actions:
            self.iface.removePluginMenu(
                self.tr(u'&Spectral Profiler'),
                action)
            self.iface.removeToolBarIcon(action)
        # remove the toolbar
        del self.toolbar

    def reopen(self):
        if self.dockwidget == None:
            self.dockwidget = SpectralProfilerDockWidget()
            self.dockwidget.parent = self
            self.iface.addDockWidget(Qt.RightDockWidgetArea, self.dockwidget)
        else:
            self.iface.addDockWidget(Qt.RightDockWidgetArea, self.dockwidget)
#--------------------------------------------------------------------------

    def run(self):
        """Run method that loads and starts the plugin"""
        if not self.pluginIsActive:
            self.pluginIsActive = True

            # dockwidget may not exist if:
            #    first run of plugin
            #    removed on close (see self.onClosePlugin method)
            if self.dockwidget == None:
                # Create the dockwidget (after translation) and keep reference
                self.dockwidget = SpectralProfilerDockWidget()
                self.dockwidget.parent = self

            # connect to provide cleanup on closing of dockwidget
            self.dockwidget.closingPlugin.connect(self.onClosePlugin)

            # show the dockwidget
            # TODO: fix to allow choice of dock location
            self.iface.addDockWidget(Qt.RightDockWidgetArea, self.dockwidget)
            self.dockwidget.show()

    # LRM : read the spectral data from the input file
    def openspc(self):

        print("SpectralProfiler : openspc : **************************************")


        fpaths = QtGui.QFileDialog.getOpenFileNames(self.dockwidget, 'Open Spectral Profiler', '*.spc')

        for fpath in fpaths:
            fname = os.path.basename(fpath)
            if fname is None or len(fname) == 0:
                return
            else:
                spectra = sp.Spectral_Profiler(fpath)
                self.spectra[fname] = spectra

                # lrm
                print("SpectralProfiler : openspc : type(spectra) = {}".format(type(spectra)))
                print("SpectralProfiler : openspc : self.spectra[fname] : {}".format(self.spectra[fname]) )
                logging.debug("SpectralProfiler : openspc : type(spectra) = %s", (type(spectra)))
                logging.debug("SpectralProfiler : openspc : fname) = %s", fname)
                logging.debug("SpectralProfiler : openspc : type(self.spectra[fname]) = %s", (type(self.spectra[fname])))



            self.draw_observations(fname)

    def draw_observations(self, fname):

        if self.v_layer is None:
            exists = False
            # Check if the layer already exists
            for layer in self.iface.legendInterface().layers():
                if layer.name() == 'sp_observations':
                    self.v_layer = layer
                    self.v_layer_provider = self.v_layer.dataProvider()
                    exists = True
                    break

            if not exists:
                self.v_layer = QgsVectorLayer("Point?crs=epsg:4326", "sp_observations", "memory")
                self.v_layer_provider = self.v_layer.dataProvider()
                self.v_layer_provider.addAttributes(FIELDS)
                self.v_layer.updateFields()

                #Set up the labels and symbology
                self.v_layer.setCustomProperty('labeling', 'pal')
                self.v_layer.setCustomProperty('labeling/fieldName', 'OBSERVATION_ID')
                self.v_layer.setCustomProperty('labeling/fontSize', '10')
                self.v_layer.setCustomProperty('labeling/placement', QgsPalLayerSettings.Line)
                self.v_layer.setCustomProperty('labeling/enabled', 'True')

                symbol = QgsMarkerSymbolV2.createSimple({'name': u'circle',
                                                         'color': u'red',
                                                         'size': u'2.5'})
                self.v_layer.rendererV2().setSymbol(symbol)

                QgsMapLayerRegistry.instance().addMapLayer(self.v_layer)

        sp = self.spectra[fname]
        latlon = sp.ancillary_data[['CENTER_LATITUDE', 'CENTER_LONGITUDE']]
        logging.debug("spectralprofiler : draw_observations : latlon = %s", latlon)

        emission_angle = sp.ancillary_data['EMISSION_ANGLE']
        logging.debug("spectralprofiler : draw_observations : emission_angle %s", emission_angle)

        

        for i in range(sp.nspectra):
            pt = QgsFeature()
            attributes = sp.ancillary_data.iloc[i].to_dict()
            for k, v in attributes.iteritems():
                attributes[k] = float(v)
            attributes['FILENAME'] = fname
            attributes['OBSERVATION_ID'] = i
            qattrs = []
            for j in range(45):
                key = FIELD_LOOKUP[j].name()
                if key in attributes.keys():

                    qattrs.append(attributes[key])
                else:
                    qattrs.append(None)
            pt.setAttributes(qattrs)

            ll = latlon.iloc[i].values
            pt.setGeometry(QgsGeometry.fromPoint(QgsPoint(ll[1], ll[0])))
            # TODO: Swap the buffer in for the geom.
            pt.geometry().buffer(500, 2)
            self.v_layer_provider.addFeatures([pt])

        self.v_layer.updateExtents()
        self.canvas.refresh()

