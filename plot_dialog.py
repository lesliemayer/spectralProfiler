import os
from PyQt4 import QtCore, QtGui, uic

# from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
# # Original
# #from matplotlib.backends.backend_qt4agg import NavigationToolbar2QTAgg as NavigationToolbar
# # L. Mayer 04/28/2017 :
# from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT as NavigationToolbar

#from matplotlib.figure import Figure

import logging

from sp_extract_plugin import SP_EXTRACT

#base = os.path.dirname(os.path.abspath(__file__))

# the dialog box for the plotting window (not the main one)
#ui_dialog, cls_dialog = uic.loadUiType(os.path.join(base, 'ui/plot_dialog.ui'))


#class PlotDialog(cls_dialog, ui_dialog):
class PlotDialog():

    #def __init__(self, *args, **kwargs):
    def __init__(self):
        #super(PlotDialog, self).__init__(*args, **kwargs)
        #self.setupUi(self)

        #self.initmpl()
        #self.inittree()

        #self.show()
        pass


    # lrm : get the spectrum data??  *** or set the spectrum data !!!!
    def set_spectra(self, spectra):

        self.data = spectra

        # for k, obs in self.data.iteritems():
        #     parent = QtGui.QStandardItem(k)
        #     for label, values in obs.iteritems():
        #         child = QtGui.QStandardItem('Observation {}'.format(label))
        #         child.setFlags(QtCore.Qt.ItemIsEnabled | QtCore.Qt.ItemIsSelectable)
        #         parent.appendRow(child)

            # obs tree - don't need
            #self.spectratree.model().appendRow(parent)


    def plot(self, emission_angle, incidence_angle, phase_angle,
             isHighlands,
             oneUmLeftShoulder, oneUmRightShoulder, twoUmLeftShoulder, twoUmRightShoulder,
             oneUmMin, oneUmMax, twoUmMin, twoUmMax, obsId):

    # dialog.plot(emission_angle, incidence_angle, phase_angle, self.isHighlands,
    #             self.oneUmLeftShoulder, self.oneUmRightShoulder, self.twoUmLeftShoulder, self.twoUmRightShoulder,
    #             self.oneUmMin, self.oneUmMax, self.twoUmMin, self.twoUmMax)

        # lrm
        # if pcorrect == 'Mare':
        #     key = 'REF2'
        # elif pcorrect == 'Highlands':
        #     key = 'REF1'
        # else:
        #     key = 'REF'  - this will make it crash : lrm

        # lrm
        key = 'REF1'  # use RAW, REF1, REF2, or QA  spectral data


        for fname, panel in self.data.iteritems():
            for k, df in panel.iteritems():

                print("plot_dialog : plot : k, df = {} {}".format(k,df))

                # lrm - this is necessary
                if key in df.columns:
                     spectra = df[key]
                else:
                     spectra = df['REF']


                print("plot_dialog : plot : spectra = {}".format(spectra))

                # *****************************************************************************
                # multiply dictionary df keys by .0001 as done in sp_extract :
                # Multiply every value in my_dict by 2
                # for key in df:
                #     df[key] *= .0001   numbers will be way off if you do this **********
                # *****************************************************************************

                df['CC'] = spectra

                logging.debug("plot_dialog : plot : oneUmLeftShoulder = %s", oneUmLeftShoulder)
                logging.debug("plot_dialog : plot : oneUmRightShoulder = %s", oneUmRightShoulder)


                # initialize 1um SP_EXTRACT
                spectraPlot = SP_EXTRACT(spectra, emission_angle, incidence_angle, phase_angle,
                                         isHighlands, oneUmLeftShoulder, oneUmRightShoulder,
                                         oneUmMin, oneUmMax, obsId, fname)
                spectraPlot.make_plots()

                # initialize 2um SP_EXTRACT
                spectraPlot = SP_EXTRACT(spectra, emission_angle, incidence_angle, phase_angle,
                                         isHighlands, twoUmLeftShoulder, twoUmRightShoulder,
                                         twoUmMin, twoUmMax, obsId, fname)
                spectraPlot.make_plots()

