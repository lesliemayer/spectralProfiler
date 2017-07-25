import os
from PyQt4 import QtCore, QtGui, uic

import logging

from sp_extract_plugin import SP_EXTRACT

#class PlotDialog(cls_dialog, ui_dialog):
class PlotDialog():

    def __init__(self):
        #super(PlotDialog, self).__init__(*args, **kwargs)
        #self.setupUi(self)

        #self.initmpl()
        #self.inittree()

        #self.show()
        pass


    def set_spectra(self, spectra):
        '''Set object data to the input spectra'''

        self.data = spectra


    def plot(self, emission_angle, incidence_angle, phase_angle,
             isHighlands,
             oneUmLeftShoulder, oneUmRightShoulder, twoUmLeftShoulder, twoUmRightShoulder,
             oneUmMin, oneUmMax, twoUmMin, twoUmMax, obsId):
        '''Plot the spectral data.
        Parameters
        ----------
        emission_angle : float
            Emission angle of the observation
        incidence_angle : float
            Incidence angle of the observation
        phase_angle : float
            Phase angle of the observation
        isHighlands : bool
            Is this Highlands albedo?
        oneUmLeftShoulder : float
            Value of 1um left shoulder continuum slope
        oneUmRightShoulder : float
            Value of 1um right shoulder continuum slope
        twoUmLeftShoulder : float
            Value of 2um left shoulder continuum slope
        twoUmRightShoulder : float
            Value of 2um right shoulder continuum slope
        oneUmMin : float
            Min wavelength to plot for 1um data
        oneUmMax : float
            Max wavelength to plot for 1um data
        twoUmMin : float
            Min wavelength to plot for 2um data
        twoUmMax : float
            Max wavelength to plot for 2um data
        obsId : integer
            Which observation number this is
        '''

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

                # Check for valid data
                #remove_bad_values(spectra)

                logging.debug("plot_dialog : plot : oneUmLeftShoulder = %s", oneUmLeftShoulder)
                logging.debug("plot_dialog : plot : oneUmRightShoulder = %s", oneUmRightShoulder)


                # initialize 1um SP_EXTRACT
                spectraPlot1um = SP_EXTRACT(spectra, emission_angle, incidence_angle, phase_angle,
                                         isHighlands, oneUmLeftShoulder, oneUmRightShoulder,
                                         oneUmMin, oneUmMax, obsId, fname)
                spectraPlot1um.make_plots()

                # initialize 2um SP_EXTRACT
                spectraPlot2um = SP_EXTRACT(spectra, emission_angle, incidence_angle, phase_angle,
                                         isHighlands, twoUmLeftShoulder, twoUmRightShoulder,
                                         twoUmMin, twoUmMax, obsId, fname)
                spectraPlot2um.make_plots()

                # calculate Band Area Ratio
                bar = spectraPlot2um.get_bandArea()/spectraPlot1um.get_bandArea()
                logging.debug("band area ratio = %s",bar)


