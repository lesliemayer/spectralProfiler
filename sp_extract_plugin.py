from __future__ import division

from struct import unpack, unpack_from

#import sys
import numpy as np
#from matplotlib.collections import LineCollection
from pylab import *
import argparse
#import os
import math
import logging  # lrm
import os  # for getting the plugin path

class SP_EXTRACT:

    def main():

        if __name__ == '__main__':
            parser = argparse.ArgumentParser(description='Spectral Profiler Reflectance Extraction Tool')

            parser.add_argument('input_data', action='store', help='The ".spc" file shipped with the SP data.')
            parser.add_argument('albedo_tab', action='store', help='The albedo table for the chosen overall reflectance (high, medium, or low).')
            parser.add_argument('-w', action='store',dest='wv_limits', default=1652, nargs=1, help='The limit wavelength to visualize to.')
            parser.add_argument('-s', action='store_true', dest='save', help='Save output to a CSV file.')
            parser.add_argument('observation', default=0,type=int, nargs='+', help='The range of observations to visualize.')
            args = parser.parse_args()

        # else:
        #     parser = argparse.ArgumentParser(description='Spectral Profiler Reflectance Extraction Tool')
        #     args = parser.parse_args()
        #     args.input_data = r'C:\Users\lrmayer\Documents\Mayer\QGIS_Plugin\Spectral_Profiler_2014\data\LeslieTest\SP_2B2_01_00896_N233_E3127.spc'
        #     args.albedo_tab = r'C:\Users\lrmayer\Documents\Mayer\QGIS_Plugin\Spectral_Profiler_2014\high_albedo_coefficients.csv'
        #     args.wv_limits = 1652
        #     args.save = True
        #     args.observation = [0]

    def __init__(self, spectra, emission_angle, incidence_angle, phase_angle, isHighlands,
                 leftShoulder, rightShoulder, plotMin, plotMax, obsID, fileName):
        # Set up logging

        # Get the directory name of the qgis plugin
        self.plugin_path = os.path.dirname(os.path.realpath(__file__))


        if (isHighlands):
        # Highlands
        #self.albedo_tab = r'C:\Users\lrmayer\Documents\Mayer\QGIS_Plugin\Spectral_Profiler_2014\high_albedo_coefficients.csv'
            #self.albedo_tab = r'C:\Users\lrmayer\.qgis2\python\plugins\SpectralProfiler\data\albedo\high_albedo_coefficients.csv'
            self.albedo_tab = self.plugin_path +  r'.\data\albedo\high_albedo_coefficients.csv'
            self.albedoTitle = 'Highlands'
        else:
            # Mare
            #self.albedo_tab = r'C:\Users\lrmayer\.qgis2\python\plugins\SpectralProfiler\data\albedo\low_albedo_coefficients.csv'
            self.albedo_tab = self.plugin_path + r'.\data\albedo\low_albedo_coefficients.csv'
            self.albedoTitle = 'Mare'

        # set plot x min, max
        self.wv_MinLimits = plotMin
        self.wv_MaxLimits = plotMax

        self.save = True
        self.observation = [0]
        self.emission_angle = np.asarray(emission_angle)
        self.incidence_angle = np.asarray(incidence_angle)
        self.phase_angle = np.asarray(phase_angle)
        self.wv_array = None

        # pass the observation spectrum data from the plugin
        logging.debug("__init__ : type(spectra) = %s", type(spectra))
        self.spectra = spectra

        # get the number of observations
        self.num_observations = len(emission_angle)
        logging.debug("sp_extract_plugin : __init__ : self.num_observations = %s", self.num_observations)


        # Set the left & right shoulders
        self.leftShoulder = leftShoulder
        self.rightShoulder = rightShoulder
        logging.debug("sp_extract_plugin : __init__ : self.leftShoulder = %s", self.leftShoulder)
        logging.debug("sp_extract_plugin : __init__ : self.rightShoulder = %s", self.rightShoulder)
        self.bandArea = None

        # number of observations
        self.obsId = obsID

        # filename
        self.fileName = fileName


    # Not using this - lrm ====================================================================================
    def openspc(self, input_data, save):

        """
        Parameters:

        input file     type .spc
                       This is the .spc file that contains the label and data.

        Returns:
        wavelength     type: ndarray
                       An array of wavelengths from all 3 detectors

        radiance       type: ndarray
                       An array of radiance values over the image.  This is binned into n observations.  In tests we have generally seen 44 observations per image.

        reflectance    type: ndarray
                       An array of reflectance values over the image.
                       This is binned into n observations.  In tests we have generally seen 44 observations per image.
        """

        label = open(input_data, 'r+b')
        for line in label:
            if "^SP_SPECTRUM_WAV" in line:
                wav_offset = int(line.split('=')[1].split(" ")[1])
            if "^SP_SPECTRUM_RAD" in line:
                rad_offset = int(line.split('=')[1].split(" ")[1])
            if "^SP_SPECTRUM_REF" in line:
                ref_offset = int(line.split('=')[1].split(" ")[1])
            if "OBJECT                               = SP_SPECTRUM_RAD" in line:
                line = label.next()
                rad_lines = int(line.split('=')[1])
            if "OBJECT                               = SP_SPECTRUM_REF" in line:
                line = label.next()
                ref_lines = int(line.split('=')[1])
            if 'NAME                         = "EMISSION_ANGLE"' in line:
                line = label.next();line = label.next(); line=label.next()
                emission_offset = int(line.split("=")[1])
            if 'NAME                         = "INCIDENCE_ANGLE"' in line:
                line = label.next();line = label.next(); line=label.next()
                incidence_offset = int(line.split("=")[1])
            if 'NAME                         = "PHASE_ANGLE"' in line:
                line = label.next();line = label.next(); line=label.next()
                phase_offset = int(line.split("=")[1])
            if 'NORMAL_SP_POINT_NUM' in line:
                num_observations = int(line.split("=")[1])
            if 'ROW_BYTES' in line:
                row_bytes = int(line.split("=")[1])
            if "^ANCILLARY_AND_SUPPLEMENT_DATA" in line:
                ancillary_offset = int(line.split("=")[1].split("<")[0])
            if "OBJECT                               = SP_SPECTRUM_QA" in line:
                #We only need ~20 lines, break before binary starts
                break

        # The Wavelength values (X axis on plot)
        label.seek(wav_offset-1) #Seek to the wavelength section
        array = np.fromstring(label.read(296*2), dtype='>H')
        wv_array = array.astype('float')
        wv_array *= 0.1

        #Radiance
        label.seek(rad_offset-1)
        array = np.fromstring(label.read(rad_lines*296*2), dtype='>H')
        rad_array = array.astype('float')
        rad_array *= 0.01
        rad_array = rad_array.reshape(rad_lines,296)
        #print rad_array

        #Reflectance
        label.seek(ref_offset-1) #Seek to the wavelength section
        array = np.fromstring(label.read(ref_lines*296*2), dtype='>H')
        ref_array = array.astype('float')
        ref_array *= 0.0001
        ref_array = ref_array.reshape(ref_lines,296)

        # lrm
        logging.debug("openspc : wv_array : %s", wv_array)
        logging.debug(" ")
        logging.debug("openspc : rad_lines : %s", rad_lines)
        logging.debug("openspc : rad_array[0,:] : %s", rad_array[0,:])
        logging.debug(" ")
        logging.debug("openspc : ref_lines : %s", ref_lines)
        logging.debug("openspc : ref_array[0,:] : %s", ref_array[0,:])

        #Parse the binary to get i, e, and phase for each observation
        angles = []
        for n in range(num_observations):
            #Emission Angle
            label.seek(ancillary_offset + (n*row_bytes-1) +  (emission_offset-1))
            emission_angle = unpack('>f', label.read(4))[0]
            #Incidence Angle
            label.seek(ancillary_offset + (n*row_bytes-1) + (incidence_offset-1))
            incidence_angle = unpack('>f', label.read(4))[0]
            #Phase Angle
            label.seek(ancillary_offset + (n*row_bytes-1) + (phase_offset-1))
            phase_angle = unpack('>f', label.read(4))[0]
            angles.append([incidence_angle, emission_angle,  phase_angle])
        angles = np.asarray(angles)

        if save == True:
            #print wv_array.shape, ref_array.T.shape
            #print np.concatenate((np.reshape(wv_array,(296,1)), ref_array.T), axis=1)
            np.savetxt('reflectance_csv.txt', np.concatenate((np.reshape(wv_array,(296,1)), ref_array.T), axis=1), delimiter=',')

        return wv_array, rad_array, ref_array, angles

        # openspc end ============================================================================================================


    def clean_data(self, array):
        """
        Parameters:

        array          type: ndarray
                       This is an array that needs to be cleaned, i.e. the data spike around 1 micron us removed and the data is clipped at 1.7 microns.

        Returns:
        cleaned array  type: ndarray
                       A 161 element array (indexed from 0 to 160).  If this is relfectance or radiance, it contains one row per observation.
        """
        try:
            #If this does not fail we have a multi-dimensional array that is either reflectance or radiance
            array.shape[1]
            mask = np.asarray(np.concatenate((np.ones(61),np.zeros(23),np.ones(212))), dtype = bool)
            array_out = array[:,mask]
            return array_out
        except:
            #We have the wavelength array
            array_out = np.delete(array, range(61,84))
            return array_out

    #def getbandnumbers(wavelengths, *args):
    def getbandnumbers(self):
        '''
        This parses the wavelength list,finds the mean wavelength closest to the
        provided wavelength, and returns the index of that value.  One (1) is added
        to the index to grab the correct band.

        Parameters
        ----------
        wavelengths: A list of wavelengths, 0 based indexing
        *args: A variable number of input wavelengths to map to bands

        Returns
        -------
        bands: A variable length list of bands.  These are in the same order they are
        provided in.  Beware that altering the order will cause unexpected results.

        '''
        bands = []  # initialize bands list

        shoulders = [self.leftShoulder, self.rightShoulder]

        # logging.debug("getbandnumbers : self.wv_array = ")
        # for ii in range(len(self.wv_array)):
        #     logging.debug("ii = %s, wv_array = %s",ii,self.wv_array[ii])

        for wv in shoulders:
            # min(iterable[,key=func]) -> value
            logging.debug("min comp : wv = %s", wv)
            bands.append(min(range(len(self.wv_array)), key=lambda i: abs(self.wv_array[i]-wv)))

        #bands = [40, 144]
        logging.debug("getbandnumbers : bands = %s", bands)

        return bands

    def parse_coefficients(self, coefficient_table):
        '''
        Parameters
        ----------

        coefficient_table     type: file path
                              The CSV file to be parsed

        Returns
        -------
        supplemental          type: list of lists
                              List of coefficients where index is the sequentially increasing wavelength.  This data is 'cleaned'.  The r_{mean} at 1003.6 is set to -999, a NoDataValue.
        '''
        d = open(coefficient_table)
        supplemental = []
        for line in d:
            line = line.split(",")
            supplemental.append([float(s) for s in line[1:]])

        return supplemental

    #def photometric_correction(self, wv, ref_vec,coefficient_table, angles, xl_fixed,c1,c2,c3):
    def photometric_correction(self, wv, ref_vec, coefficient_table, xl_fixed, c1, c2, c3):


        '''
        TODO: Docs here
        This function performs the photometric correction.
        '''
        # incidence_angle = angles[:,0]
        # emission_angle = angles[:,1]
        # phase_angle = angles[:,2]



        def _phg(g, phase_angle):
            '''This function allows positive and neg. g to be passed in'''
            phg = (1.0-g**2) / (1.0+g**2-2.0*g*np.cos(np.radians(phase_angle))**(1.5))
            return phg

        #The ref_array runs to the detector limit, but the coefficient table truncates at 1652.1, we therefore only correct the wavelengths that we know the coefficents for.
        #Column  = ref_array[:,wv]
        b_naught = coefficient_table[wv][0]
        h = coefficient_table[wv][1]
        c = coefficient_table[wv][2]
        g = coefficient_table[wv][3]

        #Compute the phase function with fixed values
        p = ((1-c)/2) * _phg(g,30) + ((1+c)/2) * _phg((-1 * g),30)
        b = b_naught / (1+(np.tan(np.radians(30/2.0))/h))
        f_fixed = (1+b)*p

        #Compute the phase function with the observation phase
        #p = (((1-c)/2) * _phg(g,phase_angle)) + (((1+c)/2)* _phg((-1 * g),phase_angle))
        p = (((1 - c) / 2) * _phg(g, self.phase_angle)) + (((1 + c) / 2) * _phg((-1 * g), self.phase_angle))
        #b = b_naught / (1+(np.tan(np.radians(phase_angle/2.0))/h))
        b = b_naught / (1 + (np.tan(np.radians(self.phase_angle / 2.0)) / h))
        f_observed = (1+b)*p

        f_ratio = f_fixed / f_observed

        #Compute the lunar lambert function
        #l = 1.0 + (c1*phase_angle) + (c2*phase_angle**2) + (c3*phase_angle**3)
        l = 1.0 + (c1 * self.phase_angle) + (c2 * self.phase_angle ** 2) + (c3 * self.phase_angle ** 3)

        #cosi = np.cos(np.radians(incidence_angle))
        cosi = np.cos(np.radians(self.incidence_angle))

        #cose = np.cos(np.radians(emission_angle))
        cose = np.cos(np.radians(self.emission_angle))


        xl_observed = 2 * l * (cosi / (cosi + cose)) + ((1-l)*cosi)
        xl_ratio = xl_fixed / xl_observed

        #Compute the photometrically corrected reflectance
        ref_vec = ref_vec * xl_ratio * f_ratio
        return ref_vec

    def continuum_correction(self, bands, ref_array, obs_id):
        y2 = ref_array[obs_id][bands[1]]
        y1 = ref_array[obs_id][bands[0]]
        wv2 = self.wv_array[bands[1]]  # lrm
        wv1 = self.wv_array[bands[0]]  # lrm

        m = (y2-y1) / (wv2 - wv1)
        b =  y1 - (m * wv1)
        y = m * self.wv_array + b  # lrm

        continuum_corrected_ref_array = ref_array[obs_id] / y
        return continuum_corrected_ref_array, y


    def make_plots(self):
        # Read in the spc file, extract necessary info, and clean the data
        # wv_array, rad_array, ref_array, angles = openspc(args.input_data, args.save)
        #elf.wv_array, rad_array, ref_array, angles = self.openspc(self.input_data, self.save)  # lrm
        self.wv_array = self.spectra.index  # lrm this is ok   lrm
        #ref_array = np.array([self.wv_array, self.spectra.values])  # lrm
        #tmp_array = self.spectra.as_matrix()
        #ref_array = np.array([(1.5, 2, 3), (4, 5, 6)])


        print("before setting ref_array **")

        # set obs 0 ref array (only need one obs)
        ref_array = np.array([(self.spectra.values)])


        logging.debug("sp_extract : make_plots :")
        logging.debug("sp_extract : make_plots : self.spectra.index = %s", self.spectra.index)
        logging.debug("sp_extract : make_plots : self.wv_array.shape) = %s", self.wv_array.shape)
        logging.debug("sp_extract : make_plots : self.wv_array = %s", self.wv_array)
        logging.debug(" ")
        logging.debug("sp_extract : make_plots : ref_array.shape) = %s", ref_array.shape)
        #logging.debug("sp_extract : make_plots : ref_array[0,:] = %s", ref_array[0,:] )
        logging.debug("sp_extract : make_plots : ref_array = %s", ref_array)

        logging.debug("make_plot : before cleaning wv_array, len(self.wv_array) = %s", len(self.wv_array))
        self.wv_array = self.clean_data(self.wv_array)  # lrm

        logging.debug("make_plot : after cleaning wv_array, len(self.wv_array) = %s", len(self.wv_array))

        ref_array = self.clean_data(ref_array)
        logging.debug("sp_extract : make_plots : After cleaning ref_array : len(ref_array[0]) = %s", len(ref_array[0]))
        logging.debug("sp_extract : make_plots : ref_array = %s", ref_array[0])

        minwv = int(self.wv_MinLimits)
        maxwv = int(self.wv_MaxLimits)

        #extent = np.where(self.wv_array <= maxwv)
        # & will give you an elementwise and (the parentheses are necessary)
        #extent = np.where((self.wv_array <= maxwv) & (self.wv_array >= minwv))
        # extent = np.where((self.wv_array <= maxwv) & (self.wv_array >= minwv) &
        #                   (input_refarray[0] >= minValid))



        #Parse the supplemental table to get photometric correction coefficients
        coefficient_table = self.parse_coefficients(self.albedo_tab)


        #Perform the photometric correction on the reflectance values
        #Compute the 'static' lunar lambert function for r(30,0,30).
        c1 = -0.019
        c2 = 0.000242
        c3 = -.00000146
        l = 1.0 + (c1*(30)) + (c2*(30**2)) + (c3*(30**3))
        cosi = math.cos(math.radians(30))
        cose = math.cos(math.radians(0))
        xl_constant = ((2*l*(cosi /(cosi + cose)))) + ((1-l)*cosi)

        #Copy the unphotometrically corrected array
        input_refarray = np.copy(ref_array)

        minValid = .0001
        extent = np.where( (self.wv_array <= maxwv) & (self.wv_array >= minwv) &
                           (input_refarray[0] >= minValid) )

        logging.debug("sp_extract : make_plots : len(extent[0]) = %s", len(extent[0]))
        logging.debug("sp_extract : make_plots : extent = %s", extent)
        logging.debug("sp_extract : make_plots  : input_refarray[0][extent] = %s", input_refarray[0][extent])

        # logging.debug("sp_extract : AFTER CHECKING VALUES OF input_refarray[0][validExtent]")
        # logging.debug("sp_extract : make_plots : len(validExtent[0]) = %s", len(validExtent[0]))
        # logging.debug("sp_extract : make_plots : validExtent = %s", validExtent)
        # logging.debug("sp_extract : make_plots : validExtent[0] = %s", validExtent[0])
        #
        #
        # logging.debug("sp_extract : make_plots  : len(input_refarray[0][extent]) = %s", len(input_refarray[0][extent]))
        # logging.debug("sp_extract : make_plots  : input_refarray[0][extent] = %s", input_refarray[0][extent])
        # logging.debug("sp_extract : make_plots  : min(input_refarray[0][extent]) = %s", min(input_refarray[0][extent]))
        # logging.debug("sp_extract : make_plots  : max(input_refarray[0][extent]) = %s", max(input_refarray[0][extent]))

        logging.debug("sp_extract : make_plots : len(coefficient_table) = %s", len(coefficient_table))

        #Perform the photometric correction
        for wv in range(len(coefficient_table)):
            # ref_array[:,wv] = self.photometric_correction(wv, ref_array[:,wv], coefficient_table,
            #                                               self.angles, xl_constant,c1,c2,c3)
            ref_array[:, wv] = self.photometric_correction(wv, ref_array[:, wv], coefficient_table,
                                                           xl_constant, c1, c2, c3)


        #Copy the photometrically corrected array
        photometrically_corrected_ref_array = np.copy(ref_array)
        continuum_slope_array = np.empty(ref_array.shape)

        #Continuum correction
        #bands = getbandnumbers(wv_array, 752.8,1547.7)
        #bands = self.getbandnumbers(self.leftShoulder, self.rightShoulder)  # lrm fix this ***
        bands = self.getbandnumbers() # lrm

        # Get the left-shoulder -> right-shoulder extent
        leftRightExtent = range(bands[0],bands[-1]+1)  # want to include band[-1] in the extent
        logging.debug("leftRightExtent : %s",leftRightExtent)

        #Continuum correct all observations
        for obs_id in range(len(ref_array)):
            ref_array[obs_id],continuum_slope_array[obs_id] = self.continuum_correction(bands, ref_array, obs_id)

        #TODO If the save flag is true we could save out all observations to CSV

        sfont = {'fontname': 'serif'}


        for obs in range(len(self.observation)):  # lrm
            #Do the plotting

            # Determine which plot this is (probably not the best way to do this!!! :
            if (self.wv_MinLimits > 1600.):
                fig = plt.figure(2, figsize=(10, 10))
                mainTitle = '2$\mu$m Observation:  '  + str(self.obsId) # + '  File: ' + self.fileName
            else:
                fig = plt.figure(1, figsize=(10, 10))
                mainTitle = '1$\mu$m Observation: ' # + str(self.obsId) # + '   File: ' + self.fileName

            # Display name of file being plotted
            fileTitle = '  File: ' + self.fileName
            figtext(.5, .95, fileTitle, fontsize=11, ha='center', **sfont)

            # Main plot title
            fig.suptitle(mainTitle, fontsize=13, **sfont)

            # Display which albedo was used
            figtext(.5, .93, 'Albedo : '+self.albedoTitle, fontsize=11, ha='center', **sfont)
            #figtext(.5, .85, 'Lorem ipsum dolor sit amet, consectetur adipiscing elit', fontsize=10, ha='center')

            # Set the horizontal spacing
            fig.subplots_adjust(hspace=0.75)

            ax1 = subplot(411)
            grid(alpha=.5)
            plot(self.wv_array[extent],input_refarray[obs][extent], linewidth=1.5)
            xlabel('Wavelength ($\eta$m)', fontsize=9)
            ax1.set_xticks(self.wv_array[extent][::4])
            ax1.set_xticklabels(self.wv_array[extent][::4], rotation=45, fontsize=8)
            ax1.set_xlim(self.wv_array[extent].min()-10, self.wv_array[extent].max()+10)
            ylabel('Reflectance', fontsize=9)
            #ax1.set_yticklabels(input_refarray[obs][extent],fontsize=8)
            title('Level 2B2 Data', fontsize=12, **sfont)

            ax2 = subplot(412)
            grid(alpha=.5)
            plot(self.wv_array[extent],photometrically_corrected_ref_array[obs][extent], linewidth=1.5)
            xlabel('Wavelength ($\eta$m)', fontsize=9)
            ax2.set_xticks(self.wv_array[extent][::4])
            ax2.set_xticklabels(self.wv_array[extent][::4], rotation=45, fontsize=8)
            ax2.set_xlim(self.wv_array[extent].min()-10, self.wv_array[extent].max()+10)
            ylabel('Reflectance', fontsize=9)
            #ax2.set_yticklabels(input_refarray[obs][extent],fontsize=8)
            title('Photometrically Corrected Data', fontsize=12, **sfont)

            ax3 = subplot(413)
            grid(alpha=.5)
            plot(self.wv_array[extent],photometrically_corrected_ref_array[obs][extent], label='Photometrically Corrected Spectrum', linewidth=1.5)
            plot(self.wv_array[extent], continuum_slope_array[obs][extent],'r--', label='Spectral Continuum', linewidth=1.5)
            xlabel('Wavelength ($\eta$m)', fontsize=9)
            ax3.set_xticks(self.wv_array[extent][::4])
            ax3.set_xticklabels(self.wv_array[extent][::4], rotation=45, fontsize=8)
            ax3.set_xlim(self.wv_array[extent].min()-10, self.wv_array[extent].max()+10)
            ylabel('Reflectance', fontsize=9)
            #ax3.set_yticklabels(input_refarray[obs][extent],fontsize=8)
            contTitle = 'Spectral Continuum  Left: ' + str(self.leftShoulder) + ' Right: ' + str(self.rightShoulder)
            #title('Continuum Slope', fontsize=12)
            title(contTitle, fontsize=12, **sfont)

            # Calculate the area under the curves
            areaContinuum = np.trapz(continuum_slope_array[obs][leftRightExtent], x=self.wv_array[leftRightExtent])
            logging.debug("Area under continuum slope curve = %s", areaContinuum)

            areaCorrected = np.trapz(photometrically_corrected_ref_array[obs][leftRightExtent], x=self.wv_array[leftRightExtent])
            logging.debug("Area under photometrically corrected curve = %s", areaCorrected)

            # Take difference of area under the 2 curves
            self.bandArea = areaContinuum - areaCorrected
            logging.debug("self.bandArea = %s", self.bandArea)


            ax4 = subplot(414)
            grid(alpha=.5)
            plot(self.wv_array[extent], ref_array[obs][extent], linewidth=1.5)
            xlabel(r'Wavelength ($\eta$m)', fontsize=9)
            ax4.set_xticks(self.wv_array[extent][::4])
            ax4.set_xticklabels(self.wv_array[extent][::4], rotation=45, fontsize=8)
            ax4.set_xlim(self.wv_array[extent].min()-10, self.wv_array[extent].max()+10)
            ylabel('Relative Reflectance', fontsize=9)
            #ax4.set_yticklabels(input_refarray[obs][extent],fontsize=8)
            title('Continuum Removed Spectrum', fontsize=12, **sfont)



            draw()

        show()

    def get_bandArea(self):
        return self.bandArea
