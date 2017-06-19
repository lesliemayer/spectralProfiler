import os
from PyQt4 import QtCore, QtGui, uic

from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
# Original
#from matplotlib.backends.backend_qt4agg import NavigationToolbar2QTAgg as NavigationToolbar
# L. Mayer 04/28/2017 :
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT as NavigationToolbar

from matplotlib.figure import Figure
import pandas as pd


from pysat.spectral import continuum
from pysat.spectral import smoothing
from pysat.spectral import analytics

base = os.path.dirname(os.path.abspath(__file__))
ui_dialog, cls_dialog = uic.loadUiType(os.path.join(base, 'ui/plot_dialog.ui'))

DEFAULT_MARKER_SIZE = 10
SELECTED_MARKER_SIZE = 15

class PlotDialog(cls_dialog, ui_dialog):

    def __init__(self, *args, **kwargs):
        super(PlotDialog, self).__init__(*args, **kwargs)
        self.setupUi(self)

        self.initmpl()
        self.inittree()
        
        #self.initgui()  # not needed - lrm
        self.show()

        self.selected_line = None
        self.annotations = []
        self.smoother_lookup = {'Box Filter': smoothing.boxcar,
                                'Gaussian': smoothing.gaussian}

    # def initgui(self):
    #
    #     # lrm
    #     # self.endpoint_slider_lower.valueChanged.connect(lambda i: self.endpoint_slider_lower_label.setText('{'
    #     #                                                                                                    '}'.format(self.wavelengths[i])))
    #
    #     # lrm
    #     # self.endpoint_slider_upper.valueChanged.connect(lambda i: self.endpoint_slider_upper_label.setText('{'
    #     #                                                                                                    '}'.format(self.wavelengths[i])))
    #
    #     # lrm
    #     # self.band_center_btn.clicked.connect(self.band_center)
    #
    #     # lrm
    #     #self.band_minima_btn.clicked.connect(self.band_minima)
    #
    #     # lrm
    #     #self.band_area_btn.clicked.connect(self.band_area)
    #
    #     # lrm
    #     #self.band_asymmetry_btn.clicked.connect(self.band_asymmetry)
    #
    #     pass  # lrm

    # lrm
    # def get_endpoints_and_spectra(self):
    #
    #     if self.selected_line:
    #         lower = float(self.endpoint_slider_lower_label.text())
    #         upper = float(self.endpoint_slider_upper_label.text())
    #
    #         l = self.selected_line
    #         xd = l.get_xdata()
    #         yd = pd.Series(l.get_ydata(), index=xd)
    #         offset = l.offset
    #         return lower, upper, yd, offset

    def band_center(self):

        # lrm
        pass
        # if self.selected_line is not None:
        #     lower, upper, yd, offset = self.get_endpoints_and_spectra()
        #     (min_idx, min_value), center_fit = analytics.band_center(yd, low_endmember=lower, high_endmember=upper)
        #     center_fit.plot(ax=self.ax, color='k', gid=4, picker=5, linewidth=4.0, alpha=0.5)
        #     self.ax.plot(min_idx, min_value, marker='*', markersize=DEFAULT_MARKER_SIZE, color='k', gid=3, picker=5)
        #     self.canvas.draw()
        #
        #     self.notes.insertPlainText('> The band center between {} and {} is {} at {}\n'.format(lower,
        #                                                                           upper,
        #                                                                               min_value,
        #                                                                               min_idx))

    def band_minima(self):

        # lrm
        pass
        # if self.selected_line is not None:
        #     lower, upper, yd, offset = self.get_endpoints_and_spectra()
        #     min_idx, min_value = analytics.band_minima(yd, low_endmember=lower, high_endmember=upper)
        #     self.ax.plot(min_idx, min_value, marker='*', markersize=DEFAULT_MARKER_SIZE, color='k', gid=3, picker=5)
        #     self.canvas.draw()
        #
        #     self.notes.insertPlainText('> The band minima between {} and {} is {} at {}\n'.format(lower,
        #                                                                                           upper,
        #                                                                               min_value,
        #                                                                               min_idx))

    def band_area(self):

        # lrm
        pass
        # if self.selected_line is not None:
        #     lower, upper, yd, offset = self.get_endpoints_and_spectra()
        #     area = analytics.band_area(yd - offset, low_endmember=lower, high_endmember=upper)
        #     fill = yd.loc[lower:upper]
        #     self.ax.fill_between(fill.index, fill.values, 1.0 + offset, facecolor='b', alpha=0.25, gid=6, picker=5)
        #     self.canvas.draw()
        #
        #     self.notes.insertPlainText('> The band area between {} and {} is {}\n'.format(lower,
        #                                                                                    upper,
        #                                                                                    area))


    # lrm
    # def band_asymmetry(self):
    #     if self.selected_line is not None:
    #         lower, upper, yd, offset = self.get_endpoints_and_spectra()
    #         self.band_area()
    #         asymmetry = analytics.band_asymmetry(yd - offset, low_endmember=lower, high_endmember=upper)
    #         self.canvas.draw()
    #
    #         self.notes.insertPlainText('> The band asymmetry between {} and {} is {}\n'.format(lower,
    #                                                                                           upper,
    #                                                                                           asymmetry))

    
    # lrm : set up the plotting window
    def initmpl(self):
        """
        Initialize the MatPlotLib Figure
        """
        self.figure = Figure()
        self.ax = self.figure.add_subplot(111)

        # Add a title - lrm
        self.ax.set_title("This is the title")
        self.ax.grid(True)  # turn on the grid
        self.canvas = FigureCanvas(self.figure)

        # Add an x label
        self.ax.set_xlabel('Wavelength (nm)')

        # Add a y label - lrm
        self.ax.set_ylabel('Reflectance')

        # lrm : add the plot
        self.mplvl.addWidget(self.canvas)

        # Add the bottom toolbar
        self.toolbar = NavigationToolbar(self.canvas, self.mplwindow, coordinates=True)
        self.mplvl.addWidget(self.toolbar)

        # lrm : select a particular spectra
        #self.figure.canvas.mpl_connect('pick_event', self.select_spectra)

        # lrm : Draw everything
        self.canvas.draw()

    def inittree(self):  # is this the observation tree view in box on the right?
        """
        Initialize the tree view
        """
        self.spectratree.setModel(QtGui.QStandardItemModel())

    # lrm : get the spectrum data??
    def set_spectra(self, spectra):

        self.data = spectra
        for k, obs in self.data.iteritems():
            parent = QtGui.QStandardItem(k)
            for label, values in obs.iteritems():
                child = QtGui.QStandardItem('Observation {}'.format(label))
                child.setFlags(QtCore.Qt.ItemIsEnabled | QtCore.Qt.ItemIsSelectable)
                parent.appendRow(child)
            self.spectratree.model().appendRow(parent)

    # lrm : not needed ----------------------------------------------------------
    # def select_spectra(self, event):
    #     def update_selection():
    #         for i, x in enumerate(self.ax.lines):
    #             if x != self.selected_line:
    #                 if x.get_gid() in [0,1]:
    #                     x.set_linewidth(1.0)
    #                     x.set_markersize(2.0)
    #                 elif x.get_gid() in [4, 5]:
    #                     x.set_linewidth(4.0)
    #                     x.set_markersize(DEFAULT_MARKER_SIZE)
    #                 elif x.get_gid() in [6]:
    #                     x.set_alpha(0.25)
    #
    #     geom = event.artist
    #     gid = geom.get_gid()
    #     self.selected_line = geom
    #     update_selection()
    #
    #     try:
    #         lineidx = self.ax.lines.index(geom)
    #     except:
    #         return
    #
    #     if event.mouseevent.button == 1:
    #         if gid in [0, 1]:
    #             self.ax.lines[lineidx].set_linewidth(2.0)
    #             self.ax.lines[lineidx].set_markersize(4.0)
    #         elif gid in [4, 5]:
    #             self.ax.lines[lineidx].set_linewidth(6.0)
    #             self.ax.lines[lineidx].set_markersize(SELECTED_MARKER_SIZE)
    #         elif gid in [6]:
    #             self.ax.lines[lineidx].set_alpha(0.75)
    #
    #     elif event.mouseevent.button == 3:
    #         self.spectra_context_menu_event(self)
    #
    #     self.canvas.draw()
    # --------------------------------------------------------------------------------

    # lrm commented out
    # def spectra_context_menu_event(self, event):
    #     menu = QtGui.QMenu()
    #     try:
    #         gid = self.selected_line.get_gid()
    #     except: return
    #
    #     if gid in [0, 5, 10]:
    #         delsp = menu.addAction('Delete Selected Spectra', self.delete_spectra)
    #         movesp = menu.addAction('Move Spectra', self.move_spectra)
    #         aspoint = menu.addAction('Show as Point', self.as_point)
    #         asline = menu.addAction('Show as Line', self.as_line)
    #         addsmooth = menu.addAction('Add Smoothed Spectra', self.smooth_spectra)
    #         ccorrect = menu.addAction('Continuum Correct', self.continuum_correct_spectra)
    #
    #     if gid in [4]:
    #         delsp = menu.addAction('Delete Parameter', self.delete_spectra)
    #
    #     menu.exec_(self.mapToGlobal(event.pos()))

    # lrm
    # def delete_spectra(self):
    #     """
    #     Remove the selected spectra from the plot.
    #     """
    #     if self.selected_line is not None:
    #         self.ax.lines.remove(self.selected_line)
    #         self.canvas.draw()

    # lrm
    # def move_spectra(self):
    #     if self.selected_line is not None:
    #         l = self.selected_line
    #         offset = float(self.offset.value())
    #         yd = l.get_ydata() + offset
    #         self.ax.plot(l.get_xdata(),
    #                      yd,
    #                      color=l.get_color(),
    #                      picker=l.get_picker(),
    #                      label=l.get_label())
    #         self.ax.lines.remove(l)
    #         self.canvas.draw()

    # lrm
    # def as_point(self):
    #     """
    #     Convert the point to a line object
    #     """
    #     if self.selected_line is not None:
    #         l = self.selected_line
    #         l.set_linestyle('')
    #         l.set_marker('o')
    #         l.set_markersize(2)
    #         self.canvas.draw()

    # def as_line(self):
    #     """
    #     Convert the line to a point object
    #     """
    #     if self.selected_line is not None:
    #         l = self.selected_line
    #         l.set_linestyle('-')
    #         l.set_marker('')
    #         self.canvas.draw()

    # def smooth_spectra(self):
    #     if self.selected_line is not None:
    #         l = self.selected_line
    #         xd = l.get_xdata()
    #         yd = pd.Series(l.get_ydata(), index=xd)
    #         color = l.get_color()
    #         smooth_func = self.smoother_lookup[self.smooth_method.currentText()]
    #         window = int(self.smooth_window.value())
    #
    #         smoothed = smooth_func(yd, window_size=window)
    #
    #         self.ax.plot(xd, smoothed, color=color, picker=5, gid=5)

    def continuum_correct_spectra(self):
        pass

    def plot(self, continuum_endpoints, correction_method, smoothing_method, smoothing_window_size, offset,
             clipping_lower, clipping_upper, pcorrect=None):

        if pcorrect == 'Mare':
            key = 'REF2'
        elif pcorrect == 'Highlands':
            key = 'REF1'
        else:
            key = 'REF'

        offset_interval = 0

        for fname, panel in self.data.iteritems():
            for k, df in panel.iteritems():
                if key in df.columns:
                    spectra = df[key]
                else:
                    spectra = df['REF']

                # lrm

                #print "plot_dialog : plot : spectra = {}".format(spectra)

                # lrm - don't do any correction
                # if correction_method != 'None':
                #     df['CC'], df['Continuum'] = continuum.continuum_correct(spectra,
                #                                                             nodes=continuum_endpoints,
                #                                                             method=correction_method.lower())
                # else:
                #    df['CC'] = spectra
                df['CC'] = spectra

                # lrm : don't do the mask
                # mask = (df.index >= clipping_lower) & (df.index <= clipping_upper)
                # spectra = df['CC'][mask]

                # lrm : don't do the offset
                # if offset:
                #     spectra += (offset + offset_interval)
                #     offset_interval += offset

                if smoothing_method != 'None':
                    smooth_func = self.smoother_lookup[smoothing_method]
                    spectra = smooth_func(spectra, window_size=smoothing_window_size)

                print("plot_dialog : plot : Just before spectra.plot")
                spectra.plot(ax=self.ax, picker=5, gid=0)

                # Custom attr to get the offset as an attribute to the line.
                setattr(self.ax.lines[-1], 'offset', (offset + offset_interval))
        self.ax.grid(True)
        self.canvas.draw()

