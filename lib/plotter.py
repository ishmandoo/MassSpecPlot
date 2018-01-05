from lib.data import Data, AuxPlots
import numpy as np
from scipy.ndimage.filters import gaussian_filter
from enum import Enum
import matplotlib.pyplot as plt
from matplotlib import animation
from matplotlib.patches import Rectangle
import os
from scipy import signal

class SpecNorm():
	GLOBAL = 0
	SCAN = 1
	MANUAL = 2
	LOCAL = 3

class Plotter:
	def __init__(self, data):
		self.data = data

	def makeAuxData(self, scan_range, aux_plot_type, aux_smoothing):
		scan_start, scan_end = scan_range

		if aux_plot_type == AuxPlots.SOURCE_CURRENT:
			data = aux_plot_type['multiplier'] * self.data.currents

		elif aux_plot_type == AuxPlots.DETECTOR_CURRENT:
			data = aux_plot_type['multiplier'] * self.data.currents

		elif aux_plot_type == AuxPlots.DETECTOR_SOURCE_RATIO:
			data = aux_plot_type['multiplier'] * self.data.intensities / self.data.currents

		elif aux_plot_type == AuxPlots.L1_VOLTAGE:
			data = aux_plot_type['multiplier'] * self.data.voltagesl1

		elif aux_plot_type == AuxPlots.L2_VOLTAGE:
			data = aux_plot_type['multiplier'] * self.data.voltagesl2

		elif aux_plot_type == AuxPlots.PRESSURE:
			data = aux_plot_type['multiplier'] * self.data.pressures

		else:
			raise ValueError('Invalid Aux Plot Type')

		smooth_data = gaussian_filter(data, aux_smoothing)
		return smooth_data

	def makeScanRanges(self, scan_range, window, step):
		scan_start, scan_end = scan_range

		return [(i, i+window) for i in range(scan_start, scan_end-window, step)]


	def findPeaks(self, masses, intensities):
		# find peaks
		#peak_indices = peakutils.indexes(intensities, thres=0.06)
		peak_indices = signal.find_peaks_cwt(intensities, np.arange(2,15))
		peak_masses = [masses[i] for i in peak_indices]
		peak_intensities = [intensities[i] for i in peak_indices]
		return peak_masses, peak_intensities

	def makeAnimation(self,  
		window = 500, 
		step = 500,
		mass_range = (None, None),
		scan_range = (None, None),
		out_name = None,
		normalization=SpecNorm.SCAN,
		markers = [],
		aux_plot_type=AuxPlots.SOURCE_CURRENT,
		aux_plot_type_2=None,
		aux_smoothing = 1,
		aux_smoothing_2 = 1,
		aux_range = (None,None),
		aux_range_2 = (None,None),
		manual_norm_range=(0,0),
		show_plot = True,
		label_peaks = True,
		local_norm_scan_range=(0,0)):

		scan_start, scan_end = scan_range
		if scan_start is None:
			scan_start = 0
		if scan_end is None:
			scan_end = self.data.scans.shape[0]

		scan_range = (scan_start, scan_end)


		mass_start, mass_end = mass_range
		if mass_start is None:
			mass_start = self.data.mass_start
		if mass_end is None:
			mass_end = self.data.mass_end

		mass_range = (mass_start, mass_end)


		ranges = self.makeScanRanges(scan_range, window, step)

		aux_data = self.data.aux_data[aux_plot_type['key']][scan_start:scan_end]

		fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10))
		plt.tight_layout(pad=2, w_pad=1.8, h_pad=5, rect=[0.05, 0, 0.85, 1])


		ln, pk, pkTxt, infoTxt = self.initSpecPlot(ax1, mass_range, normalization, manual_norm_range, local_norm_scan_range)
		aux, rect = self.initAuxPlot(ax2, aux_data, aux_range, scan_range, markers, aux_plot_type, "black")

		ax3 = None
		aux_2 = None
		if aux_plot_type_2:
			ax3 = ax2.twinx()
			aux_data_2 = self.data.aux_data[aux_plot_type['key']][scan_start:scan_end]
			aux_2, rect_2 = self.initAuxPlot(ax3, aux_data_2, aux_range_2, scan_range, markers, aux_plot_type_2, 'blue')
			rect_2.set_alpha(0)

		# a function to update the plot for each frame
		def update(frame):
			self.updateSpecPlot(ax1, mass_range, ranges[frame], ln, pk, pkTxt, label_peaks, normalization)
			self.updateAuxPlot(ax2, aux, infoTxt, rect, aux_data, ranges[frame], scan_start)
			if aux_plot_type_2:
				self.updateAuxPlot(ax3, aux_2, infoTxt, rect, aux_data_2, source_currents, detector_currents, ratios, L1_voltages, L2_voltages, ranges[frame], scan_start)
			return [ln, pk, aux, rect] + pkTxt




		ani = animation.FuncAnimation(fig, update, frames=range(len(ranges)), blit=False)
		#Writer = writers['ffmpeg']
		#FFWriter = Writer(fps=15)
		#writer = animation.AVConvFileWriter()
		if not out_name == None:
			ani.save(os.path.join(self.data.path, out_name))


		# make an alert noise
		print('\a')

		if show_plot:
			plt.show()

	def initSpecPlot(self, ax, mass_range, normalization, manual_norm_range, local_norm_scan_range):
		mass_start, mass_end = mass_range

		# the size of the pool of text objects for displaying peak masses
		peak_text_pool_size = 100

		ln, = ax.plot([], [], 'k-', animated=False, lw=1)
		pk, = ax.plot([], [], 'r.', animated=False)
		pkTxt = [ax.text(0, 0, '', verticalalignment='bottom', horizontalalignment='center', color='red', rotation=45) for i in range(peak_text_pool_size)]
		infoTxt = ax.text(1.02, 0.95, "", verticalalignment='top', horizontalalignment='left', transform=ax.transAxes)


		ax.set_xlim(mass_start - 0.05 * mass_end, 1.05 * mass_end)
		if normalization == SpecNorm.GLOBAL:
			ax.set_ylim(0, 1.2 * np.max(self.data.scans))

		elif normalization == SpecNorm.MANUAL:
			ax.set_ylim(manual_norm_range[0], manual_norm_range[1])

		elif normalization == SpecNorm.LOCAL:
			scan_start, scan_end = local_norm_scan_range
			ax.set_ylim(0, 1.2 * np.max(self.data.scans[local_norm_scan_range[0]:local_norm_scan_range[1],mass_start:mass_end]))

		ax.set_xlabel("m/q (amu/e)")
		ax.set_ylabel("intensity")
		ax.get_yaxis().set_ticks([])

		return ln, pk, pkTxt, infoTxt

	def initAuxPlot(self, ax, aux_data, manual_range, scan_range, markers, aux_plot_type, color):
		scan_start, scan_end = scan_range
		manual_start, manual_end = manual_range
		times = self.data.times

		aux, = ax.plot([], [], color=color, animated=False, lw=1)
		ax.tick_params(axis='y', colors=color)
		ax.xaxis.label.set_color(color)

		y_lim_top = manual_end
		y_lim_bot = manual_start

		if y_lim_top is None:
			y_lim_top = 1.1 * np.percentile(aux_data, 99.5)

		if y_lim_bot is None:
			y_lim_bot = np.percentile(aux_data, 0.5) - 0.05 * y_lim_top

		# make the highlight rectangle showing the scan window
		rect = Rectangle((0,y_lim_bot),0, 2*y_lim_top, color=(0,0,0), alpha=0.2)
		ax.add_artist(rect)

		for marker in markers:
			try:
				if scan_start < marker.scan < scan_end:
					time = times[marker.scan - scan_start]
					ax.plot([time, time], [y_lim_bot, y_lim_top], color = "red", dashes=[5,5], lw=1)
					ax.text(time, y_lim_top, marker.title, verticalalignment='bottom', horizontalalignment='left', color='red', rotation=45)
			except IndexError:
				warn("Marker out of plot range")

		# set axis limits
		ax.set_xlim(times[0] - 0.05 * times[-1], 1.05 * times[-1])
		ax.set_ylim(y_lim_bot, y_lim_top)

		ax.set_xlabel("time (s)")
		ax.set_ylabel(aux_plot_type['x'], color=color)

		return aux, rect

	def findPeaks(self, scan):
		# find peaks
		peak_indices = signal.find_peaks_cwt(scan, np.arange(2,15))
		peak_intensities = [scan[i] for i in peak_indices]
		return peak_indices, peak_intensities


	def updateSpecPlot(self, ax, mass_range, window_range, spec, peak, peak_text_pool, label_peaks, normalization):
		window_start, window_end = window_range
		mass_start, mass_end = mass_range

		spectrum = np.sum(self.data.scans[window_start:window_end, mass_start:mass_end],axis=0)

		peak_masses, peak_intensities = self.findPeaks(spectrum)

		# find max peak height for use in setting axis limits
		if len(peak_intensities) > 0 and max(peak_intensities) > 0:
			max_intensity = max(peak_intensities)
		else:
			max_intensity = np.max(spectrum)

		if max_intensity == 0:
			max_intensity = 1.
		print(max_intensity)

		# set the y limits for the spectrum
		if normalization == SpecNorm.SCAN:
			ax.set_ylim(0, 1.2 * max_intensity)

		# loop over the peak text pool and update as many of them as necessary
		# set the rest to blank
		if label_peaks:
			for i in range(len(peak_text_pool)):
				if i < len(peak_masses):
					peak_text_pool[i].set_text(str(peak_masses[i]))
					peak_text_pool[i].set_x(peak_masses[i] + mass_start)
					peak_text_pool[i].set_y(peak_intensities[i] + max_intensity * 0.01)

				else:
					peak_text_pool[i].set_text('')
					peak_text_pool[i].set_x(0)
					peak_text_pool[i].set_y(0)


		spec.set_data(range(mass_start, mass_end), spectrum)
		if label_peaks:
			peak.set_data(peak_masses + mass_start, peak_intensities)

	def updateAuxPlot(self, ax, aux, info_text, highlight_rect, data, window_range, scan_start):
		# unpack the frame window start and end
		window_start, window_end = window_range
		times = self.data.times

		# find the window start and end frames relative to the animation start scan
		# use these for indexing into currents or times
		window_start_relative = window_start - scan_start
		window_end_relative = window_end - scan_start

		# find the start and end time for the window
		window_start_time = times[window_start_relative]
		window_end_time = times[window_end_relative]

		# update the highlighted window
		highlight_rect.set_x(window_start_time)
		highlight_rect.set_width(window_end_time - window_start_time)

		aux.set_data(times[:window_end_relative],data[:window_end_relative])
		info_text.set_text("scans %d - %d\n%.0fs - %.0fs\nsource %.2f nA\ndetector %.2e\nratio %.2e\nL1 %d V\nL2 %d V" % (
			window_start, 
			window_end, 
			window_start_time, 
			window_end_time, 
			np.mean(self.data.currents[window_start:window_end]), 
			np.mean(self.data.intensities[window_start:window_end]), 
			np.mean(self.data.ratios[window_start:window_end]), 
			np.mean(self.data.voltagesl1[window_start:window_end]),
			np.mean(self.data.voltagesl2[window_start:window_end])
			))