import numpy as np
import math
import os
import csv
from scipy.io import netcdf
from scipy.ndimage.filters import gaussian_filter
from collections import defaultdict
from scipy import signal
import matplotlib.pyplot as plt
from matplotlib import animation
from matplotlib.patches import Rectangle
from bisect import bisect_left
import glob
from warnings import warn
import re
import time
#plt.rcParams["animation.convert_path"] = u"C:\\Program Files\\ImageMagick-7.0.7-Q16\\magick.exe"
#plt.rcParams["animation.convert_path"] = u"magick"



class SpecNorm:
	GLOBAL = 0
	SCAN = 1
	MANUAL = 2
	LOCAL = 3

class AuxPlots:
	SOURCE_CURRENT = {'x':'current (nA)', 'multiplier': 1000000000.}
	DETECTOR_CURRENT = {'x':'detector counts (kcps)', 'multiplier': 1.}
	DETECTOR_SOURCE_RATIO = {'x':'ratio (kcps/nA)', 'multiplier': .000000001}
	L1_VOLTAGE = {'x':'L1 voltage (V)', 'multiplier': 1000.}
	L2_VOLTAGE = {'x':'L2 voltage (V)', 'multiplier': 1000.}
	PRESSURE = {'x':'pressure (kPa)', 'multiplier': 1.}

class Scan:
	def __init__(self, n, t, intensity, i, vl1, vl2, p, intensity_cutoff=None):
		self.n = n
		self.masses = defaultdict(np.uint64)
		self.t = t
		self.i = i
		self.vl1 = vl1
		self.vl2 = vl2
		self.p = p
		self.intensity = intensity
		if not intensity_cutoff is None:
			sum([masses[mass] for mass in masses if mass > intensity_cutoff])
	def __add__(self, other):
		new_scan = Scan(self.n, self.t, self.i, self.vl1, self.vl2, self.p, self.intensity)
		new_scan.masses = self.masses.copy()
		for mass in other.masses:
			new_scan.masses[mass] += other.masses[mass]
		return new_scan
	def __radd__(self, other):
		if other == 0:
			return self
		return self.__add__(other)
	def __repr__(self):
		return str(self.masses)

class Marker:
	def __init__(self, title, scan):
		self.title = title
		self.scan = scan

class Spectrum:
	def __init__(self, scans = []):
		self.scans = scans
		self.path = None
		self.mass_start = None
		self.mass_end = None
		self.s0 = 0

	def load(self, path, specPath=None, auxPath=None, logPath=None):
		lines = None
		self.path = path
		if specPath is None:
			files = glob.glob(os.path.join(path,"*.cd*"))
			if len(files) > 0:
				specPath = files[0]

		if auxPath is None:
			files = glob.glob(os.path.join(path,"*.tsv"))
			if len(files) > 0:
				auxPath = files[0]

		if logPath is None:
			files = glob.glob(os.path.join(path,"*.log"))
			if len(files) > 0:
				logPath = files[0]

		if specPath is None:
			raise IOError("Missing spectrum file")

		if auxPath is None:
			raise IOError("Missing aux file")

		if logPath is None:
			warn("Missing aux file")

		print("loading\ncdf -> %s\naux -> %s\nlog -> %s"%(specPath, auxPath, logPath))


		if not logPath is None:
			with open(logPath, 'r') as f:
				log = f.read()
				pattern = re.compile("[0-1][0-9]:[0-5][0-9]:[0-5][0-9]")
				self.s0 = pattern.search(log).group(0).split(":")[-1]

		with open(auxPath,'r') as tsv:
		    lines = [line.strip().split('\t') for line in tsv]

		with netcdf.netcdf_file(specPath, 'r', mmap=False) as f:
			v = f.variables

			scan_start_indices = v['scan_index']
			scan_start_times = v['scan_acquisition_time']
			scan_total_intensities = v['total_intensity']

			self.mass_start = int(v['mass_range_min'][0])
			self.mass_end = int(v['mass_range_max'][0])


			masses = v['mass_values']
			intensities = v['intensity_values']

			#(i, vl1, vl2, p, y0, m0, d0, h0, m0, s0) = lines[0]
			#h0, m0, s0 = float(h0), float(m0), float(s0)

			time_stamp = f.experiment_date_time_stamp
			h0, m0, s0 = time_stamp[8:10], time_stamp[10:12], self.s0

			times = [toSeconds((h,m,s),(h0,m0,s0)) for (i, vl1, vl2, p, y, m, d, h, m, s) in lines]
			current = [float(i) for (i, vl1, vl2, p, y, m, d, h, m, s) in lines]
			voltagel1 = [float(vl1) for (i, vl1, vl2, p, y, m, d, h, m, s) in lines]
			voltagel2 = [float(vl2) for (i, vl1, vl2, p, y, m, d, h, m, s) in lines]
			pressure =[float(p) for (i, vl1, vl2, p, y, m, d, h, m, s) in lines]

			# loop over the scans
			for i in range(scan_start_indices.shape[0]-1):
				# get the starting time for the scan
				t = scan_start_times[i]

				# find the index of the line in the aux data file closest to this time
				j = bisect_left(times, t)
				if j == len(times):
					j = len(times) - 1

				# create a new Scan object for each
				new_scan = Scan(i, t, scan_total_intensities[i], current[j], voltagel1[j], voltagel2[j], pressure[j])

				# loop over the masses in the scan
				for point in range(scan_start_indices[i], scan_start_indices[i+1]):
					mass = math.floor(masses[point])
					intensity = intensities[point]
					# if the mass is non-zero, add it to the dictionary, combining if necessary
					if not mass == 0:
						new_scan.masses[mass] += intensity

				self.scans.append(new_scan)


	def currentPlot(self, scan_range):

		scan_start, scan_end = scan_range
		i = [scan.i for scan in self.scans[scan_start:scan_end]]
		t = [scan.t for scan in self.scans[scan_start:scan_end]]
		plt.plot(t, i)
		plt.show()

	def makeScanRanges(self, scan_range, window, step):
		scan_start, scan_end = scan_range

		return [(i, i+window) for i in range(scan_start, scan_end-window, step)]

	def makeSpectrum(self, window_range):
		window_start, window_end = window_range
		return sum(self.scans[window_start:window_end]).masses

	def makeSpectra(self, mass_range, scan_ranges):
		# compile the overall spectrum for each mass window

		mass_start, mass_end = mass_range
		spec_list = [self.makeSpectrum(window_start, window_end) for window_start, window_end in scan_ranges]

		# generate the list of masses
		masses = range(mass_start, mass_end)

		# generate the list of intensities for each mass window
		intensities_list = [[spec[m] for m in masses] for spec in spec_list]

		return masses, intensities_list

	def makeAuxData(self, scan_range, aux_plot_type, aux_smoothing):
		scan_start, scan_end = scan_range

		if aux_plot_type == AuxPlots.SOURCE_CURRENT:
			data = np.array([scan.i * aux_plot_type['multiplier'] for scan in self.scans[scan_start:scan_end]])

		elif aux_plot_type == AuxPlots.DETECTOR_CURRENT:
			data = np.array([scan.intensity * aux_plot_type['multiplier'] for scan in self.scans[scan_start:scan_end]])

		elif aux_plot_type == AuxPlots.DETECTOR_SOURCE_RATIO:
			data = np.array([(scan.intensity / scan.i) * aux_plot_type['multiplier'] for scan in self.scans[scan_start:scan_end]])

		elif aux_plot_type == AuxPlots.L1_VOLTAGE:
			data = np.array([scan.vl1 * aux_plot_type['multiplier'] for scan in self.scans[scan_start:scan_end]])

		elif aux_plot_type == AuxPlots.L2_VOLTAGE:
			data = np.array([scan.vl2 * aux_plot_type['multiplier'] for scan in self.scans[scan_start:scan_end]])

		elif aux_plot_type == AuxPlots.PRESSURE:
			data = np.array([scan.p * aux_plot_type['multiplier'] for scan in self.scans[scan_start:scan_end]])

		else:
			raise ValueError('Invalid Aux Plot Type')

		smooth_data = gaussian_filter(data, aux_smoothing)
		times = np.array([scan.t for scan in self.scans[scan_start:scan_end]])
		return times, smooth_data

	def findPeaks(self, masses, intensities):
		# find peaks
		#peak_indices = peakutils.indexes(intensities, thres=0.06)
		peak_indices = signal.find_peaks_cwt(intensities, np.arange(2,15))
		peak_masses = [masses[i] for i in peak_indices]
		peak_intensities = [intensities[i] for i in peak_indices]
		return peak_masses, peak_intensities

	def updateSpecPlot(self, ax, spec, peak, masses, intensities, peak_text_pool, normalization):
		peak_masses, peak_intensities = self.findPeaks(masses, intensities)

		# find max peak height for use in setting axis limits
		if len(peak_intensities) > 0 and max(peak_intensities) > 0:
			max_intensity = max(peak_intensities)
		else:
			 max_intensity = max(intensities)

		if max_intensity == 0:
			max_intensity = 1.

		# set the y limits for the spectrum
		if normalization == SpecNorm.SCAN:
			ax.set_ylim(0, 1.2 * max_intensity)

		# loop over the peak text pool and update as many of them as necessary
		# set the rest to blank
		for i in range(len(peak_text_pool)):
			if i < len(peak_masses):
				peak_text_pool[i].set_text(str(peak_masses[i]))
				peak_text_pool[i].set_x(peak_masses[i])
				peak_text_pool[i].set_y(peak_intensities[i] + max_intensity * 0.01)

			else:
				peak_text_pool[i].set_text('')
				peak_text_pool[i].set_x(0)
				peak_text_pool[i].set_y(0)


		spec.set_data(masses, intensities)
		peak.set_data(peak_masses, peak_intensities)

	def updateAuxPlot(self, ax, aux, info_text, highlight_rect, times, data, source_currents, detector_currents, ratios, L1_voltages, L2_voltages, scan_range, scan_start):
		# unpack the frame window start and end
		window_start, window_end = scan_range

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
		info_text.set_text("scans %d - %d\n%.0fs - %.0fs\nsource %.2f nA\ndetector %.0f\nratio %.2f\nL1 %.2f V\nL2 %.2f V" % (
			window_start, 
			window_end, 
			window_start_time, 
			window_end_time, 
			source_currents[window_start_relative], 
			detector_currents[window_start_relative], 
			ratios[window_start_relative], 
			L1_voltages[window_start_relative], 
			L2_voltages[window_start_relative]
			))

	def initSpecPlot(self, ax, intensities_list, mass_range, normalization, manual_norm_range, local_norm_scan_range):
		mass_start, mass_end = mass_range

		# the size of the pool of text objects for displaying peak masses
		peak_text_pool_size = 100

		ln, = ax.plot([], [], 'k-', animated=False, lw=1)
		pk, = ax.plot([], [], 'r.', animated=False)
		pkTxt = [ax.text(0, 0, '', verticalalignment='bottom', horizontalalignment='center', color='red', rotation=45) for i in range(peak_text_pool_size)]
		infoTxt = ax.text(1.02, 0.95, "", verticalalignment='top', horizontalalignment='left', transform=ax.transAxes)


		ax.set_xlim(mass_start - 0.05 * mass_end, 1.05 * mass_end)
		if normalization == SpecNorm.GLOBAL:
			ax.set_ylim(0, 1.2 * max([max(intensities) for intensities in intensities_list]))

		elif normalization == SpecNorm.MANUAL:
			ax.set_ylim(manual_norm_range[0], manual_norm_range[1])

		elif normalization == SpecNorm.LOCAL:
			scan_start, scan_end = local_norm_scan_range
			ax.set_ylim(0, 1.2 * max([max(intensities) for intensities in intensities_list[local_norm_scan_range[0]:local_norm_scan_range[1]] ]))

		ax.set_xlabel("m/q (amu/e)")
		ax.set_ylabel("intensity")
		ax.get_yaxis().set_ticks([])

		return ln, pk, pkTxt, infoTxt

	def initAuxPlot(self, ax, times, aux_data, manual_range, scan_range, markers, aux_plot_type, color):
		scan_start, scan_end = scan_range
		manual_start, manual_end = manual_range
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
		local_norm_scan_range=(0,0)):
		# unpack the beginning and end of the mass ranges and scan ranges
		scan_start, scan_end = scan_range
		if scan_start is None:
			scan_start = 0
		if scan_end is None:
			scan_end = self.scans[-1].n

		scan_range = (scan_start, scan_end)


		mass_start, mass_end = mass_range
		if mass_start is None:
			mass_start = self.mass_start
		if mass_end is None:
			mass_end = self.mass_end

		mass_range = (mass_start, mass_end)


		ranges = self.makeScanRanges(scan_range, window, step)
		#print("making animation with scan ranges %s"%str(ranges))


		# generate the list of intensities for each mass window
		t = time.time()
		masses, intensities_list = self.makeSpectra(mass_range, ranges)
		print(time.time()-t)

		#peaks_indices_list = [self.findPeakIndicies(intensities), for intensities in intensities_list]

		# make numpy arrays with the current for each scan and the time for each scan

		times, aux_data = self.makeAuxData(scan_range, aux_plot_type, aux_smoothing)


		_, source_currents = self.makeAuxData(scan_range, AuxPlots.SOURCE_CURRENT, aux_smoothing)
		_, detector_currents = self.makeAuxData(scan_range, AuxPlots.DETECTOR_SOURCE_RATIO, aux_smoothing)
		_, ratios = self.makeAuxData(scan_range, AuxPlots.DETECTOR_CURRENT, aux_smoothing)
		_, L1_voltages = self.makeAuxData(scan_range, AuxPlots.L1_VOLTAGE, aux_smoothing)
		_, L2_voltages = self.makeAuxData(scan_range, AuxPlots.L2_VOLTAGE, aux_smoothing)

		# make the figure, upper and lower plots, and text object pool
		fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10))
		plt.tight_layout(pad=2, w_pad=1.8, h_pad=5, rect=[0.05, 0, 0.85, 1])


		ln, pk, pkTxt, infoTxt = self.initSpecPlot(ax1, intensities_list, mass_range, normalization, manual_norm_range, local_norm_scan_range)
		aux, rect = self.initAuxPlot(ax2, times, aux_data, aux_range, scan_range, markers, aux_plot_type, "black")

		ax3 = None
		aux_2 = None
		if aux_plot_type_2:
			ax3 = ax2.twinx()
			_, aux_data_2 = self.makeAuxData(scan_range, aux_plot_type_2, aux_smoothing_2)
			aux_2, rect_2 = self.initAuxPlot(ax3, times, aux_data_2, aux_range_2, scan_range, markers, aux_plot_type_2, 'blue')
			rect_2.set_alpha(0)

		# a function to update the plot for each frame
		def update(frame):
			self.updateSpecPlot(ax1, ln, pk, masses, intensities_list[frame], pkTxt, normalization)
			self.updateAuxPlot(ax2, aux, infoTxt, rect, times, aux_data, source_currents, detector_currents, ratios, L1_voltages, L2_voltages, ranges[frame], scan_start)
			if aux_plot_type_2:
				self.updateAuxPlot(ax3, aux_2, infoTxt, rect, times, aux_data_2, source_currents, detector_currents, ratios, L1_voltages, L2_voltages, ranges[frame], scan_start)
			return [ln, pk, aux, rect] + pkTxt




		ani = animation.FuncAnimation(fig, update, frames=range(len(ranges)), blit=False)
		#Writer = writers['ffmpeg']
		#FFWriter = Writer(fps=15)
		#writer = animation.AVConvFileWriter()
		if not out_name == None:
			ani.save(os.path.join(self.path, out_name))


		# make an alert noise
		print('\a')

		if show_plot:
			plt.show()


	def makeAuxPlot(self,
		markers = [],
		aux_plot_type=AuxPlots.SOURCE_CURRENT,
		aux_plot_type_2=None,
		aux_smoothing = 1,
		aux_smoothing_2 = 1,
		aux_range = (None,None),
		aux_range_2 = (None,None),
		out_name = None
		):
		scan_start, scan_end = scan_range
		if scan_start is None:
			scan_start = 0
		if scan_end is None:
			scan_end = self.scans[-1].n

		scan_range = (scan_start, scan_end)

		fig, ax = plt.subplots()		

		aux, rect = self.initAuxPlot(ax, times, aux_data, aux_range, scan_range, markers, aux_plot_type, "black")

		ax2 = None
		aux_2 = None
		if aux_plot_type_2:
			ax2 = ax.twinx()
			_, aux_data_2 = self.makeAuxData(scan_range, aux_plot_type_2, aux_smoothing_2)
			aux_2, rect_2 = self.initAuxPlot(ax2, times, aux_data_2, aux_range_2, scan_range, markers, aux_plot_type_2, 'blue')


		self.updateAuxPlot(ax, aux, infoTxt, rect, times, aux_data, source_currents, detector_currents, ratios, L1_voltages, L2_voltages, ranges[frame], scan_start)
		if aux_plot_type_2:
			self.updateAuxPlot(ax2, aux_2, infoTxt, rect, times, aux_data_2, source_currents, detector_currents, ratios, L1_voltages, L2_voltages, ranges[frame], scan_start)



def toSeconds(hms, hms0):
	h, m, s = hms
	h0, m0, s0 = hms0
	return 3600*(float(h)-float(h0)) + 60*(float(m)-float(m0)) + float(s) - float(s0)

# TODO: Change animation to use imageio?
# TODO: Change indices to always start at scan 0, with blank entries if necessary
