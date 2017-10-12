from scipy.io import netcdf 
from collections import defaultdict 
import math 
import matplotlib.pyplot as plt 
from copy import copy 
from scipy import signal
from matplotlib import animation
import numpy as np 
import csv  
from bisect import bisect_left
from matplotlib.patches import Rectangle



class SpecNorm:
	GLOBAL = 0
	SCAN = 1

class AuxPlots:
	SOURCE_CURRENT = {'x':'current (nA)', 'multiplier': 1000000000.}
	DETECTOR_CURRENT = {'x':'detector counts (kcps)', 'multiplier': 1.}
	L1_VOLTAGE = {'x':'voltage (V)', 'multiplier': 1.}
	L2_VOLTAGE = {'x':'voltage (V)', 'multiplier': 1.}
	PRESSURE = {'x':'pressure (kPa)', 'multiplier': 1.}

class Scan:
	def __init__(self, t, intensity, i, vl1, vl2, p):
		self.masses = defaultdict(long)
		self.t = t
		self.i = i
		self.vl1 = vl1
		self.vl2 = vl2
		self.p = p
		self.intensity = intensity
	def __add__(self, other):
		new_scan = Scan(self.t, self.i, self.vl1, self.vl2, self.p, self.intensity)
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

	def load(self, specPath, auxPath):
		lines = None
		with open(auxPath,'r') as tsv:
		    lines = [line.strip().split('\t') for line in tsv]
		    
		(i, vl1, vl2, p, y0, m0, d0, h0, m0, s0) = lines[0]
		h0, m0, s0 = float(h0), float(m0), float(s0)

		times = [toSeconds((h,m,s),(h0,m0,s0)) for (i, vl1, vl2, p, y, m, d, h, m, s) in lines]
		current = [float(i) for (i, vl1, vl2, p, y, m, d, h, m, s) in lines]
		voltagel1 = [float(vl1) for (i, vl1, vl2, p, y, m, d, h, m, s) in lines]
		voltagel2 = [float(vl2) for (i, vl1, vl2, p, y, m, d, h, m, s) in lines]
		pressure =[float(p) for (i, vl1, vl2, p, y, m, d, h, m, s) in lines]

		with netcdf.netcdf_file(specPath, 'r') as f:
			v = copy(f.variables)

			scan_start_indices = copy(v['scan_index'])
			scan_start_times = copy(v['scan_acquisition_time'])
			scan_total_intensities = copy(v['total_intensity'])

			masses = copy(v['mass_values'])
			intensities = copy(v['intensity_values'])

			# loop over the scans
			for i in range(scan_start_indices.shape[0]-1):
				# get the starting time for the scan
				t = scan_start_times[i]

				# find the index of the line in the aux data file closest to this time
				j = bisect_left(times, t)
				if j == len(times):
					j = len(times) - 1

				# create a new Scan object for each
				new_scan = Scan(t, scan_total_intensities[i], current[j], voltagel1[j], voltagel2[j], pressure[j])

				# loop over the masses in the scan
				for point in range(scan_start_indices[i], scan_start_indices[i+1]):
					mass = math.floor(masses[point])
					intensity = intensities[point]
					# if the mass is non-zero, add it to the dictionary, combining if necessary
					if not mass == 0:
						new_scan.masses[mass] += intensity

				self.scans.append(new_scan)

	def plot(self):
		spec = sum(self.scans).masses
		masses = spec.keys()
		intensities = spec.values()

		peak_masses, peak_intensities = self.findPeaks(masses, intensities)

		plt.plot(masses, intensities)		
		plt.scatter(peak_masses, peak_intensities)
		plt.show()

	def currentPlot(self, scan_range):

		scan_start, scan_end = scan_range
		i = [scan.i for scan in self.scans[scan_start:scan_end]]
		t = [scan.t for scan in self.scans[scan_start:scan_end]]
		plt.plot(t, i)
		plt.show()

	def makeScanRanges(self, scan_range, window, step):
		scan_start, scan_end = scan_range

		return [(i, i+window) for i in range(scan_start, scan_end-window, step)]

	def makeSpectra(self, mass_range, scan_ranges):
		mass_start, mass_end = mass_range

		# compile the overall spectrum for each mass window
		spec_list = [sum(self.scans[window_start:window_end]).masses for window_start, window_end in scan_ranges]

		# generate the list of masses
		masses = range(mass_start, mass_end)

		# generate the list of intensities for each mass window
		intensities_list = [[spec[m] for m in masses] for spec in spec_list]

		return masses, intensities_list

	def makeAuxData(self, scan_range, aux_plot_type):
		scan_start, scan_end = scan_range

		if aux_plot_type == AuxPlots.SOURCE_CURRENT:
			data = np.array([scan.i * aux_plot_type['multiplier'] for scan in self.scans[scan_start:scan_end]])
		elif aux_plot_type == AuxPlots.DETECTOR_CURRENT:
			data = np.array([scan.intensity * aux_plot_type['multiplier'] for scan in self.scans[scan_start:scan_end]])
		elif aux_plot_type == AuxPlots.L1_VOLTAGE:
			data = np.array([scan.vl1 * aux_plot_type['multiplier'] for scan in self.scans[scan_start:scan_end]])
		elif aux_plot_type == AuxPlots.L2_VOLTAGE:
			data = np.array([scan.vl2 * aux_plot_type['multiplier'] for scan in self.scans[scan_start:scan_end]])
		elif aux_plot_type == AuxPlots.PRESSURE:
			data = np.array([scan.p * aux_plot_type['multiplier'] for scan in self.scans[scan_start:scan_end]])
		else:
			raise ValueError('Invalid Aux Plot Type')
		times = np.array([scan.t for scan in self.scans[scan_start:scan_end]])
		return times, data

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
		max_intensity = max(peak_intensities)

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

	def updateAuxPlot(self, ax, aux, info_text, highlight_rect, times, data, currents, scan_range, scan_start):
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
		info_text.set_text("scans %d - %d\n%.0fs - %.0fs\n%.2f nA" % (window_start, window_end, window_start_time, window_end_time, currents[window_start_relative]))

	def initSpecPlot(self, ax, intensities_list, mass_range, normalization):
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

		ax.set_xlabel("m/q (amu/e)")
		ax.set_ylabel("intensity")
		ax.get_yaxis().set_ticks([])

		return ln, pk, pkTxt, infoTxt

	def initAuxPlot(self, ax, times, aux_data, scan_range, markers, aux_plot_type, color):
		scan_start, scan_end = scan_range
		aux, = ax.plot([], [], color=color, animated=False, lw=1)
		ax.tick_params(axis='y', colors=color)
		ax.xaxis.label.set_color(color)

		y_lim_top = 1.1 * np.percentile(aux_data, 99.5)
		y_lim_bot = np.percentile(aux_data, 0.5) - 0.05 * y_lim_top

		# make the highlight rectangle showing the scan window
		rect = Rectangle((0,y_lim_bot),0, 2*y_lim_top, color=(0,0,0), alpha=0.2)
		ax.add_artist(rect)

		for marker in markers:
			time = times[marker.scan - scan_start]
			ax.plot([time, time], [y_lim_bot, y_lim_top], color = "red", dashes=[5,5], lw=1)
			ax.text(time, y_lim_top, marker.title, verticalalignment='bottom', horizontalalignment='left', color='red', rotation=45)

		# set axis limits
		ax.set_xlim(times[0] - 0.05 * times[-1], 1.05 * times[-1])
		ax.set_ylim(y_lim_bot, y_lim_top)

		ax.set_xlabel("time (s)")
		ax.set_ylabel(aux_plot_type['x'], color=color)

		return aux, rect

	def makeAnimation(self, mass_range, scan_range, window, step, normalization=SpecNorm.SCAN, markers = [], aux_plot_type=AuxPlots.SOURCE_CURRENT, aux_plot_type_2=None):
		# unpack the beginning and end of the mass ranges and scan ranges
		scan_start, scan_end = scan_range

		ranges = self.makeScanRanges(scan_range, window, step)


		# generate the list of intensities for each mass window
		masses, intensities_list = self.makeSpectra(mass_range, ranges)

		#peaks_indices_list = [self.findPeakIndicies(intensities), for intensities in intensities_list]

		# make numpy arrays with the current for each scan and the time for each scan
		times, aux_data = self.makeAuxData(scan_range, aux_plot_type)		

		if aux_plot_type == AuxPlots.SOURCE_CURRENT:
			currents = aux_data
		else:
			_, currents = self.makeAuxData(scan_range, AuxPlots.SOURCE_CURRENT)

		# make the figure, upper and lower plots, and text object pool
		fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10))
		plt.tight_layout(pad=2, w_pad=1.8, h_pad=5, rect=[0.05, 0, 0.85, 1])


		ln, pk, pkTxt, infoTxt = self.initSpecPlot(ax1, intensities_list, mass_range, normalization)
		aux, rect = self.initAuxPlot(ax2, times, aux_data, scan_range, markers, aux_plot_type, "black")

		ax3 = None
		aux_2 = None
		if aux_plot_type_2:
			ax3 = ax2.twinx()
			_, aux_data_2 = self.makeAuxData(scan_range, aux_plot_type_2)	
			aux_2, rect_2 = self.initAuxPlot(ax3, times, aux_data_2, scan_range, markers, aux_plot_type_2, 'blue')
			rect_2.set_alpha(0)
		
		# a function to update the plot for each frame
		def update(frame):
			self.updateSpecPlot(ax1, ln, pk, masses, intensities_list[frame], pkTxt, normalization)
			self.updateAuxPlot(ax2, aux, infoTxt, rect, times, aux_data, currents, ranges[frame], scan_start)
			if aux_plot_type_2:
				self.updateAuxPlot(ax3, aux_2, infoTxt, rect, times, aux_data_2, currents, ranges[frame], scan_start)
			return [ln, pk, aux, rect] + pkTxt


		

		ani = animation.FuncAnimation(fig, update, frames=range(len(ranges)), blit=False)
		#Writer = writers['ffmpeg']
		#FFWriter = Writer(fps=15)
		writer = animation.AVConvFileWriter()
		ani.save('test3.gif')
		plt.show()


def toSeconds(hms, hms0):
	h, m, s = hms
	h0, m0, s0 = hms0
	return 3600*(float(h)-float(h0)) + 60*(float(m)-float(m0)) + float(s) - float(s0)



s = Spectrum()
#s.load("C:\\Users\\ishma\\Dropbox (SteinLab)\\spectra\\(ID     186nm) 2016-07-26-tip15\\mass spec\\FEB 23 2017 16-51.cd","C:\\Users\\ishma\\Dropbox (SteinLab)\\spectra\\(ID     186nm) 2016-07-26-tip15\\ivpt\\ivpt1.tsv")
s.load("../(ID     186nm) 2016-07-26-tip15/mass spec/FEB 23 2017 16-51.cd","../(ID     186nm) 2016-07-26-tip15/ivpt/ivpt1.tsv")

s.makeAnimation((50,1000), (50,5000), 500, 50, normalization=SpecNorm.SCAN, markers=[Marker("10 mM", 1600)], aux_plot_type=AuxPlots.DETECTOR_CURRENT, aux_plot_type_2=AuxPlots.SOURCE_CURRENT)

# TODO: Change animation to use imageio?
# TODO: Change indices to always start at scan 0, with blank entries if necessary