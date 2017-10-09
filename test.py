from scipy.io import netcdf 
from collections import defaultdict 
import math 
import matplotlib.pyplot as plt 
from copy import copy 
import peakutils 
from scipy import signal
from matplotlib.animation import FuncAnimation, ArtistAnimation 
import numpy as np 
import csv  
from bisect import bisect_left
from matplotlib.patches import Rectangle



class SpecNorm:
	GLOBAL = 0
	SCAN = 1

class AuxPlots:
	SOURCE_CURRENT = 0
	DETECTOR_CURRENT = 1
	L1_VOLTAGE = 2
	L2_VOLTAGE = 3
	PRESSURE = 4

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

		peaks_indices = peakutils.indexes(intensities, thres=0.06)

		peak_masses = [masses[i] for i in peaks_indices]
		peak_intensities = [intensities[i] for i in peaks_indices]


		plt.plot(masses, intensities)		
		plt.scatter(peak_masses, peak_intensities)
		plt.show()

	def currentPlot(self, scan_range):

		scan_start, scan_end = scan_range
		i = [scan.i for scan in self.scans[scan_start:scan_end]]
		t = [scan.t for scan in self.scans[scan_start:scan_end]]
		plt.plot(t, i)
		plt.show()

	def animation(self, mass_range, scan_range, window, step, normalization=SpecNorm.SCAN, markers = []):
		# unpack the beginning and end of the mass ranges and scan ranges
		mass_start, mass_end = mass_range
		scan_start, scan_end = scan_range

		ranges = [(i, i+window) for i in range(scan_start, scan_end-window, step)]

		# make numpy arrats with the current for each scan and the time for each scan
		currents = np.array([scan.i for scan in self.scans[scan_start:scan_end]])
		times = np.array([scan.t for scan in self.scans[scan_start:scan_end]])

		# the size of the pool of text objects for displaying peak masses
		peak_text_pool_size = 100

		# make the figure, upper and lower plots, and text object pool
		fig, (ax1, ax2) = plt.subplots(2, 1)
		plt.tight_layout(pad=1.8, w_pad=1.8, h_pad=1.8, rect=[0, 0, 0.8, 1])
		ln, = ax1.plot([], [], 'k-', animated=False)
		pk, = ax1.plot([], [], 'r.', animated=False)
		pkTxt = [ax1.text(0, 0, '', verticalalignment='bottom', horizontalalignment='center', color='red') for i in range(peak_text_pool_size)]
		infoTxt = ax1.text(1.02, 0.95, "", verticalalignment='top', horizontalalignment='left', transform=ax1.transAxes)

		aux, = ax2.plot([], [], 'k-', animated=False)

		# make the highlight rectangle showing the scan window
		rect = Rectangle((0,0),1, 1.5 * max(currents), color=(0,0,0), alpha=0.2)
		ax2.add_artist(rect)

		y_lim_bot = 0.9 * np.percentile(currents, 0.5)
		y_lim_top = 1.1 * np.percentile(currents, 99.5)

		for marker in markers:
			time = times[marker.scan - scan_start]
			ax2.plot([time, time], [y_lim_bot, y_lim_top], 'b-', dashes=[5,5], lw=1)
			ax2.text(time, y_lim_top, marker.title, verticalalignment='bottom', horizontalalignment='left', color='blue')

		# set axis limits
		ax1.set_xlim(mass_start - 0.05 * mass_end, 1.05 * mass_end)
		ax2.set_xlim(times[0] - 0.05 * times[-1], 1.05 * times[-1])
		ax2.set_ylim(y_lim_bot, y_lim_top)

		ax1.set_xlabel("m/q (amu/e)")
		ax1.set_ylabel("intensity")
		ax1.get_yaxis().set_ticks([])

		ax2.set_xlabel("time (s)")
		ax2.set_ylabel("current (nA)")

		# compile the overall spectrum for each mass window
		spec_list = [sum(self.scans[window_start:window_end]).masses for window_start, window_end in ranges]

		# generate the list of masses
		masses = range(mass_start, mass_end)

		# generate the list of intensities for each mass window
		intensities_list = [[spec[m] for m in masses] for spec in spec_list]

		if normalization == SpecNorm.GLOBAL:
			ax1.set_ylim(0, 1.1 * max([max(intensities) for intensities in intensities_list]))



		# a function to update the plot for each frame
		def update(frame):
			# unpack the frame window start and end
			window_start, window_end = ranges[frame]

			# find the window start and end frames relative to the animation start scan
			# use these for indexing into currents or times
			window_start_relative = window_start - scan_start
			window_end_relative = window_end - scan_start

			intensities = intensities_list[frame]

			# find peaks
			#peak_indices = peakutils.indexes(intensities, thres=0.06)
			peak_indices = signal.find_peaks_cwt(intensities, np.arange(2,15))
			peak_masses = [masses[i] for i in peak_indices]
			peak_intensities = [intensities[i] for i in peak_indices]

			# find max peak height for use in setting axis limits
			max_intensity = max(peak_intensities)

			# loop over the peak text pool and update as many of them as necessary
			# set the rest to blank
			for i in range(peak_text_pool_size):
				if i < len(peak_masses):
					pkTxt[i].set_text(str(peak_masses[i]))
					pkTxt[i].set_x(peak_masses[i])
					pkTxt[i].set_y(peak_intensities[i] + max_intensity * 0.03)

				else:
					pkTxt[i].set_text('')
					pkTxt[i].set_x(0)
					pkTxt[i].set_y(0)

			# find the start and end time for the window
			window_start_time = times[window_start_relative]
			window_end_time = times[window_end_relative]-times[window_start_relative]

			# update the highlighted window
			rect.set_x(window_start_time)
			rect.set_width(window_end_time)

			# set the y limits for the spectrum
			if normalization == SpecNorm.SCAN:
				ax1.set_ylim(0, 1.1 * max_intensity)

			infoTxt.set_text("scans %d - %d\n%.0fs - %.0fs\n%.2f nA" % (window_start, window_end, window_start_time, window_end_time, currents[window_start_relative]))

			ln.set_data(masses, intensities)
			pk.set_data(peak_masses, peak_intensities)
			aux.set_data(times[:window_end],currents[:window_end])
			return [ln, pk, aux, rect] + pkTxt


		

		ani = FuncAnimation(fig, update, frames=range(len(ranges)), blit=False)
		plt.show()


def toSeconds(hms, hms0):
	h, m, s = hms
	h0, m0, s0 = hms0
	return 3600*(float(h)-float(h0)) + 60*(float(m)-float(m0)) + float(s) - float(s0)

def findPeakIndices(data, width):
	smoothedData

s = Spectrum()
#s.load("C:\Users\ishma\Dropbox (SteinLab)\spectra\(ID     186nm) 2016-07-26-tip15\mass spec\FEB 23 2017 16-51.cd","C:\Users\ishma\Dropbox (SteinLab)\spectra\(ID     186nm) 2016-07-26-tip15\ivpt\ivpt1.tsv")
s.load("../(ID     186nm) 2016-07-26-tip15/mass spec/FEB 23 2017 16-51.cd","../(ID     186nm) 2016-07-26-tip15/ivpt/ivpt1.tsv")

s.animation((50,1000), (50,5000), 500, 200, normalization=SpecNorm.GLOBAL, markers=[Marker("Test", 300)])

#s.currentPlot((0,20000))

#for scan in s.scans:
#	print scan.masses

'''
a = Scan()
a.masses[1] = 1
a.masses[2] = 2

b = Scan()
b.masses[2] = 3
b.masses[3] = 4

print a
print b
print a+b
print a
'''