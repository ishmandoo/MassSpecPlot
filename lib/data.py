import numpy as np
import math
import os
import csv
from scipy.io import netcdf
from bisect import bisect_left
import glob
from warnings import warn
from enum import Enum
import re

class AuxPlots():
	SOURCE_CURRENT = {'key': 0, 'x':'current (nA)', 'multiplier': 1000000000.}
	DETECTOR_CURRENT = {'key': 1, 'x':'detector counts (kcps)', 'multiplier': 1.}
	DETECTOR_SOURCE_RATIO = {'key': 2, 'x':'ratio (kcps/nA)', 'multiplier': .000000001}
	L1_VOLTAGE = {'key': 3, 'x':'L1 voltage (V)', 'multiplier': 1000.}
	L2_VOLTAGE = {'key': 4, 'x':'L2 voltage (V)', 'multiplier': 1000.}
	PRESSURE = {'key': 5, 'x':'pressure (kPa)', 'multiplier': 1.}

class Data:
	def __init__(self):
		self.scans = None
		self.times = None
		self.currents = None
		self.voltagesl1 = None
		self.voltagesl2 = None
		self.pressures = None
		self.intensities = None
		self.ratios = None

		self.aux_data = None

		self.path = None
		self.mass_start = None
		self.mass_end = None
		self.s0 = 0

	def toSeconds(self, hms, hms0):
		h, m, s = hms
		h0, m0, s0 = hms0
		return 3600*(float(h)-float(h0)) + 60*(float(m)-float(m0)) + float(s) - float(s0)

	def load(self, path, specPath=None, auxPath=None, logPath=None, onblast_cutoff=None):
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

			# data from the cdf file
			masses = v['mass_values']
			intensities = v['intensity_values']

			scan_start_indices = v['scan_index']
			scan_start_times = v['scan_acquisition_time']
			scan_total_intensities = v['total_intensity']


			self.mass_start = int(v['mass_range_min'][0])
			self.mass_end = int(v['mass_range_max'][0])

			self.scans = np.zeros((scan_start_indices.shape[0], self.mass_end))

			time_stamp = f.experiment_date_time_stamp

			h0, m0, s0 = time_stamp[8:10], time_stamp[10:12], self.s0

			# data from the aux file
			aux_times = [self.toSeconds((h,m,s),(h0,m0,s0)) for (i, vl1, vl2, p, y, m, d, h, m, s) in lines]
			currents = [float(i) for (i, vl1, vl2, p, y, m, d, h, m, s) in lines]
			voltagesl1 = [float(vl1) for (i, vl1, vl2, p, y, m, d, h, m, s) in lines]
			voltagesl2 = [float(vl2) for (i, vl1, vl2, p, y, m, d, h, m, s) in lines]
			pressures =[float(p) for (i, vl1, vl2, p, y, m, d, h, m, s) in lines]

			self.times = scan_start_times

			# find the indices of the line in the aux data file closest to this time
			aux_indices = [max(bisect_left(aux_times, t), len(aux_times)-1) for t in scan_start_times]
			self.currents = np.array([currents[i] for i in aux_indices]) * AuxPlots.SOURCE_CURRENT['multiplier']
			self.voltagesl1 = np.array([voltagesl1[i] for i in aux_indices]) * AuxPlots.L1_VOLTAGE['multiplier']
			self.voltagesl2 = np.array([voltagesl2[i] for i in aux_indices]) * AuxPlots.L2_VOLTAGE['multiplier']
			self.pressures = np.array([pressures[i] for i in aux_indices]) * AuxPlots.PRESSURE['multiplier']

			# loop over the scans
			for i in range(scan_start_indices.shape[0]-1):


				# loop over the masses in the scan
				for point in range(scan_start_indices[i], scan_start_indices[i+1]):
					mass = math.floor(masses[point])
					intensity = intensities[point]
					# if the mass is non-zero, add it to the dictionary, combining if necessary
					self.scans[i, mass] += intensity

			if onblast_cutoff is None:
				self.intensities = np.array(scan_total_intensities.data)
			else:
				self.intensities = np.sum(self.scans[:,onblast_cutoff:], axis=1)
			self.ratios = self.intensities * AuxPlots.DETECTOR_SOURCE_RATIO['multiplier'] / self.currents
			self.intensities *= AuxPlots.DETECTOR_CURRENT['multiplier']

			self.aux_data = [
				self.currents,
				self.intensities,
				self.ratios,
				self.voltagesl1,
				self.voltagesl2,
				self.pressures
			]