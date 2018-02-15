from lib.msplot import *
import easygui
from lib.data import Data, AuxPlots
from lib.plotter import Plotter
import numpy as np

if __name__ == '__main__':
	path = easygui.diropenbox()
	d = Data()

	d.load(path)

	p = Plotter(d)
	print(d.aux_data)


# if __name__ == '__main__':
# 	path = easygui.diropenbox()

# 	s = Spectrum()

# 	s.load(path, onblast_cutoff=150)
# 	#s.load("C:\\Users\\ishma\\Dropbox (SteinLab)\\spectra\\(ID     186nm) 2016-07-26-tip15")
# 	#s.load("../(ID     186nm) 2016-07-26-tip15/")
	
	p.makeAnimation(
		window = 500,
		step = 20,
		scan_range = (None, None),
		mass_range = (None, None),
		out_name = "test.mp4",
		normalization=SpecNorm.SCAN,
		label_peaks=True,
		aux_plot_type=AuxPlots.DETECTOR_CURRENT,
		aux_plot_type_2=AuxPlots.DETECTOR_SOURCE_RATIO,
		aux_smoothing=10
		)