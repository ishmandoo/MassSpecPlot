from lib.msplot import *
import easygui

if __name__ == '__main__':
	path = easygui.diropenbox()

	s = Spectrum()

	s.load(path)
	#s.load("C:\\Users\\ishma\\Dropbox (SteinLab)\\spectra\\(ID     186nm) 2016-07-26-tip15")
	#s.load("../(ID     186nm) 2016-07-26-tip15/")

	s.makeAnimation(
		window = 1000,
		step = 1000,
		scan_range = (None, None),
		mass_range = (None, None),
		out_name = "test2.mp4",
		normalization=SpecNorm.SCAN,
		label_peaks=True,
		aux_plot_type=AuxPlots.EFFICIENCY,
		#aux_plot_type_2=AuxPlots.DETECTOR_SOURCE_RATIO,
		aux_smoothing=100,
		font_size=14
		)