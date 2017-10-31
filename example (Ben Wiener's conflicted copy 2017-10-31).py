from lib.msplot import *
import easygui

path = easygui.diropenbox()

s = Spectrum()


s.load(path)
#s.load("C:\\Users\\ishma\\Dropbox (SteinLab)\\spectra\\(ID     186nm) 2016-07-26-tip15")
#s.load("../(ID     186nm) 2016-07-26-tip15/")

#s.makeAnimation("test.mp4", (50,1000), (50,5000), 500, 50, normalization=SpecNorm.SCAN, markers=[Marker("10 mM", 1600)], aux_plot_type=AuxPlots.DETECTOR_CURRENT, aux_plot_type_2=AuxPlots.SOURCE_CURRENT)
s.makeAnimation((50,1000), (50,-1), 500, 500,
	out_name = "test.mp4",
	normalization=SpecNorm.SCAN,
	markers=[Marker("10 mM", 1600)],
	aux_plot_type=AuxPlots.DETECTOR_SOURCE_RATIO,
	aux_plot_type_2=AuxPlots.DETECTOR_CURRENT,
	aux_smoothing=10,
	aux_smoothing_2=10,
	aux_range=(None,1.25e8),
	aux_range_2=(1000000000,None)
	)
