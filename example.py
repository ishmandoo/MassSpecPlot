from lib.msplot import *

s = Spectrum()
#s.load("C:\\Users\\ishma\\Dropbox (SteinLab)\\spectra\\(ID     186nm) 2016-07-26-tip15\\mass spec\\FEB 23 2017 16-51.cd","C:\\Users\\ishma\\Dropbox (SteinLab)\\spectra\\(ID     186nm) 2016-07-26-tip15\\ivpt\\ivpt1.tsv")
s.load("../(ID     186nm) 2016-07-26-tip15/mass spec/FEB 23 2017 16-51.cd","../(ID     186nm) 2016-07-26-tip15/ivpt/ivpt1.tsv")

#s.makeAnimation("test.mp4", (50,1000), (50,5000), 500, 50, normalization=SpecNorm.SCAN, markers=[Marker("10 mM", 1600)], aux_plot_type=AuxPlots.DETECTOR_CURRENT, aux_plot_type_2=AuxPlots.SOURCE_CURRENT)
s.makeAnimation("test.mp4", (50,1000), (50,5000), 500, 50, normalization=SpecNorm.SCAN, markers=[Marker("10 mM", 1600)], aux_plot_type=AuxPlots.DETECTOR_SOURCE_RATIO, aux_plot_type_2=AuxPlots.DETECTOR_CURRENT)
