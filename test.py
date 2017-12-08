from lib.msplot import *
import easygui

if __name__ == '__main__':
	path = easygui.diropenbox()

	s = Spectrum()

	s.load(path)
	#s.load("C:\\Users\\ishma\\Dropbox (SteinLab)\\spectra\\(ID     186nm) 2016-07-26-tip15")
	#s.load("../(ID     186nm) 2016-07-26-tip15/")

	s.plotSpectrum(out_name="test.png",label_peaks=False)