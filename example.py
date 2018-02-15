from lib.msplot import *
import easygui

if __name__ == '__main__':
	path = easygui.diropenbox()

	s = Spectrum()

	s.load(path)
	#s.load("C:\\Users\\ishma\\Dropbox (SteinLab)\\spectra\\(ID     186nm) 2016-07-26-tip15")
	#s.load("../(ID     186nm) 2016-07-26-tip15/")

	s.makeAnimation(
		window = 500,
		step = 500,
		scan_range = (None, None),
		mass_range = (None, None),
		out_name = "test.mp4",
		normalization=SpecNorm.SCAN,
		label_peaks=True,
		markers=[Marker("10 mM", 11000)],
		aux_plot_type=AuxPlots.L1_VOLTAGE,
		aux_plot_type_2=AuxPlots.L2_VOLTAGE,
		aux_smoothing=3,
		aux_smoothing_2=3,
		spec_smoothing=2,
		aux_range=(None,None),
		aux_range_2=(None,None),
		font_size=12,
		frame_rate=10
		)


'''
	optional parameters for makeAnimation:
	--------------------------------------
	window -- the number of scans in a scan window
	       -- default is 500
	step -- the number of scans the window moves with each step
	     -- default is 500
	scan_range -- a tuple giving the starting and ending scans in the animation
	           -- a None in the first position starts at scan 0
	           -- a None in the second ends at the last scan
	           -- default is (None, None)
	mass_range -- a tuple giving the starting and ending mass
			   -- a None in the first position starts at the first mass
	           -- a None in the second ends at the last mass
	           -- default is (None, None)
	out_name -- the name of the output file
			 -- default is None, giving no output file
	show_plot -- a boolean which controls whether or not a plot shows up after processing
			  -- default is True
	normalization -- sets the normalization of the spectra, see below for options
	              -- default is SpecNorm.GLOBAL
	label_peaks -- a boolean indicating whether peaks should be labeled in the spectrum plot
				-- default is True
	markers -- a list of markers in the aux plot, see below for details
			-- default is an empty list, no markers
	aux_plot_type -- the type of the aux plot left axis, see below for options
				  -- default is AuxPlots.SOURCE_CURRENT
	aux_smoothing -- the width of the gaussian smoothing on the aux plot left axis
				  -- default is no smoothing
	spec_smoothing -- the width of the gaussian smoothing on the aux plot left axis
				  -- default is no smoothing
	aux_range -- a tuple giving the bounds of the y-axis of the aux plot left axis
			  -- a None in either positions gives an automatic bound
			  -- default is (None, None)	
	aux_plot_type_2 -- the type of the aux plot right axis, see below for options
					-- None gives no right axis
				    -- default is None
	aux_smoothing_2 -- the width of the gaussian smoothing on the aux plot right axis
				    -- default is 1
	aux_range_2 -- a tuple giving the bounds of the y-axis of the aux plot right axis
			    -- a None in either position gives an automatic bound
			    -- default is (None, None)
	manual_norm_range -- a tuple giving the bounds of the y-axis of the spectrum plot when using SpecNorm.MANUAL
					  -- default is (0,0), which breaks the plot
	local_norm_scan_range -- a tuple giving the range of windows to use for the normalization
						  -- the y-axis of the spectra will be set using the highest peak in the window range
						  -- default is (0,0), which breaks the plot
	frame_rate -- a number setting the frame rate in frames per second
			   -- default is 3 fps
	font_size -- a number setting the font size of the titles
			  -- other font sizes are set relative to this
			  -- default is 12


	AuxPlots options:
	-----------------
	AuxPlots.SOURCE_CURRENT -- Keithley Current
	AuxPlots.DETECTOR_CURRENT -- Extrel Detector Current
	AuxPlots.DETECTOR_SOURCE_RATIO -- DETECTOR_CURRENT/SOURCE_CURRENT
	AuxPlots.L1_VOLTAGE -- voltage on L1
	AuxPlots.L2_VOLTAGE -- voltage on L2
	AuxPlots.PRESSURE -- chamber pressure

	SpecNorm options:
	-----------------
	SpecNorm.GLOBAL -- normalized to the highest peak in the whole animation
	SpecNorm.SCAN  -- normalized to the highest peak in each window
	SpecNorm.MANUAL -- normalized to the level specified in the optional paramter manual_norm_range
	SpecNorm.LOCAL -- normalized to the highest peak in the window range specified in local_norm_scan_range

	markers:
	--------
	create a marker option using the marker constructor Marker(name, scan_number)
'''