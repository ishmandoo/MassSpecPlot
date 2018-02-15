class Animation:
	def __init__(self,  
		window = 500, 
		step = 500,
		mass_range = (None, None),
		scan_range = (None, None),
		out_name = None,
		normalization=SpecNorm.SCAN,
		markers = [],
		aux_plot_type=AuxPlots.SOURCE_CURRENT,
		aux_plot_type_2=None,
		aux_smoothing = None,
		aux_smoothing_2 = None,
		spec_smoothing = None,
		aux_range = (None,None),
		aux_range_2 = (None,None),
		manual_norm_range=(0,0),
		show_plot = True,
		label_peaks = True,
		font_size = 12,
		frame_rate = 3,
		local_norm_scan_range=(0,0)):
		self.window = window
		self.step = step
		self.out_name = out_name
		self.normalization = normalization
		self.markers = markers
		self.aux_plot_type = aux_plot_type
		self.aux_plot_type_2 = aux_plot_type_2
		self.aux_smoothing = aux_smoothing
		self.aux_smoothing_2 = aux_smoothing_2
		self.aux_range = aux_range
		self.aux_range_2 = aux_range_2
		self.manual_norm_range = manual_norm_range
		self.show_plot = show_plot
		self.label_peaks = label_peaks
		self.font_size = font_size
		self.frame_rate = frame_rate
		self.local_norm_scan_range = local_norm_scan_range

	def makeAuxPlot