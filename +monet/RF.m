%{
monet.RF (computed) # my newest table
-> rf.Sync
-> rf.Trace
-> monet.RFMethod
-----
nbins              : smallint                      # temporal bins
bin_size           : float                         # (ms) temporal bin size
degrees_x          : float                         # degrees along x
degrees_y          : float                         # degrees along y
rf_map             : longblob                      # receptive field map

-----
# add additional attributes
%}

classdef RF < dj.Relvar & dj.AutoPopulate

	properties
		popRel = rf.Sync*rf.Segment*monet.RFMethod  & psy.MovingNoise & rf.Trace
	end

	methods(Access=protected)

		function makeTuples(self, key)
			self.insert(key)
		end
	end

end