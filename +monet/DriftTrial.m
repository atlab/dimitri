%{
monet.DriftTrial (computed) # noise drift trials
-> monet.DriftTrialSet
trial:  smallint   # trial index
-----
direction:  float   # (degrees) direction of drift
onset: double  # (s) onset time in rf.Sync times
offset: double # (s) offset time in rf.Sync times
%}

classdef DriftTrial < dj.Relvar


	methods(Access=protected)

		function makeTuples(self, key)
			self.insert(key)
		end
	end

end