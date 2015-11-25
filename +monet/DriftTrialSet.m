%{
monet.DriftTrialSet (computed) # all drift trials for this scan
-> rf.Sync
-----
%}

classdef DriftTrialSet < dj.Relvar & dj.AutoPopulate

	properties
		popRel = rf.Sync & psy.MovingNoise
	end

	methods(Access=protected)

		function makeTuples(self, key)
		%!!! compute missing fields for key here
			self.insert(key)
		end
	end

end