%{
monet.DriftResponse (computed) # calcium responses to drift trials
-> rf.Trace
-> monet.DriftTrialSet
-----
response: float  # averaged response
%}

classdef DriftResponse < dj.Relvar & dj.AutoPopulate

	properties
		popRel  = rf.Sync & monet.DriftTrialSet & rf.Trace
	end

	methods(Access=protected)

		function makeTuples(self, key)
		%!!! compute missing fields for key here
			self.insert(key)
		end
	end

end