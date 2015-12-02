%{
monet.DriftResponseSet (computed) # my newest table
-> monet.DriftTrialSet
-----
%}

classdef DriftResponseSet < dj.Relvar & dj.AutoPopulate

	properties
		popRel = monet.DriftTrialSet & monet.DriftTrial & rf.Trace
	end

	methods(Access=protected)

		function makeTuples(self, key)
			self.insert(key)
            makeTuples(monet.DriftResponse, key)
		end
	end

end