%{
monet.OriDesign (computed) # design matrix for orientation tuning
-> monet.DriftTrialSet
-----
design_matrix :longblob  #  design matrix for drift responses
%}

classdef OriDesign < dj.Relvar & dj.AutoPopulate

	properties
		popRel  = monet.DriftTrialSet
	end

	methods(Access=protected)

		function makeTuples(self, key)
		%!!! compute missing fields for key here
			self.insert(key)
		end
	end

end