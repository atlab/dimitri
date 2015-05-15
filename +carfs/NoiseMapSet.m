%{
carfs.NoiseMapSet (computed) # RFs computed from gaussian noise
-> carfs.TraceSet
-----
-> vis2p.VisStims
%}

classdef NoiseMapSet < dj.Relvar & dj.AutoPopulate

	properties
        popRel  = pro(carfs.TraceSet) * vis2p.VisStims & 'exp_type="MoviesExperiment"'
	end

	methods(Access=protected)

		function makeTuples(self, key)
			self.insert(key)
            makeTuples(carfs.NoiseMap,rmfield(key,'stim_idx'))
		end
	end

end