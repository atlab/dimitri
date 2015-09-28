%{
carfs.NoiseMapSet (computed) # RFs computed from gaussian noise
-> carfs.TraceSet
---
-> vis2p.VisStims
bin_size=0.1                : float                         # (s) temporal bin size
%}

classdef NoiseMapSet < dj.Relvar & dj.AutoPopulate

	properties
        popRel  = pro(carfs.TraceSet) & (vis2p.VisStims & 'exp_type="MoviesExperiment"')
	end

	methods(Access=protected)

		function makeTuples(self, key)
            tuple = fetch(pro(carfs.TraceSet) * vis2p.VisStims & 'exp_type="MoviesExperiment"' & key);
            assert(length(tuple)==1, 'only one movie experiment allowed')
			binSize = 0.1;
            self.insert(setfield(tuple, 'bin_size', binSize)) %#ok<SFLD>
            makeTuples(carfs.NoiseMap, key, binSize)
		end
	end

end
