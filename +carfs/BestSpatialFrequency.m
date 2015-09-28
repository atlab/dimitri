%{
carfs.BestSpatialFrequency (computed) # spatial frquency with highest response for tuned cells only
-> carfs.Trace
-----
best_spatial_freq  : float 
%}

classdef BestSpatialFrequency < dj.Relvar & dj.AutoPopulate

	properties
		popRel = pro(pro(carfs.TraceSet), carfs.VonMises, 'count(distinct spatial_freq) -> n_freqs') & 'n_freqs > 1'
	end

	methods(Access=protected)

		function makeTuples(self, key)
            data = fetch(carfs.VonMises & key & 'von_p_value<0.05', 'von_amp1', 'spatial_freq');
            [amp, masknums, spatialFreqs] = dj.struct.tabulate(data,  'von_amp1', 'masknum', 'spatial_freq');
            [~,j] = nanmax(amp, [], 2);
            for i=1:length(masknums)
                tuple = key;
                tuple.masknum = masknums(i);
                tuple.best_spatial_freq = spatialFreqs(j(i));                
                self.insert(tuple)
            end
		end
	end

end