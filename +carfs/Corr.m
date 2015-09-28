%{
carfs.Corr (computed) # my newest table
-> carfs.TraceSet
-----
masknums    : longblob  # just in case the order is changed
corr_matrix : longblob  # partial correlation coefficient 
%}

classdef Corr < dj.Relvar & dj.AutoPopulate

	properties
		popRel = pro(carfs.TraceSet) & carfs.SubfieldSet & carfs.VonMises
	end

	methods(Access=protected)

		function makeTuples(self, key)
            
            % fetch traces
            disp 'fetching traces...'
            traceTimes = fetch1(carfs.TraceSet & key,'trace_times');
            [traces, traceKeys] = fetchn(carfs.Trace & key, 'trace');
            traces = [traces{:}];
            nTraces = size(traces,2);
            
            % downsample by dfactor
            dfactor = 3;
            k = hamming(dfactor*2+1); 
            k=k/sum(k);
            traces = conv2(traces, k, 'valid');
            traceTimes = conv2(traceTimes(:), k, 'valid');
            traces = traces(1:dfactor:end,:);
            traceTimes = traceTimes(1:dfactor:end);
            
            disp 'loading movies...'
            stim = fetch1(vis2p.VisStims & key & 'exp_type="MoviesExperiment"', 'stim_file');
            stimTimes = [stim.events.syncedTimes];
            stimTimes = stimTimes(ismember([stim.events.types], ...
                find(ismember(stim.eventTypes,{'showStimulus','endStimulus'}))))/1000;
            ix = traceTimes > min(stimTimes) & traceTimes < max(stimTimes);
            traces = traces(ix, :);
            traceTimes = traceTimes(ix);
            
            
            % compute correlations
            traces = reshape(traces,[1 1 size(traces)]);
            hyperParameters = {exp(-7:.05:-1) exp(-5:.05:-1)};
            [hyperParameters, bestDelta, visited, losses] = ...
                    cove.crossEstimateHyper(traces, 1, 'lv-glasso', hyperParameters);
            
            [R,M,V,extras] = cove.estimate(traces, 1, 'lv-glasso', hyperParameters);
            
            key.corr_matrix = -corrcov(inv(R));
            key.masknums = [traceKeys.masknum];
            self.insert(key)
		end
	end

end