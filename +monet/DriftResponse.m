%{
monet.DriftResponse (computed) # calcium responses to drift trials
-> monet.DriftResponseSet
-> monet.DriftTrial
-> rf.Trace
-----
response: float  # averaged response
%}

classdef DriftResponse < dj.Relvar

	methods

		function makeTuples(self, key)
            disp 'preparing traces...'
            frameTimes = fetch1(rf.Sync & key, 'frame_times');
            dt = median(diff(frameTimes));
            [traces, traceKeys] = fetchn(rf.Trace & key, 'ca_trace');
            
            % high-pass filter
            traces = [traces{:}];
            cutoff = 0.05;
            k = hamming(round(1/dt/cutoff)*2+1);
            k = k/sum(k);
            m = mean(traces);
            traces = bsxfun(@rdivide,traces-ne7.dsp.convmirr(double(traces),k),m);
            
            disp 'deconvolving...'
            for i=1:size(traces,2)
                traces(:,i) = fast_oopsi(double(traces(:,i)),struct('dt',dt),struct('lambda',0.3));
            end

            disp 'snipping...'
            latency = 0.03;
            for trialKey = fetch(monet.DriftTrial & key)'
                [onset, offset] = fetch1(monet.DriftTrial & trialKey, 'onset', 'offset');
                responses = num2cell(mean(traces(frameTimes > onset+latency & frameTimes < offset+latency,:)));
                tuples = dj.struct.join(traceKeys, trialKey);
                [tuples.response] = deal(responses{:});
                self.insert(tuples)
            end
            
            disp done
		end
	end

end