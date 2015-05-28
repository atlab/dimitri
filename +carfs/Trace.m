%{
carfs.Trace (computed) # deconvolved traces
-> carfs.TraceSet
-> vis2p.MaskTraces
---
mask_type                   : enum('neuron','red','astrocyte') # fluorescent label
trace                       : longblob                      # deconvovled trace
%}

classdef Trace < dj.Relvar
    methods
        
        function makeTuples(self, key, frameTimes, traceTimes)
            disp 'fetching raw traces...'
            tuples = ...
                fetch(carfs.TraceSet*pro(vis2p.MaskTraces,'mask_type')*pro(vis2p.MaskTracesRaw,'calcium_trace') & key, ...
                'masknum', 'mask_type', 'calcium_trace');
            dt = median(diff(frameTimes));
            
            disp 'high-pass filtration...'
            cutoff = 0.05;
            k = hamming(round(1/dt/cutoff)*2+1);
            k = k/sum(k);
            offset = 495420;   % This is specific to Manolis' setup
            traces = [tuples.calcium_trace]+offset;
            mn = ne7.dsp.convmirr(double(traces),k);
            traces = (traces-mn)./mn;
            
            disp 'projecting out first principal component...'
            [U,D,V] = svds(double(traces),1);
            traces = traces - U*D*V';
            
            disp 'downsampling...'
            k = hamming(ceil(median(diff(traceTimes))/dt)*2+1);
            k = k/sum(k);
            traces = ne7.dsp.convmirr(double(traces),k);
            traces = interp1(frameTimes,traces,traceTimes);
            dt = diff(traceTimes([1 2]));
            
            disp 'deconvolving...'
            for i=1:size(traces,2)
                traces(:,i) = single(fast_oopsi(traces(:,i),struct('dt',dt),struct('lambda',0.3)));
            end
            
            disp 'inserting...'
            for i=1:length(tuples)
                tuple = tuples(i);
                tuple.trace = traces(:,i);
                self.insert(rmfield(tuple,'calcium_trace'))
            end
            fprintf done\n
        end
    end
    
end
