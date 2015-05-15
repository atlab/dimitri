%{
carfs.TraceSet (computed) # deconvolved traces
-> vis2p.MaskGroup
-> vis2p.MaskGroupRaw
-----
-> vis2p.VisStims
trace_times   : longblob  #  (s) frame times taken from vis2p.VisStims
%}

classdef TraceSet < dj.Relvar & dj.AutoPopulate
    
    properties
        popRel  = vis2p.MaskGroup * pro(vis2p.MaskGroupRaw);
        %            & pro(vis2p.VisStims & 'exp_type="NatImExperiment"')
    end
    
    methods(Access=protected)
        
        function makeTuples(self, key)
            stimKeys = fetch(vis2p.VisStims & key, 'frame_timestamps');
            assert(all(arrayfun(@(a) all(a.frame_timestamps==stimKeys(1).frame_timestamps),stimKeys)), ...
                'discrepancy of synchronization')
            tuple = key;
            tuple.stim_idx = stimKeys(1).stim_idx;
            
            frameTimes = fetch1(vis2p.MaskGroupRaw & key, 'frame_timestamps')/1000;
            dt = 1/20;  % Hz
            tuple.trace_times = frameTimes(1):dt:frameTimes(end);
            
            self.insert(tuple)
            makeTuples(carfs.Trace,key,frameTimes,tuple.trace_times)
        end
    end
    
end