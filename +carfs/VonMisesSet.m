%{
carfs.VonMisesSet (computed) # my newest table
-> carfs.TraceSet
-----
-> vis2p.VisStims
nshuffles : int  # shuffles for resampling test
%}

classdef VonMisesSet < dj.Relvar & dj.AutoPopulate
    
    properties
        popRel = pro(carfs.TraceSet) & (vis2p.VisStims & 'exp_type="MultDimExperiment"')
    end
    
    methods(Access=protected)     
        function makeTuples(self, key)
            tuple = dj.struct.join(key, ...
                fetch(vis2p.VisStims & key & 'exp_type="MultDimExperiment"'));
            tuple.nshuffles = 10000;
            self.insert(tuple)
            
            makeTuples(carfs.VonMises, tuple)
        end
    end    
end