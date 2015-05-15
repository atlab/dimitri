%{
carfs.ImageMapSet (computed) # receptive fields computed from natural images
-> carfs.TraceSet
-----
->vis2p.VisStims
%}

classdef ImageMapSet < dj.Relvar & dj.AutoPopulate
    
    properties
        popRel  = pro(carfs.TraceSet) * vis2p.VisStims & 'exp_type="NatImExperiment"'
    end
    
    methods(Access=protected)
        
        function makeTuples(self, key)
            self.insert(key)
            makeTuples(carfs.ImageMap, rmfield(key,'stim_idx'))
        end
    end
    
end