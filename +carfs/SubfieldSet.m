%{
carfs.SubfieldSet (computed) # my newest table
-> carfs.NoiseMap
-----
%}

classdef SubfieldSet < dj.Relvar & dj.AutoPopulate
    
    properties
        popRel = carfs.NoiseMap & pro(pro(carfs.TraceSet)*carfs.VonMises)
        
        subordinates = {carfs.Subfield}
    end
    
    methods(Access=protected)
        
        function makeTuples(self, key)
            self.insert(key)
            makeTuples(carfs.Subfield, key)
        end
    end
    
end