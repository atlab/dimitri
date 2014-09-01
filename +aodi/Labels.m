%{
aodi.Labels (computed) # cell labels found in site
-> aodi.TrialTraces
celltype   : varchar(255)   # cell types found in site
-----
%}

classdef Labels < dj.Relvar & dj.AutoPopulate
    
    properties
        popRel = aodi.TrialTraces
    end
    
    methods(Access=protected)
        
        function makeTuples(self, key)
            for cellType = unique(fetch1(aodi.TrialTraces & key, 'celltypes'))'           
                self.insert(setfield(key,'celltype',cellType{1})) %#ok<SFLD>
            end
        end
    end
    
end