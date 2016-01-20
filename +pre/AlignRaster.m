%{
pre.AlignRaster (computed) # raster alignment for bidirectional tuning
-> pre.ScanCheck
-----
raster_phase             : float                         # shift of odd vs even raster lines
%}

classdef AlignRaster < dj.Relvar & dj.AutoPopulate
    
    properties
        popRel = pre.ScanCheck
    end
    
    methods(Access=protected)
        
        function makeTuples(self, key)
            [template, bidirectional, fill_fraction] = fetch1(pre.ScanCheck*pre.ScanInfo & key, ...
                'template', 'bidirectional', 'fill_fraction');
            if bidirectional
                key.raster_phase = ne7.ip.computeRasterCorrection(template, fill_fraction);
            end
            self.insert(key)
        end
    end
    
end