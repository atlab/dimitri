%{
aodi.Filter (lookup) # specification of calcium trace filters
filter_id       : smallint               # 
---
algorithm                   : enum('dF/F','fast_oopsi')     # algorithm specification
lo_cutoff                   : float                         # (Hz) low-frequency cutoff
lambda                      : float                         # expected firing rate
subtract_pc=0               : tinyint                       # number of princ comps to subtract
requested_binsize=100       : smallint                      # (ms) bin duration
%}

classdef Filter < dj.Relvar
    methods
        function fill(self)
            self.inserti({
                1    'dF/F'        .05   .3  0   100
                2    'fast_oopsi'  .05   .3  0   100
                4    'fast_oopsi'  .05   .3  1   100
                6    'fast_oopsi'  .05   .3  24  100
                
                11    'fast_oopsi' .05   .3  1  150
                })
        end
    end
end
