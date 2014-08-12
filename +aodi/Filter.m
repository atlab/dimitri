%{
aodi.Filter (lookup) # specification of calcium trace filters
filter_id  : smallint
-----
algorithm : enum('dF/F','fast_oopsi')    #  algorithm specification
lo_cutoff : float      # (Hz) low-frequency cutoff
lambda    : float      # expected firing rate
%}

classdef Filter < dj.Relvar
    methods
        function fill(self)
            self.insert({
                1    'dF/F'        .05   .3
                2    'fast_oopsi'  .05   .3
                })
        end
    end
end