%{
carfs.SpaceTime (lookup) # spatiotemporal frequency selection, populated by dependent tables as needed.
spatial_freq    : decimal(5,5)           # cycles/degree (-1 = marginalize)
temp_freq       : decimal(4,2)           # Hz (-1 = marginalize)
---
%}

classdef SpaceTime < dj.Relvar
end
