%{
carfs.CorrPair (computed) # my newest table
-> carfs.Corr
role : enum('first','second')  # order in pair does not matter
-----
-> carfs.Trace
%}

classdef CorrPair < dj.Relvar
end