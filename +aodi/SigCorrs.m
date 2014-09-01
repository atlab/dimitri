%{
aodi.SigCorrs (computed) # my newest table
-> aodi.TrialTraces
-----
sig_corrs  :  longblob   #  signal correlation matrix
%}

classdef SigCorrs < dj.Relvar & dj.AutoPopulate
    
    properties
        popRel = aodi.TrialTraces
    end
    
    methods(Access=protected)
        
        function makeTuples(self, key)
            [X, evokedBins] = fetch1(aodi.TrialTraces & key, 'trace_segments', 'evoked_bins');
            X = X(1:evokedBins,:,:,:);
            X = nansum(X,3);
            X = reshape(X,[],size(X,4));
            key.sig_corrs = corrcov(nancov(X));
            self.insert(key)
        end
    end
end