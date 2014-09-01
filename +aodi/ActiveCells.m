%{
aodi.ActiveCells (computed) # selection of cells for noise correlation analysis
-> aodi.TrialTraces
-----
ncells  : smallint
selection : blob
%}

classdef ActiveCells < dj.Relvar & dj.AutoPopulate
    properties
        popRel = aodi.TrialTraces & 'min_trials>200'
    end
    
    methods(Access = protected)
        function makeTuples(self, key)
            [X,evokedBins] = fetch1(aodi.TrialTraces & key, 'trace_segments','evoked_bins');
            X = X(1:evokedBins,:,:,:);  % consider only the evoked bins
            [nBins, nConds, nTrials, nCells] = size(X);
            
            % selection criterion:
            % The activity in the least active continuous quater of the trials, ordered
            % sequentially, must exceed 0.01 of activity in most active quarter
            % and the total activity must exceed 0.01 of the median
            % activity across the population.
            
            X = reshape(X,[],nCells);
            X = X(~any(isnan(X),2),:);
            [j,i] = meshgrid(1:size(X,2), ceil(linspace(eps,4,size(X,1))));
            g = accumarray({i(:),j(:)},X(:));
            m = median(g(:));
            key.selection = all(g>0.01*m) & all(bsxfun(@gt,g,0.01*max(g,[],1)));
            key.ncells = sum(key.selection);
            self.insert(key)
        end
    end
    
end