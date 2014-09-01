%{
aodi.Quality (computed) # study of signal quality
-> vis2p.MaskGroup
-> vis2p.VisStims
-----
# add additional attributes
%}

classdef Quality < dj.Relvar & dj.AutoPopulate

	properties
        popRel  = vis2p.MaskGroup * vis2p.VisStims & aodi.CovMatrix
    end

	methods(Access=protected)

		function makeTuples(self, key)
            disp 'fetching raw traces...'
            [tuple.cellnums, tuple.celltypes, x, y, z, traces] = ...
                fetchn(vis2p.MaskTraces*vis2p.MaskCells & key, ...
                'masknum', 'mask_type', 'img_x', 'img_y', 'img_z', 'calcium_trace');
            tuple.ntraces = numel(tuple.cellnums);
            traces = double([traces{:}]);
            
            disp 'parsing stimulus...'
            
            times = fetch1(vis2p.VisStims & key, 'frame_timestamps');
            times = times/1000;
            dt = median(diff(times));
            
            M = mean(traces);
            traces = bsxfun(@minus,traces,M);  % zero-mean
            
            disp 'subtracting largest principal components'
            ncomps = 16;
            [U,D,V] = svds(traces,ncomps);
            ix = zeros(1,ncomps);
            for iComp=1:ncomps
                [pxx,f] = pwelch(U(:,iComp),[],[],[],1/dt);
                ix(iComp) = sum(pxx(f>0.2 & f<5))/sum(pxx);  % fraction of energy in 0.1 - 10 Hz
            end
            stem(ix)
            ix = ix < 0.5;  % at least a fifth of energy must be in frequency range of calcium signals
            traces = traces - U(:,ix)*D(ix,ix)*V(:,ix)';
            
            % self.insert(key)
		end
	end

end