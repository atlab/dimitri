%{
aodi.TrialTraces (computed) # my newest table
-> vis2p.MaskGroup
-> vis2p.MaskGroupRaw
-> vis2p.VisStims
-> aodi.Filter
---
cellnums                    : blob                          # cell selection
celltypes                   : blob                          # cell types: neuron, PV, SST, VIP, etc
ntraces                     : smallint                      # number neurons
cell_xyz                    : blob                          # cell positions
evoked_bins                 : tinyint                       # number of evoked bins -- the remaining bins are spontaneous
latency                     : smallint                      # (ms) presumed screen-to-V1 latency
binsize                     : smallint                      # (ms) bin duration
ntrials                     : blob                          # the number of trials for each direction
min_trials                  : smallint                      # min number of repeats in each condition
ndirs                       : tinyint                       # number of condition
directions                  : blob                          # directions of motion
trace_segments              : longblob                      # trace segements: nBins x nDirs x nTrials x nCells
%}

classdef TrialTraces < dj.Relvar & dj.AutoPopulate
    
    properties
        popRel  = vis2p.MaskGroup * pro(vis2p.MaskGroupRaw) * vis2p.VisStims * aodi.Filter & 'algorithm="fast_oopsi"' & 'filter_id in (4,11)'
    end
    
    methods(Access=protected)
        
        function makeTuples(self, key)
            tuple = key;
            
            disp 'fetching raw traces...'
            [tuple.cellnums, tuple.celltypes, x, y, z, traces] = ...
                fetchn(pro(vis2p.MaskTraces,'mask_type')*pro(vis2p.MaskCells,'img_x','img_y','img_z')*pro(vis2p.MaskTracesRaw,'calcium_trace') & key, ...
                'masknum', 'mask_type', 'img_x', 'img_y', 'img_z', 'calcium_trace');
            tuple.ntraces = numel(tuple.cellnums);
            traces = double([traces{:}]);
            
            % compute xyz (microns)
            volumeDims = [200 200 100]; % x,y,z microns
            
            depth = fetchn(vis2p.Depth & key, 'depth');
            if isempty(depth)
                depth = fetch1(vis2p.Scans*vis2p.Depth & key, 'z-surfz->depth');
            end
            [xsize,ysize,zsize] = fetch1(vis2p.Movies & key, 'xsize','ysize','zsize');
            x =  x/xsize*volumeDims(1);
            y =  y/ysize*volumeDims(2);
            z = -z/zsize*volumeDims(3) + depth;
            tuple.cell_xyz = [x y z];
            
            
            disp 'parsing stimulus...'
            filt = fetch(aodi.Filter & key,'*');
            times = fetch1(vis2p.MaskGroupRaw & key, 'frame_timestamps');
            times = times/1000;
            dt = median(diff(times));
            
            tuple.binsize = filt.requested_binsize;  % ms
            tuple.latency = 30;  % ms   visual latency
            times = times - tuple.latency/1000;
            
            stimFile = fetch1(vis2p.VisStims & key, 'stim_filename');
            stim = load(getLocalPath(stimFile));
            stim = stim.stim;
            directions = [stim.params.conditions.orientation];
            if numel(unique(diff(sort(directions))))>1
                warning 'nonuniform directions'
            end
            onsetType = find(strcmp(stim.eventTypes,'showSubStimulus'));
            clearScreenType = find(strcmp(stim.eventTypes,'clearScreen'));
            assert(~isempty(onsetType) && ~isempty(clearScreenType))
            nTrials = numel(stim.params.trials);
            assert(numel(stim.events)==nTrials)
            trialNum = 0;
            presentations = [];
            for iTrial = 1:nTrials
                interest = ismember(stim.events(iTrial).types, [onsetType,clearScreenType]);
                trialTimes = stim.events(iTrial).syncedTimes(interest)/1000;
                types = stim.events(iTrial).types(interest);
                assert(types(end)~=onsetType, 'last event cannot be an onset')
                iEvents = find(types==onsetType);
                onsets = trialTimes(iEvents);
                durations = trialTimes(iEvents+1) - onsets;
                for i = 1:length(onsets)
                    trialNum = trialNum + 1;
                    rec.trial_num = trialNum;
                    rec.direction = directions(stim.params.trials(iTrial).conditions(i));
                    rec.onset = onsets(i);
                    rec.duration = durations(i);
                    presentations = [presentations; rec]; %#ok<AGROW>
                end
            end
            duration = median([presentations.duration]);
            assert(mean(abs(1-[presentations.duration]/duration)>0.05)<0.99, ...
                'presentations do not have uniform durations')
            intervals = diff(sort([presentations.onset]));
            interval = median(intervals);
            assert(mean(intervals < 0.95*interval)<0.99, ...
                'short intertrial intervals exist')
            
            trialBins = round(interval/(tuple.binsize/1000));
            tuple.evoked_bins = min(trialBins,round((duration+0.2)/(tuple.binsize/1000)));
            
            
            M = mean(traces);
            traces = bsxfun(@minus,traces,M);  % zero-mean
            if filt.subtract_pc
                disp 'subtracting largest principal components'
                [U,D,V] = svds(traces,filt.subtract_pc);
                ix = true(1,filt.subtract_pc);
                if filt.subtract_pc > 1
                    disp 'Select components interactively in ix'
                    keyboard                
                end
                traces = traces - U(:,ix)*D(ix,ix)*V(:,ix)';
            end
            
            disp 'prefiltering (high-pass)...'
            k = hamming(2*floor(1/filt.lo_cutoff/dt)+1);
            k = k/sum(k);
            traces = traces - convmirr(traces,k);
            
            disp 'downsampling to 20 Hz ...'
            cutoff = 20;  % Hz
            d = round(1/cutoff/dt);
            k = hamming(2*d+1);
            k = k/sum(k); 
            traces = convmirr(traces,k);
            traces = traces(ceil(d/2):d:end,:);
            times = times(ceil(d/2):d:end);
            dt = mean(diff(times));
            traces = bsxfun(@rdivide, traces, std(traces));
            
            
            switch filt.algorithm
                case 'fast_oopsi'
                    disp 'deconvolving traces...'
                    deconv = @(x) fast_oopsi(x,struct('dt',dt), struct('lambda',filt.lambda));
                    for i=1:size(traces,2)
                        traces(:,i) = deconv(traces(:,i));
                    end
                otherwise
                    error('algorithm "%s" has not been implemented', filt.algorithm)
            end
            
            disp 'binnning traces'
            % low-pass filter for binning
            k = hamming(2*floor(tuple.binsize/1000/dt)+1);
            k = k/sum(k);
            traces = convmirr(traces,k);
            snippets = cell(size(presentations));
            times = times(1:size(traces,1));  % sometimes times has one extra sample
            % resample signals relative to the onset
            for i=1:numel(presentations)
                snippets{i} = interp1q(times',traces,presentations(i).onset + (0:trialBins-1)'*tuple.binsize/1000);
            end
            [presentations.snippet] = deal(snippets{:});
            [snip,directions] = dj.struct.tabulate(presentations,'snippet','direction');
            tuple.directions = directions;
            tuple.ndirs = length(directions);
            tuple.ntrials = sum(~cellfun(@(x) isempty(x) || any(isnan(x(:))),snip),2);  % trials per condition
            tuple.min_trials = min(tuple.ntrials);
            
            snip(cellfun(@isempty,snip)) = {nan(trialBins,tuple.ntraces)};
            tuple.trace_segments = permute((cell2mat(permute(snip,[3,4,1,2]))),[1 3 4 2]);
            
            self.insert(tuple)
        end
    end
    
end
