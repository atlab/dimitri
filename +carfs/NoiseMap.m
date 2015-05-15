%{
carfs.NoiseMap (computed) # receptive fields computed from gaussian noise
-> carfs.NoiseMapSet
-> carfs.Trace
-----
noise_map : longblob  # receptive field for each trace
%}

classdef NoiseMap < dj.Relvar
    
    
    methods
        
        function dump(self)
            for key = self.fetch'
                disp(key)
                map = fetch1(self & key, 'noise_map');
                if any(isnan(map(:)))
                    disp 'found nans'
                else
                    mx = max(abs(map(:)));
                    map = round(map/mx*31.5 + 32.5);
                    cmap = ne7.vis.doppler;

                    for i=1:size(map,3)
                        im = reshape(cmap(map(:,:,i),:),[size(map,1) size(map,2) 3]);
                        f = sprintf('~/Desktop/carfs/%d/%d-%d_%02d.png',...
                            key.mouse_id, key.scan_idx, key.masknum, i);
                        imwrite(im,f)
                    end
                end
            end
        end
        
        
        function makeTuples(self, key)
            
            % temporal binning
            nBins = 12;
            binSize = 0.04;  %s
            
            disp 'fetching traces...'
            traceTimes = fetch1(carfs.TraceSet & key,'trace_times');
            dt = median(diff(traceTimes));
            [traces, traceKeys] = fetchn(carfs.Trace & key, 'trace');
            traces = [traces{:}];
            nTraces = size(traces,2);
            
            disp 'loading movies...'
            stimFile = fetch1(vis2p.VisStims & key & 'exp_type="MoviesExperiment"', 'stim_filename');
            stim = load(getLocalPath(stimFile));
            stim = stim.stim;
            
            moviePath = getLocalPath('~/Desktop/ManolisNoiseMovies');
            nTrials = length(stim.events);
            fps = 30;  
            
            % resample to RF times, aligned on onset
            k = hamming(ceil(binSize/dt)*2+1);  % kernel for downsampling calcium traces
            k = k/sum(k);            
            traces = ne7.dsp.convmirr(traces,k);
            maps = num2cell(zeros(nTraces,1));

            for iTrial = 1:nTrials
                fprintf .
                iMovie = stim.params.conditions(stim.params.trials(iTrial).condition).movieNumber;
                movie = read(VideoReader(fullfile(moviePath,sprintf('noiseMap%d.avi',iMovie)))); %#ok<TNMLP>
                movie = double(squeeze(movie(:,:,1,:)))/127-1;
                onset = stim.events(iTrial).syncedTimes(ismember(stim.eventTypes(stim.events(iTrial).types),'showStimulus'))/1000;
                offset = stim.events(iTrial).syncedTimes(ismember(stim.eventTypes(stim.events(iTrial).types),'endStimulus'))/1000;
                offset = min(offset,onset+(size(movie,3)-1)/fps);
                assert(isscalar(onset) && isscalar(offset))
                assert(onset > traceTimes(1), 'onsets should not start before trace')
                
                t0 = max(traceTimes(1),onset+(nBins-1)*binSize)+0.05;  % start time for traces
                t1 = t0-(nBins-1)*binSize;                             % start time for movie
                t2 = min(traceTimes(end),offset);                      % end time for both
                
                traceSnippets = interp1(traceTimes,traces,t0:binSize:t2);
                stimTimes = onset+(0:size(movie,3)-1)/fps;
                movie = permute(interp1(stimTimes',permute(movie,[3 1 2]),(t1:binSize:t2)','linear'),[2 3 1]);
                
                assert(~any(isnan(movie(:))), 'nans found in interpolated movie. Check timing.')

                % compute receptive field maps by spike-triggered averaging
                sz = size(movie);
                movie = fliplr(reshape(movie,sz(1)*sz(2),sz(3)));  % flipped to make work with convolution
                for iTrace=1:nTraces
                    maps{iTrace} = maps{iTrace} + reshape(conv2(movie,traceSnippets(:,iTrace)','valid'),sz(1),sz(2),[])/nTrials;
                end
            end
            fprintf \n
            
            disp inserting..
            for iTrace = 1:nTraces
                tuple = traceKeys(iTrace);
                tuple.noise_map = single(maps{iTrace});
                self.insert(tuple)
            end
            
        end
    end
    
end