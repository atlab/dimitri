%{
carfs.ImageMap (computed) # receptive fields computed from natural img
-> carfs.ImageMapSet
->carfs.Trace
-----
image_map : longblob  # receptive field
%}

classdef ImageMap < dj.Relvar
    
    methods
        
        function makeTuples(self, key)
            
            disp 'fetching traces...'
            times = fetch1(carfs.TraceSet & key,'trace_times');
            [traces, traceKeys] = fetchn(carfs.Trace & key, 'trace');
            traces = [traces{:}];
            nTraces = size(traces,2);
            
            disp 'loading image database...'
            img = load(getLocalPath('/lab/users/Manolis/Matlab/stimulation/stimuli/vanhateren_iml/imageDatabase'));
            sources = {img.imageDatabase.filename};
            img = cat(3,img.imageDatabase.image);
            
            disp 'averaging responses...'
            stimFile = fetch1(vis2p.VisStims & key & 'exp_type="NatImExperiment"', 'stim_filename');
            stim = load(getLocalPath(stimFile));
            stim = stim.stim;
            nTrials = length(stim.events);
            imgIdx = nan(nTrials,1);
            onsets = nan(nTrials,1);
            offsets = nan(nTrials,1);
            for iTrial = 1:nTrials
                imgIdx(iTrial) = find(strcmp(sources,stim.params.conditions(stim.params.trials(iTrial).condition).source));
                onset = stim.events(iTrial).syncedTimes(ismember(stim.eventTypes(stim.events(iTrial).types),'showStimulus'))/1000;
                offset = stim.events(iTrial).syncedTimes(ismember(stim.eventTypes(stim.events(iTrial).types),'endStimulus'))/1000;
                assert(isscalar(onset) && isscalar(offset))
                onsets(iTrial) = onset + 0.15;  % 50 ms latency
                offsets(iTrial) = onset + 0.4;
            end
            cumTraces = cumsum(traces);
            responses =  interp1(times',cumTraces,offsets) - interp1(times',cumTraces,onsets);
            
            disp 'whitening images...'
            img = double(img(:,:,imgIdx));  % select the right images
            sz = size(img);
            img = bsxfun(@minus  ,img,reshape(mean(reshape(img,[],sz(3))),1,1,sz(3)));  % normalize by mean
            img = reshape(pinv(reshape(img,sz(1)*sz(2),nTrials))',sz);         % whiten sequence
            
            disp 'computing RFs...'
            R = reshape(reshape(img,prod(sz(1:2)),nTrials)*responses,sz(1),sz(2),nTraces);
            
            disp 'inserting'
            tuples = dj.struct.join(key,traceKeys);
            for i=1:length(tuples)
                tuple = tuples(i);
                tuple.image_map = single(R(:,:,i));
                self.insert(tuple)
            end
        end
    end
    
end