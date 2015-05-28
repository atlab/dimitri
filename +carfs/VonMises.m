%{
carfs.VonMises (computed) # my newest table
-> carfs.VonMisesSet
-> carfs.Trace
-> carfs.SpatialFreq
---
von_p_value                 : float                         # von mises shuffle p-values
von_r2                      : float                         # R-squared of response
von_pref                    : float                         # von mises preferred direction
von_base                    : float                         # von mises bases
von_amp1                    : float                         # von mises preferred amplitude
von_amp2                    : float                         # von mises anti-preferred amplitude
von_sharp                   : float                         # von mises sharpness
von_osi                     : float                         # orientation selectivity index (pref-ortho)/(pref+ortho)
von_dsi                     : float                         # direction selectivity index (pref-anti)/(pref+anti)
%}

classdef VonMises < dj.Relvar
    methods
        
        function makeTuples(self, key)
            
            tuple = key;
            
            disp 'fetching processed traces...'
            [tuple.cellnums, tuple.celltypes, traces, traceKeys] = ...
                fetchn(carfs.Trace & key, 'masknum', 'mask_type', 'trace');
            nTraces = length(traceKeys);
            traces = double([traces{:}]);
            times = fetch1(carfs.TraceSet & rmfield(key,'stim_idx'), 'trace_times');            
            assert(size(traces,1)==length(times), 'synchronization error')
                        
            visual_latency = 30;  % ms
            times = times - visual_latency/1000;
            
            disp 'parsing stimulus...'
            stimFile = fetch1(vis2p.VisStims & key, 'stim_filename');
            stim = load(getLocalPath(stimFile));
            stim = stim.stim;
            directions = [stim.params.conditions.orientation];
            spatialFreqs = [stim.params.conditions.spatialFreq];
            if numel(unique(diff(unique(directions))))>1
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
                    condIx = stim.params.trials(iTrial).conditions(i);
                    rec.direction = directions(condIx);
                    rec.spatialFreq = spatialFreqs(condIx);
                    rec.onset = onsets(i);
                    rec.duration = durations(i);
                    rec.responses = sum(traces(times>rec.onset & times < rec.onset + rec.duration, :));
                    presentations = [presentations; rec]; %#ok<AGROW>
                end
            end
            
            disp 'tabulating...'
            duration = median([presentations.duration]);
            [responses, directions, spatialFreqs] = dj.struct.tabulate(presentations, 'responses', 'direction', 'spatialFreq');
            
            for spatialIndex = 1:length(spatialFreqs)
                r = cellfun(@(x) fill(x,nTraces), squeeze(responses(:,spatialIndex,:)),'uni',false);
                r = permute(cell2mat(r),[3 1 2]);
                
                disp 'computing von Mises tuning...'
                [von, r2, p] = ne7.rf.VonMises2.computeSignificance(r, key.nshuffles);
                for iTrace = 1:nTraces
                    tuple = traceKeys(iTrace);
                    tuple.spatial_freq = spatialFreqs(spatialIndex);
                    tuple.von_r2 = r2(iTrace);
                    tuple.von_base = von.w(iTrace,1);
                    tuple.von_amp1 = von.w(iTrace,2);
                    tuple.von_amp2 = von.w(iTrace,3);
                    tuple.von_sharp= von.w(iTrace,4);
                    tuple.von_pref = von.w(iTrace,5);
                    tuple.von_p_value = p(iTrace);
                    
                    pref = tuple.von_base + tuple.von_amp1 + tuple.von_amp2*exp(-2*tuple.von_sharp);
                    anti = tuple.von_base + tuple.von_amp2 + tuple.von_amp1*exp(-2*tuple.von_sharp);
                    orth = tuple.von_base + (tuple.von_amp1 + tuple.von_amp2)*exp(-tuple.von_sharp);
                    
                    tuple.von_osi = (pref-orth)/(pref+orth);
                    tuple.von_dsi = (pref-anti)/(pref+anti);
                    self.insert(tuple)
                end
            end
            
            
        end
    end
    
end


function ret = fill(x,nTraces)
if isempty(x)
    ret = nan(1,1,nTraces);
else
    ret = reshape(x,[1,1,nTraces]);
end
end
