%{
carfs.GratingResponseSet (computed) # calcium responses for each trace and grating trials
-> carfs.TraceSet
-> vis2p.VisStims
-> carfs.SpaceTime
-----
%}

classdef GratingResponseSet < dj.Relvar & dj.AutoPopulate
    
    properties
        popRel = pro(carfs.TraceSet) * (vis2p.VisStims & 'exp_type="MultDimExperiment"')
    end
    
    methods(Access=protected)
        
        function makeTuples(self, key)
                        
            disp 'fetching processed traces...'
            [traces, traceKeys] = fetchn(carfs.Trace & key, 'trace');
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
            tempFreqs = [stim.params.conditions.speed];
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
                    rec.tempFreq = tempFreqs(condIx);
                    rec.onset = onsets(i);
                    rec.duration = durations(i);
                    ix = find(times>rec.onset & times < rec.onset + rec.duration);
                    if ix>=1
                        rec.responses = sum(traces(ix,:),1);
                        if length(rec.responses) ~= nTraces
                            error 'wrong number of traces'
                        end
                        presentations = [presentations; rec]; %#ok<AGROW>
                    end
                end
            end
            
            disp 'tabulating...'
            duration = median([presentations.duration]);
            [responses, directions, spatialFreqs, tempFreqs] = ...
                dj.struct.tabulate(presentations, 'responses', 'direction', 'spatialFreq', 'tempFreq');
            
            disp 'inserting...'
            for spatialIndex = 1:length(spatialFreqs)
                for tempIndex = 1:length(tempFreqs)
                    k = struct(...
                        'spatial_freq', spatialFreqs(spatialIndex), ...
                        'temp_freq', tempFreqs(tempIndex));
                    inserti(carfs.SpaceTime, k)
                    k = dj.struct.join(key, k);
                    self.inserti(dj.struct.join(k,key));
                    r = cell2mat(...
                        cellfun(@(x) fill(x,nTraces), squeeze(responses(:,spatialIndex,tempIndex,:)),'uni',false));
                    for iTrace = 1:length(traceKeys)
                        insert(carfs.GratingResponses, ...
                            setfield(dj.struct.join(k, traceKeys(iTrace)), ...
                            'spike_responses', single(r(:,:,iTrace)))) %#ok<SFLD>                            
                    end
                end
            end
            disp done
        end
    end
    
end


function ret = fill(x,nTraces)
if isempty(x)
    ret = nan(1,1,nTraces);
else
    if numel(x) ~= nTraces
        error 'wrong dimensions'
    else
        ret = reshape(x,[1,1,nTraces]);
    end
    
end
end
