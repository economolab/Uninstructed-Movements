function [params,obj] = processData(obj,params,prbnum)

%% STANDARD ROUTINES
% find trials to use (only need to do this once)
if prbnum==1 || ~isfield(params,'trialid')
    params.trialid = findTrials(obj, params.condition);
    disp(' ')
    disp('--Trials Found')
    disp(' ')
end

% find clusters to use
params.cluid{prbnum} = findClusters({obj.clu{prbnum}(:).quality}', params.quality);
disp(' ')
disp('--Clusters Found')
disp(' ')

% warp spikes post go cue according to median lick duration for each lick
if params.timeWarp
    [obj,obj.lick] = timeWarp(obj,params,prbnum);
    disp(' ')
    disp('--Time Warp Finished')
    disp(' ')
    
end

% align spikes in every cluster to an event
obj = alignSpikes(obj,params,prbnum);
disp(' ')
disp('--Spikes Aligned')
disp(' ')

% get trial avg psth and single trial data
obj = getSeq(obj,params,prbnum);
disp(' ')
disp('--PSTHs and Single Trial Data Done')
disp(' ')


% remove low fr clusters
[obj, params.cluid{prbnum}] = removeLowFRClusters(obj,params.cluid{prbnum},params.lowFR,prbnum);
disp(' ')
disp('--Removed Low Firing Rate Clusters')
disp(' ')

% get mean fr and std dev for each cluster/trial type during presample (baseline fr)
[obj.presampleFR{prbnum}, obj.presampleSigma{prbnum}] = baselineFR(obj,params,prbnum);
disp(' ')
disp('--Presample Statisitics Calculated')
disp(' ')


end 










