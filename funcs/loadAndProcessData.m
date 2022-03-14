function [meta,params,obj] = loadAndProcessData(meta,params)

%% LOAD DATA

dat = load(fullfile(meta.datapth,meta.datafn));
obj = dat.obj;


%% STANDARD ROUTINES

% find trials to use
params.trialid = findTrials(obj, params.condition);

% find clusters to use
params.cluid = findClusters({obj.clu{params.probe}(:).quality}', params.quality);

% warp spikes post go cue according to median lick duration for each lick
if params.timeWarp
    [obj,obj.lick] = timeWarp(obj,params);
end

% align spikes in every cluster to an event
obj = alignSpikes(obj,params);

% get trial avg psth and single trial data
obj = getSeq(obj,params);

% remove low fr clusters
[obj, params.cluid] = removeLowFRClusters(obj,params.cluid);

% get mean fr and std dev for each cluster/trial type during presample (baseline fr)
[obj.presampleFR, obj.presampleSigma] = baselineFR(obj,params);


end 
