function [obj,params] = processDataForFA(obj,params)

% get trials to use based on params.condition defined in main script
params.trialid = findTrials(obj, params.condition);
% get clusters to use based on params.quality defined in main script
params.cluid = findClusters({obj.clu{params.probe}(:).quality}', params.quality);
% align all spike times to params.alignEvent
obj = alignSpikesForFA(obj,params);
% bin single trials and compute PSTHs per condition (PSTHs not used at the
% moment)
obj = getSeqForFA(obj,params);
% remove clusters that have a mean firing rate across all time and trials
% less than params.lowFR
[obj, params.cluid] = removeLowFRClustersForFA(obj,params.cluid,params.lowFR);
% compute presample mean firing rate and standard deviation for each
% cluster (good to have handy for many analyses, including z-scoring)
[obj.presampleFR, obj.presampleSigma] = baselineFRForFA(obj,params);


end