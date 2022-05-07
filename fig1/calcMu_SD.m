function [mu, sd] = calcMu_SD(obj,params,cond,epoch,alignEvent,sessix)

% which trials to use for each condition used for finding the mode
trials = getTrialsForModeID(obj,cond);

% find time in each trial corresponding to epoch
epochix = nan(obj.bp.Ntrials,2);
for trix = 1:obj.bp.Ntrials
    epochix(trix,:) = findedges(obj.time,obj.bp,params.dt,epoch,trix,alignEvent); % (idx1,idx2)
end

epochMean = getEpochMean(obj,epochix,trials,params,sessix);

[mu,sd] = getEpochStats(epochMean,params,trials,sessix);

end 