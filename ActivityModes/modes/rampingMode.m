function rampingmode = rampingMode(obj,meta,cond,epoch,alignEvent)
% ramping mode: in correct trials
%       (hit_presample - hit_delay) / sqrt(sum(sd for each tt ^2));

% which trials to use for each condition used for finding the mode
trials = getTrialsForModeID(obj,cond);

% find time in each trial corresponding to epoch
sampepochix = nan(obj.bp.Ntrials,2);
delayepochix = nan(obj.bp.Ntrials,2);
for trix = 1:obj.bp.Ntrials
    sampepochix(trix,:)  = findedges(obj.time,obj.bp,meta.dt,epoch{1},trix,alignEvent); % (idx1,idx2)
    delayepochix(trix,:) = findedges(obj.time,obj.bp,meta.dt,epoch{2},trix,alignEvent); % (idx1,idx2)
end

sampEpochMean  = getEpochMean(obj,sampepochix,trials,meta);
delayEpochMean = getEpochMean(obj,delayepochix,trials,meta);

[sampmu,sampsd] = getEpochStats(sampEpochMean,meta,trials);
[delaymu,delaysd] = getEpochStats(delayEpochMean,meta,trials);
sd =  [sampsd delaysd]; % (clu,cond)

rampingmode = (sampmu - delaymu) ./ sqrt(sum(sd.^2,2));
rampingmode(isnan(rampingmode)) = 0;
rampingmode = rampingmode./sum(abs(rampingmode)); % (ncells,1)


end % outcomeMode