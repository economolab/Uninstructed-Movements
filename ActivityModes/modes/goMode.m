function gomode = goMode(obj,meta,cond,epoch,alignEvent)
% go mode: 0.1 sec before or after go cue
%       (hit_after - hit_before) / sqrt(sum(sd for each tt ^2));

% which trials to use for each condition used for finding the mode
trials = getTrialsForModeID(obj,cond);

% find time in each trial corresponding to epoch
postepochix = nan(obj.bp.Ntrials,2);
preepochix = nan(obj.bp.Ntrials,2);
for trix = 1:obj.bp.Ntrials
    postepochix(trix,:)  = findedges(obj.time,obj.bp,meta.dt,epoch{1},trix,alignEvent); % (idx1,idx2)
    preepochix(trix,:) = findedges(obj.time,obj.bp,meta.dt,epoch{2},trix,alignEvent); % (idx1,idx2)
end

postEpochMean  = getEpochMean(obj,postepochix,trials,meta);
preEpochMean = getEpochMean(obj,preepochix,trials,meta);

[postmu,postsd] = getEpochStats(postEpochMean,meta,trials);
[premu,presd] = getEpochStats(preEpochMean,meta,trials);
sd =  [postsd presd]; % (clu,cond)

gomode = (postmu - premu) ./ sqrt(sum(sd.^2,2));
gomode(isnan(gomode)) = 0;
gomode = gomode./sum(abs(gomode)); % (ncells,1)


end % outcomeMode