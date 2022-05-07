function cd = calcCD(obj,meta,cond,epoch,alignEvent)
% 2. choice mode: defined during delay period
%       ((hitR  - hitL)) / sqrt(sum(sd for each tt ^2));

% which trials to use for each condition used for finding the mode
trials = getTrialsForModeID(obj,cond);

% find time in each trial corresponding to epoch
epochix = nan(obj.bp.Ntrials,2);
for trix = 1:obj.bp.Ntrials
    epochix(trix,:) = findedges(obj.time,obj.bp,meta.dt,epoch,trix,alignEvent); % (idx1,idx2)
end

epochMean = getEpochMean(obj,epochix,trials,meta);

[mu,sd] = getEpochStats(epochMean,meta,trials);

% calculate mode according to definition
cd = ((mu(:,1)-mu(:,2)))./ sqrt(sum(sd.^2,2));
cd(isnan(cd)) = 0;
cd = cd./sum(abs(cd)); % (ncells,1)

end % choiceMode