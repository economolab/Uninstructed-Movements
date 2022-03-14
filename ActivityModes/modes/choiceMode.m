function choicemode = choiceMode(obj,meta,cond,epoch,alignEvent)
% 2. choice mode: defined during delay period
%       ((hitR - missR) + (missL - hitL)) / sqrt(sum(sd for each tt ^2));

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
choicemode = ((mu(:,1)-mu(:,3)) + (mu(:,4)-mu(:,2)))./ sqrt(sum(sd.^2,2));
choicemode(isnan(choicemode)) = 0;
choicemode = choicemode./sum(abs(choicemode)); % (ncells,1)

end % choiceMode