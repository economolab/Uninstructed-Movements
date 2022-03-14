function stimmode = stimulusMode(obj,params,cond,epoch,alignEvent)
% 1. stimulus mode: defined during stimulus (sample) period
%       ((hitR - missL) + (missR - hitL)) / sqrt(sum(sd for each tt ^2));

% which trials to use for each condition used for finding the mode
trials = getTrialsForModeID(obj,cond);

% find time in each trial corresponding to epoch
epochix = nan(obj.bp.Ntrials,2);
for trix = 1:obj.bp.Ntrials
    epochix(trix,:) = findedges(obj.time,obj.bp,params.dt,epoch,trix,alignEvent); % (idx1,idx2)
end

epochMean = getEpochMean(obj,epochix,trials,params);

[mu,sd] = getEpochStats(epochMean,params,trials);


% calculate mode according to definition
stimmode = ((mu(:,1)-mu(:,4)) + (mu(:,3)-mu(:,2)))./ sqrt(sum(sd.^2,2));
stimmode(isnan(stimmode)) = 0;
stimmode = stimmode./sum(abs(stimmode)); % (ncells,1)

end % stimulusMode





