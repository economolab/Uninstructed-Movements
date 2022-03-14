function outcomemode = outcomeModeMulti(objs,meta,cond,epoch,alignEvent)
% outcome mode: defined during response epoch (0 to 1.3 s rel go cue)
%       ((hitR - missR) + (hitL - missL)) / sqrt(sum(sd for each tt ^2));

mu = cell(numel(objs),1);
sd = cell(numel(objs),1);
for i = 1:numel(objs)
    % which trials to use for each condition used for finding the mode
    trials = getTrialsForModeID(objs{i},cond);
    
    % find time in each trial corresponding to epoch
    epochix = nan(objs{i}.bp.Ntrials,2);
    for trix = 1:objs{i}.bp.Ntrials
        epochix(trix,:) = findedges(objs{i}.time,objs{i}.bp,meta(i).dt,epoch,trix,alignEvent); % (idx1,idx2)
    end
    
    epochMean = getEpochMean(objs{i},epochix,trials,meta(i));
    
    [mu{i},sd{i}] = getEpochStats(epochMean,meta(i),trials);
    
end

mu = cell2mat(mu);
sd = cell2mat(sd);

outcomemode = ((mu(:,1)-mu(:,3)) + (mu(:,2)-mu(:,4)))./ sqrt(sum(sd.^2,2));
outcomemode(isnan(outcomemode)) = 0;
outcomemode = outcomemode./sum(abs(outcomemode)); % (ncells,1)


end % outcomeMode