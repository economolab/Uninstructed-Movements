function actionmode = actionModeMulti(objs,meta,cond,epoch,alignEvent)
% action mode: defined during mvmt init (0.1 to 0.3 s rel to go cue)
%       (hitR - hitL) / sqrt(sum(sd for each tt ^2));

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

% calculate mode according to definition
actionmode = (mu(:,1)-mu(:,2))./ sqrt(sum(sd.^2,2));
actionmode(isnan(actionmode)) = 0;
actionmode = actionmode./sum(abs(actionmode)); % (ncells,1)


end % actionMode