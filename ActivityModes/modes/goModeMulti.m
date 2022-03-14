function gomode = goModeMulti(objs,meta,cond,epoch,alignEvent)
% go mode: 0.1 sec before or after go cue
%       (hit_after - hit_before) / sqrt(sum(sd for each tt ^2));

postmu = cell(numel(objs),1);
premu = cell(numel(objs),1);
sd = cell(numel(objs),1);
for i = 1:numel(objs)
    % which trials to use for each condition used for finding the mode
    trials = getTrialsForModeID(objs{i},cond);
    
    % find time in each trial corresponding to epoch
    postepochix = nan(objs{i}.bp.Ntrials,2);
    preepochix = nan(objs{i}.bp.Ntrials,2);
    for trix = 1:objs{i}.bp.Ntrials
        postepochix(trix,:)  = findedges(objs{i}.time,objs{i}.bp,meta(i).dt,epoch{1},trix,alignEvent); % (idx1,idx2)
        preepochix(trix,:)   = findedges(objs{i}.time,objs{i}.bp,meta(i).dt,epoch{2},trix,alignEvent); % (idx1,idx2)
    end
    
    postEpochMean  = getEpochMean(objs{i},postepochix,trials,meta(i));
    preEpochMean   = getEpochMean(objs{i},preepochix,trials,meta(i));
    
    [postmu{i},postsd] = getEpochStats(postEpochMean,meta(i),trials);
    [premu{i},presd]   = getEpochStats(preEpochMean,meta(i),trials);
    sd{i} =  [postsd presd]; % (clu,cond)
    
end

postmu = cell2mat(postmu);
premu = cell2mat(premu);
sd = cell2mat(sd);


gomode = (postmu - premu) ./ sqrt(sum(sd.^2,2));
gomode(isnan(gomode)) = 0;
gomode = gomode./sum(abs(gomode)); % (ncells,1)


end % outcomeMode