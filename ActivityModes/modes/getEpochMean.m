function epochMean = getEpochMean(obj,epochix,trials,params)
% slice data by trials and epochix
% compute the mean during the epoch for each trial, cluster
nTrials = max(sum(trials.ix));
nTime = max(diff(epochix'));
psth = zeros(nTime,numel(params.cluid),nTrials,trials.N); % (time,clu,trials,cond)
epochMean = zeros(numel(params.cluid),nTrials,trials.N);
for cluix = 1:numel(params.cluid) % for each cluster
    for cnd = 1:trials.N % for each cond
        cndtrid = find(trials.ix(:,cnd));
        for trix = 1:numel(cndtrid) % for every trial in cnd
            trid = cndtrid(trix);
            e1 = epochix(trid,1);
            e2 = epochix(trid,2);
            
            psth(1:diff(epochix(trid,:))+1,cluix,trix,cnd) = ...
                obj.trialdat(e1:e2,cluix,trid);
            
            % calculate the avg firing rate during the epoch for trial
            epochMean(cluix,trix,cnd) = nanmean(psth(:,cluix,trix,cnd),1);
        end
    end
end

end % getEpochMean