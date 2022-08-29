function psth = getEpochPsth(obj,epochix,trials,meta)
% slice data by trials and epochix
nTrials = max(sum(trials.ix));
nTime = max(diff(epochix'));
psth = nan(nTime,numel(meta.cluid),nTrials,trials.N); % (time,clu,trials,cond)
for cluix = 1:numel(meta.cluid) % for each cluster
    for cnd = 1:trials.N % for each cond
        cndtrid = find(trials.ix(:,cnd));
        for trix = 1:numel(cndtrid) % for every trial in cnd
            trid = cndtrid(trix);
            e1 = epochix(trid,1);
            e2 = epochix(trid,2);
            
            psth(1:diff(epochix(trid,:))+1,cluix,trix,cnd) = ...
                obj.trialpsth(e1:e2,cluix,trid);

        end
    end
end

end % getEpochPsth