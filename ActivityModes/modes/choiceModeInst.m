function choicemodeinst = choiceModeInst(obj,meta,cond,epoch,alignEvent)
% 2. choice mode: defined during delay period
%       ((hitR - missR) + (missL - hitL)) / sqrt(sum(sd for each tt ^2));

% which trials to use for each condition used for finding the mode
trials = getTrialsForModeID(obj,cond);

% find time in each trial corresponding to epoch
epochix = nan(obj.bp.Ntrials,2);
for trix = 1:obj.bp.Ntrials
    epochix(trix,:) = findedges(obj.time,obj.bp,meta.dt,epoch,trix,alignEvent); % (idx1,idx2)
end

epochPsth = getEpochPsth(obj,epochix,trials,meta); % (time,clu,trials,cond)

epochPsth = squeeze(nanmean(epochPsth,3));

for t = 1:size(epochPsth,1)
    choicemodeinst(t,:) = (epochPsth(t,:,1) - epochPsth(t,:,3)) + (epochPsth(t,:,4) - epochPsth(t,:,2));
end
choicemodeinst(isnan(choicemodeinst)) = 0;
choicemodeinst = choicemodeinst ./ sum(choicemodeinst);
choicemodeinst = pca(choicemodeinst - mean(choicemodeinst),'NumComponents',3);


end % choiceModeInst