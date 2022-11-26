function [corr,lags] = findCrossCorrelations(orthproj,kin,kinfns,numlags)

nTrials = size(orthproj.jawPos,2);
for i = 1:numel(kinfns)                 % For each mode...
    currproj = orthproj.(kinfns{i});
    currfeat = kin.(kinfns{i});
    R = NaN((2*numlags)+1,nTrials);   % Store correlations for each trial (numlags+1 x numTrials)
    for ii = 1:nTrials                % For each trial...
        [R(:,ii),lags] = xcorr(currproj(:,ii),currfeat(:,ii),numlags);  % Find the cross-correlation between the projection and the feature position
    end
    corr.(kinfns{i}) = mean(R,2,'omitnan');     % Average the cross-correlation values across trials
end

end