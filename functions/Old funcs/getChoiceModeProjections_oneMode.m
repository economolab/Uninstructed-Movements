function [projection] = getChoiceModeProjections_oneMode(obj,allModes,smooth,conditions)
projection.ALLTrix = cell(1,numel(conditions));     % Cell array: 1 x # projection conditions
projection.NOEarly = cell(1,numel(conditions));
avgPSTH_noEarlyMove = cell(1,numel(conditions));

% psth = obj.psth - mean(obj.psth,1);     % Mean center the PSTHs
psth = obj.psth;     % Mean center the PSTHs
for j = 1:numel(conditions)             % For each condition...
    cond = conditions(j);
    if ~isempty(allModes.choice_mode)
        projection.ALLTrix{j} = mySmooth(psth(:,:,cond)*allModes.choice_mode,smooth);          % Project the trials from curr condition onto normal mode
    end
end

% Re-create PSTH without early move trials and then re-do for loop
for j = 1:numel(conditions)                                         % For each condition...
    avgPSTH_noEarlyMove{j} = mean(obj.trialpsth_noEarly{j},3);      % Find new trial-averaged PSTH
end
psth = cat(3,avgPSTH_noEarlyMove{1},avgPSTH_noEarlyMove{2});        % Make into (time x clusters x conditions)
psth = psth - mean(psth,1);                                         % Mean center

for j = 1:numel(conditions)             % For each condition...
    cond = conditions(j);
    if ~isempty(allModes.choice_mode)
        projection.NOEarly{j} = mySmooth(psth(:,:,cond)*allModes.choice_mode,smooth);          % Project the trials from curr condition onto normal mode
    end
end