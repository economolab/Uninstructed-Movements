% Function for finding the direction (R/L) selectivity for all cells in a session 
% Right - left trial-averaged PSTH; then normalize all the selectivity
% values to the max selectivity value in the session

% INPUTS: obj = the current data object you want to look at 
% psth_to_use = if you want to find direction selectivity for all
% trials, enter 'obj.psth'. Early move trials: 'obj.trialpsth_Early' or
% early move trials removed: 'obj.trialpsth_noEarly'

% OUTPUT: (time x num cells) struct array.  With normalized direction
% selectivity for each cell over the course of a trial (averaged for the
% whole session)

function selectNorm = findDirectionSelectivity(obj,psth_to_use,smooth)

if strcmp(psth_to_use, 'obj.psth')
    Rpsth = mySmooth(obj.psth(:,:,1),smooth);
    Lpsth = mySmooth(obj.psth(:,:,2),smooth);
elseif strcmp(psth_to_use, 'obj.trialpsth_Early')
    Rpsth = mySmooth(mean(obj.trialpsth_Early{1},3),smooth);
    Lpsth = mySmooth(mean(obj.trialpsth_Early{2},3),smooth);
elseif strcmp(psth_to_use, 'obj.trialpsth_noEarly')
    Rpsth = mySmooth(mean(obj.trialpsth_noEarly{1},3),smooth);
    Lpsth = mySmooth(mean(obj.trialpsth_noEarly{2},3),smooth);
end

DirSelect = Rpsth-Lpsth;                        % Find R-L selectivity for all cells
maxSelect =  max(abs(DirSelect),[],'all');      % Find max directional selectivity
selectNorm = DirSelect./maxSelect;              % Normalize all selectivity values to max selectivity

end