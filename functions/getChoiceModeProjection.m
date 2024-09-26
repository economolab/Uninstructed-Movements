% Function for projecting trials onto choice modes identified with and without early movements

% INPUTS: allModes = variable that contains regular choice mode
% allModes_NoMove = variable that contains choice mode identified without early
% movement trials
% conditions = conditions that you want to project onto these modes

% OUTPUTS = allModes and allModes_NoMove will have a field called 'latentChoice'
% that is a [1 x c] cell array where 'c' is the number of conditions
% In each cell array will contain the projection of the trial-averaged PSTH
% for that condition onto the corresponding mode

function latentChoice = getChoiceModeProjection(obj,choice_mode,smooth,conditions)
latentChoice = cell(1,numel(conditions));     % Cell array: 1 x # projection conditions

%psth = obj.psth - mean(obj.psth,1);     % Mean center the PSTHs
psth = obj.psth;

for j = 1:numel(conditions)             % For each condition...
    cond = conditions(j);
    if ~isempty(choice_mode)
        latentChoice{j} = mySmooth(psth(:,:,cond)*choice_mode,smooth);          % Project the trials from curr condition onto normal mode 
    end
end

end %getModeProjections