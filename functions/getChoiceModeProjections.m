% Function for projecting trials onto choice modes identified with and without early movements 

% INPUTS: rez = variable that contains regular choice mode 
% removeEarly = variable that contains choice mode identified without early
% movement trials
% conditions = conditions that you want to project onto these modes

% OUTPUTS = rez and removeEarly will have a field called 'latentChoice'
% that is a [1 x c] cell array where 'c' is the number of conditions
% In each cell array will contain the projection of the trial-averaged PSTH
% for that condition onto the corresponding mode

function [rez,removeEarly] = getChoiceModeProjections(obj,rez,removeEarly,smooth,conditions)

rez.latentChoice = cell(1,numel(conditions));     % Cell array: 1 x # projection conditions
removeEarly.latentChoice = cell(1,numel(conditions));

    psth = obj.psth - mean(obj.psth,1);     % Mean center the PSTHs
    for j = 1:numel(conditions)             % For each condition...
        cond = conditions(j);
        if ~isempty(rez.choice_mode)
            rez.latentChoice{j} = mySmooth(psth(:,:,cond)*rez.choice_mode,smooth);          % Project the trials from curr condition onto normal mode
            removeEarly.latentChoice{j} = mySmooth(psth(:,:,cond)*removeEarly.choice_mode,smooth);  % Project the trials from curr condition onto mode w/o early move
        end
    end

end %getModeProjections