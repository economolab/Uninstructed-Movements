% Function for finding single trial projections onto a specific activity
% mode

% INPUT: cd = the activity mode that you want to project onto (i.e.
% rez.choice_mode)
% OUTPUT: lat_single = [time x nTrials] struct array 

function latent = getTrialLatents(obj,cd,conditions,met)
smooth = 15;
latent = cell(1,numel(conditions));
for cond = 1:numel(conditions)
    nTrials = numel(met.trialid{cond});
    temp = nan(size(obj.trialdat, 1), nTrials);  % time x num trials
    for j = 1:nTrials                % For each trial...
        trix = met.trialid{cond}(j);
        ts = obj.trialdat(:,:,trix);                  % Get the PSTHs for each neuron
        temp(:, j) = mySmooth(ts*cd,smooth);           % Multiply PSTHs by coding dimension and smooth it
    end
    latent{cond} = temp;
end
end    %getTrialLatents