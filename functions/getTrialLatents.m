% Function for finding single trial projections onto a specific activity
% mode

% INPUT: cd = the activity mode that you want to project onto (i.e.
% rez.choice_mode)
% OUTPUT: lat_single = [time x nTrials] struct array 

function lat_single = getTrialLatents(obj,cd)
smooth = 15;
             
latent = nan(size(obj.trialpsth, 1), size(obj.trialpsth, 3));  % time x num trials
for j = 1:size(obj.trialpsth, 3)                % For each trial...
    ts = obj.trialpsth(:,:,j);                  % Get the PSTHs for each neuron
    latent(:, j) = mySmooth(ts*cd,smooth);      % Multiply PSTHs by coding dimension and smooth it
end
lat_single = latent;
end    %getTrialLatents