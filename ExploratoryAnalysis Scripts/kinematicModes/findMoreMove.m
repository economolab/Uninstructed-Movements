% Identify the condition in which the animal moves more
function [modecond,modekin] = findMoreMove(conditions,kin,kinfeat,params,obj,times)
% Pre-allocate values to be stored for each condition
jv = cell(1,numel(conditions));
maxmin = NaN(1,numel(conditions));
stdeviation = NaN(1,numel(conditions));
avg = NaN(1,numel(conditions));
cnt = 0;

% Define time-points over which you want to find the average kinematic
% measure
e1 = find(obj.time>times(1),1,'first');
e2 = find(obj.time>times(2),1,'first');

kinix = find(strcmp(kin.featLeg,kinfeat));
% For each condition, find the average kinematic measure during the
% specified time epoch
for cc = 1:numel(conditions)
    cond = conditions(cc);
    condtrix = params.trialid{cond};
    % Grab all of the kinematic feature values for this condition
    jv{cc} = kin.dat_std(:,condtrix,kinix);

    % Find the average during the specified time epoch on each trial
    jv_avg{cc} = mean(jv{cc}(e1:e2,:),1);

    maxmin(cc) = max(jv_avg{cc})-min(jv_avg{cc});  % Find the range of kin measure values in the current condition
    stdeviation(cc) = std(jv_avg{cc});         % Find the stdev of kin measure values in the current condition
    avg(cc) = mean(jv_avg{cc});                % Find the average of kin measure values in the current condition
end
touse = NaN(1,3);
% Find the condition with the largest range, stdev, and average kinematic
% measure 
touse(1) = find(max(maxmin)); touse(2) = find(max(stdeviation)); touse(3) = find(max(avg));
if sum(touse)==3                            % If the largest range, stdev, and avg occurs in the first condition...
    modecond = 1;                           % Use condition 1 to find the kinematic modes
elseif sum(touse)==6                        % If condition 2 has the larger distribution of kinematic measure...
    modecond = 2;                           % Use condition 2 to find the kinematic modes
else                                        % If there are discrepancies between which condition has the widest distribution...
    modecond = touse(3);                    % Just choose the condition with the largest average kinematic measure   
end
modekin = jv{modecond};
end