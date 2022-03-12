% Function for finding trial indices where the tongue is moving for more
% than 15% of the delay period

% Tongue is moving when it is visible for at least 30 frames during the
% delay period

function obj = findEarlyTongue(obj)

currtraj = obj.traj{1};     % Get side-view video trajectories from current data object
earlyTongueix = [];
for j = 1:size(currtraj,2)            % Go through all trials...
    tonguex = currtraj(j).ts(:,1,1);                % Get x-values for tongue trajectories for current trial
    goCue = obj.bp.ev.goCue(j);                     % Find time of goCue for current trial
    preGoix = find(currtraj(j).frameTimes < goCue); % Find frame indices that occur before the goCue
    earlytongue = find(~isnan(tonguex(preGoix)));   % In pre-goCue frames, find frames that have tongue visible
    for e = 1:numel(earlytongue)
        ix = earlytongue(e);
        if (ix+30 < numel(tonguex))
            if isnan(tonguex(ix+30))                    % For each index that the tongue is visible, check to see that it is visible for 30 frames
                earlytongue(e) = 0;
            end
        end
    end

    moveThresh = 0.15*length(preGoix);                  % Find amount of frames that is 15% of the pre-Go cue frame times 

    if sum(earlytongue)>moveThresh                     % If the tongue is visible for more than 15% of the frames before the go-Cue...
        earlyTongueix = [earlyTongueix,j];              % Store the trial number
    end
    t = currtraj(j).frameTimes - 0.5;
    %figure(i);
    %         plot(t,10*j+tonguex); hold on;
end
obj.earlyTongueix = earlyTongueix';             % Store the trial numbers for current object that have early tongue movements
end  %findEarlyTongue