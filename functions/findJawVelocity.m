% Find jaw velocity at all points during the trial
% INPUT:
% edges = time during trial that you want to find jaw velocity for
% obj = current data object
% OUTPUT = [ntrials x length of trial] array of jaw velocity values (for a single
% trial)

function jaw = findJawVelocity(edges, obj,conditions, met)
traj = obj.traj{1};                             % Get the video data
jaw = cell(1,numel(conditions));

for cond = 1:numel(conditions)
    nTrials = numel(met.trialid{cond});
    tempjaw = nan(numel(edges), nTrials);    % (time x num trials in curr condition)

    for i = 1:nTrials                        % For every trial in the condition
        trix = met.trialid{cond}(i);                      
        if isnan(traj(trix).NdroppedFrames )                       % If the video data from this trial isn't good, skip it
            continue;
        end

        if ~isnan(traj(trix).frameTimes)                           % If the video data from this trial is good...
            ts = mySmooth(traj(trix).ts(:, 2, 2), 21);                                               % Side-view, up and down position of the jaw, smoothed
            tsinterp = interp1(traj(trix).frameTimes-0.5-mode(obj.bp.ev.goCue), ts, edges);          % Linear interpolation of jaw position to keep number of time points consistent across trials
            basederiv = median(diff(tsinterp),'omitnan');                                         % Find the median jaw velocity (aka baseline)
        end
        %Find the difference between the jaw velocity and the
        %baseline jaw velocity
        tempjaw(2:end, i) = abs(diff(tsinterp)-basederiv);      % Values > 0 = jaw is moving
    end

    jaw{cond} = tempjaw;
end

end  % findJawVelocity