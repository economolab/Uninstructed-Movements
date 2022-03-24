function [xvel, yvel] = findVelocity(edges, obj, conditions, met, view, feat)
traj = obj.traj{view};                             % Get the video data
xvel = cell(1,numel(conditions));
yvel = cell(1,numel(conditions));

for cond = 1:numel(conditions)
    nTrials = numel(met.trialid{cond});
    tempx = nan(numel(edges), nTrials);    % (time x num trials in curr condition)
    tempy = nan(numel(edges), nTrials);    % (time x num trials in curr condition)

    for i = 1:nTrials                        % For every trial in the condition
        trix = met.trialid{cond}(i);                      
        if isnan(traj(trix).NdroppedFrames )                       % If the video data from this trial isn't good, skip it
            continue;
        end

        if ~isnan(traj(trix).frameTimes)                           % If the video data from this trial is good...
            ts = mySmooth(traj(trix).ts(:, 1:2, feat), 21);                                               % Side-view, up and down position of the jaw, smoothed
            tsinterp = interp1(traj(trix).frameTimes-0.5-mode(obj.bp.ev.goCue), ts, edges);          % Linear interpolation of jaw position to keep number of time points consistent across trials
            basederiv = median(diff(tsinterp),'omitnan');                                         % Find the median jaw velocity (aka baseline)
        end
        %Find the difference between the jaw velocity and the
        %baseline feature velocity
        tempx(2:end, i) = abs(diff(tsinterp(:, 1))-basederiv(1));      
        tempy(2:end, i) = abs(diff(tsinterp(:, 2))-basederiv(2));      
    end

    xvel{cond} = tempx;
    yvel{cond} = tempy;
end

end  % findVelocity