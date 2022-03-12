% Find jaw velocity at all points during the trial
% OUTPUT = [1 x length of trial] array of jaw velocity values (for a single
% trial)

function jaw = findJawVelocity(edges, obj)

jaw = nan(numel(edges), size(obj.trialpsth, 3));

traj = obj.traj{1};
nTrials = size(obj.trialpsth, 3);
for i = 1:nTrials                % For every trial...
    if isnan(traj(i).NdroppedFrames )
        continue;
    end
    if ~isnan(traj(i).frameTimes)
        ts = mySmooth(traj(i).ts(:, 2, 2), 21);
        tsinterp = interp1(traj(i).frameTimes-0.5-mode(obj.bp.ev.goCue), ts, edges);   %Linear interpolation of jaw position to keep number of time points consistent across trials
        basederiv = median(diff(tsinterp),'omitnan');                                         %Find the median jaw velocity (aka baseline)
    end
    %Find the difference between the jaw velocity and the
    %baseline jaw velocity
    jaw(2:end, i) = abs(diff(tsinterp)-basederiv);      % Values > 0 = jaw is moving
end
end  % findJawVelocity