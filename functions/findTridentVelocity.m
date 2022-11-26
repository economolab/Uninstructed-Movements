% Find jaw velocity at all points during the trial
% INPUT:
% edges = time during trial that you want to find jaw velocity for
% obj = current data object
% OUTPUT = [ntrials x length of trial] array of jaw velocity values (for a single
% trial)

function trident = findTridentVelocity(edges, obj,conditions, met,toplot)
traj = obj.traj{1};                             % Get the video data
trident = cell(1,numel(conditions));

for cond = 1:numel(conditions)
    nTrials = numel(met.trialid{cond});
    temptri = nan(numel(edges), nTrials);    % (time x num trials in curr condition)
    derivthresh = 0.5;
    for i = 1:nTrials                        % For every trial in the condition
        trix = met.trialid{cond}(i);
        if isfield(traj,'NdroppedFrames')
            if isnan(traj(trix).NdroppedFrames )                       % If the video data from this trial isn't good, skip it
                continue;
            end

            if ~isnan(traj(trix).frameTimes)                               % If the video data from this trial is good...
                ts = mySmooth(traj(trix).ts(:, 2, 5), 21);                                               % Side-view, up and down position of the trident, smoothed
                tsinterp = interp1(traj(trix).frameTimes-0.5-mode(obj.bp.ev.goCue), ts, edges);          % Linear interpolation of trident position to keep number of time points consistent across trials
                basederiv = median(diff(tsinterp),'omitnan');                                         % Find the median jaw velocity (aka baseline)
            end

        else
            ts = mySmooth(traj(trix).ts(:, 2, 5), 21);                                               % Side-view, up and down position of the trident, smoothed
            nFrames = length(ts);
            trialLen = nFrames*(1/400);
            frameTimes = 0:(1/400):trialLen;
            tsinterp = interp1(frameTimes(2:end)-obj.bp.ev.goCue(trix), ts, edges);                           % Linear interpolation of trident position to keep number of time points consistent across trials
            basederiv = median(diff(tsinterp),'omitnan');                                            % Find the median jaw velocity (aka baseline)
        end

        %Find the difference between the trident velocity and the
        %baseline trident velocity
        if strcmp(toplot,'prob')
            temptri(2:end, i) = abs(diff(tsinterp)-basederiv)>derivthresh;      % Values > 0 = jaw is moving
        elseif strcmp(toplot,'vel')
            temptri(2:end, i) = abs(diff(tsinterp)-basederiv);              % Values > 0 = jaw is moving
        end
    end
    trident{cond} = temptri;
end

end  % findTridentVelocity