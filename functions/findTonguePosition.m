function [xpos, ypos] = findTonguePosition(edges, obj, conditions, met, view, feat,params)

traj = obj.traj{view};                             % Get the video data
xpos = cell(1,numel(conditions));
ypos = cell(1,numel(conditions));

for cond = 1:numel(conditions)
    c = conditions{cond};
    nTrials = numel(met.trialid{c});
    tempx = nan(numel(edges), nTrials);    % (time x num trials in curr condition)
    tempy = nan(numel(edges), nTrials);    % (time x num trials in curr condition)

    for i = 1:nTrials                        % For every trial in the condition
        trix = met.trialid{c}(i);
        if isnan(traj(trix).NdroppedFrames )                       % If the video data from this trial isn't good, skip it
            continue;
        end

        if ~isnan(traj(trix).frameTimes)                           % If the video data from this trial is good...
                        
            ts = traj(trix).ts(:, 1:2, feat);

            tsinterp = interp1(traj(trix).frameTimes-0.5- obj.bp.ev.(params.alignEvent)(trix), ts, edges);               % Linear interpolation of jaw position to keep number of time points consistent across trials
        end
        tempx(:, i) = abs(tsinterp(:, 1));
        tempy(:, i) = abs(tsinterp(:, 2));
    end

    xpos{cond} = tempx;
    ypos{cond} = tempy;
end

end  % findPosition