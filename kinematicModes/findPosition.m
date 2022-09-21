function [xpos, ypos] = findPosition(edges, obj, conditions, met, view, feat,params)

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
            if strcmp(traj(trix).featNames{feat},'tongue')         % Don't smooth the DLC trajectory if the feature is the tongue
                ts = traj(trix).ts(:, 1:2, feat);
            else
                ts = mySmooth(traj(trix).ts(:, 1:2, feat), 21);                                               % Specified view, x and y position
            end
            tsinterp = interp1(traj(trix).frameTimes-0.5 - obj.bp.ev.(params.alignEvent)(trix), ts, edges);               % Linear interpolation of jaw position to keep number of time points consistent across trials
            basederiv = median(tsinterp(1:100, :),'omitnan');                                             % Find the median jaw velocity (aka baseline)
        end
        %Find the difference between the feat velocity and the
        %baseline feat velocity (IF NOT TONGUE)
        if strcmp(traj(trix).featNames{feat},'tongue')% || strcmp(traj(trix).featNames{feat},'nose')
            tempx(:, i) = abs(tsinterp(:, 1));      
            tempy(:, i) = abs(tsinterp(:, 2));
        else
            tempx(:, i) = abs(tsinterp(:, 1)-basederiv(1));      
            tempy(:, i) = abs(tsinterp(:, 2)-basederiv(2));  
        end
    end

    xpos{cond} = tempx;
    ypos{cond} = tempy;
end

end  % findPosition