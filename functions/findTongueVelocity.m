function [xvel, yvel] = findTongueVelocity(edges, obj, conditions, met, view, feat)
traj = obj.traj{view};                             % Get the video data
xvel = cell(1,numel(conditions));
yvel = cell(1,numel(conditions));

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
            if strcmp(traj(trix).featNames{feat},'tongue')
                ts = traj(trix).ts(:, 1:2, feat);                                               % Side-view, up and down position of the jaw, smoothed
            else
                ts = mySmooth(traj(trix).ts(:, 1:2, feat), 21);
            end
            tsinterp = interp1(traj(trix).frameTimes-0.5-mode(obj.bp.ev.goCue), ts, edges);          % Linear interpolation of jaw position to keep number of time points consistent across trials
            basederiv = median(diff(tsinterp),'omitnan');                                            % Find the median jaw velocity (aka baseline)
        end
        
        %Find the difference between the feat velocity and the
        %baseline feature velocity (NOT FOR TONGUE)
        if strcmp(traj(trix).featNames{feat},'tongue') %|| strcmp(traj(trix).featNames{feat},'nose')
            tempx(2:end, i) = abs(diff(tsinterp(:, 1)));
            tempy(2:end, i) = abs(diff(tsinterp(:, 2)));
        else
            tempx(2:end, i) = abs(diff(tsinterp(:, 1))-basederiv(1));
            tempy(2:end, i) = abs(diff(tsinterp(:, 2))-basederiv(2));
        end
    
    end

    xvel{cond} = tempx;
    yvel{cond} = tempy;
end

end  % findVelocity