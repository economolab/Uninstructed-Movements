function [xvel, yvel] = findVelocity(taxis, obj, trialnums, view, feat, alignEv)


traj = obj.traj{view};

featix = findDLCFeatIndex(obj.traj,view,{feat});

xvel = nan(numel(taxis),numel(trialnums));
yvel = nan(numel(taxis),numel(trialnums));

for i = 1:numel(trialnums)                        % For each trial
    trix = trialnums(i);
    if isnan(traj(trix).NdroppedFrames )                       % If the video data from this trial isn't good, skip it
        continue;
    end
    
    if ~isnan(traj(trix).frameTimes)                           % If the video data from this trial is good...
        ts = mySmooth(traj(trix).ts(:, 1:2, featix), 21);                                               % get time series for x,y position of feat and current view
        tsinterp = interp1(traj(trix).frameTimes-0.5-obj.bp.ev.(alignEv)(trix), ts, taxis);               % Linear interpolation of x,y position to keep number of time points consistent across trials
        basederiv = median(diff(tsinterp),'omitnan');                                            % Find the median  velocity (aka baseline)
    end
    
    %Find the difference between the feat velocity and the
    %baseline feature velocity (NOT FOR TONGUE)
    if contains(feat,'tongue') %|| strcmp(traj(trix).featNames{feat},'nose')
        xvel(:,i) = myDiff(tsinterp(:, 1),1/400);
        yvel(:,i) = myDiff(tsinterp(:, 2),1/400);
    else
        xvel(:,i) = myDiff(tsinterp(:, 1),1/400)-basederiv(1);
        yvel(:,i) = myDiff(tsinterp(:, 2),1/400)-basederiv(2);
    end
    
    
end

end  % findVelocity





