function [xvel, yvel] = findVelocity(taxis, obj, trialnums, view, feat, alignEv, warp)


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
        if contains(feat,'tongue')
            ts = traj(trix).ts(:, 1:2, featix);
        else
            ts = mySmooth(traj(trix).ts(:, 1:2, featix), 21);                                               % get time series for x,y position of feat and current view
        end

        if ~warp
            try
                tsinterp = interp1(traj(trix).frameTimes - 0.5 - obj.bp.ev.(alignEv)(trix), ts, taxis);               % Linear interpolation of x,y position to keep number of time points consistent across trials
            catch
                tsinterp = interp1(traj(trix).frameTimes(1:end-1)  - 0.5 - obj.bp.ev.(alignEv)(trix), ts, taxis);
            end
        else
            try
                tsinterp = interp1(traj(trix).frameTimes_warped - obj.bp.ev.(alignEv)(trix), ts, taxis); % for behav only data objects, frame times aren't exact, so someimtes have an extra index in there
            catch
                tsinterp = interp1(traj(trix).frameTimes_warped(1:end-1) - obj.bp.ev.(alignEv)(trix), ts, taxis);
            end
        end
        basederiv = median(diff(tsinterp),'omitnan');                                            % Find the median  velocity (aka baseline)
    end

    %Find the difference between the feat velocity and the
    %baseline feature velocity (NOT FOR TONGUE)
    if contains(feat,'tongue') %|| strcmp(traj(trix).featNames{feat},'nose')
        xvel(:,i) = gradient(tsinterp(:,1));
        yvel(:,i) = gradient(tsinterp(:,2));
    else
        xvel(:,i) = gradient(tsinterp(:, 1))-basederiv(1);
        yvel(:,i) = gradient(tsinterp(:, 2))-basederiv(2);
    end

    %         fill missing values for all features except tongue
    if contains(feat,'tongue')
        continue
    else
        xvel(:,i) = fillmissing(xvel(:,i),'nearest');
        yvel(:,i) = fillmissing(yvel(:,i),'nearest');
    end


end

% 'a'

end  % findVelocity





