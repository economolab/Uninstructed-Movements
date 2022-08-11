function [xpos, ypos] = findPosition(taxis, obj, trialnums, view, feat, alignEv, warp)

traj = obj.traj{view};

featix = findDLCFeatIndex(obj.traj,view,{feat});

xpos = nan(numel(taxis),numel(trialnums));
ypos = nan(numel(taxis),numel(trialnums));

for i = 1:numel(trialnums)                        % For each trial
    trix = trialnums(i);

    if isnan(traj(trix).NdroppedFrames )                       % If the video data from this trial isn't good, skip it
        continue;
    end


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

    xpos(:,i) = tsinterp(:,1);
    ypos(:,i) = tsinterp(:,2);



    %         fill missing values for all features except tongue
    if contains(feat,'tongue')
        continue
    else
        xpos(:,i) = fillmissing(xpos(:,i),'nearest');
        ypos(:,i) = fillmissing(ypos(:,i),'nearest');
    end

end

end  % findPosition
