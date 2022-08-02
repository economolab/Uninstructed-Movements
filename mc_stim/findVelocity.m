function [xvel, yvel] = findVelocity(taxis, obj, trialnums, view, feat, alignEv)


traj = obj.traj{view};

featix = findDLCFeatIndex(obj.traj,view,{feat});

xvel = nan(numel(taxis),numel(trialnums));
yvel = nan(numel(taxis),numel(trialnums));

for i = 1:numel(trialnums)                        % For each trial
    trix = trialnums(i);
    try
        if isnan(traj(trix).NdroppedFrames )                       % If the video data from this trial isn't good, skip it
            continue;
        end
    catch
        % hack fix for new data objs that don't have NdrppedFrames or frameTimes
    end

    if ~isfield(traj(trix),'frameTimes')
        traj(trix).frameTimes = (1:size(traj(trix).ts,1)) ./ 400;
    end

    if ~isnan(traj(trix).frameTimes)                           % If the video data from this trial is good...
        if contains(feat,'tongue')
            ts = traj(trix).ts(:, 1:2, featix);
%             ts = fillmissing(ts,'nearest');
        else
            ts = mySmooth(traj(trix).ts(:, 1:2, featix), 21);                                               % get time series for x,y position of feat and current view
            %                     tssmooth = fillmissing(tssmooth,'nearest');
        end

        try
            tsinterp = interp1(traj(trix).frameTimes-obj.bp.ev.(alignEv)(trix), ts, taxis);               % Linear interpolation of x,y position to keep number of time points consistent across trials
        catch
            tsinterp = interp1(traj(trix).frameTimes(1:end-1)-obj.bp.ev.(alignEv)(trix), ts, taxis); % for behav only data objects, frame times aren't exact, so someimtes have an extra index in there
        end
        %         ts = mySmooth(traj(trix).ts(:, 1:2, featix), 21);                                               % get time series for x,y position of feat and current view
        basederiv = median(diff(tsinterp),'omitnan');                                            % Find the median  velocity (aka baseline)
    end

    %Find the difference between the feat velocity and the
    %baseline feature velocity (NOT FOR TONGUE)
    if contains(feat,'tongue') %|| strcmp(traj(trix).featNames{feat},'nose')
%         xvel(:,i) = myDiff(tsinterp(:, 1),1/400);
%         yvel(:,i) = myDiff(tsinterp(:, 2),1/400);
        xvel(:,i) = gradient(tsinterp(:,1));
        yvel(:,i) = gradient(tsinterp(:,2));
    else
%         xvel(:,i) = myDiff(tsinterp(:, 1),1/400)-basederiv(1);
%         yvel(:,i) = myDiff(tsinterp(:, 2),1/400)-basederiv(2);
        xvel(:,i) = gradient(tsinterp(:, 1))-basederiv(1);
        yvel(:,i) = gradient(tsinterp(:, 2))-basederiv(2);
    end

    %         fill missing values for all features except tongue
    if contains(feat,'tongue')
        continue
%         xvel(:,i) = fillmissing(xvel(:,i),'constant',nanmin(xvel(:,i)));
%         yvel(:,i) = fillmissing(yvel(:,i),'constant',nanmin(xvel(:,i)));
    else
        xvel(:,i) = fillmissing(xvel(:,i),'nearest');
        yvel(:,i) = fillmissing(yvel(:,i),'nearest');
    end


end

% 'a'

end  % findVelocity





