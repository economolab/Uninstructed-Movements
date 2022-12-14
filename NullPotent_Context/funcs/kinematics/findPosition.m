function [xpos, ypos] = findPosition(taxis, obj, nTrials, view, feat, alignEv)

traj = obj.traj{view};

featix = findDLCFeatIndex(obj.traj,view,{feat});

xpos = nan(numel(taxis),nTrials);
ypos = nan(numel(taxis),nTrials);

for i = 1:nTrials % for each trial

    trix = i;

    % check if video data for trial is good, skip if not
    if isfield(traj(trix),'NdroppedFrames')
        if isnan(traj(trix).NdroppedFrames )
            continue;
        end
    end

    % check if frameTimes exists, if not, create it
    if ~isfield(traj(trix),'frameTimes')
        traj(trix).frameTimes = (1:size(traj(trix).ts,1)) ./ 400;
    end

    % get feat trajs for valid trials
    if ~isnan(traj(trix).frameTimes)
        ts = traj(trix).ts(:,1:2,featix);
        % if tongue, don't smooth
        if ~contains(feat,'tongue')
            ts = mySmooth(ts, 21);
        end
        % linear interpolation of frameTimes to taxis
        tsinterp = interp1(traj(trix).frameTimes-0.5-obj.bp.ev.(alignEv)(trix), ts, taxis);
    end

    xpos(:,i) = tsinterp(:,1);
    ypos(:,i) = tsinterp(:,2);



    % fill missing values for all features except tongue
    if ~contains(feat,'tongue') || ~contains(feat,'paw')
        xpos(:,i) = fillmissing(xpos(:,i),'nearest');
        ypos(:,i) = fillmissing(ypos(:,i),'nearest');
    end

end


end  % findPosition
