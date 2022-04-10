function [xpos, ypos] = findPosition(taxis, obj, trialnums, view, feat, alignEv)

traj = obj.traj{view};

featix = findDLCFeatIndex(obj.traj,view,{feat});

xpos = nan(numel(taxis),numel(trialnums));
ypos = nan(numel(taxis),numel(trialnums));

for i = 1:numel(trialnums)                        % For each trial
    trix = trialnums(i);
    if isnan(traj(trix).NdroppedFrames )                       % If the video data from this trial isn't good, skip it
        continue;
    end
    
    if ~isnan(traj(trix).frameTimes)                           % If the video data from this trial is good...
        ts = mySmooth(traj(trix).ts(:, 1:2, featix), 21);                                               % get time series for x,y position of feat and current view
        tsinterp = interp1(traj(trix).frameTimes-0.5-obj.bp.ev.(alignEv)(trix), ts, taxis);               % Linear interpolation of x,y position to keep number of time points consistent across trials
        xpos(:,i) = tsinterp(:,1);
        ypos(:,i) = tsinterp(:,2);
    end
end

end  % findPosition

% figure; 
% subplot(1,2,1)
% imagesc(xpos')
% subplot(1,2,2)
% imagesc(ypos')