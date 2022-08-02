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
            tsinterp = interp1(traj(trix).frameTimes - 0.5 -obj.bp.ev.(alignEv)(trix), ts, taxis);               % Linear interpolation of x,y position to keep number of time points consistent across trials
        catch
            tsinterp = interp1(traj(trix).frameTimes(1:end-1)  - 0.5 - obj.bp.ev.(alignEv)(trix), ts, taxis); % for behav only data objects, frame times aren't exact, so someimtes have an extra index in there
        end
        %         tsinterpsmooth = interp1(traj(trix).frameTimes-0.5-obj.bp.ev.(alignEv)(trix), tssmooth, taxis);               % Linear interpolation of x,y position to keep number of time points consistent across trials
        %
        %         figure(1); hold on
        %         plot(taxis,tsinterp,'LineWidth',2);
        %         plot(taxis,tsinterpsmooth,'LineWidth',2)
        %         xlim([-2.5 2.5])
        %         pause
        %         clf

        xpos(:,i) = tsinterp(:,1);
        ypos(:,i) = tsinterp(:,2);



        %         fill missing values for all features except tongue
        if contains(feat,'tongue')
            continue
%             xpos(:,i) = fillmissing(xpos(:,i),'constant',nanmin(xpos(:,i)));
%             ypos(:,i) = fillmissing(ypos(:,i),'constant',nanmin(xpos(:,i)));
        else
            xpos(:,i) = fillmissing(xpos(:,i),'nearest');
            ypos(:,i) = fillmissing(ypos(:,i),'nearest');
        end

    end
end

end  % findPosition
