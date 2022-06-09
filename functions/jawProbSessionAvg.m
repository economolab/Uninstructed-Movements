% Function for finding the average probability of jaw movement at all time
% points in a trial

% INPUT to 'jawkinAnimalAvg' function: obj, met for current session
% conditions: the trial conditions that you want to look at

% OUTPUT: jawprob: cell array (1 x num conditions)
% Each cell will contain the probability of jaw movement for one specified
% condition

function [jawprob,jawstd] = jawProbSessionAvg(obj,met,conditions,edges,params)
traj = obj.traj{1};                             % Get the video data
if strcmp(traj(1).featNames{2},'jaw')
    featnum = 2;
elseif strcmp(traj(1).featNames{4},'jaw')
    featnum = 4;
end
jawprob = cell(1,numel(conditions));
jawstd = cell(1,numel(conditions));
for c = 1:numel(conditions)
    cond = conditions{c};
    derivthresh = 0.5;  %If the velocity of the jaw crosses this threshold, the jaw is considered to be moving
    nTrials = length(met.trialid{cond});
    temptri = nan(numel(edges), nTrials);    % (time x num trials in curr condition)
    for i = 1:nTrials                        % For every trial in the condition
        trix = met.trialid{cond}(i);
        if isfield(traj,'NdroppedFrames')
            if isnan(traj(trix).NdroppedFrames )                       % If the video data from this trial isn't good, skip it
                continue;
            end

            if ~isnan(traj(trix).frameTimes)                               % If the video data from this trial is good...
                ts = mySmooth(traj(trix).ts(:, 2, featnum), 21);                                               % Side-view, up and down position of the jaw, smoothed
                tsinterp = interp1(traj(trix).frameTimes-0.5-obj.bp.ev.(params.alignEvent)(trix), ts, edges);          % Linear interpolation of jaw position to keep number of time points consistent across trials
                basederiv = median(diff(tsinterp),'omitnan');                                         % Find the median jaw velocity (aka baseline)
            end

        else
            ts = mySmooth(traj(trix).ts(:, 2, featnum), 21);                                               % Side-view, up and down position of the jaw, smoothed
            nFrames = length(ts);
            frameTimes = (1:nFrames)./400;
            tsinterp = interp1(frameTimes-obj.bp.ev.(params.alignEvent)(trix), ts, edges);                           % Linear interpolation of jaw position to keep number of time points consistent across trials
            basederiv = median(diff(tsinterp),'omitnan');                                            % Find the median jaw velocity (aka baseline)
        end

        %Find the difference between the jaw velocity and the
        %baseline jaw velocity
        temptri(2:end, i) = abs(diff(tsinterp)-basederiv)>derivthresh;      % Values > 0 = jaw is moving

    end

    jawstd{c} = std(temptri,0,2,'omitnan');
    jawdat = mean(temptri,2,'omitnan');         % Take the average probability of jaw movement across trials
    jawdat = medfilt1(jawdat,10);           % Smooth the avg probability of jaw movement

    jawprob{c} = jawdat;                % Store the avg prob of jaw movement for the current condition
end
%end % jawProbSessionAvg