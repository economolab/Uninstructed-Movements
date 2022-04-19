% Function for finding the average probability of jaw movement at all time
% points in a trial

% INPUT to 'jawkinAnimalAvg' function: obj, met for current session
% conditions: the trial conditions that you want to look at 

% OUTPUT: jawprob: cell array (1 x num conditions)
% Each cell will contain the probability of jaw movement for one specified
% condition

function tridentprob = tridentProbSessionAvg(obj,met,conditions)

tridentprob = cell(1,numel(conditions));
for i = 1:numel(conditions)
    cond = conditions{i};
    trialsToUse = met.trialid{cond};

    derivthresh = 0.3;  %If the velocity of the jaw crosses this threshold, the jaw is considered to be moving
    %edges = 0:0.005:5.5;
    edges = met.tmin:met.dt:met.tmax;
    Ntrials = length(trialsToUse);
    traj = obj.traj{1};

    if isfield(traj,'frameTimes')
        good_ix = NaN(Ntrials, 1);
        %Find the trials that have numerical values for frameTimes (aka no NaNs)
        for ii = 1:Ntrials
            ix = trialsToUse(ii);
            if ~isnan(traj(ix).frameTimes)
                good_ix(ii,1) = ix;       %If the trial is okay, note the trial number
            end
        end

        %Get trial numbers of usable trials
        use_trials = find(~isnan(good_ix))';
        trialsToUse = trialsToUse(use_trials);
    end
    
    trident = NaN(numel(edges), length(trialsToUse));

    for ii = 1:length(trialsToUse)
        q = trialsToUse(ii);
        ts = mySmooth(traj(q).ts(:, 2, 5), 21);

        if isfield(traj,'frameTimes')
            tsinterp = interp1(traj(q).frameTimes-0.5, ts, edges);   %Linear interpolation of jaw position to keep number of time points consistent across trials
            basederiv = median(diff(tsinterp),"omitnan");                   %Find the median jaw velocity (aka baseline)
        else
            nFrames = length(ts);
            trialLen = nFrames*(1/400);
            frameTimes = 0:(1/400):trialLen;
            tsinterp = interp1(frameTimes(2:end)-obj.bp.ev.goCue(q), ts, edges);
            basederiv = median(diff(tsinterp),"omitnan");                   %Find the median jaw velocity (aka baseline)
        end
        %Find when the difference between the jaw velocity and the
        %baseline jaw
        %velocity is above a given threshold (when is jaw moving?)
        trident(2:end, ii) = abs(diff(tsinterp)-basederiv)>derivthresh;% | abs(tsinterp(2:end)-basepos)>posthresh;

    end    
    
    tridentdat = mean(trident,2,'omitnan');         % Take the average probability of jaw movement across trials
    tridentdat = medfilt1(tridentdat,10);           % Smooth the avg probability of jaw movement 
    tridentprob{i} = tridentdat;                % Store the avg prob of jaw movement for the current condition 
end
%end tridentProbSessionAvg