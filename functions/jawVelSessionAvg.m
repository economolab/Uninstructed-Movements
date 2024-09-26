% Function for finding the average probability of jaw movement at all time
% points in a trial

% INPUT to 'jawkinAnimalAvg' function: obj, met for current session
% conditions: the trial conditions that you want to look at 

% OUTPUT: jawprob: cell array (1 x num conditions)
% Each cell will contain the probability of jaw movement for one specified
% condition

function [jawprob,jawstd] = jawVelSessionAvg(obj,met,conditions)

jawprob = cell(1,numel(conditions));
jawstd = cell(1,numel(conditions));
for i = 1:numel(conditions)
    cond = conditions{i};
    trialsToUse = met.trialid{cond};

    derivthresh = 0.3;  %If the velocity of the jaw crosses this threshold, the jaw is considered to be moving
    edges = 0:0.005:5.5;
    Ntrials = length(trialsToUse);

    good_ix = NaN(Ntrials, 1);

    traj = obj.traj{1};

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

    jaw = NaN(numel(edges), length(trialsToUse));

    for ii = 1:length(trialsToUse)
        q = trialsToUse(ii);
        
        ts = mySmooth(traj(q).ts(:, 2, 2), 21);
        tsinterp = interp1(traj(q).frameTimes-0.5, ts, edges);          %Linear interpolation of jaw position to keep number of time points consistent across trials
        basederiv = median(diff(tsinterp),"omitnan");                   %Find the median jaw velocity (aka baseline)

        %Find when the difference between the jaw velocity and the
        %baseline jaw
        %velocity is above a given threshold (when is jaw moving?)
        jaw(2:end, ii) = abs(diff(tsinterp)-basederiv);%>derivthresh;% | abs(tsinterp(2:end)-basepos)>posthresh;

    end    
    
jawstd{i} = std(jaw,0,2,'omitnan');
jawdat = mean(jaw,2,'omitnan');         % Take the average probability of jaw movement across trials
jawdat = medfilt1(jawdat,10);           % Smooth the avg probability of jaw movement

jawprob{i} = jawdat;                % Store the avg prob of jaw movement for the current condition

%     figure; imagesc(jaw)
%     figure; imagesc(jaw(:, nojawmovix))
end
%end % jawProbSessionAvg