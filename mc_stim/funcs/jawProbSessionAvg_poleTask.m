% Function for finding the average probability of jaw movement at all time
% points in a trial

% INPUT to 'jawkinAnimalAvg' function: obj, met for current session
% conditions: the trial conditions that you want to look at 

% OUTPUT: jawprob: cell array (1 x num conditions)
% Each cell will contain the probability of jaw movement for one specified
% condition

function jawprob = jawProbSessionAvg_poleTask(traj,trials,featix,thresh)

derivthresh = thresh;


conditions = fieldnames(trials);
jawprob = cell(1,numel(conditions));
for i = 1:numel(conditions)
    cond = conditions{i};
    trix = trials.(cond);

%     derivthresh = 0.3;  %If the velocity of the jaw crosses this threshold, the jaw is considered to be moving

    % find number of time points across all trials
    alltrials = [trials.R ; trials.L];
    for tt = 1:numel(alltrials)
        n(tt) = size(traj(alltrials(tt)).ts,1);
    end
    skip = [];
    for tt = 1:numel(n)
        if n(tt) ~= mode(n)
            skip = [skip ; alltrials(tt)];
        end
    end
    
    
    jaw = nan(n(10), length(trix));
%     jawvel = jaw;

    for ii = 1:length(trix)

        q = trix(ii);
        
        if ismember(q,skip)
            continue
        end
        
        ts = mySmooth(traj(q).ts(:, 2, featix), 21);
        basederiv = median(diff(ts),"omitnan");                   %Find the median jaw velocity (aka baseline)

        %Find when the difference between the jaw velocity and the
        %baseline jaw
        %velocity is above a given threshold (when is jaw moving?)
%         jawvel(2:end,ii) = abs(diff(ts)-basederiv);
        jaw(2:end,ii) = abs(diff(ts)-basederiv)>derivthresh;% | abs(tsinterp(2:end)-basepos)>posthresh;

    end
    
%     nojawmovix = ~ismember(trialsToUse, find(obj.earlyMoveix));
%     jawdat_noearly = mean(jaw(:, nojawmovix), 2, 'omitnan');
%     jawdat_noearly = medfilt1(jawdat_noearly, 10);
    
    
    jawdat = mean(jaw,2,'omitnan');         % Take the average probability of jaw movement across trials
    jawdat = medfilt1(jawdat,10);           % Smooth the avg probability of jaw movement 
    jawprob{i} = jawdat;                % Store the avg prob of jaw movement for the current condition
    
%     figure; imagesc(jaw)
%     figure; imagesc(jaw(:, nojawmovix))
end
%end % jawProbSessionAvg