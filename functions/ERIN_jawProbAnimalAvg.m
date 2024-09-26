% Function for finding the average probability of jaw movement at all time
% points in a trial

% INPUT to 'jawkinAnimalAvg' function: obj, met for current session
% conditions: the trial conditions that you want to look at

% OUTPUT: jawprob: cell array (1 x num conditions)
% Each cell will contain the probability of jaw movement for one specified
% condition

function alljaw = ERIN_jawProbAnimalAvg(objs,meta,conditions)
jawdat = cell(numel(objs),numel(conditions));
jawvelocity = cell(numel(objs),numel(conditions));

for oo = 1:numel(objs)
    obj = objs{oo};
    met = meta(oo);

    jawprob = cell(1,numel(conditions));
    for i = 1:numel(conditions)
        cond = conditions{i};
        trialsToUse = met.trialid{cond};

        derivthresh = 0.15;  %If the velocity of the jaw crosses this threshold, the jaw is considered to be moving
        edges = 0:0.005:5.5;
        Ntrials = length(trialsToUse);

        traj = obj.traj{1};

        jaw = NaN(1401, Ntrials);
        jawvel = NaN(1401, Ntrials);

        for ii = 1:Ntrials
            q = ii;
            ts = mySmooth(traj(q).ts(:, 2, 2), 21);
            ts = ts(1:1401,:);
            %tsinterp = interp1(traj(q).frameTimes-0.5, ts, edges);   %Linear interpolation of jaw position to keep number of time points consistent across trials
            basederiv = median(diff(ts),"omitnan");                   %Find the median jaw velocity (aka baseline)

            %Find when the difference between the jaw velocity and the
            %baseline jaw
            %velocity is above a given threshold (when is jaw moving?)
            jaw(2:end, ii) = abs(diff(ts)-basederiv)>derivthresh;% | abs(tsinterp(2:end)-basepos)>posthresh;
            jawvel(2:end,ii) = abs(diff(ts)-basederiv);
        end     % end one condition for one session
        jawdat{oo,i} = jaw; 
        jawvelocity{oo,i} = jawvel; 

    end 
end  

alljaw = cell(1,numel(conditions));
for i = 1:numel(conditions)
    temp = [];
    for oo = 1:numel(objs)
        temp = [temp,jawvelocity{oo,i}];
    end
%     temp = mean(temp,2,'omitnan');
%     temp = medfilt1(temp,10);
    alljaw{i} = temp;
end

end  % end ERIN_jawProbAnimalAvg