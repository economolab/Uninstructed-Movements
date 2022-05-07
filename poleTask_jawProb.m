
clear,clc,close all

addpath(genpath(pwd))

%%


meta.pth = '/Users/Munib/Documents/Economo-Lab/data/poleTask';
meta.anm = 'anm395709'; % anm395710  anm395709   anm395709
meta.date = '180111';   % 180111     180111      180111

[obj,traj] = loadPoleTaskData(meta); % traj is just the side cam

%% 

% get left and right hits, not early, not stim, trials

% 
% hit = obj.yes_no_multi_pole_sensoryInputobj_hit_history;
% R = obj.yes_no_multi_pole_sensoryInputobj_correct_R_history;
% L = obj.yes_no_multi_pole_sensoryInputobj_correct_L_history;
% stim = obj.yes_no_multi_pole_sensoryInputobj_stim_delay_in_ms_history(1:end-1);
% 
% trials.R = hit&R&~stim;
% trials.L = hit&L&~stim;
% 
% trials.R = find(trials.R(1:numel(traj)));
% trials.L = find(trials.L(1:numel(traj)));


% right trial is pole is posterior and jaw is not visible
% left trial otherwise (this is probably not right since there may be stim, nr, and error trials in here)
R = false(numel(traj),1);
L = false(numel(traj),1);
for i = 1:numel(traj)
    if any(isnan(traj(i).ts(:,2,4)))
        R(i) = true;
    else
        L(i) = true;
    end
end
trials.R = find(R);
trials.L = find(L);


%%

thresh = 0.2;

smth = 21;



jawix = 4;
jawprob = jawProbSessionAvg_poleTask(traj,trials,jawix,thresh);

time = (1:numel(jawprob{1})) ./ 400;
cut = 15;
time = time(cut:end);

jawR = mySmooth(jawprob{1}(cut:end),smth);
jawL = mySmooth(jawprob{2}(cut:end),smth);
figure; plot(time,jawR,'b'); hold on; plot(time,jawL,'r')
title([meta.anm ' | ' meta.date]);
xlabel('Time')
ylabel('Prob Jaw Move')
ax = gca;
ax.FontSize = 20;

tridentix = 5;
tridentprob = jawProbSessionAvg_poleTask(traj,trials,tridentix,thresh);
tridentR = mySmooth(tridentprob{1}(cut:end),smth);
tridentL = mySmooth(tridentprob{2}(cut:end),smth);
figure; plot(time,tridentR,'b'); hold on; plot(time,tridentL,'r')
title([meta.anm ' | ' meta.date]);
ylabel('Prob Trident Move')
xlabel('Time')
ax = gca;
ax.FontSize = 20;


%%

figure; clf; hold on
for i = 1:numel(traj)
    plot(traj(i).ts(:,2,4) + 10*(i-1))
end
hold off

figure; clf; hold on
for i = 1:numel(trials.R)
    plot(traj(trials.R(i)).ts(:,2,4) + 10*(i-1))
end
hold off

figure; clf; hold on
for i = 1:numel(trials.L)
    plot(traj(trials.L(i)).ts(:,2,4) + 10*(i-1))
end
hold off

%% Helper Functions

function [obj,traj] = loadPoleTaskData(meta)

contents = dir(meta.pth);
fns = {contents.name};

toMatchObj = {'data',meta.anm,meta.date};
toMatchTraj = {'Traj',meta.anm,meta.date};

[fn,~] = patternMatchCellArray(fns',toMatchObj','all');

dat = load(fullfile(meta.pth,fn{1}));
obj = dat.saved;
clear dat

[fn,~] = patternMatchCellArray(fns',toMatchTraj','all');

dat = load(fullfile(meta.pth,fn{1}));
traj = dat.traj;
clear dat

end


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
end % jawProbSessionAvg