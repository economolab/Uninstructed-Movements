
clear,clc,close all

addpath(genpath(pwd))

%%


meta.pth = '/Users/Munib/Documents/Economo-Lab/data/poleTask';
meta.anm = 'anm395710'; % anm395710  anm395709
meta.date = '180111';   % 180111     180111

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