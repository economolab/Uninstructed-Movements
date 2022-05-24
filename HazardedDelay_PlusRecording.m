% Script for looking at Hazarded Delay Recording + Video
clear; clc; close all;
%%
addpath(genpath('C:\Users\Jackie\Documents\Grad School\Economo Lab\Code\ActivityModes'));
addpath(genpath('C:\Users\Jackie\Documents\Grad School\Economo Lab\Code\Data-Loading-Scripts'));
addpath(genpath('C:\Users\Jackie\Documents\Grad School\Economo Lab\Code\Uninstructed-Movements'));
addpath(genpath('C:\Users\Jackie\Documents\Grad School\Economo Lab\Code\Utils'));

% addpath(genpath('C:\Code\ActivityModes'));
% addpath(genpath('C:\Code\Uninstructed Movements\Uninstructed-Movements\DataLoadingScripts'));
% addpath(genpath('C:\Code\Uninstructed-Movements'));
% addpath(genpath('C:\Code\Utils'));
% addpath(genpath('C:\Code\DataLoadingScripts'));


% Saving params
outputdir = 'C:\Users\Jackie\Documents\Grad School\Economo Lab\Figures\Uninstructed Movements';
toSave = 'no';
%% SET RUN PARAMS

params.alignEvent          = 'delay';   % delay,goCue, or firstLick
params.lowFR               = 1; % remove clusters firing less than this val
params.dt = 0.05;
params.jawMeasure          = 'sideJaw'; % sideJaw or Trident

% set conditions to use for projections
params.condition(1) = {'R&hit&~stim.enable&autowater.nums==2&~early'}; % right hits, no stim, aw off, no early response
params.condition(2) = {'L&hit&~stim.enable&autowater.nums==2&~early'}; % left hits, no stim, aw of, no early response


% set conditions used for finding the modes
aw = '2'; % 1-on, 2-off
stim = '0'; % 0-off
params.modecondition(1) = {['R&hit&autowater.nums==' aw '&stim.num==' stim '&~early']};     % R hits, 2afc, stim on/off, not early
params.modecondition(2) = {['L&hit&autowater.nums==' aw '&stim.num==' stim '&~early']};     % L hits, 2afc, stim on/off, not early
params.modecondition(3) = {['R&miss&autowater.nums==' aw '&stim.num==' stim '&~early']};    % R miss, 2afc, stim on/off, not early
params.modecondition(4) = {['L&miss&autowater.nums==' aw '&stim.num==' stim '&~early']};    % L miss, 2afc, stim on/off, not early
params.modecondition(5) = {['hit&autowater.nums==' aw '&stim.num==' stim '&~early']};       % All hits, 2afc, stim on/off, not early
params.modecondition(6) = {['hit&autowater.nums==1&stim.num==' stim '&~early']};        % All hits, aw on, stim on/off, not early


params.delay(1) = 0.3000;
params.delay(2) = 0.6000;
params.delay(3) = 1.2000;
params.delay(4) = 1.8000;
params.delay(5) = 2.4000;
%% SET METADATA FROM ALL RELEVANT SESSIONS/ANIMALS
meta = [];
meta = loadJEB11_ALMVideo(meta);
meta = loadJEB12_ALMVideo(meta);
meta(end).tmin = -2.5; % (s) relative to params.alignEvent
meta(end).tmax = 3;  % (s) relative to params.alignEvent
meta(end).dt = 0.005;
taxis = meta(end).tmin:meta(end).dt:meta(end).tmax;   % get time-axis with 0 as time of event you aligned to
taxis = taxis(1:end-1);
%% PREPROCESS DATA
objs = loadObjs(meta);

for i = 1:numel(meta)
    obj = objs{i};
    obj.condition = params.condition;
    % get trials and clusters to use
    meta(i).trialid = findTrials(obj, obj.condition);   % Get which trials pertain to the behavioral conditions you are looking at
    cluQuality = {obj.clu{meta(i).probe}(:).quality}';  % Get clusters that are of the qualities that you specified
    meta(i).cluid = findClusters(cluQuality, meta(i).quality);
    % align data
    obj = alignSpikes(obj,meta(i),params);              % Align the spike times to the event that you specified
    % get trial avg psth, single trial data, and single trial data grouped
    % by condition (aka R 2AFC, R AW, etc.)
    obj = getPSTHs(obj,meta(i));
    objs{i} = obj;
end
%% Remove unwanted sessions

% remove sessions with less than 40 trials of rhit and lhit each (same as
% hidehiko ppn paper)
use = false(size(objs));
for i = 1:numel(use)
    met = meta(i);
    check1 = numel(met.trialid{1}) > 40;
    check2 = numel(met.trialid{2}) > 40;
    if check1 && check2
        use(i) = true;
    end
end

meta = meta(use);
objs = objs(use);
%% Sanity check--plot prob of jaw movement averaged across trials with the same delay length
obj = objs{2};
met = meta(2);

delaylen = obj.bp.ev.goCue - obj.bp.ev.delay;       % Find the delay length for all trials
conditions = {1,2};
% 
met = getDelayTrialID(met,conditions,delaylen);     % Group the trials in each condition based on their delay length
numdelay = NaN(1,length(params.delay));
for g=1:length(params.delay)
    numdelay(g) = length(find(met.del_trialid{1}==g))+length(find(met.del_trialid{2}==g));
end
figure();
bar(numdelay)
ylabel('Num trials')
xlabel('Delay length (s)')
xticklabels(params.delay)
title('Distribution of delay lengths for example session')

% % Find the probability of jaw [Jaw] movement at all time points in the session for trials of
% % specific conditions
% jaw_by_cond = findJawVelocity(taxis, obj,conditions,met,'prob');    % (1 x conditions cell array)
% % Each cell: (time x trials in that condition)
% 
% % Find average jaw velocity for each delay length
% jawvel.left = cell(1,length(params.delay));         % (1 x number of delay lengths)
% jawvel.right = cell(1,length(params.delay));
% for g = 1:length(params.delay)                  % For each delay length...
%     gix = find(met.del_trialid{1}==g);              % Get the trial IDs in the first condition that have the current delay length
%     tempjaw = nanmean(jaw_by_cond{1}(:,gix),2);     % Find avg jaw velocity for first condition trials with that delay
%     jawvel.right{g} = medfilt1(tempjaw,10);         % Apply median filter
% 
%     gix = find(met.del_trialid{2}==g);              % Same thing for second condition
%     tempjaw = nanmean(jaw_by_cond{2}(:,gix),2);
%     jawvel.left{g} = medfilt1(tempjaw,10);
% end
% 
% % Plot probability of jaw[Jaw] movement for each delay length
% figure();
% colors = {[0 0 1],[1 0 0]};
% plotJawProb_HazardDel(taxis, jawvel,params)
%% Load all of the data

for gg = 1:length(meta)
    figure();
    sesh = gg;
    obj = objs{sesh};     
    met = meta(sesh);

    anm = obj.pth.anm;                  % Animal name
    date = obj.pth.dt;                  % Session date
    probenum = string(met.probe);       % Which probe was used

    %%%% FIND JAW VEL %%%%
    % Find the prob of jaw movement for all points in the trial; aligned to
    % delay
    conditions = {1,2};
    colors = {[0 0 1],[1 0 0]};
    if strcmp(params.jawMeasure,'sideJaw')
        plotJawProb_SessAvg(obj,met,conditions,colors,taxis,'no',params)
    elseif strcmp(params.jawMeasure,'Trident')
        plotTridentVel_SessAvg(obj,met,conditions,colors,'no')
    end
    xlim([-1.4 0.3])
    xline(0,'LineStyle','--')
    xline(-1.3,'LineStyle','--')
    legend('Right','Left','Location','best')
    sesstitle = strcat('Prob of jaw movement for',{' '},anm,{' '},date);  % Name/title for session
    title(sesstitle,'FontSize',13)
    
    %%%% JAW MOVEMENT SELECTIVITY %%%%
    
    edges = met.tmin:met.dt:met.tmax;
    [jawprob,~] = jawProbSessionAvg(obj,met,conditions,edges,params);
    jawselectivity = jawprob{1}-jawprob{2};
    figure();
    plot(taxis,jawselectivity(2:end),'LineWidth',2,'Color','black')
    xlim([-1.4 0.3])
    xline(0,'LineStyle','--')
    xline(-1.3,'LineStyle','--')
    ylabel('Selectivity (probability of jaw movement)')
    xlabel('Time since delay onset (s)')
    sesstitle = strcat('Selectivity in jaw movements for',{' '},anm,{' '},date);  % Name/title for session
    title(sesstitle,'FontSize',13)

    %%%% FIND CHOICE MODE %%%%
    rez.time = objs{1}.time;
    rez.condition = objs{1}.condition;
    rez.alignEvent = params.alignEvent;

    % Find CDchoice (coding dimension during delay period)
    cond{1} = params.modecondition{1};
    cond{2} = params.modecondition{2};
    epoch = 'fulldelay';
    choice_mode = choiceMode(obj,met,cond,epoch,rez.alignEvent,'no');


    % Get the projection of specified conditions onto the choice mode
    smooth = 51;
    conditions = [1,2];
    latentChoice = getChoiceModeProjection(obj,choice_mode,smooth,conditions);

    % Plot
    figure();
    colors = {[0 0 1],[1 0 0]};
    for c = 1:numel(conditions)
        temp = conditions(c);
        col = colors{c};
        plot(taxis,latentChoice{temp},'Color',col,'LineWidth',3)
        hold on;
    end
    col = 'black';
    xline(0,'LineStyle','--')
    xline(-1.3,'LineStyle','--')
    xlim([-1.4 0.3])
    xlabel('Time since delay onset')
    ylabel('a.u.')
    %set(gca, 'YDir','reverse')
    legend('Right','Left','Location','best')
    sesstitle = strcat('Projection onto CDchoice for',{' '},anm,{' '},date);  % Name/title for session
    title(sesstitle,'FontSize',13)

    %%%% FIND NEURAL SELECTIVITY %%%%
end

