% Script for quantifying jaw movements on hazarded delay 2AFC
%%
clear; clc; close all;

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

% Which method you want to use to identify early movement trials:
% 'motionEnergy' or 'DeepLabCut'
params.alignEvent          = 'goCue';   % goCue or firstLick
params.dt = 0.05;

% set conditions to use for projections
params.condition(1) = {'R&hit&~stim.enable&~early'}; % right hits, no stim, aw off
params.condition(2) = {'L&hit&~stim.enable&~early'}; % left hits, no stim, aw off

% Delay period length that you want to warp all delay lengths to
params.desiredDelay = 0.9000;       

% Different delay period lengths that were used in hazarded delay
params.delay(1) = 0.3000;
params.delay(2) = 0.6000;
params.delay(3) = 1.2000;
params.delay(4) = 1.8000;
% params.delay(5) = 2.4000;
%%  LOAD META DATA
meta = [];
meta = loadJEB11_BehavVid(meta);
meta = loadJEB12_BehavVid(meta);
meta(end).tmin = -2.5; % (s) relative to params.alignEvent
meta(end).tmax = 3;  % (s) relative to params.alignEvent
meta(end).dt = 0.005;

taxis = meta(end).tmin:meta(end).dt:meta(end).tmax;   % get time-axis with 0 as time of event you aligned to
taxis = taxis(1:end-1);
%% PREPROCESS DATA
objs = loadBehavVid(meta);
for i = 1:numel(meta)
    obj = objs{i};
    meta(i).trialid = findTrials(obj, params.condition);   % Get which trials pertain to the behavioral conditions you are looking at
end
%%  PLOT SUMMARY FIGURES FOR EACH SESSION
%   Figure 1 = prob of jaw movement for each delay duration (averaged
%   throughout the session)
for gg = 1:numel(meta)
    obj = objs{gg};
    met = meta(gg);

    anm = met.anm;
    date = met.date;
    sesstitle = strcat(anm,date);  % Name/title for session

    delaylen = obj.bp.ev.goCue - obj.bp.ev.delay;       % Find the delay length for all trials
    conditions = {1,2};

    met = getDelayTrialID(met,conditions,delaylen);     % Group the trials in each condition based on their delay length

    % Find the probability of jaw [Jaw] movement at all time points in the session for trials of
    % specific conditions
    jaw_by_cond = findJawVelocity(taxis, obj,conditions,met,'prob');    % (1 x conditions cell array)
    % Each cell: (time x trials in that condition)

    % Find average jaw velocity for each delay length
    jawvel.left = cell(1,length(params.delay));         % (1 x number of delay lengths)
    jawvel.right = cell(1,length(params.delay));
    for g = 1:length(params.delay)                  % For each delay length...
        gix = find(met.del_trialid{1}==g);              % Get the trial IDs in the first condition that have the current delay length
        tempjaw = nanmean(jaw_by_cond{1}(:,gix),2);     % Find avg jaw velocity for first condition trials with that delay
        jawvel.right{g} = medfilt1(tempjaw,10);         % Apply median filter

        gix = find(met.del_trialid{2}==g);              % Same thing for second condition
        tempjaw = nanmean(jaw_by_cond{2}(:,gix),2);
        jawvel.left{g} = medfilt1(tempjaw,10);
    end

    % Plot probability of jaw[Jaw] movement for each delay length
    figure();
    colors = {[0 0 1],[1 0 0]};
    %     plotJawProb_SessAvg(obj,met,conditions,colors)
    %     xline(0,'LineStyle','--')
    %     xline(-0.9,'LineStyle','-.')
    %     xline(-2.2,'LineStyle','-.')
    plotJawProb_HazardDel(taxis, jawvel,params)
    sgtitle(sesstitle)
end

%% PLOT SUMMARY FIGURES FOR EACH SESSION
%   Heatmaps of jaw velocity on each trial
for gg = 1:numel(meta)
    obj = objs{gg};
    met = meta(gg);

    anm = met.anm;
    date = met.date;
    sesstitle = strcat(anm,date);  % Name/title for session

    delaylen = obj.bp.ev.goCue - obj.bp.ev.delay;       % Find the delay length for all trials
    conditions = {1,2};

    met = getDelayTrialID(met,conditions,delaylen);     % Group the trials in each condition based on their delay length

    % Find jaw velocity for each condition
    vel_by_cond = findJawVelocity(taxis, obj,conditions,met,'vel');

    % Combining jaw vel information across conditions into a single struct
    val1 = obj.bp.ev.goCue(met.trialid{1});
    [go1, ix1] = sort(val1, 'descend');          % Sort right trials in descending order by delay length
    vel_by_cond{1} = vel_by_cond{1}(:,ix1);      % Sort jaw vel on right trials in that same order

    val2 = obj.bp.ev.goCue(met.trialid{2});     % Same for left trials
    [go2, ix2] = sort(val2, 'descend');
    vel_by_cond{2} = vel_by_cond{2}(:,ix2);

    % Plot heatmap of jaw velocity; separated by delay length
    plotHazardedDelHeatmap (met,vel_by_cond, delaylen,taxis,ix1,ix2)
    sgtitle(sesstitle)
end
%% Plot summary figures for each animal
% Probability of jaw movement averaged across all sessions for an animal,
% separated by delay length

% Find average for each animal (across all of it's sessions)
[jaw_allAnm,nAnimals,uc,sessbyAnm] = findHazardedJaw_Multi(objs,meta,conditions,taxis,params);

for a = 1:nAnimals
    clear jawvel
    figure();
    jawvel.right = jaw_allAnm.right(a,:);
    jawvel.left = jaw_allAnm.left(a,:);
    plotJawProb_HazardDel(taxis, jawvel,params)
    anm = uc{a};
    nSess = num2str(sessbyAnm{a});
    sesstitle = strcat('Average across  ',' ',nSess,' sessions for',' ',anm);  % Name/title for session
    sgtitle(sesstitle)
end
%% Plot summary figure for average across all animals

AvgAll.right = cell(1,length(params.delay));
AvgAll.left = cell(1,length(params.delay));
for p = 1:length(params.delay)
    temp.right = [];
    temp.left = [];
    for a = 1:nAnimals
        temp.right = [temp.right,jaw_allAnm.right{a,p}];
        temp.left = [temp.left,jaw_allAnm.left{a,p}];
    end
    AvgAll.right{p} = mean(temp.right,2,'omitnan');
    AvgAll.left{p} = mean(temp.left,2,'omitnan');
end

figure();
plotJawProb_HazardDel(taxis, AvgAll,params)
sgtitle('Avg across both animals')




