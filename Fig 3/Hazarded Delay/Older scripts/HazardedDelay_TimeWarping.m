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

% Different delay period lengths that were used in hazarded delay
params.delay(1) = 0.3000;
params.delay(2) = 0.6000;
params.delay(3) = 1.2000;
params.delay(4) = 1.8000;
% params.delay(5) = 2.4000;

% Parameters for warping the delay period to 0.9 s
warp.desiredDelay = 1.4500;                     % Desired start time for delay period within the trial
warp.desiredGoCue = warp.desiredDelay + 0.9;    % Desired end time for delay period within the trial
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

for gg = 1:numel(meta)
    obj = objs{gg};
    met = meta(gg);

    anm = met.anm;
    date = met.date;
    sesstitle = strcat(anm,date);  % Name/title for session

    delaylen = obj.bp.ev.goCue - obj.bp.ev.delay;       % Find the delay length for all trials
    conditions = {1,2};

    met = getDelayTrialID(met,conditions,delaylen);     % Group the trials in each condition based on their delay length

    [obj] = warpDelayPeriod(obj,met,conditions,warp,'behaviorOnly');
    
    % Find the probability of jaw movement at all time points in the session for trials of
    % specific conditions
    jaw_by_cond = findWarpedJawVelocity(taxis, obj,conditions,met,'prob');    % (1 x conditions cell array)
                                                                        % Each cell: (time x trials in that condition)
    % Get the average prob of jaw movement across warped trials 
    meanjaw = cell(1,numel(conditions));
    for c = 1:numel(conditions)
        temp = jaw_by_cond{c};
        tempjaw = mean(temp,2,'omitnan');
        meanjaw{c} = medfilt1(tempjaw,10);
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

