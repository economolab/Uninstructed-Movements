% Script for quantifying jaw movements on hazarded delay 2AFC
%%
clear; clc; close all;

addpath(genpath('C:\Users\Jackie\Documents\Grad School\Economo Lab\Code\ActivityModes\funcs'));
addpath(genpath('C:\Users\Jackie\Documents\Grad School\Economo Lab\Code\Data-Loading-Scripts'));
addpath(genpath('C:\Users\Jackie\Documents\Grad School\Economo Lab\Code\Uninstructed-Movements\functions'));
addpath(genpath('C:\Users\Jackie\Documents\Grad School\Economo Lab\Code\Uninstructed-Movements\utils'));
addpath(genpath('C:\Users\Jackie\Documents\Grad School\Economo Lab\Code\Utils'));

% Saving params
outputdir = 'C:\Users\Jackie\Documents\Grad School\Economo Lab\Figures\Uninstructed Movements';
toSave = 'no';
%% SET RUN PARAMS

% Which method you want to use to identify early movement trials:
% 'motionEnergy' or 'DeepLabCut'
params.alignEvent          = 'delay';   % goCue or firstLick
params.dt = 0.05;

% set conditions to use for projections
params.condition(1) = {'R&hit&~stim.enable&~early'};  % All hits, no stim, aw off
params.condition(end+1) = {'L&hit&~stim.enable&~early'};  % All hits, no stim, aw off

% Delay period length that you want to warp all delay lengths to
params.desiredDelay = 0.9000;       

% Different delay period lengths that were used in hazarded delay
params.delay(1) = 0.3000;
params.delay(2) = 0.6000;
params.delay(3) = 1.2000;
params.delay(4) = 1.8000;
params.delay(5) = 2.4000;
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
    disp(['Loading Hazarded Delay behavior obj for session ' num2str(i) ' out of ' num2str(numel(meta))])
    meta(i).trialid = findTrials(obj, params.condition);   % Get which trials pertain to the behavioral conditions you are looking at
end
%% LOAD DATA

% ----------------------------------------------
% -- Neural Data --
% obj (struct array) - one entry per session
% params (struct array) - one entry per session
% ----------------------------------------------
% [obj,params] = loadSessionData(meta,params);

%% For each session--Trials separated by delay length
plotIndiv = 'no';
% Get the trialIDs corresponding to each delay length
% Find the PSTH for R and L trials of each delay length
% Find the jaw velocities for R and L trials of each delay length
for gg = 1:length(meta)
    sesh = gg;
    met = meta(sesh);

    del(sesh).delaylen = objs{sesh}.bp.ev.goCue - objs{sesh}.bp.ev.delay;       % Find the delay length for all trials
    conditions = {1,2};
    del(sesh).del_trialid = getDelayTrix(met,conditions,params,del(sesh));     % Group the trials in each condition based on their delay length
    %del(sesh).delPSTH = getPSTHbyDel(params,del(sesh),objs{sesh});             % Get avg PSTH for each delay length

    % Find the probability of jaw [Jaw] movement at all time points in the session for trials of
    % specific conditions
    del(sesh).taxis = taxis;
    jaw_by_cond = findJawVelocity(del(sesh).taxis, objs{sesh},conditions,'vel',params,met);    % (1 x conditions cell array)
    % Each cell: (time x trials in that condition)

    % Find average jaw velocity for each delay length
    jawvel.left = cell(1,length(params.delay));         % (1 x number of delay lengths)
    jawvel.right = cell(1,length(params.delay));
    for g = 1:length(params.delay)                  % For each delay length...
        gix = find(del(sesh).del_trialid{1}==g);              % Get the trial IDs in the first condition that have the current delay length
        tempjaw = nanmean(jaw_by_cond{1}(:,gix),2);     % Find avg jaw velocity for first condition trials with that delay
        jawvel.right{g} = medfilt1(tempjaw,10);         % Apply median filter

        gix = find(del(sesh).del_trialid{2}==g);              % Same thing for second condition
        tempjaw = nanmean(jaw_by_cond{2}(:,gix),2);
        jawvel.left{g} = medfilt1(tempjaw,10);
    end
    del(sesh).jawvel = jawvel;
end
%% Re-organize data
% Will end up with:
% jawVel = (1 x nDelLengths).  Each DelLength will have a (time x conditions) array of jaw velocities
for sesh = 1:length(meta)
    ev.sample = objs{sesh}.bp.ev.sample;
    ev.delay = objs{sesh}.bp.ev.delay;
    ev.goCue = objs{sesh}.bp.ev.goCue;
    ev.(params.alignEvent) = objs{sesh}.bp.ev.(params.alignEvent);

    del(sesh).alignEvent = params.alignEvent;

    tempj = cell(1,length(params.delay));
    tempsel = cell(1,length(params.delay));
    for d = 1:length(params.delay)
        tempj{d} = cat(2,del(sesh).jawvel.right{d},del(sesh).jawvel.left{d});
        tempsel{d} = tempj{d}(:,1)-tempj{d}(:,2);
    end
    del(sesh).jawvel = tempj; del(sesh).selectivity = tempsel;
end
%% Avg selectivity in jaw vel across sessions
clear temp    
del2use = 3;
temp.sel = NaN(length(del(1).taxis),length(meta));   % time x sessions
for sesh = 1:length(meta)
    temp.sel(:,sesh) = mySmooth(del(sesh).selectivity{del2use},51);
end
Avg.selectivity = mean(temp.sel,2,'omitnan');  Std.selectivity = std(temp.sel,0,2,'omitnan');

clear temp
%% Plot the selectivity in jaw vel across sessions
figure();
haz_taxis = del(1).taxis + 0.5;
col = [0 0 0];
nSessions = length(meta);
alpha = 0.3;

toplot = Avg.selectivity;
err = (1.96*(Std.selectivity/nSessions));
ax = gca;
shadedErrorBar(taxis, toplot, err ,{'Color',col,'LineWidth',2}, alpha, ax)

title('Difference in jaw velocity')
xlim([-1.5 1.25])
xline(0,'LineStyle','-.','LineWidth',1.15,'Color',[0 0 0])
xline(-1.3,'LineStyle','-.','LineWidth',1,'Color',[0 0 0])
xlabel(['Time from ' params(1).alignEvent ' (s)'])

