clear,clc,close all

% add paths for data loading scripts, all fig funcs, and utils
utilspth = 'C:\Users\munib\Documents\Economo-Lab\code\uninstructedMovements_v2';
addpath(genpath(fullfile(utilspth,'DataLoadingScripts')));
addpath(genpath(fullfile(utilspth,'funcs')));
addpath(genpath(fullfile(utilspth,'utils')));


% add paths for figure specific functions
addpath(genpath(pwd))

%% PARAMETERS
params.alignEvent          = 'goCue'; % 'jawOnset' 'goCue'  'moveOnset'  'firstLick'  'lastLick'

% time warping only operates on neural data for now.
% TODO: time warp for video and bpod data
params.timeWarp            = 0;  % piecewise linear time warping - each lick duration on each trial gets warped to median lick duration for that lick across trials
params.nLicks              = 20; % number of post go cue licks to calculate median lick duration for and warp individual trials to

params.lowFR               = 1; % remove clusters with firing rates across all trials less than this val

% set conditions to calculate PSTHs for
params.condition(1) = {'L&~stim.enable&~autowater'};             % right hits, no stim, aw off
params.condition(end+1) = {'R&~stim.enable&~autowater'};             % left hits, no stim, aw off

params.tmin = -2.5;
params.tmax = 2.5;
params.dt = 1/100;

% smooth with causal gaussian kernel
params.smooth = 15;

% cluster qualities to use
params.quality = {'all'}; % accepts any cell array of strings - special character 'all' returns clusters of any quality


params.traj_features = {{'tongue','left_tongue','right_tongue','jaw','trident','nose'},...
    {'top_tongue','topleft_tongue','bottom_tongue','bottomleft_tongue','jaw','top_paw','bottom_paw','top_nostril','bottom_nostril'}};

params.feat_varToExplain = 80; % num factors for dim reduction of video features should explain this much variance


params.advance_movement = 0.0;


%% SPECIFY DATA TO LOAD

datapth = '/Users/Munib/Documents/Economo-Lab/data/';

meta = [];

% meta = loadJEB6_ALMVideo(meta,datapth);
% meta = loadJEB7_ALMVideo(meta,datapth);
% meta = loadEKH1_ALMVideo(meta,datapth);
% meta = loadEKH3_ALMVideo(meta,datapth);
meta = loadJGR2_ALMVideo(meta,datapth);
% meta = loadJGR3_ALMVideo(meta,datapth);
% meta = loadJEB15_ALMVideo(meta,datapth);


params.probe = {meta.probe}; % put probe numbers into params, one entry for element in meta, just so i don't have to change code i've already written


%% LOAD DATA

meta = meta(1);
params.probe = params.probe(1);
[obj,params] = loadSessionData(meta,params);


%%
close all

clrs = getColors();

cond = [2,3];
lw = 1.5;
ms = 10;

goCue = 0;
sample = mode(obj.bp.ev.sample) - mode(obj.bp.ev.goCue);
delay = mode(obj.bp.ev.delay) - mode(obj.bp.ev.goCue);

nTrials = min(numel(params.trialid{1}),numel(params.trialid{2}));

% left trials
f = figure; hold on;
f.Position = [520   480   314   499];
trialOffset = 1;
for trix = 1:nTrials
    check1 = 0;
    check2 = 0;

    trial = params.trialid{1}(trix);
    lickL =  obj.bp.ev.lickL{trial} - obj.bp.ev.goCue(trial);
    lickL(lickL > 2) = [];
    lickR =  obj.bp.ev.lickR{trial} - obj.bp.ev.goCue(trial);
    lickR(lickR > 2) = [];

    plot([sample sample], trialOffset+[-0.5 0.5], 'k:', 'LineWidth', lw);
    plot([delay delay], trialOffset+[-0.5 0.5], 'k:', 'LineWidth', lw);
    plot([goCue goCue], trialOffset+[-0.5 0.5], 'k:', 'LineWidth', lw);

    if ~isempty(lickL)
        plot(lickL, trialOffset*ones(size(lickL)), '.', 'Color', clrs.lhit, 'MarkerSize',ms);
        check1 = 1;
    end

    if ~isempty(lickR)
        plot(lickR, trialOffset*ones(size(lickR)), '.', 'Color', clrs.rhit, 'MarkerSize',ms);
        check2 = 1;
    end

    if obj.bp.hit(trial)
        fill([2.4 2.65 2.65 2.4], [trialOffset-0.5 trialOffset-0.5 trialOffset+0.5 trialOffset+0.5], [150 150 150]./255,'EdgeColor','none')
    elseif obj.bp.miss(trial)
        fill([2.4 2.65 2.65 2.4], [trialOffset-0.5 trialOffset-0.5 trialOffset+0.5 trialOffset+0.5], [0 0 0]./255,'EdgeColor','none')
    elseif obj.bp.no(trial)
%         fill([2.4 2.65 2.65 2.4], [trialOffset trialOffset trialOffset+1 trialOffset+1], [1 1 1]./255,'EdgeColor','k')
    end

    if check1 || check2
        trialOffset = trialOffset + 1;
    end



end
xlim([-2.5 2.7]);
xlabel('Time (s) from go cue')
ylabel('Trials')
ax = gca;
ax.FontSize = 12;


% right trials
f = figure; hold on;
f.Position = [520   480   314   499];
trialOffset = 1;
for trix = 1:nTrials
    check1 = 0;
    check2 = 0;

    trial = params.trialid{2}(trix);
    lickL =  obj.bp.ev.lickL{trial} - obj.bp.ev.goCue(trial);
    lickL(lickL > 2) = [];
    lickR =  obj.bp.ev.lickR{trial} - obj.bp.ev.goCue(trial);
    lickR(lickR > 2) = [];

    plot([sample sample], trialOffset+[-0.5 0.5], 'k:', 'LineWidth', lw);
    plot([delay delay], trialOffset+[-0.5 0.5], 'k:', 'LineWidth', lw);
    plot([goCue goCue], trialOffset+[-0.5 0.5], 'k:', 'LineWidth', lw);

    if ~isempty(lickL)
        plot(lickL, trialOffset*ones(size(lickL)), '.', 'Color', clrs.lhit, 'MarkerSize',ms);
        check1 = 1;
    end

    if ~isempty(lickR)
        plot(lickR, trialOffset*ones(size(lickR)), '.', 'Color', clrs.rhit, 'MarkerSize',ms);
        check2 = 1;
    end

    if obj.bp.hit(trial)
        fill([2.4 2.65 2.65 2.4], [trialOffset-0.5 trialOffset-0.5 trialOffset+0.5 trialOffset+0.5], [150 150 150]./255,'EdgeColor','none')
    elseif obj.bp.miss(trial)
        fill([2.4 2.65 2.65 2.4], [trialOffset-0.5 trialOffset-0.5 trialOffset+0.5 trialOffset+0.5], [0 0 0]./255,'EdgeColor','none')
    elseif obj.bp.no(trial)
%         fill([2.4 2.65 2.65 2.4], [trialOffset trialOffset trialOffset+1 trialOffset+1], [1 1 1]./255,'EdgeColor','k')
    end

    if check1 || check2
        trialOffset = trialOffset + 1;
    end



end
xlim([-2.5 2.7]);
xlabel('Time (s) from go cue')
ylabel('Trials')
ax = gca;
ax.FontSize = 12;












