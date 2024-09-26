% DECODING CDContext FROM ALL KINEMATIC FEATURES
clear,clc,close all

% add paths for data loading scripts, all fig funcs, and utils
utilspth = 'C:\Users\Jackie\Documents\Grad School\Economo Lab\Code\Uninstructed-Movements\NullPotent_Context';
addpath(genpath(fullfile(utilspth,'DataLoadingScripts')));
addpath(genpath(fullfile(utilspth,'funcs')));
addpath(genpath(fullfile(utilspth,'utils')));
addpath(genpath(fullfile(utilspth,'fig1')));
addpath(genpath(fullfile(utilspth,'Context_funcs')));
otherpth = 'C:\Users\Jackie\Documents\Grad School\Economo Lab\Code\Uninstructed-Movements';
addpath(genpath(fullfile(otherpth,'Decoding Analysis')));

%% PARAMETERS
params.alignEvent          = 'firstLick'; % 'jawOnset' 'goCue'  'moveOnset'  'firstLick'  'lastLick'

% time warping only operates on neural data for now.
% TODO: time warp for video and bpod data
params.timeWarp            = 0;  % piecewise linear time warping - each lick duration on each trial gets warped to median lick duration for that lick across trials
params.nLicks              = 20; % number of post go cue licks to calculate median lick duration for and warp individual trials to

params.lowFR               = 1; % remove clusters with firing rates across all trials less than this val

% set conditions to calculate PSTHs for
params.condition(1)     = {'(hit|miss|no)'};                             % all trials
params.condition(end+1) = {'hit&~stim.enable&~autowater&~early'};               % 2AFC hits, no stim
params.condition(end+1) = {'hit&~stim.enable&autowater&~early'};                % AW hits, no stim
params.condition(end+1) = {'miss&~stim.enable&~autowater&~early'};              % 2AFC miss, no stim, aw off
params.condition(end+1) = {'miss&~stim.enable&autowater&~early'};               % AW miss, no stim

params.tmin = -3;
params.tmax = 2.5;
params.dt = 1/100;

% smooth with causal gaussian kernel
params.smooth = 15;

% cluster qualities to use
params.quality = {'all'}; % accepts any cell array of strings - special character 'all' returns clusters of any quality

params.feat_varToExplain = 80; % num factors for dim reduction of video features should explain this much variance

params.N_varToExplain = 80; % keep num dims that explains this much variance in neural data (when doing n/p)

params.advance_movement = 0;

params.bctype = 'reflect'; % options are : reflect  zeropad  none
%% SPECIFY DATA TO LOAD

datapth = 'C:\Users\Jackie\Documents\Grad School\Economo Lab';

meta = [];

% --- ALM ---
meta = loadJEB6_ALMVideo(meta,datapth);
meta = loadJEB7_ALMVideo(meta,datapth);
meta = loadEKH1_ALMVideo(meta,datapth);
meta = loadEKH3_ALMVideo(meta,datapth);
meta = loadJGR2_ALMVideo(meta,datapth);
meta = loadJGR3_ALMVideo(meta,datapth);

params.probe = {meta.probe}; % put probe numbers into params, one entry for element in meta, just so i don't have to change code i've already written
%% LOAD DATA

% ----------------------------------------------
% -- Neural Data --
% obj (struct array) - one entry per session
% params (struct array) - one entry per session
% ----------------------------------------------
[obj,params] = loadSessionData(meta,params);
%% Find trials in which the animal switches between contexts
for sessix = 1:numel(meta)
    bp = obj(sessix).bp;
    [toAW_ix, toAFC_ix] = findSwitchTrials(bp);
    % toAW_ix = the first trial in an AW block; toAFC_ix = the first trial in a 2AFC block
    obj(sessix).toAW_ix = toAW_ix;  obj(sessix).toAFC_ix = toAFC_ix;
end
%% Find avg context mode across trials in a given session on switch trials
clear temp

sessix = 1;

nBufferTrials = 10;                              % Number of trials that you want to look at before and after switches
% Time period over which you want to average a cell's FR
trialstart = median(obj(1).bp.ev.bitStart)-median(obj(1).bp.ev.(params(1).alignEvent));
start = find(obj(1).time>trialstart,1,'first');
samp = median(obj(1).bp.ev.sample)-median(obj(1).bp.ev.(params(1).alignEvent));
stop = find(obj(1).time<samp,1,'last');
nCells = size(obj(sessix).trialdat,2);

for c = 1:nCells
    cellPSTH = obj(sessix).trialdat(:,c,:);         % (time x trials) for this cell

    % Find avg CDContext on 2AFC --> AWswitch trials (during the presample period)
    switchtype = 'toAW_ix';
    toAW = findCDCont_SwitchAligned(nBufferTrials, obj, sessix, cellPSTH,switchtype,start,stop);

    % Find avg CDContext on AW --> 2AFC switch trials (during the presample period)
    switchtype = 'toAFC_ix';
    to2AFC = findCDCont_SwitchAligned(nBufferTrials, obj, sessix, cellPSTH,switchtype,start,stop);

    
    tRange = -nBufferTrials:nBufferTrials;           % Number of trials that you want to plot
    probenum = meta(sessix).probe;
    cellqual = obj(sessix).clu{probenum}(c).quality;
    plotSwitchAlignedFR(toAW, to2AFC, tRange,nBufferTrials,c,cellqual)
    pause
end
%% FUNCTIONS

function stdCD = getStd(toAW,to2AFC)

stdCD.toAW = std(toAW,0,1,'omitnan');     % Get standard deviation of switch aligned presamp CDContext across all trials
stdCD.to2AFC = std(to2AFC,0,1,'omitnan');
end

function blah = normalizeCDCont(blah)
maxblah = max(abs(blah.toAW_CDCont)); blah.toAW_CDCont = blah.toAW_CDCont./maxblah;
maxblah = max(abs(blah.to2AFC_CDCont)); blah.to2AFC_CDCont = blah.to2AFC_CDCont./maxblah;
end

% Plotting
function plotSwitchAlignedFR(toAW, to2AFC, tRange,nBufferTrials,cellnum, cellqual)
subplot(1,2,1)
toplot = mean(toAW,1,'omitnan');
plot(tRange,toplot,'LineWidth',2)
hold on;
ax = gca;
xline(0,'k--')
title('2AFC to AW switches')
xlabel('Trials to context switch')
xlim([-nBufferTrials, nBufferTrials]);
ylabel('Firing rate')
hold off;

subplot(1,2,2)
toplot = mean(to2AFC,1,'omitnan');
plot(tRange,toplot,'LineWidth',2)
hold on;
ax = gca;
xline(0,'k--')
title('AW to 2AFC switches')
xlabel('Trials to context switch')
xlim([-nBufferTrials, nBufferTrials]);
ylabel('Firing rate')
hold off;

sgtitle(['Cell #' num2str(cellnum) '; Quality = ' cellqual])
end