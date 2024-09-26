% Aligning presample MOVE-CDCont to switches in context blocks
clear,clc,close all

% add paths for data loading scripts, all fig funcs, and utils
utilspth = 'C:\Users\Jackie\Documents\Grad School\Economo Lab\Code\Munib Uninstruct Move\uninstructedMovements_v2';
addpath(genpath(fullfile(utilspth,'DataLoadingScripts')));
addpath(genpath(fullfile(utilspth,'funcs')));
addpath(genpath(fullfile(utilspth,'utils')));
addpath(genpath(fullfile(utilspth,'fig1')));
figpth = 'C:\Users\Jackie\Documents\Grad School\Economo Lab\Code\Uninstructed-Movements\Fig 3';
addpath(genpath(fullfile(figpth,'funcs')));

%% PARAMETERS
params.alignEvent          = 'goCue'; % 'jawOnset' 'goCue'  'moveOnset'  'firstLick'  'lastLick'

% time warping only operates on neural data for now.
% TODO: time warp for video and bpod data
params.timeWarp            = 0;  % piecewise linear time warping - each lick duration on each trial gets warped to median lick duration for that lick across trials
params.nLicks              = 20; % number of post go cue licks to calculate median lick duration for and warp individual trials to

params.lowFR               = 1; % remove clusters with firing rates across all trials less than this val

% set conditions to calculate PSTHs for
params.condition(1)     = {'(hit|miss|no)'};                             % all trials
params.condition(end+1) = {'~no&~stim.enable&~autowater&~early'};               % 2AFC response trials, no stim
params.condition(end+1) = {'~no&~stim.enable&autowater&~early'};                % AW hits response trials, no stim

params.tmin = -3;
params.tmax = 2.5;
params.dt = 1/100;

% smooth with causal gaussian kernel
params.smooth = 15;

% cluster qualities to use
params.quality = {'all'}; % accepts any cell array of strings - special character 'all' returns clusters of any quality


params.traj_features = {{'tongue','left_tongue','right_tongue','jaw','trident','nose'},...
    {'top_tongue','topleft_tongue','bottom_tongue','bottomleft_tongue','jaw','top_nostril','bottom_nostril'}};

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

% ------------------------------------------
% -- Motion Energy --
% me (struct array) - one entry per session
% ------------------------------------------
for sessix = 1:numel(meta)
    me(sessix) = loadMotionEnergy(obj(sessix), meta(sessix), params(sessix), datapth);
end

%% Load kinematic data
nSessions = numel(meta);
for sessix = 1:numel(meta)
    message = strcat('----Getting kinematic data for session',{' '},num2str(sessix), {' '},'out of',{' '},num2str(nSessions),'----');
    disp(message)
    kin(sessix) = getKinematics(obj(sessix), me(sessix), params(sessix));
end
%% Calculate all CDs and find single trial projections
clearvars -except obj meta params me sav kin

disp('----Calculating MoveCDContext Mode----')
cond2use = [2,3];
cond2proj = [2,3];
regr = getMovement_CDContext(obj,params,cond2use,cond2proj,kin,'standardize');

disp('----Projecting single trials of movement onto MoveCDContext----')
cd = 'context';
regr = getMove_SingleTrialProjs(regr,obj,cd,kin);
%% Find trials in which the animal switches between contexts
for sessix = 1:numel(meta)
    bp = obj(sessix).bp;
    [toAW_ix, toAFC_ix] = findSwitchTrials(bp);
    % toAW_ix = the first trial in an AW block; toAFC_ix = the first trial in a 2AFC block
    obj(sessix).toAW_ix = toAW_ix;  obj(sessix).toAFC_ix = toAFC_ix;
end
%% Find avg context mode across trials in a given session on switch trials
clear temp
nBufferTrials = 10;                              % Number of trials that you want to look at before and after switches
% Time period over which you want to average CDContext
trialstart = median(obj(1).bp.ev.bitStart)-median(obj(1).bp.ev.(params(1).alignEvent));
start = find(obj(1).time>trialstart,1,'first');
samp = median(obj(1).bp.ev.sample)-median(obj(1).bp.ev.(params(1).alignEvent));
stop = find(obj(1).time<samp,1,'last');

for sessix = 1:numel(meta)
    CDCont = regr(sessix).singleMoveCDProj;
    
    % Find avg CDContext on 2AFC --> AWswitch trials (during the presample period)
    switchtype = 'toAW_ix';
    contswitch(sessix).toAW_CDCont = findCDCont_SwitchAligned(nBufferTrials, obj, sessix, CDCont,switchtype,start,stop);

    % Find avg CDContext on AW --> 2AFC switch trials (during the presample period)
    switchtype = 'toAFC_ix';
    contswitch(sessix).to2AFC_CDCont = findCDCont_SwitchAligned(nBufferTrials, obj, sessix, CDCont,switchtype,start,stop);
end
%%
% Normalize the CDContext values to the max of the absolute values (for
% each session--i.e. each session will have CDCont values between -1 and 1)
for sessix = 1:length(meta)
    blah = contswitch(sessix);
    blah = normalizeCDCont(blah);
    contswitch(sessix) = blah;
end
%% Plot switch-aligned CDContext
% Concatenate all switch-aligned CDContexts from across sessions
toAW = []; to2AFC = [];
for sessix = 1:length(meta)
    toAW = [toAW; contswitch(sessix).toAW_CDCont];
    to2AFC = [to2AFC; contswitch(sessix).to2AFC_CDCont];
end
col = [0.35 0.35 0.35];
nSessions = length(meta);
tRange = -nBufferTrials:nBufferTrials;           % Number of trials that you want to plot
stdCD = getStd(toAW,to2AFC);                     % Get the standard deviation of CDCont across all trials (from all sessions)
alpha = 0.2;                                     % Opacity of confidence intervals

plotSwitchAlignedCDCont(toAW, to2AFC,stdCD, nSessions, tRange, alpha, col,nBufferTrials)

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
function plotSwitchAlignedCDCont(toAW, to2AFC,stdCD, nSessions, tRange, alpha, col,nBufferTrials)
figure();
subplot(1,2,1)
toplot = mean(toAW,1,'omitnan');
plot(tRange,toplot)
hold on;
ax = gca;
err = 1.96*(stdCD.toAW/nSessions);
shadedErrorBar(tRange, toplot, err ,{'Color',col,'LineWidth',2}, alpha, ax)
xline(0,'k--')
title('2AFC to AW switches')
xlabel('Trials to context switch')
xlim([-nBufferTrials, nBufferTrials]);
ylabel('Normalized MOVE-CDContext proj (a.u.)')

subplot(1,2,2)
toplot = mean(to2AFC,1,'omitnan');
plot(tRange,toplot)
hold on;
ax = gca;
err = 1.96*(stdCD.toAW/nSessions);
shadedErrorBar(tRange, toplot, err ,{'Color',col,'LineWidth',2}, alpha, ax)
xline(0,'k--')
title('AW to 2AFC switches')
xlabel('Trials to context switch')
xlim([-nBufferTrials, nBufferTrials]);
ylabel('Normalized MOVE-CDContext proj (a.u.)')
end

