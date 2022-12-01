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
    {'top_tongue','topleft_tongue','bottom_tongue','bottomleft_tongue','jaw','top_paw','bottom_paw','top_nostril','bottom_nostril'}};

params.feat_varToExplain = 80; % num factors for dim reduction of video features should explain this much variance

params.N_varToExplain = 80; % keep num dims that explains this much variance in neural data (when doing n/p)

params.advance_movement = 0;
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
%% Calculate MOVE-CDCont and find single trial projections
clearvars -except obj meta params me sav kin

disp('----Calculating MoveCDContext Mode----')
cond2use = [2,3];
cond2proj = [2,3];
regr = getMovement_CDContext(obj,params,cond2use,cond2proj,kin);

disp('----Projecting single trials of movement onto MoveCDContext----')
cd = 'contextPresamp';
regr = getMove_SingleTrialProjs(regr,obj,cd,kin);
%% Calculate CDCont and find single trial projections

disp('----Calculating CDContext Mode----')
cond2use = [2,3];
cond2proj = [2,3];
neur = getCDContext(obj,params,cond2use,cond2proj);

disp('----Projecting single trials of movement onto CDContext----')
cd = 'contextPresamp';
neur = getSingleTrialProjs(neur,obj,cd);
%% Find avg CDCont and MOVE-CDCont during the presample period on each trial
% Time period over which you want to average CDContext
start = find(obj(1).time>-3,1,'first');
samp = median(obj(1).bp.ev.sample)-median(obj(1).bp.ev.(params(1).alignEvent));
stop = find(obj(1).time<samp,1,'last');

for sessix = 1:length(meta)
    regr(sessix).presampAvg = mean(regr(sessix).singleMoveCDProj(start:stop,:),1,'omitnan');
    neur(sessix).presampAvg = mean(neur(sessix).singleProj(start:stop,:),1,'omitnan');
end
%% Make a scatter plot for each session of CDCont vs MOVE-CDCont during presample period on each trial
for sessix = 1:length(meta)
    figure();
    scatter(regr(sessix).presampAvg, neur(sessix).presampAvg)
    xlabel('MOVE-CDContext')
    ylabel('Neural CDContext')
end