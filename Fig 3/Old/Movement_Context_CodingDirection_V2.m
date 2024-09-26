% Calculating a Movement-CDContext (instead of a CD being a weight assigned
% to each neuron, its a weight assigned to each kinematic feature)
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
% If you want the standardized kineamtic data to be used in finding the Move-CDCont, set the last input to be 'standardize'
regr = getMovement_CDContext(obj,params,cond2use,cond2proj,kin,'standardize');

disp('----Projecting single trials of movement onto MoveCDContext----')
cd = 'context';
regr = getMove_SingleTrialProjs(regr,obj,cd,kin);

%% For single sessions, show all trials projected onto MoveCDContext, and condition-averaged
colors = {'black','magenta'};
alpha = 0.2;
for sessix = 1:length(meta)
    figure();
    for c = 1:length(cond2use)
        condtrix = params(sessix).trialid{cond2use(c)};
        temp = regr(sessix).singleMoveCDProj(:,condtrix);
        %plot(obj(1).time,temp,'Color',colors{c}); hold on;
        regr(sessix).condavgProj(:,c) = mean(temp,2,'omitnan');
        temperr = 1.96*(std(temp,0,2,'omitnan')/length(condtrix));
        plot(obj(1).time,mean(temp,2,'omitnan'),'Color',colors{c},'LineWidth',2); hold on
        ax = gca;
        shadedErrorBar(obj(1).time, mean(temp,2,'omitnan'), temperr ,{'Color',colors{c},'LineWidth',2}, alpha, ax)
    end
end
%% Get condition-averaged MoveCDContext across all sessions
% Take avg across sessions 
sessavg = NaN(length(obj(1).time),length(cond2use));        % (time x conditions)
sessstd = NaN(length(obj(1).time),length(cond2use));
for c = 1:length(cond2use)
    temp = [];
    for sessix = 1:length(meta)
        temp = [temp,regr(sessix).condavgProj(:,c)];
    end
    sessavg(:,c) = mean(temp,2,'omitnan');
    sessstd(:,c) = std(temp,0,2,'omitnan');
end
%% Plot condition-averaged MoveCDContext across all sessions
samp = median(obj(1).bp.ev.sample)-median(obj(1).bp.ev.(params(1).alignEvent));
delay = median(obj(1).bp.ev.delay)-median(obj(1).bp.ev.(params(1).alignEvent));
trialstart = median(obj(1).bp.ev.bitStart)-median(obj(1).bp.ev.(params(1).alignEvent));
% Set plot params and plot
colors = {'black','magenta'};
alpha = 0.2;
figure();

for c = 1:length(cond2use)
    temperr = 1.96*(sessstd(:,c)/length(meta));
    plot(obj(1).time,sessavg(:,c),'Color',colors{c},'LineWidth',2); hold on
    ax = gca;
    shadedErrorBar(obj(1).time, sessavg(:,c), temperr ,{'Color',colors{c},'LineWidth',2}, alpha, ax)
end
xline(0,'k','Linestyle','--')
xline(delay,'k','Linestyle','-.')
xline(samp,'k','Linestyle','-.')
xlabel('Time from go cue (s)')
ylabel('Projection onto Movement-CDContext (a.u.)')
xlim([trialstart 2.5])
