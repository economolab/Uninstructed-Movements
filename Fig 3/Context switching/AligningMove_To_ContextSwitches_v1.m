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
params.condition(end+1) = {'hit&~stim.enable&~autowater&~early'};               % 2AFC hits, no stim
params.condition(end+1) = {'hit&~stim.enable&autowater&~early'};                % AW hits, no stim
params.condition(end+1) = {'miss&~stim.enable&~autowater&~early'};              % 2AFC miss, no stim, aw off
params.condition(end+1) = {'miss&~stim.enable&autowater&~early'};               % AW miss, no stim

params.tmin = -3;
params.tmax = 2.5;
params.dt = (1/100)*3;

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
%% Find trials in which the animal switches between contexts
for sessix = 1:numel(meta)
    bp = obj(sessix).bp;
    [toAW_ix, toAFC_ix] = findSwitchTrials(bp);
    % toAW_ix = the first trial in an AW block; toAFC_ix = the first trial in a 2AFC block
    obj(sessix).toAW_ix = toAW_ix;  obj(sessix).toAFC_ix = toAFC_ix;
end
%% Do PCA on all move features
clear temp

Ndims = 10; % number of principal components to look at
for sessix = 1:length(meta)
    % coeff is the principal components. Should be of size (Nfeatures,Ndims). Each column
    % is a principal component, and each element of each column is the weight associated with a feature.
    [coeff, score] = doMovePCA(kin(sessix).dat, Ndims);


    % project each condition onto the principal components
    latents = zeros(size(kin(sessix).dat,1),Ndims,size(kin(sessix).dat,2));           % (time,nDims,numTrials)
    sm = 200;                                                               % How much smoothing you want to do

    % The latents are the same as the 'scores' given by the PCA function
    % But these are smoothed and separated in the third dimension by trial
    for i = 1:size(kin(sessix).dat,2)                                       % For all trials...
        for j = 1:Ndims                                                     % For all principal components...
            temp = squeeze(kin(sessix).dat(:,i,:));
            temp(isnan(temp)) = 0;
            latents(:,j,i) = mySmooth(temp*coeff(:,j),sm);   % Assign the 'latents' to be trial's movement features weighted by their coefficients for the given PC
        end                                                                 % And smooth the latent
    end
    moveswitch(sessix).pcaLatents = latents;
end

%% Find avg movement for specified feature across trials in a given session on switch trials
nBufferTrials = 10;                              % Number of trials that you want to look at before and after switches
% Time period over which you want to average Feature Movement
start = find(obj(1).time>-3,1,'first');
stop = find(obj(1).time<0,1,'last');

pc2use = 1;
for sessix = 1:numel(meta)
    Move = moveswitch(sessix).pcaLatents(:,pc2use,:);
    
    % Find avg movement PCA on 2AFC --> AWswitch trials (during the presample period)
    switchtype = 'toAW_ix';
    moveswitch(sessix).toAW_Move = findCDCont_SwitchAligned(nBufferTrials, obj, sessix, Move,switchtype,start,stop);

    % Find avg movement PCA on AW --> 2AFC switch trials (during the presample period)
    switchtype = 'toAFC_ix';
    moveswitch(sessix).to2AFC_Move = findCDCont_SwitchAligned(nBufferTrials, obj, sessix, Move,switchtype,start,stop);
end
%%
% Normalize the Movement values to the max of the absolute values (for
% each session--i.e. each session will have Movement values between -1 and 1)
for sessix = 1:length(meta)
    blah = moveswitch(sessix);
    blah = normalizeCDCont(blah);
    moveswitch(sessix) = blah;
end
%% Plot switch-aligned CDContext
% Concatenate all switch-aligned CDContexts from across sessions
toAW = []; to2AFC = [];
for sessix = 1:length(meta)
    toAW = [toAW; moveswitch(sessix).toAW_Move];
    to2AFC = [to2AFC; moveswitch(sessix).to2AFC_Move];
end

tRange = -nBufferTrials:nBufferTrials;           % Number of trials that you want to plot
stdCD = getStd(toAW,to2AFC);                     % Get the standard deviation of CDCont across all trials (from all sessions)
alpha = 0.2;                                     % Opacity of confidence intervals
col = [0.5 0.5 0.5];

plotSwitchAlignedMove(toAW, to2AFC,stdCD, nSessions, tRange, alpha, col,nBufferTrials,feat)

%% FUNCTIONS

function stdCD = getStd(toAW,to2AFC)

stdCD.toAW = std(toAW,0,1,'omitnan');     % Get standard deviation of switch aligned presamp CDContext across all trials
stdCD.to2AFC = std(to2AFC,0,1,'omitnan');
end

function blah = normalizeCDCont(blah)
maxblah = max(abs(blah.toAW_Move)); blah.toAW_Move = blah.toAW_Move./maxblah;
maxblah = max(abs(blah.to2AFC_Move)); blah.to2AFC_Move = blah.to2AFC_Move./maxblah;
end

% Plotting
function plotSwitchAlignedMove(toAW, to2AFC,stdCD, nSessions, tRange, alpha, col,nBufferTrials,feat)
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
ylabel(['Normalized' feat '(a.u.)'])

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
ylabel(['Normalized' feat '(a.u.)'])
end