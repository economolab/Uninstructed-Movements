% Predicting context from random null and potent dimensions
% -------------------------------------------------------------------------------------
% Using all 2AFC and all AW trials to find the Null and Potent Spaces
% -------------------------------------------------------------------------------------
clear,clc,close all

% add paths
utilspth = 'C:\Users\Jackie\Documents\Grad School\Economo Lab\Code\Munib Uninstruct Move\uninstructedMovements_v2';
addpath(genpath(fullfile(utilspth,'DataLoadingScripts')));
addpath(genpath(fullfile(utilspth,'funcs')));
addpath(genpath(fullfile(utilspth,'utils')));
addpath(genpath(fullfile(utilspth,'fig3')));
figpth = 'C:\Users\Jackie\Documents\Grad School\Economo Lab\Code\Uninstructed-Movements\Fig 4';
addpath(genpath(fullfile(figpth,'funcs')));
addpath(genpath(fullfile(figpth,'Context switching')));
figpth = 'C:\Users\Jackie\Documents\Grad School\Economo Lab\Code\Uninstructed-Movements\Fig 7';
addpath(genpath(fullfile(figpth,'funcs')));
addpath(genpath(fullfile(figpth,'Context_funcs')));
figpth = 'C:\Users\Jackie\Documents\Grad School\Economo Lab\Code\Uninstructed-Movements\Fig 6';
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

params.condition(end+1) = {'R&hit&~stim.enable&~autowater&~early'};             % right hits, no stim, aw off
params.condition(end+1) = {'L&hit&~stim.enable&~autowater&~early'};             % left hits, no stim, aw off
params.condition(end+1) = {'R&miss&~stim.enable&~autowater'};                   % error right, no stim, aw off
params.condition(end+1) = {'L&miss&~stim.enable&~autowater'};                   % error left, no stim, aw off
params.condition(end+1) = {'~early&~stim.enable&~autowater'};                          %  no stim, aw off

params.condition(end+1) = {'R&hit&~stim.enable&autowater&~early'};             % right hits, no stim, aw on
params.condition(end+1) = {'L&hit&~stim.enable&autowater&~early'};             % left hits, no stim, aw on
params.condition(end+1) = {'R&miss&~stim.enable&autowater'};                   % error right, no stim, aw on
params.condition(end+1) = {'L&miss&~stim.enable&autowater'};                   % error left, no stim, aw on
params.condition(end+1) = {'~early&~stim.enable&autowater'};                          %  no stim, aw on

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
%% Null and Potent Space

clearvars -except obj meta params me sav

% -----------------------------------------------------------------------
% -- Curate Input Data --
% zscore single trial neural data (time*trials,neurons), for all trials
% -- Calculate null and potent spaces --
% null space from quiet time points
% potent space from moving time points
% -- Calculate coding directions from null and potent spaces --
% -----------------------------------------------------------------------

for sessix = 1:numel(meta)
    % -- input data
    trialdat_zscored = zscore_singleTrialNeuralData(obj(sessix).trialdat, obj(sessix));

    % -- Calculate the null and potent spaces for each session
    cond2use = [2 3 7 8];    % All 2AFC hit trials, all AW hit trials (NUMBERING ACCORDING TO PARAMS.CONDITION)
    nullalltime = 0;      % use all time points to estimate null space if 1
    cond2proj = 2:11;     % (NUMBERING ACCORDING TO PARAMS.CONDITION)
    rez(sessix) = singleTrial_elsayed_np(trialdat_zscored, obj(sessix), me(sessix), params(sessix), cond2use, cond2proj,nullalltime);
end
%% Build a support vector machine to classify context 
% Wu et al (Axel, Shadlen), 2020: Support vector machine to classify trial-type from activity of single neurons.  
% We used the firing rates in the 500 ms time window before animal's first lick. The decoding capability of each area was 
% estimated by using varying numbers of randomly selected neurons that are recorded simultaneously in a session. 
% The classifier was trained on randomly selected 90% of the trials in each session and then tested on the remaining 10% of the trials.  
% The training/testing was repeated 50 times for every given number of neurons and for all the sessions that may be included. 
% The correct rates from the 50 repetitions were then averaged. 
clearvars -except obj meta params me sav rez

classparams.randomized = 'random';                          % Whether you want the dimensions to be randomly selected or to be picked in order
classparams.nfolds = 4;                                                 % # of folds for cross-validation in SVM classifier
classparams.nIterations = 50;                                           % # of iterations that you randomly select Train trials
classparams.nDims = 12;                                                 % Total number of dimensions that will be utilized as an input to the classifier
conds2class = [6,11];                                       % Specify the conditions that you want to be included in the classifier (Rn: ~early AW and ~early 2AFC)
for sessix = 1:length(meta)
    % Times to look at the Null/Potent PSTHs for each dimension
    trialstart = median(obj(1).bp.ev.bitStart)-median(obj(1).bp.ev.(params(1).alignEvent));
    params(sessix).times.start = find(obj(1).time>trialstart,1,'first');
    samp = median(obj(1).bp.ev.sample)-median(obj(1).bp.ev.(params(1).alignEvent));
    params(sessix).times.stop = find(obj(1).time<samp,1,'last');

    %%% Use SVM to classify trial-type (2AFC or AW) based on pre-sample full neural population projected onto 
    %%% randomly selected #s of null and potent dimensions
    % Null space
    spacename = 'Null';
    percentCorr = doContextClassifier(classparams,params(sessix), obj(sessix), rez(sessix),spacename,conds2class);
    acc(sessix).Null = percentCorr;

    % Potent space
    spacename = 'Potent';
    percentCorr = doContextClassifier(classparams,params(sessix), obj(sessix), rez(sessix),spacename,conds2class);
    acc(sessix).Potent = percentCorr;
end

%% Concatenate the SVM accuracies for all sessions
clearvars -except obj meta params me sav acc classparams rez

all.Null = []; all.Potent = [];
for sessix = 1:length(meta)
    all.Null = [all.Null; acc(sessix).Null];
    all.Potent = [all.Potent; acc(sessix).Potent];
end
%% Plot
alpha = 0.2;
figure();
ax = gca;

% Plot avg accuracy of classifier using Null space dimensions
col = [0 1 0];
temperr = 1.96*(std(all.Null,0,1)/sqrt(length(meta)));
shadedErrorBar(1:classparams.nDims,mean(all.Null,1,'omitnan'),temperr,{'Color',col,'LineWidth',2}, alpha, ax);
hold on;

% Plot avg accuracy of classifier using Potent space dimensions
col = [1 0 1];
temperr = 1.96*(std(all.Potent,0,1)/sqrt(length(meta)));
shadedErrorBar(1:classparams.nDims,mean(all.Potent,1,'omitnan'),temperr,{'Color',col,'LineWidth',2}, alpha, ax);
xlim([1 classparams.nDims])
ylabel('Classifier accuracy (%)')
xlabel('# dimensions used')
%%
sessix = 4;
null = rez(sessix).N_null_psth;
potent = rez(sessix).N_potent_psth;
cond1 = 5;
cond2 = 10;
cond3 = 6;
cond4 = 7;
plot3(mySmooth(potent(55:80,1,cond1),51),mySmooth(potent(55:80,2,cond1),51),mySmooth(null(55:80,1,cond1),51),'r','LineWidth', 2)
hold on;
plot3(mySmooth(potent(55:80,1,cond2),51),mySmooth(potent(55:80,2,cond2),51),mySmooth(null(55:80,1,cond2),51),'b','LineWidth', 2)

% plot3(mySmooth(potent(55:80,1,cond3),51),mySmooth(potent(55:80,2,cond3),51),mySmooth(null(55:80,1,cond3),51),'g','LineWidth', 2)
% hold on;
% plot3(mySmooth(potent(55:80,1,cond4),51),mySmooth(potent(55:80,2,cond4),51),mySmooth(null(55:80,1,cond4),51),'k','LineWidth', 2)
xlabel('Potent dim 1')
ylabel('Potent dim 2')
zlabel('Null dim 1')



