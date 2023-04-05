% DECODING CDlate FROM ALL KINEMATIC FEATURES
clear,clc,close all

% add paths for data loading scripts, all fig funcs, and utils
utilspth = 'C:\Users\Jackie\Documents\Grad School\Economo Lab\Code\Uninstructed-Movements\NullPotent_Context';
addpath(genpath(fullfile(utilspth,'DataLoadingScripts')));
addpath(genpath(fullfile(utilspth,'funcs')));
addpath(genpath(fullfile(utilspth,'utils')));
% rmpath(genpath(fullfile(utilspth,'fig1/')));
% rmpath(genpath(fullfile(utilspth,'mc_stim/')));

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
params.condition(1)     = {'(hit|miss|no)'};                             % all trials
params.condition(end+1) = {'R&hit&~stim.enable&~autowater'};             % R 2AFC hits, no stim
params.condition(end+1) = {'L&hit&~stim.enable&~autowater'};             % L 2AFC hits, no stim
params.condition(end+1) = {'R&miss&~stim.enable&~autowater'};            % R error 2AFC, no stim, aw off
params.condition(end+1) = {'L&miss&~stim.enable&~autowater'};            % L error 2AFC, no stim

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
meta = loadJEB14_ALMVideo(meta,datapth);
meta = loadJEB15_ALMVideo(meta,datapth);

% --- M1TJ ---
% meta = loadJEB14_M1TJVideo(meta,datapth);

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
%% Calculate all CDs and find single trial projections
clearvars -except obj meta params me sav

disp('----Calculating coding dimensions----')
cond2use = [2,3];
cond2proj = [2,3];
rez = getCodingDimensions(obj,params,cond2use,cond2proj);

disp('----Projecting single trials onto CDlate----')
cd = 'late';
rez = getSingleTrialProjs(rez,obj,cd);
%% Load kinematic data
nSessions = numel(meta);
for sessix = 1:numel(meta)
    message = strcat('Getting kinematic data for session',{' '},num2str(sessix), {' '},'out of',{' '},num2str(nSessions));
    disp(message)
    kin(sessix) = getKinematics(obj(sessix), me(sessix), params(sessix));
end
%% Organize the input data for the regression (each column is a single regressor (i.e. kinematic feature))
nSessions = numel(meta);

trainFraction = 0.5;
cond2use = [2,3];
regr = getTrainTestSplit(params,trainFraction,cond2use);

windowlength = 0.2;             % in seconds
for sessix = 1:length(meta)
    message = strcat('Creating regression matrices for session',{' '},num2str(sessix), {' '},'out of',{' '},num2str(nSessions)); 
    disp(message)
    
    regr(sessix).regressorMatrix = organizeRegressorMatrix(obj(sessix),params(sessix),kin(sessix),windowlength,regr(sessix));
    regr(sessix).predictorMatrix = organizePredictorMatrix(obj(sessix),params(sessix),windowlength,regr(sessix),rez(sessix));
end

%% Fit the regression model with the training data
for sessix = 1:3%length(meta)
    message = strcat('Fitting regression model for session',{' '},num2str(sessix), {' '},'out of',{' '},num2str(nSessions)); 
    disp(message)
    
    y = regr(sessix).predictorMatrix;
    X = regr(sessix).regressorMatrix;
    regr(sessix).Model = fitrlinear(X,y);
end
%% Predict CDlate on testing trials from the model

%%

% ------ Set up data to do the regression------%
% x1 = Weight;        % Regressor 1
% x2 = Horsepower;    % Regressor 2
% y = MPG;            % Thing you are trying to predict
% 
% X = [ones(size(x1)) x1 x2 x1.*x2];  % Concatenate all regressors into one matrix

% -------Fit the model-------------------------%
% b = regress(y,X);   % Find the weights for each regressor


% ------Predict values from model--------------%
% x1fit = min(x1):100:max(x1);        % Values of x1 (Regressor 1) that you did not measure and that you want to predict
% x2fit = min(x2):10:max(x2);         % Values of x2 (Regressor 2) that you did not measure and that you want to predict
% YFIT = b(1) + b(2)*X1FIT + b(3)*X2FIT + b(4)*X1FIT.*X2FIT;   % Model prediction of Y (based on your regressor weights, fit for X1Fit and X2Fit)
             

