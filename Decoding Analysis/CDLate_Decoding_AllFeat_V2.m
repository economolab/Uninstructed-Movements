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

params.tmin = -2.5;
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
clearvars -except obj meta params me sav kin

disp('----Calculating coding dimensions----')
cond2use = [2,3];
cond2proj = [2,3];
regr = getCodingDimensions(obj,params,cond2use,cond2proj);

disp('----Projecting single trials onto CDlate----')
cd = 'late';
regr = getSingleTrialProjs(regr,obj,cd);
%% Load kinematic data
nSessions = numel(meta);
for sessix = 1:numel(meta)
    message = strcat('----Getting kinematic data for session',{' '},num2str(sessix), {' '},'out of',{' '},num2str(nSessions),'----');
    disp(message)
    kin(sessix) = getKinematics(obj(sessix), me(sessix), params(sessix));
end
%% Predict CDLate from DLC features

clearvars -except datapth kin me meta obj params regr nSessions

% params
rez.nFolds = 4;                                     % number of iterations (bootstrap)

rez.binSize = 30;                                   % Bin size to decode over (in milliseconds)
difTime = params(1).dt*1000;                        % Convert time-difference from s to ms
rez.dt = floor(rez.binSize / difTime);              % How many samples confer your bin size
rez.tm = obj(1).time(1:rez.dt:numel(obj(1).time));  % Create a new time axis over which to decode (based on the bin size that you want)
rez.numT = numel(rez.tm);                           % Number of time-points

rez.train = 0.5;                                      % fraction of trials to use for training (1-train for testing)

% match number of right and left hits, and right and left misses
cond2use = 2:5;
hitcond = [1 2];
misscond = [3 4];

% True Values of CDlate for each session
trueVals.Rhit = cell(nSessions,1);
trueVals.Lhit = cell(nSessions,1);

% Model prediction of CDlate for each session, predicted by each feature
modelpred.Rhit = cell(nSessions,1);
modelpred.Lhit = cell(nSessions,1);


mask = true(numel(kin(1).featLeg),1);
rez.featix = find(mask);                    % Indices in the feature legend that correspond to specified body part


% Do the decoding for each session
for sessix = 1:numel(obj)
    disp(['Decoding session ' num2str(sessix) ' / ' num2str(numel(obj))])

    % Getting the proper number of trials
    [trials,trials_hit] = EqualizeTrialNums(params(sessix),cond2use,hitcond,misscond);

    % Organize predictor and regressors from the proper trials and
    % kinematic features
    [Y,X] = getPredictorsRegressors(trials,regr(sessix),kin(sessix),rez);

    % Train/Test Split
    [trials,in] = TrainTestSplit(trials,Y,X,rez);

    % Decoding
    pred = DLC_CD_Decoder(in,rez);
    pred = pred';

    % Divide true data and model predictions into R and L hit trials
    trials.RHit.TestIX = ismember(trials.test,trials_hit{1});       % Logical array for all of the test trials indicating whether they were a right hit or not
    trials.RHit.Test = trials.test(trials.RHit.TestIX);             % Get the trial numbers that are a R hit and were used in the test
    trials.LHit.TestIX = ismember(trials.test,trials_hit{2});       % Logical array for all of the test trials indicating whether they were a left hit or not
    trials.LHit.Test = trials.test(trials.LHit.TestIX);             % Get the trial numbers that are a L hit and were used in the test

    % Label test data
    trueVals.Rhit{sessix} = in.test.y(:,trials.RHit.TestIX);
    trueVals.Lhit{sessix} = in.test.y(:,trials.LHit.TestIX);
    modelpred.Rhit{sessix} = pred(:,trials.RHit.TestIX);
    modelpred.Lhit{sessix} = pred(:,trials.LHit.TestIX);



    % %     % shuffle labels for a 'null' distribution
    % %
    % %
    % %     Y = randsample(Y,numel(Y));
    % %
    % %     % train/test split
    % %
    % %     in.train.y = Y(trials.trainidx);
    % %     in.test.y  = Y(trials.testidx);
    % %
    % %     acc_shuffled(:,sessix) = DLC_ChoiceDecoder(in,rez,trials);
end
%% Get R^2 values between true CDlate and prediction on each trial
R2 = correlateTrue_PredictedCD(trueVals, modelpred);
%% Plot histogram of R2 values
nBins = 8;
histogram(R2,nBins,'FaceColor',[0.25 0.25 0.25])
xlabel('R^2 value')
ylabel('Num sessions')
%% Plot an example session of CDlate prediction vs true value
sessix = 1;                                                     % Which session you want to plot for
sesstitle = strcat(meta(sessix).anm, {' '},meta(sessix).date);

% Calculate averages and standard deviation for true CD and predicted CD
% for this session
[avgCD,stdCD] = getAvgStd(trueVals,modelpred);

% Get confidence intervals
[upperci, lowerci] = getConfInt(meta, avgCD, stdCD);

% Plot average CDlate projections and predictions for the session
plotExampleCDLatePrediction(obj,rez,avgCD,upperci,lowerci,sesstitle)
%% FUNCTIONS
function [avgCD,stdCD] = getAvgStd(trueVals,modelpred)

avgCD.Rhit.true = mean(trueVals.Rhit{sessix},2,'omitnan');      % Get average true CDlate for R and L hits for this session
avgCD.Lhit.true = mean(trueVals.Lhit{sessix},2,'omitnan');
stdCD.Rhit.true = std(trueVals.Rhit{sessix},0,2,'omitnan');     % Get standard deviation of true CDlate for R and L hits
stdCD.Lhit.true = std(trueVals.Lhit{sessix},0,2,'omitnan');

avgCD.Rhit.pred = mean(modelpred.Rhit{sessix},2,'omitnan');     % Get average predicted CDlate for R and L hits for this session
avgCD.Lhit.pred = mean(modelpred.Lhit{sessix},2,'omitnan');
stdCD.Rhit.pred = std(modelpred.Rhit{sessix},0,2,'omitnan');    % Get stdev of predicted CDlate for R and L hits for this session
stdCD.Lhit.pred = std(modelpred.Lhit{sessix},0,2,'omitnan');
end