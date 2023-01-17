% DECODING CDlate FROM ALL KINEMATIC FEATURES
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
params.condition(end+1) = {'R&hit&~stim.enable&~autowater&~early'};             % R 2AFC hits, no stim
params.condition(end+1) = {'L&hit&~stim.enable&~autowater&~early'};             % L 2AFC hits, no stim
params.condition(end+1) = {'R&miss&~stim.enable&~autowater&~early'};            % R error 2AFC, no stim, aw off
params.condition(end+1) = {'L&miss&~stim.enable&~autowater&~early'};            % L error 2AFC, no stim

params.tmin = -2.5;
params.tmax = 2.5;
params.dt = (1/100)*3;

% smooth with causal gaussian kernel
params.smooth = 15;

% cluster qualities to use
params.quality = {'all'}; % accepts any cell array of strings - special character 'all' returns clusters of any quality


params.traj_features = {{'tongue','left_tongue','right_tongue','jaw','trident','nose'},...
    {'top_tongue','topleft_tongue','bottom_tongue','bottomleft_tongue','jaw','top_nostril','bottom_nostril'}};

params.feat_varToExplain = 80; % num factors for dim reduction of video features should explain this much variance

params.N_varToExplain = 80; % keep num dims that explains this much variance in neural data (when doing n/p)

params.advance_movement = 0;
%% SPECIFY DATA TO LOAD

datapth = 'C:\Users\Jackie\Documents\Grad School\Economo Lab';

meta = [];

% --- ALM ---
% meta = loadJEB6_ALMVideo(meta,datapth);
% meta = loadJEB7_ALMVideo(meta,datapth);
% meta = loadEKH1_ALMVideo(meta,datapth);
% meta = loadEKH3_ALMVideo(meta,datapth);
% meta = loadJGR2_ALMVideo(meta,datapth);
% meta = loadJGR3_ALMVideo(meta,datapth);
meta = loadJEB13_ALMVideo(meta,datapth);
% meta = loadJEB14_ALMVideo(meta,datapth);
% meta = loadJEB15_ALMVideo(meta,datapth);

params.probe = {meta.probe}; % put probe numbers into params, one entry for element in meta, just so i don't have to change code i've already written

%% LOAD DATA

% ----------------------------------------------
% -- Neural Data --
% obj (struct array) - one entry per session
% params (struct array) - one entry per session
% ----------------------------------------------
[obj,params] = loadSessionData(meta,params);
%%
% ------------------------------------------
% -- Motion Energy --
% me (struct array) - one entry per session
% ------------------------------------------
for sessix = 1:numel(meta)
    disp(['Loading ME for session ' num2str(sessix)])
    me(sessix) = loadMotionEnergy(obj(sessix), meta(sessix), params(sessix), datapth);
end
%% Calculate all CDs and find single trial projections
clearvars -except obj meta params me sav kin

disp('----Calculating coding dimensions----')
cond2use = [2,3];
cond2proj = [2,3];
orthogonalize = 'non-orthog';                                       % Set to orthogonalize if you want the projections to be onto the orthogonalized CDs
regr = getCodingDimensions_2afc(obj,params,cond2use,cond2proj,orthogonalize);

disp('----Projecting single trials onto CDlate----')
cd = 'late';
regr = getSingleTrialProjs(regr,obj,cd,orthogonalize);
%% Load kinematic data
nSessions = numel(meta);
for sessix = 1:numel(meta)
    message = strcat('----Getting kinematic data for session',{' '},num2str(sessix), {' '},'out of',{' '},num2str(nSessions),'----');
    disp(message)
    kin(sessix) = getKinematics(obj(sessix), me(sessix), params(sessix));
end
%% Predict CDLate from DLC features

clearvars -except datapth kin me meta obj params regr nSessions
clc;

% params
rez.nFolds = 4;                                     % number of iterations (bootstrap)

rez.binSize = 30;                                   % Bin size to decode over (in milliseconds)
difTime = params(1).dt*1000;                        % Convert time-difference from s to ms
rez.dt = floor(rez.binSize / difTime);              % How many samples confer your bin size
rez.tm = obj(1).time(1:rez.dt:numel(obj(1).time));  % Create a new time axis over which to decode (based on the bin size that you want)
rez.numT = numel(rez.tm);                           % Number of time-points

rez.train = 1;                                      % fraction of trials to use for training (1-train for testing)

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

    % Organize predictor and regressors from the proper trials and kinematic features
    % If you want to use standardized kinematic data, specify last function input as 'standardize'
    [Y,X] = getPredictorsRegressors(trials,regr(sessix),kin(sessix),rez,'standardize');

    % Train/Test Split
    [trials,in] = TrainTestSplit(trials,Y,X,rez);

    % Decoding
    pred = DLC_CD_Decoder(in,rez);
    pred = mySmooth((pred'),31);

    % Divide true data and model predictions into R and L hit trials
    if rez.train==1                                                     % If all data is used to train the model (Cross-validated, the true data is just all of it i.e. the training data)
        trials.RHit.TrainIX = ismember(trials.train,trials_hit{1});       
        trials.RHit.Train = trials.train(trials.RHit.TrainIX);             
        trials.LHit.TrainIX = ismember(trials.train,trials_hit{2});       
        trials.LHit.Train = trials.train(trials.LHit.TrainIX);
    else                                                                % If model is not cross-validated and just doing a train/test split
        trials.RHit.TestIX = ismember(trials.test,trials_hit{1});       % Logical array for all of the test trials indicating whether they were a right hit or not
        trials.RHit.Test = trials.test(trials.RHit.TestIX);             % Get the trial numbers that are a R hit and were used in the test
        trials.LHit.TestIX = ismember(trials.test,trials_hit{2});       % Logical array for all of the test trials indicating whether they were a left hit or not
        trials.LHit.Test = trials.test(trials.LHit.TestIX);             % Get the trial numbers that are a L hit and were used in the test
    end

    % Label test data as true data and model prediction
    if rez.train==1                                                     % If you are using the whole data set as the training set (if cross-validating), then the test data is just the train data
        trueVals.Rhit{sessix} = in.train.y(:,trials.RHit.TrainIX);
        trueVals.Lhit{sessix} = in.train.y(:,trials.LHit.TrainIX);
        modelpred.Rhit{sessix} = pred(:,trials.RHit.TrainIX);
        modelpred.Lhit{sessix} = pred(:,trials.LHit.TrainIX);
    else                                                                % If you are using a train/test split, then the assign the test data as the data not used for training
        trueVals.Rhit{sessix} = in.test.y(:,trials.RHit.TestIX);
        trueVals.Lhit{sessix} = in.test.y(:,trials.LHit.TestIX);
        modelpred.Rhit{sessix} = pred(:,trials.RHit.TestIX);
        modelpred.Lhit{sessix} = pred(:,trials.LHit.TestIX);
    end
end
%% Get R2 value for all sessions during delay
JEB13_delR2 = [];
% Each dot = an average value of CDlate during the delay period 
start = find(obj(1).time>-0.9,1,'first');
stop = find(obj(1).time<-0.05,1,'last');
for sessix = 1:length(meta)
    tempR = mean(trueVals.Rhit{sessix}((start+1):(stop+1),:),1,'omitnan');        % For each trial, get the average CDlate during the delay period
    tempL = mean(trueVals.Lhit{sessix}((start+1):(stop+1),:),1,'omitnan');
    truedat = [tempR, tempL];

    tempR = mean(modelpred.Rhit{sessix}(start:stop,:),1,'omitnan');        % For each trial, get the average predicted CDlate during the delay period
    tempL = mean(modelpred.Lhit{sessix}(start:stop,:),1,'omitnan');
    modeldat = [tempR, tempL];

    R2 = corrcoef(truedat,modeldat);
    R2 = R2(2);
    JEB13_delR2 = [JEB13_delR2, R2];
end
%% Plot bar plot to show average R2 values across sessions
load('C:\Users\Jackie\Documents\Grad School\Economo Lab\Code\Uninstructed-Movements\Fig 3\Decoding Analysis\delR2Vals_NoJEB13.mat')
delR2_ALL = [nonJEB13.delR2, JEB13_delR2];                              % Group the R2 values for JEB13 and all other animals
anmNames = nonJEB13.anmNames;
JEB13_anmNames = {'JEB13','JEB13','JEB13','JEB13','JEB13','JEB13'};
anmNames_all = [anmNames,JEB13_anmNames];                               % Animal names for each session
nSessions = numel(anmNames_all);

uniqueAnm = unique(anmNames_all);

exsess = 22;                                                            % The index of the session that you want to be highlighted
markerSize = 60;
figure();
bar(mean(delR2_ALL),'FaceColor',[0.75 0.75 0.75]); hold on;         % Plot the average R2 value across all sessions
for sessix = 1:nSessions
    curranm = anmNames_all{sessix};                 % Get the name of the animal for this session
    switch curranm                                  % Switch the marker shape depending on which animal is being plotted
        case uniqueAnm{1}
            shape = 'o';
        case uniqueAnm{2}
            shape = '<';
        case uniqueAnm{3}
            shape = '^';
        case uniqueAnm{4}
            shape = 'v';
        case uniqueAnm{5}
            shape = '>';
        case uniqueAnm{6}
            shape = 'square';
        case uniqueAnm{7}
            shape = 'diamond';
        case uniqueAnm{8}
            shape = 'hexagram';
        case uniqueAnm{9}
            shape = 'pentagram';
    end
    scatter(1,delR2_ALL(sessix),markerSize,'filled',shape,'black'); hold on;
end
scatter(1,delR2_ALL(exsess),markerSize,'filled','pentagram','cyan','MarkerEdgeColor','black')
legend([' ',anmNames_all])
ylim([0.4 1])
ax = gca;
ax.FontSize = 16;