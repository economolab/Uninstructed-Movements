% DECODING CDContext FROM ALL KINEMATIC FEATURES
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
%% Calculate all CDs and find single trial projections
clearvars -except obj meta params me sav kin

disp('----Calculating Context Mode----')
cond2use = [2,3];
cond2proj = [2,3];
regr = getCDContext(obj,params,cond2use,cond2proj);

disp('----Projecting single trials onto CDContext----')
cd = 'contextPresamp';
regr = getSingleTrialProjs(regr,obj,cd);
%% Load kinematic data
nSessions = numel(meta);
for sessix = 1:numel(meta)
    message = strcat('----Getting kinematic data for session',{' '},num2str(sessix), {' '},'out of',{' '},num2str(nSessions),'----');
    disp(message)
    kin(sessix) = getKinematics(obj(sessix), me(sessix), params(sessix));
end
%% Predict CDContext from DLC features

clearvars -except datapth kin me meta obj params regr nSessions

% params
rez.nFolds = 4;                                     % number of iterations (bootstrap)

rez.binSize = 30;                                   % Bin size to decode over (in milliseconds)
difTime = params(1).dt*1000;                        % Convert time-difference from s to ms
rez.dt = floor(rez.binSize / difTime);              % How many samples confer your bin size
rez.tm = obj(1).time(1:rez.dt:numel(obj(1).time));  % Create a new time axis over which to decode (based on the bin size that you want)
rez.numT = numel(rez.tm);                           % Number of time-points

rez.train = 1;                                      % fraction of trials to use for training (1-train for testing)

% match number of 2AFC and AW hits, and 2AFC and AW misses
cond2use = 2:5;                                     % Which of the 'params.conditions' you are including in the decoding 
hitcond = [1 2];                                    % Out of the cond2use conditions, specify which are hits and which are misses
misscond = [3 4];

% True Values of CDContext for each session
trueVals.AFChit = cell(nSessions,1);
trueVals.FWhit = cell(nSessions,1);

% Model prediction of CDContext for each session, predicted by each feature
modelpred.AFChit = cell(nSessions,1);
modelpred.FWhit = cell(nSessions,1);

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
    pred = mySmooth((pred'),31);

     % Divide true data and model predictions into R and L hit trials
    if rez.train==1                                                     % If all data is used to train the model (Cross-validated, the true data is just all of it i.e. the training data)
        trials.AFCHit.TrainIX = ismember(trials.train,trials_hit{1});       
        trials.AFCHit.Train = trials.train(trials.AFCHit.TrainIX);             
        trials.FWHit.TrainIX = ismember(trials.train,trials_hit{2});       
        trials.FWHit.Train = trials.train(trials.FWHit.TrainIX);
    else                                                                % If model is not cross-validated and just doing a train/test split
        trials.AFCHit.TestIX = ismember(trials.test,trials_hit{1});       % Logical array for all of the test trials indicating whether they were a right hit or not
        trials.AFCHit.Test = trials.test(trials.AFCHit.TestIX);             % Get the trial numbers that are a 2AFC hit and were used in the test
        trials.FWHit.TestIX = ismember(trials.test,trials_hit{2});       % Logical array for all of the test trials indicating whether they were a left hit or not
        trials.FWHit.Test = trials.test(trials.FWHit.TestIX);             % Get the trial numbers that are a AW hit and were used in the test
    end

   
    % Label test data as true data and model prediction
    if rez.train==1                                                          % If you are using the whole data set as the training set (if cross-validating), then the test data is just the train data
        trueVals.AFChit{sessix} = in.train.y(:,trials.AFCHit.TrainIX);
        trueVals.FWhit{sessix} = in.train.y(:,trials.FWHit.TrainIX);
        modelpred.AFChit{sessix} = pred(:,trials.AFCHit.TrainIX);
        modelpred.FWhit{sessix} = pred(:,trials.FWHit.TrainIX);
    else                                                                     % If you are using a train/test split, then the assign the test data as the data not used for training
        trueVals.AFChit{sessix} = in.test.y(:,trials.AFCHit.TestIX);
        trueVals.FWhit{sessix} = in.test.y(:,trials.FWHit.TestIX);
        modelpred.AFChit{sessix} = pred(:,trials.AFCHit.TestIX);
        modelpred.FWhit{sessix} = pred(:,trials.FWHit.TestIX);
    end
end
%% Plot a scatter plot for a single session of true CDContext and predicted CDContext for each trial
delR2 = [];
% Each dot = an average value of CDContext during the presample period 
trialstart = median(obj(1).bp.ev.bitStart)-median(obj(1).bp.ev.(params(1).alignEvent));
start = find(obj(1).time>trialstart,1,'first');
samp = median(obj(1).bp.ev.sample)-median(obj(1).bp.ev.(params(1).alignEvent));
stop = find(obj(1).time<samp,1,'last');
for sessix = 1:length(meta)
    tempR = mean(trueVals.AFChit{sessix}((start+1):(stop+1),:),1,'omitnan');        % For each trial, get the average CDlate during the delay period
    tempL = mean(trueVals.FWhit{sessix}((start+1):(stop+1),:),1,'omitnan');
    truedat = [tempR, tempL];

    tempR = mean(modelpred.AFChit{sessix}(start:stop,:),1,'omitnan');        % For each trial, get the average predicted CDlate during the delay period
    tempL = mean(modelpred.FWhit{sessix}(start:stop,:),1,'omitnan');
    modeldat = [tempR, tempL];

    R2 = corrcoef(truedat,modeldat);
    R2 = R2(2);
    delR2 = [delR2, R2];
    R2 = num2str(R2);
    figure();
    coeff = polyfit(truedat,modeldat,1);                   % Find the line of best fit
    sesstitle = strcat(meta(sessix).anm, {' '},meta(sessix).date);
    scatter(truedat,modeldat,'filled','black')
    hline = refline(coeff);
    hline.LineStyle = '--';
    hline.Color = 'k';
    xlabel('True data')
    ylabel('Model prediction')
    legend('data',['R^2 = ' R2],'Location','best')
    title(sesstitle)
end
%% Plot bar plot to show average R2 values across sessions
figure();
bar(mean(delR2),'FaceColor',[0.75 0.75 0.75]); hold on;
scatter(1,delR2,'filled')
%% For each trial, find the correlation between the true CDcontext and the predicted CDcontext (in the given time range
% Get R^2 values between true CDlate and prediction on each trial
%Specify time-points over which you want to find the correlation 
% stop = find(obj(1).time<0.05,1,'last');
% timeTrue = 2:stop;
% timePred = 1:(stop-1);

% R2 = correlateTrue_PredictedCD(trueVals, modelpred,meta,timeTrue,timePred);
% %% Plot histogram of R2 values
% % nBins = 5;
% % figure();
% % histogram(R2,nBins,'FaceColor',[0.25 0.25 0.25])
% % xlabel('R^2 value')
% % ylabel('Num sessions')
%% Plot an example session of CDlate prediction vs true value
for sessix = 1:length(meta)
    sesstitle = strcat(meta(sessix).anm, {' '},meta(sessix).date);

    % Plot all true and predicted CDContext projections for 2AFC trials
%     figure();
%     subplot(2,1,1)
%     plot(obj(1).time,trueVals.AFChit{sessix},'Color',[0 0 1]); hold on;
%     plot(rez.tm(1:end-1),modelpred.AFChit{sessix},'Color',[0.5 0.5 1]);
% 
%     % Plot all true and predicted CDContext projections for 2AW trials 
%     subplot(2,1,2)
%     plot(obj(1).time,trueVals.FWhit{sessix},'Color',[1 0 0]); hold on;
%     plot(rez.tm(1:end-1),modelpred.FWhit{sessix},'Color',[1 0.5 0.5]);

    % Calculate averages and standard deviation for true CD and predicted CD
    % for this session
    [avgCD,stdCD] = getAvgStd(trueVals,modelpred,sessix);

    % Get confidence intervals
    [upperci, lowerci] = getConfInt(meta, avgCD, stdCD);

    % Plot average CDlate projections and predictions for the session
    plotExampleCDContextPrediction(obj,rez,avgCD,upperci,lowerci,sesstitle)
end
%% FUNCTIONS


function [avgCD,stdCD] = getAvgStd(trueVals,modelpred,sessix)

avgCD.AFChit.true = mean(trueVals.AFChit{sessix},2,'omitnan');      % meta(sessix)aGet average true CDlate for R and L hits for this session
avgCD.FWhit.true = mean(trueVals.FWhit{sessix},2,'omitnan');
stdCD.AFChit.true = std(trueVals.AFChit{sessix},0,2,'omitnan');     % Get standard deviation of true CDlate for R and L hits
stdCD.FWhit.true = std(trueVals.FWhit{sessix},0,2,'omitnan');

avgCD.AFChit.pred = mean(modelpred.AFChit{sessix},2,'omitnan');     % Get average predicted CDlate for R and L hits for this session
avgCD.FWhit.pred = mean(modelpred.FWhit{sessix},2,'omitnan');
stdCD.AFChit.pred = std(modelpred.AFChit{sessix},0,2,'omitnan');    % Get stdev of predicted CDlate for R and L hits for this session
stdCD.FWhit.pred = std(modelpred.FWhit{sessix},0,2,'omitnan');
end

function [upperci, lowerci] = getConfInt(meta, avgCD, stdCD)
nSessions = length(meta);

upperci.AFC.true = avgCD.AFChit.true(2:end)+1.96*(stdCD.AFChit.true(2:end)/nSessions);  % Find the upper 95% confidence interval for each condition
lowerci.AFC.true = avgCD.AFChit.true(2:end)-1.96*(stdCD.AFChit.true(2:end)/nSessions);  % Find lower 95% condifence interval for each condition
upperci.FW.true = avgCD.FWhit.true(2:end)+1.96*(stdCD.FWhit.true(2:end)/nSessions);  
lowerci.FW.true = avgCD.FWhit.true(2:end)-1.96*(stdCD.FWhit.true(2:end)/nSessions);  

upperci.AFC.pred = avgCD.AFChit.pred(2:end)+1.96*(stdCD.AFChit.pred(2:end)/nSessions);  
lowerci.AFC.pred = avgCD.AFChit.pred(2:end)-1.96*(stdCD.AFChit.pred(2:end)/nSessions);  
upperci.FW.pred = avgCD.FWhit.pred(2:end)+1.96*(stdCD.FWhit.pred(2:end)/nSessions);  
lowerci.FW.pred = avgCD.FWhit.pred(2:end)-1.96*(stdCD.FWhit.pred(2:end)/nSessions);  
end