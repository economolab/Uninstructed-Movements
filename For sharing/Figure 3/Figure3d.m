%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 3d -- Using video data to predict projections onto CDchoice on
% single trials
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear,clc,close all
%%
whichcomp = 'LabPC';                                                % LabPC or Laptop

% Base path for code depending on laptop or lab PC
if strcmp(whichcomp,'LabPC')
    basepth = 'C:\Code';
elseif strcmp(whichcomp,'Laptop')
    basepth = 'C:\Users\Jackie\Documents\Grad School\Economo Lab\Code';
end

% add paths
utilspth = [basepth '\Munib Uninstruct Move\uninstructedMovements_v2'];
addpath(genpath(fullfile(utilspth,'DataLoadingScripts')));
addpath(genpath(fullfile(utilspth,'funcs')));
addpath(genpath(fullfile(utilspth,'utils')));
addpath(genpath(fullfile(utilspth,'fig1')));
figpth = [basepth  '\Uninstructed-Movements\Fig 2'];
addpath(genpath(fullfile(figpth,'funcs')));
%% PARAMETERS
params.alignEvent          = 'goCue';

params.timeWarp            = 0;  % piecewise linear time warping - each lick duration on each trial gets warped to median lick duration for that lick across trials
params.nLicks              = 20; % number of post go cue licks to calculate median lick duration for and warp individual trials to

params.lowFR               = 1; % remove clusters with firing rates across all trials less than this val

% set conditions to calculate PSTHs for
params.condition(1)     = {'(hit|miss|no)'};                                    % (1) all trials
params.condition(end+1) = {'R&hit&~stim.enable&~autowater&~early'};             % (2) right DR hits, no stim, aw off
params.condition(end+1) = {'L&hit&~stim.enable&~autowater&~early'};             % (3) left DR hits, no stim, aw off
params.condition(end+1) = {'R&miss&~stim.enable&~autowater&~early'};            % (4) error right,DR, no stim, aw off
params.condition(end+1) = {'L&miss&~stim.enable&~autowater&~early'};            % (5) error left,DR, no stim, aw off
params.condition(end+1) = {'R&no&~stim.enable&~autowater&~early'};              % (6) no response right,DR, no stim, aw off
params.condition(end+1) = {'L&no&~stim.enable&~autowater&~early'};              % (7) no response left,DR, no stim, aw off
params.condition(end+1) = {'hit&~stim.enable&~autowater&~early'};               % (8) all DR hits, no stim, aw off

% Parameters to set up time-axis
params.tmin = -2.5;             % Min time-point (in s) with respect to alignEvent
params.tmax = 2.5;              % Max time-point (in s) with respect to alignEvent
params.dt = 1/100;              % Time window (s)

% smooth PSTHs with causal gaussian kernel
params.smooth = 15;

% cluster (sorted unit) qualities to use
params.quality = {'all'}; % accepts any cell array of strings - special character 'all' returns clusters of any quality

% Kinematic features to load
params.traj_features = {{'tongue','left_tongue','right_tongue','jaw','trident','nose'},...
    {'top_tongue','topleft_tongue','bottom_tongue','bottomleft_tongue','jaw','top_nostril','bottom_nostril'}};

params.feat_varToExplain = 80; % num factors for dim reduction of video features should explain this much variance

params.N_varToExplain = 80; % keep num dims that explains this much variance in neural data (when doing n/p)

params.advance_movement = 0;
params.bctype = 'reflect'; % options are : reflect  zeropad  none
%% SPECIFY DATA TO LOAD

if strcmp(whichcomp,'LabPC')
    datapth = 'C:\Users\Jackie Birnbaum\Documents\Data';
elseif strcmp(whichcomp,'Laptop')
    datapth = 'C:\Users\Jackie\Documents\Grad School\Economo Lab';
end

meta = [];

% Load all sessions
meta = loadJEB6_ALMVideo(meta,datapth);
meta = loadJEB7_ALMVideo(meta,datapth);
meta = loadEKH1_ALMVideo(meta,datapth);
meta = loadEKH3_ALMVideo(meta,datapth);
meta = loadJGR2_ALMVideo(meta,datapth);
meta = loadJGR3_ALMVideo(meta,datapth);
meta = loadJEB13_ALMVideo(meta,datapth);
meta = loadJEB14_ALMVideo(meta,datapth);
meta = loadJEB15_ALMVideo(meta,datapth);
meta = loadJEB19_ALMVideo(meta,datapth);

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

% Set parameters for calculating coding dimensions
cond2use = [2 3];           % Conditions used for calculating CDchoice (corresponding to PARAMS.CONDITION)
inclramp = 'yes';
rampcond = 8;               % Condition used for calculating CDramp
cond2proj = 2:7;            % Condition-averaged PSTHs to project onto CDchoice
% right hits, left hits, right miss, left miss, right no, left no (corresponding to null/potent psths in rez)
cond2use_trialdat = [2 3];  % for calculating selectivity explained in full neural pop
disp('----Calculating coding dimensions----')
regr = getCodingDimensions_2afc(obj,params,cond2use,cond2proj);


% Project single trials onto CD
cd = 'late';                % 'late' for projecting onto CDchoice
smooth = 60;                % How much to smooth single trial projections by
disp('----Projecting single trials onto CDchoice----')
regr = getSingleTrialProjs(regr,obj,cd,smooth);
%% Load kinematic data
nSessions = numel(meta);
for sessix = 1:numel(meta)
    message = strcat('----Getting kinematic data for session',{' '},num2str(sessix), {' '},'out of',{' '},num2str(nSessions),'----');
    disp(message)
    kin(sessix) = getKinematics(obj(sessix), me(sessix), params(sessix));
end
%% SET DECODING PARAMETERS %%

% input data = neural data (time*trials,neurons)
% output data = kin data   (time*trials,kin feats)

par.pre=6; % time bins prior to output used for decoding
par.post=0; % time bins after output used for decoding
par.dt = params(1).dt; % moving time bin
par.pre_s = par.pre .* params(1).dt; % dt, pre_s, and post_s just here to know how much time you're using
par.post_s = par.post .* params(1).dt;

trialstart = mode(obj(1).bp.ev.bitStart)-mode(obj(1).bp.ev.(params(1).alignEvent)+0.4);     % Timing of when you want to start the decoding
start = find(obj(1).time>trialstart,1,'first');
go = mode(obj(1).bp.ev.goCue)-mode(obj(1).bp.ev.(params(1).alignEvent));                    % Timing of go cue
stop = find(obj(1).time>go,1,'first');

par.timerange = start:stop;

% cross validation folds
par.nFolds = 4;

% data sets
par.train = 0.6; % fraction of trials to be used for training
par.test = 1 - par.train;   % fraction of trials to be used for testing

% features to use to decode
par.feats = {'motion','nos','jaw'};         % Use all kinematic features measured of motion energy, nose, and jaw
temp = cellfun(@(x) patternMatchCellArray(kin(1).featLeg,{x},'all') , par.feats,'UniformOutput',false);
par.feats = cat(1, temp{:});

% trials to be used in the decoding
par.cond2use = [2 3]; % R and L hits, ~autowater (with reference to params.conditions)

% regularization parameters
par.regularize = 1; % 1 = ridge regression
par.lambdas = logspace(-3,3,20); % regularization parameters to search over

%% DECODING -- Predict CDchoice from DeepLabCut features

close all

for sessix = 1:numel(meta)          % For all sessions...
    disp(['Decoding for session ' ' ' num2str(sessix) ' / ' num2str(numel(meta))])

    % Organize data into predictors and regressors
    [X,Y,par] = preparePredictorsRegressors_v2(par, sessix, kin, regr,params);

    % Search over regularization parameters for the best model
    for ilambda = 1:numel(par.lambdas)
        lambda = par.lambdas(ilambda);
        cvmdl{ilambda} = fitrlinear(X.train,Y.train,'Learner','svm','KFold',par.nFolds,'Regularization','ridge','Lambda',lambda);
        loss(ilambda) = cvmdl{ilambda}.kfoldLoss;
    end

    % Pick the model with the best lambda
    [~,bestmodel] = min(loss);
    mdl = cvmdl{bestmodel};
    % can predict test data using best cv mdl
    for i = 1:par.nFolds
        mdl_loss(i) = mdl.Trained{i}.loss(X.train,Y.train);
    end
    [~,bestmodel] = min(mdl_loss);
    testmdl = mdl.Trained{bestmodel};

    pred = testmdl.predict(X.test);   % use left-out test data for the prediction

    % Save the predictor coefficients from the testmdl
    loadings(:,sessix) = testmdl.Beta;  % Average the coefficients for each predictor term across folds; save these for each time point

    % Save original input data and predicted data
    y = reshape(Y.test,Y.size.test(1),Y.size.test(2)); % original input data (standardized)
    yhat = reshape(pred,Y.size.test(1),Y.size.test(2)); % prediction

    % find which trials were Rhit and which were Lhit in test set
    Rhit_trials = ismember(par.trials.test,par.trials.Rhit);
    Lhit_trials = ismember(par.trials.test,par.trials.Lhit);

    % neural data, split into trial type
    trueVals.Rhit{sessix} = y(:,Rhit_trials);
    trueVals.Lhit{sessix} = y(:,Lhit_trials);

    % predicted neural data, split into ground truth
    modelpred.Rhit{sessix} = yhat(:,Rhit_trials);
    modelpred.Lhit{sessix} = yhat(:,Lhit_trials);
end

disp('---FINISHED DECODING FOR ALL SESSIONS---')
clearvars -except datapth kin me meta obj params regr nSessions exsess modelpred trueVals par loadings
%% Baseline subtract CDchoice
% Times that you want to use to baseline subtract CDchoice
trialstart = mode(obj(1).bp.ev.bitStart)-mode(obj(1).bp.ev.(params(1).alignEvent));
start = find(obj(1).time>trialstart,1,'first');
samp = mode(obj(1).bp.ev.sample)-mode(obj(1).bp.ev.(params(1).alignEvent));
stop = find(obj(1).time<samp,1,'last');

cond2use = {'Rhit','Lhit'};
for sessix = 1:length(meta)             % For each session...
    for c = 1:length(cond2use)          % For each condition...
        cond = cond2use{c};
        curr = modelpred.(cond){sessix};    % Get the predicted CDchoice from that condition
        presampcurr = curr(start:stop,:);   % Get pre-sample values of predicted CDchoice [time x trials]
        presampcurr = mean(presampcurr,1,'omitnan');    % Average across time [1 x trials]
        presampcurr = mean(presampcurr,'omitnan');      % Average across trials

        modelpred.(cond){sessix} = modelpred.(cond){sessix}-presampcurr;    % Baseline subtract

        curr = trueVals.(cond){sessix};     % Do the same thing for ground truth CDchoice values
        presampcurr = curr(start:stop,:);
        presampcurr = mean(presampcurr,1,'omitnan');
        presampcurr = mean(presampcurr,'omitnan');

        trueVals.(cond){sessix} = trueVals.(cond){sessix}-presampcurr;
    end
end

%% Plot example of motion energy and CDChoice, averaged across right and left trials
% colors = getColors();
%
% exsess = 9;
% nSessions = numel(meta);
% cond2plot = [2 3];
% MEix = find(strcmp(kin(1).featLeg,'motion_energy'));
%
% times.trialstart = median(obj(1).bp.ev.bitStart)-median(obj(1).bp.ev.(params(1).alignEvent));
% times.startix = find(obj(1).time>times.trialstart,1,'first');
% times.samp = median(obj(1).bp.ev.sample)-median(obj(1).bp.ev.(params(1).alignEvent));
% times.stopix = find(obj(1).time<times.samp,1,'last');
%
% cols = {colors.rhit,colors.lhit};
% alph = 0.2;
% numtrix2plot = 30;
% for sessix =  1:numel(meta)
%     figure();
%     plotCondAvgMEandCD(cond2plot, sessix, params, obj, meta, kin, trueVals, alph, cols, MEix,times,par,'invert')
%
%     figure();
%     plotExampleChoiceSelME(numtrix2plot, cond2plot, sessix, params, obj, meta, kin, MEix,times)
% end
%% Make heatmaps for a single session showing CDchoice across trials and predicted CDchoice
sm = 20;        % How much to smooth single trial projections by

load('C:\Code\Uninstructed-Movements\LeftRightDiverging_Colormap.mat')

% Times that you want to use to sort CDChoice
del = mode(obj(1).bp.ev.delay)-mode(obj(1).bp.ev.(params(1).alignEvent));
start = find(obj(1).time>del,1,'first');
resp = mode(obj(1).bp.ev.goCue)-mode(obj(1).bp.ev.(params(1).alignEvent))-0.05;
stop = find(obj(1).time<resp,1,'last');

sess2use = [4,6,19,21];         % Session to use

cond2plot = {'Lhit','Rhit'};    % Which conditions to plot
for sessix = sess2use                                                                 % For each session...
    figure();
    cnt = 0;
    tempTrue = [];
    tempPred = [];
    l1 = size(trueVals.(cond2plot{1}){sessix},2)+0.5;   % Where you will draw horizontal line to delineate right and left
    % Combine the true values for CDchoice and the model predicted
    % values across conditions
    % tempTrue = (time x [num left trials + num right trials])
    for c = 1:length(cond2plot)                                                 % For left and right trials...
        cond = cond2plot{c};
        currTrue = trueVals.(cond){sessix};                                     % Get the true single trial CDchoice projections for that condition and session
        [~,sortix] = sort(mean(currTrue(start:stop,:),1,'omitnan'),'descend');  % Sort the true projections by average magnitude during the delay period
        tempTrue = [tempTrue,currTrue(:,sortix)];
        currPred = modelpred.(cond){sessix};                                    % Get the model predicted single trial CDchoice projections
        tempPred = [tempPred,currPred(:,sortix)];
    end
    nTrials = size(tempTrue,2);                                                 % Total number of trials that are being plotted
    % Plot ground truth CDchoice projections in left subplot
    ax1 = subplot(1,2,1);                                                       % Plot true CDchoice data on left subplot
    imagesc(obj(sessix).time(par.timerange),1:nTrials,-1*tempTrue'); hold on    % Heatmap of true data (sorted left trials will be on top, then a white line, then sorted right trials)
    line([obj(sessix).time(1),obj(sessix).time(end)],[l1,l1],'Color','black','LineStyle','--')

    % Plot CDchoice projections predicted from video in right subplot
    ax2 = subplot(1,2,2);
    imagesc(obj(sessix).time(par.timerange),1:nTrials,mySmooth(-1*tempPred,sm)'); hold on

    line([obj(sessix).time(1),obj(sessix).time(end)],[l1,l1],'Color','black','LineStyle','--')
    title(ax1,'CDchoice - neural data')
    colorbar(ax1)
    clim(ax1,[-4 4])
    colormap(LeftRightDiverging_Colormap)
    xlabel(ax1,'Time from go cue (s)')
    xline(ax1,0,'k--','LineWidth',1)
    xline(ax1,-0.9,'k--','LineWidth',1)
    xline(ax1,-2.2,'k--','LineWidth',1)
    xlim(ax1,[-2.5 0])
    set(gca,'TickDir','out');

    title(ax2,'Video prediction')
    xlabel(ax2,'Time from go cue (s)')
    xline(ax2,0,'k--','LineWidth',1)
    xline(ax2,-0.9,'k--','LineWidth',1)
    xline(ax2,-2.2,'k--','LineWidth',1)
    colorbar(ax2)
    colormap(LeftRightDiverging_Colormap)
    xlim(ax2,[-2.5 0])
    clim(ax2,[-1.5 1.5])
    set(gca,'TickDir','out');

    sgtitle(['Example session:  ' meta(sessix).anm ' ' meta(sessix).date])
end

%% Example plots by session for relating predicted and true CDchoice
del = mode(obj(1).bp.ev.delay)-mode(obj(1).bp.ev.(params(1).alignEvent));
start = find(obj(1).time>del,1,'first');
resp = mode(obj(1).bp.ev.goCue)-mode(obj(1).bp.ev.(params(1).alignEvent))-0.05;
stop = find(obj(1).time<resp,1,'last');

delR2_ALL = [];

plotexample = 'no';

if strcmp(plotexample,'yes')
    plotrange = exsess;
else
    plotrange = 1:length(meta);
end


for sessix = plotrange
    %%% Plot a scatter plot for a single session of true CDlate and predicted CDlate for each trial
    %%% Each dot = an average value of CDlate during the delay period
    figure();
    subplot(1,2,2)
    tempR2 = Scatter_ModelPred_TrueCDchoice(trueVals, modelpred, sessix, start, stop,meta,'invert');
    % Save R2 value for that session
    delR2_ALL = [delR2_ALL, tempR2];

    % Calculate averages and standard deviation for true CD and predicted CD  this session
    [avgCD,stdCD] = getAvgStd(trueVals,modelpred,sessix);

    colors = getColors();
    alph  = 0.2;

    %%% Plot an example session of CDlate prediction vs true value
    subplot(1,2,1)
    plotExampleCDchoice_Pred(colors, obj, par, meta, avgCD, stdCD, sessix, trueVals,alph, tempR2,'invert');

    %     if strcmp(plotexample,'no')
    %         close all
    %     end
end
%% Plot bar plot to show average R2 values across sessions
exsess = 3;

colors = getColors();
delR2_ALL = abs(delR2_ALL);

nSessions =length(meta);

% The index of the session that you want to be highlighted
markerSize = 60;
figure();
b = bar(mean(delR2_ALL),'FaceColor',colors.afc); hold on;                   % Plot the average R2 value across all sessions
ix2plot = 1:nSessions;
ix2plot(exsess) = [];
scatter(ones(nSessions-1,1),delR2_ALL(ix2plot),markerSize,'filled','o','MarkerFaceColor',...,
    'k','XJitter','randn','XJitterWidth',0.25); hold on;
scatter(1,delR2_ALL(exsess),markerSize,'o','MarkerEdgeColor',colors.afc)
errorbar(b.XEndPoints,mean(delR2_ALL,'omitnan'),std(delR2_ALL,'omitnan'),'LineStyle','none','Color','k','LineWidth',1)
ylim([0 0.8])
set(gca,'TickDir','out');
ax = gca;
ax.FontSize = 16;
title(['CDChoice: Ex session = ' meta(exsess).anm meta(exsess).date])
%% Print summary statistics
disp('---Summary statistics for CDChoice prediction---')
disp(['Average R2 across all sessions (n = ' num2str(length(meta)) ' ) = ' num2str(mean(delR2_ALL))])
disp(['Standard deviation across all sessions = ' num2str(std(delR2_ALL))])
if par.regularize==1
    regtype = 'Ridge';
else
    regtype = 'Linear; No regularization';
end
disp(['Regularization type: ' regtype])
disp(['Train percentage: ' num2str(100*par.train) ' %'])
t = datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss Z');
disp(t)

%% FUNCTIONS
function [avgCD,stdCD] = getAvgStd(trueVals,modelpred,sessix)

avgCD.Rhit.true = mean(mySmooth(trueVals.Rhit{sessix},81),2,'omitnan');      % Get average true CDlate for R and L hits for this session
avgCD.Lhit.true = mean(mySmooth(trueVals.Lhit{sessix},81),2,'omitnan');
stdCD.Rhit.true = std(mySmooth(trueVals.Rhit{sessix},81),0,2,'omitnan');     % Get standard deviation of true CDlate for R and L hits
stdCD.Lhit.true = std(mySmooth(trueVals.Lhit{sessix},81),0,2,'omitnan');

modelpred.Rhit{sessix} = fillmissing(modelpred.Rhit{sessix},"nearest");
modelpred.Lhit{sessix} = fillmissing(modelpred.Lhit{sessix},"nearest");
infix = find(isinf(modelpred.Rhit{sessix})); modelpred.Rhit{sessix}(infix) = 0;
infix = find(isinf(modelpred.Lhit{sessix})); modelpred.Lhit{sessix}(infix) = 0;
avgCD.Rhit.pred = mean(mySmooth(modelpred.Rhit{sessix},81),2,'omitnan');     % Get average predicted CDlate for R and L hits for this session
avgCD.Lhit.pred = mean(mySmooth(modelpred.Lhit{sessix},81),2,'omitnan');
stdCD.Rhit.pred = std(mySmooth(modelpred.Rhit{sessix},81),0,2,'omitnan');    % Get stdev of predicted CDlate for R and L hits for this session
stdCD.Lhit.pred = std(mySmooth(modelpred.Lhit{sessix},81),0,2,'omitnan');
end