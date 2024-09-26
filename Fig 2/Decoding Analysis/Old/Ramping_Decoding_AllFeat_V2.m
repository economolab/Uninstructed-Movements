% DECODING CDlate FROM ALL KINEMATIC FEATURES
clear,clc,close all

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
addpath(genpath(fullfile(utilspth,'figNP')));
figpth = [basepth  '\Uninstructed-Movements\Fig 2'];
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
params.condition(end+1) = {'R&miss&~stim.enable&~autowater&~early'};            % error right, no stim, aw off
params.condition(end+1) = {'L&miss&~stim.enable&~autowater&~early'};            % error left, no stim, aw off
params.condition(end+1) = {'R&no&~stim.enable&~autowater&~early'};              % no right, no stim, aw off
params.condition(end+1) = {'L&no&~stim.enable&~autowater&~early'};              % no left, no stim, aw off
params.condition(end+1) = {'hit&~stim.enable&~autowater&~early'};               % all hits, no stim, aw off


params.tmin = -2.5;
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

if strcmp(whichcomp,'LabPC')
    datapth = 'C:\Users\Jackie Birnbaum\Documents\Data';
elseif strcmp(whichcomp,'Laptop')
    datapth = 'C:\Users\Jackie\Documents\Grad School\Economo Lab';
end

meta = [];

% --- ALM ---
meta = loadJEB13_ALMVideo(meta,datapth);
meta = loadJEB6_ALMVideo(meta,datapth);
meta = loadJEB7_ALMVideo(meta,datapth);
meta = loadEKH1_ALMVideo(meta,datapth);
meta = loadEKH3_ALMVideo(meta,datapth);
meta = loadJGR2_ALMVideo(meta,datapth);
meta = loadJGR3_ALMVideo(meta,datapth);
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

disp('----Calculating coding dimensions----')
cond2use = [2 3 6 7]; % right hits, left hits (corresponding to PARAMS.CONDITION)
inclramp = 'yes';
rampcond = 8;
cond2proj = 2:7;  % right hits, left hits, right miss, left miss, right no, left no (corresponding to null/potent psths in rez)
cond2use_trialdat = [2 3]; % for calculating selectivity explained in full neural pop
regr = getCodingDimensions_2afc(obj,params,cond2use,rampcond,cond2proj, inclramp);

disp('----Projecting single trials onto CDlate----')
cd = 'ramping';
regr = getSingleTrialProjs(regr,obj,cd);
%% Load kinematic data
nSessions = numel(meta);
for sessix = 1:numel(meta)
    message = strcat('----Getting kinematic data for session',{' '},num2str(sessix), {' '},'out of',{' '},num2str(nSessions),'----');
    disp(message)
    kin(sessix) = getKinematics(obj(sessix), me(sessix), params(sessix));
end
%% Predict CDRamping from DLC features
clearvars -except datapth kin me meta obj params regr nSessions exsess
clc;

% input data = neural data (time*trials,neurons)
% output data = kin data   (time*trials,kin feats)
par.pre=5; % time bins prior to output used for decoding
par.post=0; % time bins after output used for decoding
par.dt = params(1).dt; % moving time bin
par.pre_s = par.pre .* params(1).dt; % dt, pre_s, and post_s just here to know how much time you're using. Set params.dt and pre/post appropriately for you analysis
par.post_s = par.post .* params(1).dt;


trialstart = mode(obj(1).bp.ev.bitStart)-mode(obj(1).bp.ev.(params(1).alignEvent)+0.4);
start = find(obj(1).time>trialstart,1,'first');
go = mode(obj(1).bp.ev.goCue)-mode(obj(1).bp.ev.(params(1).alignEvent));
stop = find(obj(1).time>go,1,'first');

par.timerange = start:stop;

% cross val folds
par.nFolds = 4;

% data sets
par.train = 1; % fraction of trials (just using cross-val here, matlab's kFoldPredict uses held out data for testing)
par.test = 1 - par.train;

% feature to use to decode
% par.feats = kin(1).featLeg;
par.feats = {'motion','nos','jaw'};
temp = cellfun(@(x) patternMatchCellArray(kin(1).featLeg,{x},'all') , par.feats,'UniformOutput',false);
par.feats = cat(1, temp{:});
% par.feats = {'tongue_angle','tongue_length','motion_energy'};

% trials
par.cond2use = [2 3];

par.regularize = 1; % if 0, linear regression. if 1, ridge regression
%% DECODING

close all

for sessix = 1:numel(meta)
    disp(['Decoding for session ' ' ' num2str(sessix) ' / ' num2str(numel(meta))])

    [X,Y] = preparePredictorsRegressors(par, sessix, kin, regr,params);

    if par.regularize
        mdl = fitrlinear(X.train,Y.train,'Learner','leastsquares','KFold',par.nFolds,'Regularization','ridge');
    else
        mdl = fitrlinear(X.train,Y.train,'Learner','leastsquares','KFold',par.nFolds);
    end

    % Save the predictor coefficients for this point in time (averaged
    % across fold iterations)
    nPredictors = size(X.train,2);                  % Number of kinematic predictors * binWidth
    loadings = NaN(nPredictors,par.nFolds);         % (# predictors x folds for CV)
    for fol = 1:par.nFolds
        loadings(:,fol) = mdl.Trained{fol}.Beta;
    end
    avgloadings(:,sessix) = mean(loadings,2,'omitnan');  % Average the coefficients for each predictor term across folds; save these for each time point

    pred = kfoldPredict(mdl);

    y = reshape(Y.train,Y.size(1),Y.size(2)); % original input data (standardized)
    yhat = reshape(pred,Y.size(1),Y.size(2)); % prediction

    %     figure()
    %     subplot(1,2,1); imagesc(y'); colorbar; subplot(1,2,2); imagesc(yhat');colorbar()
    %     figure()
    %     plot(yhat)
    %     pause

    cnt = 1;
    for c = 1:length(par.cond2use)
        if c==1
            cond = 'Rhit';
        else
            cond = 'Lhit';
        end
        ncondtrix = length(params(sessix).trialid{par.cond2use(c)});
        ixrange = cnt:cnt+ncondtrix-1;
        trueVals.(cond){sessix} = y(:,ixrange);
        modelpred.(cond){sessix} = yhat(:,ixrange);
        cnt = ncondtrix+1;
    end
end

disp('---FINISHED DECODING FOR ALL SESSIONS---')
%% Baseline subtract CDRamping
% Times that you want to use to baseline normalize CDTrialType
trialstart = mode(obj(1).bp.ev.bitStart)-mode(obj(1).bp.ev.(params(1).alignEvent));
start = find(obj(1).time>trialstart,1,'first');
samp = mode(obj(1).bp.ev.sample)-mode(obj(1).bp.ev.(params(1).alignEvent));
stop = find(obj(1).time<samp,1,'last');

cond2use = 'Rhit';
for sessix = 1:length(meta)
    curr = modelpred.(cond2use){sessix};
    presampcurr = curr(start:stop,:);
    presampcurr = mean(presampcurr,1,'omitnan');
    presampcurr = mean(presampcurr,'omitnan');

    modelpred.Rhit{sessix} = modelpred.Rhit{sessix}-presampcurr;
    modelpred.Lhit{sessix} = modelpred.Lhit{sessix}-presampcurr;

    curr = trueVals.(cond2use){sessix};
    presampcurr = curr(start:stop,:);
    presampcurr = mean(presampcurr,1,'omitnan');
    presampcurr = mean(presampcurr,'omitnan');

    trueVals.Rhit{sessix} = trueVals.Rhit{sessix}-presampcurr;
    trueVals.Lhit{sessix} = trueVals.Lhit{sessix}-presampcurr;
end

%% Plot example of motion energy and ramping mode, averaged across right and left trials
exsess = 19;
colors = getColors();

nSessions = numel(meta);
cond2plot = [2 3];
MEix = find(strcmp(kin(1).featLeg,'motion_energy'));

times.trialstart = median(obj(1).bp.ev.bitStart)-median(obj(1).bp.ev.(params(1).alignEvent));
times.startix = find(obj(1).time>times.trialstart,1,'first');
times.samp = median(obj(1).bp.ev.sample)-median(obj(1).bp.ev.(params(1).alignEvent));
times.stopix = find(obj(1).time<times.samp,1,'last');

cols = {colors.rhit,colors.lhit};
alph = 0.2;
sm = 21;
for sessix = exsess %1:numel(meta)
    condME = [];
    condRamp = [];
    for c = 1:length(cond2plot)
        cond = cond2plot(c);
        condtrix = params(sessix).trialid{cond};
        temp = squeeze(kin(sessix).dat(:,condtrix,MEix));
        condME = [condME,temp];

         if c==1
            dir = 'Rhit';
        else
            dir = 'Lhit';
        end
        temp = trueVals.(dir){sessix};
        condRamp = [condRamp,temp];
    end
    % Baseline subtract ME where baseline is presample
    %     presampME = mean(condME(times.startix:times.stopix,:),1,'omitnan');
    %     presampME = mean(presampME);
    %     condME = condME-presampME;
    % Baseline subtract ME where baseline is 1st percentile
    pctME = prctile(condME,1);
    condME = condME-pctME;
    nTrials = size(condME,2);

    ax1 = subplot(2,1,1);
    toplot = mean(mySmooth(condME,31),2,'omitnan');
    err = 1.96*(std(mySmooth(condME,31),0,2,'omitnan')./sqrt(nTrials));
    ax = gca;
    shadedErrorBar(obj(sessix).time,mySmooth(toplot,sm),err,{'Color',[0.2 0.2 0.2],'LineWidth',2},alph,ax); hold on;

    ax2 = subplot(2,1,2);
    toplot = mean(mySmooth(condRamp(par.timerange,:),21),2,'omitnan');
    err = 1.96*(std(mySmooth(condRamp(par.timerange,:),21),0,2,'omitnan')./sqrt(nTrials));
    ax = gca;
    shadedErrorBar(obj(sessix).time(par.timerange),toplot,err,{'Color',[0.2 0.2 0.2],'LineWidth',2},alph,ax); hold on;
end
xlabel(ax1,'Time from go cue (s)')
ylabel(ax1,'Motion energy (a.u.)')
xline(ax1,0,'k--','LineWidth',1)
xline(ax1,-0.9,'k--','LineWidth',1)
xline(ax1,-2.2,'k--','LineWidth',1)
xlim(ax1,[-2.3 0])
set(ax1,'TickDir','out'); 

xlabel(ax2,'Time from go cue (s)')
ylabel(ax2,'CDRamping (a.u.)')
xline(ax2,0,'k--','LineWidth',1)
xline(ax2,-0.9,'k--','LineWidth',1)
xline(ax2,-2.2,'k--','LineWidth',1)
xlim(ax2,[-2.3 0])
set(ax2,'TickDir','out'); 
%% Make heatmaps for a single session showing CDTrialType across trials and predicted CDTrialType
plotheatmap = 'yes';
sm = 15;
invertCD = 'invert';                    % 'Invert' or 'no' for whether or not you want to flip the sign of the CD projection

load('C:\Code\Uninstructed-Movements\LeftRightDiverging_Colormap.mat')

if strcmp(plotheatmap,'yes')
    % Times that you want to use to sort CDTrialType
    del = mode(obj(1).bp.ev.delay)-mode(obj(1).bp.ev.(params(1).alignEvent));
    start = find(obj(1).time>del,1,'first');
    resp = mode(obj(1).bp.ev.goCue)-mode(obj(1).bp.ev.(params(1).alignEvent))-0.05;
    stop = find(obj(1).time<resp,1,'last');

    goodsess = [4,6,19,21];

    cond2plot = {'Lhit','Rhit'};
    for sessix = exsess                                                                  % For each session...
        figure();
        cnt = 0;
        tempTrue = [];
        tempPred = [];
        % Combine the true values for CDRamping and the model predicted
        % values across conditions
        % tempTrue = (time x [num left trials + num right trials])
        for c = 1:length(cond2plot)                                                 % For left and right trials...
            cond = cond2plot{c};
            currTrue = trueVals.(cond){sessix};                                     % Get the true single trial CDRamping projections for that condition and session
            tempTrue = [tempTrue,currTrue];
            currPred = modelpred.(cond){sessix};                                    % Get the model predicted single trial CDRamping projections
            tempPred = [tempPred,currPred];
        end
        [~,sortix] = sort(mean(tempTrue(start:stop,:),1,'omitnan'),'descend');      % Sort the true projections by average magnitude during the delay period
        True2plot = mySmooth(tempTrue(:,sortix),31);
        Pred2plot = tempPred(:,sortix);                                             % Sort the model predictions in the same order

        nTrials = size(True2plot,2);                                                 % Total number of trials that are being plotted
        ax1 = subplot(1,2,1);                                                       % Plot true CDRamping data on left subplot
        imagesc(obj(sessix).time(par.timerange),1:nTrials,True2plot'); hold on                      % Heatmap of true data (sorted left trials will be on top, then a white line, then sorted right trials)

        ax2 = subplot(1,2,2);
        imagesc(obj(sessix).time(par.timerange),1:nTrials,mySmooth(Pred2plot,21)'); hold on
        title(ax1,'CDRamping - data')
        colorbar(ax1)
        colormap(LeftRightDiverging_Colormap)
        clim(ax1,[-3 3])
        xlabel(ax1,'Time from go cue (s)')
        xline(ax1,0,'k--','LineWidth',1)
        xline(ax1,-0.9,'k--','LineWidth',1)
        xline(ax1,-2.2,'k--','LineWidth',1)
        xlim(ax1,[-2.5 0])

        title(ax2,'Model prediction')
        xlabel(ax2,'Time from go cue (s)')
        xline(ax2,0,'k--','LineWidth',1)
        xline(ax2,-0.9,'k--','LineWidth',1)
        xline(ax2,-2.2,'k--','LineWidth',1)
        xlim(ax2,[-2.5 0])

        colorbar(ax2)
        colormap(LeftRightDiverging_Colormap)
        clim(ax2,[-2.5 2.5])

        sgtitle(['Example session:  ' meta(sessix).anm ' ' meta(sessix).date])
    end
end
%% Example plots by session for relating predicted and true CDTrialType
del = mode(obj(1).bp.ev.delay)-mode(obj(1).bp.ev.(params(1).alignEvent));
start = find(obj(1).time>del,1,'first');
resp = mode(obj(1).bp.ev.goCue)-mode(obj(1).bp.ev.(params(1).alignEvent))-0.05;
stop = find(obj(1).time<resp,1,'last');

delR2_ALL = [];

plotexample = 'yes';

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
    tempR2 = Scatter_ModelPred_TrueCDTrialType(trueVals, modelpred, sessix, start, stop,meta,'no');
    % Save R2 value for that session
    delR2_ALL = [delR2_ALL, tempR2];

    % Calculate averages and standard deviation for true CD and predicted CD  this session
    [avgCD,stdCD] = getAvgStd(trueVals,modelpred,sessix);

    colors = getColors();
    alph  = 0.2;

    %%% Plot an example session of CDlate prediction vs true value
    subplot(1,2,1)
    plotExampleCDRamp_Pred(colors, obj, par, meta, avgCD, stdCD, sessix, trueVals,alph, tempR2,'no');

    %     if strcmp(plotexample,'no')
    %         close all
    %     end
end
%% Plot bar plot to show average R2 values across sessions
colors = getColors();
delR2_ALL = abs(delR2_ALL);
nSessions = length(meta);

% The index of the session that you want to be highlighted
markerSize = 60;
figure();
bar(mean(delR2_ALL),'FaceColor',colors.afc); hold on;                   % Plot the average R2 value across all sessions
for sessix = 1:nSessions
 
    scatter(1,delR2_ALL(sessix),markerSize,'filled','o','MarkerFaceColor',[0.65 0.65 0.65]); hold on;
end
scatter(1,delR2_ALL(exsess),markerSize,'filled','o','cyan','MarkerEdgeColor','black')
ylim([0 1])
ax = gca;
ax.FontSize = 16;
title(['Ex session = ' meta(exsess).anm meta(exsess).date])
%% Print summary statistics 
disp('---Summary statistics for CDRamping prediction---')
disp(['Average R2 across all sessions (n = ' num2str(length(meta)) ' ) = ' num2str(mean(delR2_ALL))])
disp(['Standard deviation across all sessions = ' num2str(std(delR2_ALL))])
if par.regularize==1
    regtype = 'Ridge';
else
    regtype = 'none';
end
disp(['Regularization type: ' regtype])
t = datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss Z');
disp(t)
%% FUNCTIONS
function [avgCD,stdCD] = getAvgStd(trueVals,modelpred,sessix)
temp = [trueVals.Rhit{sessix},trueVals.Lhit{sessix}];
avgCD.allhits.true = mean(mySmooth(temp,81),2,'omitnan');
stdCD.allhits.true = mean(mySmooth(temp,81),2,'omitnan');
avgCD.allhits.nTrials = size(temp,2);

avgCD.Rhit.true = mean(mySmooth(trueVals.Rhit{sessix},81),2,'omitnan');      % Get average true CDlate for R and L hits for this session
avgCD.Lhit.true = mean(mySmooth(trueVals.Lhit{sessix},81),2,'omitnan');
stdCD.Rhit.true = std(mySmooth(trueVals.Rhit{sessix},81),0,2,'omitnan');     % Get standard deviation of true CDlate for R and L hits
stdCD.Lhit.true = std(mySmooth(trueVals.Lhit{sessix},81),0,2,'omitnan');

modelpred.Rhit{sessix} = fillmissing(modelpred.Rhit{sessix},"nearest");
modelpred.Lhit{sessix} = fillmissing(modelpred.Lhit{sessix},"nearest");
infix = find(isinf(modelpred.Rhit{sessix})); modelpred.Rhit{sessix}(infix) = 0;
infix = find(isinf(modelpred.Lhit{sessix})); modelpred.Lhit{sessix}(infix) = 0;

temp = [modelpred.Rhit{sessix},modelpred.Lhit{sessix}];
avgCD.allhits.pred = mean(mySmooth(temp,81),2,'omitnan');
stdCD.allhits.pred = mean(mySmooth(temp,81),2,'omitnan');

avgCD.Rhit.pred = mean(mySmooth(modelpred.Rhit{sessix},81),2,'omitnan');     % Get average predicted CDlate for R and L hits for this session
avgCD.Lhit.pred = mean(mySmooth(modelpred.Lhit{sessix},81),2,'omitnan');
stdCD.Rhit.pred = std(mySmooth(modelpred.Rhit{sessix},81),0,2,'omitnan');    % Get stdev of predicted CDlate for R and L hits for this session
stdCD.Lhit.pred = std(mySmooth(modelpred.Lhit{sessix},81),0,2,'omitnan');
end