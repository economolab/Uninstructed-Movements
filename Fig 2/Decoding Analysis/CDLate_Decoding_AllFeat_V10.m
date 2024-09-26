% DECODING CDlate FROM ALL KINEMATIC FEATURES (with ridge regression,
% regularization; train/test split, cross-validation)
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
addpath(genpath(fullfile(utilspth,'fig1')));
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
params.condition(1)     = {'(hit|miss|no)'};                             % (1) all trials
params.condition(end+1) = {'R&hit&~stim.enable&~autowater&~early'};             % (2) right hits, no stim, aw off
params.condition(end+1) = {'L&hit&~stim.enable&~autowater&~early'};             % (3) left hits, no stim, aw off
params.condition(end+1) = {'R&miss&~stim.enable&~autowater&~early'};            % (4) error right, no stim, aw off
params.condition(end+1) = {'L&miss&~stim.enable&~autowater&~early'};            % (5) error left, no stim, aw off
params.condition(end+1) = {'R&no&~stim.enable&~autowater&~early'};              % (6) no right, no stim, aw off
params.condition(end+1) = {'L&no&~stim.enable&~autowater&~early'};              % (7) no left, no stim, aw off
params.condition(end+1) = {'hit&~stim.enable&~autowater&~early'};               % (8) all hits, no stim, aw off


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

disp('----Calculating coding dimensions----')
cond2use = [2 3]; % right hits, left hits (corresponding to PARAMS.CONDITION)
inclramp = 'yes';
rampcond = 8;
cond2proj = 2:7;  % right hits, left hits, right miss, left miss, right no, left no (corresponding to null/potent psths in rez)
cond2use_trialdat = [2 3]; % for calculating selectivity explained in full neural pop
regr = getCodingDimensions_2afc(obj,params,cond2use,cond2proj);

disp('----Projecting single trials onto CDlate----')
cd = 'late';
smooth = 60;
regr = getSingleTrialProjs(regr,obj,cd,smooth);
%% Load kinematic data
nSessions = numel(meta);
for sessix = 1:numel(meta)
    message = strcat('----Getting kinematic data for session',{' '},num2str(sessix), {' '},'out of',{' '},num2str(nSessions),'----');
    disp(message)
    kin(sessix) = getKinematics(obj(sessix), me(sessix), params(sessix));
end
%% Sanity check -- why do averages of single trial projections onto CDChoice look so much smaller than the projection of the cond avg PSTHs?
% figure()
% for sessix = 1:length(meta)
%     subplot(1,3,1)
%     plot(obj(sessix).time,regr(sessix).cd_proj(:,1),'r'); hold on;
%     plot(obj(sessix).time,regr(sessix).cd_proj(:,2),'b'); hold off;
% 
%     subplot(1,3,2)
%     temptrix = params(sessix).trialid{2};
%     plot(obj(sessix).time,regr(sessix).singleProj(:,temptrix),'r'); hold on;
% 
%     temptrix = params(sessix).trialid{3};
%     plot(obj(sessix).time,regr(sessix).singleProj(:,temptrix),'b'); hold off;
% 
%     subplot(1,3,3)
%     temptrix = params(sessix).trialid{2};
%     plot(obj(sessix).time,mean(regr(sessix).singleProj(:,temptrix),2,'omitnan'),'r'); hold on;
% 
%     temptrix = params(sessix).trialid{3};
%     plot(obj(sessix).time,mean(regr(sessix).singleProj(:,temptrix),2,'omitnan'),'b'); hold off;
% 
%     sgtitle([meta(sessix).anm ' ; ' meta(sessix).date])
%     pause
% end
%% Predict CDTrialType from DLC features
%%% DECODING PARAMETERS %%%

% input data = neural data (time*trials,neurons)
% output data = kin data   (time*trials,kin feats)

par.pre=6; % time bins prior to output used for decoding
par.post=0; % time bins after output used for decoding
par.dt = params(1).dt; % moving time bin
par.pre_s = par.pre .* params(1).dt; % dt, pre_s, and post_s just here to know how much time you're using. Set params.dt and pre/post appropriately for you analysis
par.post_s = par.post .* params(1).dt;

trialstart = mode(obj(1).bp.ev.bitStart)-mode(obj(1).bp.ev.(params(1).alignEvent)+0.4);
start = find(obj(1).time>trialstart,1,'first');
go = mode(obj(1).bp.ev.goCue)-mode(obj(1).bp.ev.(params(1).alignEvent));
stop = find(obj(1).time>go,1,'first');

par.timerange = start:stop;

% these parameters above are important for decoding accuracy. for actual
% analyses (like a final analysis to be included in a publication),
% you should vary these parameters only if you have a validation
% set that you won't test on until these parameters are set. Otherwise,
% there's a high risk of overfitting

% cross val folds
par.nFolds = 4;

% data sets
par.train = 0.6; % fraction of trials 
par.test = 1 - par.train;

% feature to use to decode
% par.feats = kin(1).featLeg;
par.feats = {'motion','nos','jaw'};
temp = cellfun(@(x) patternMatchCellArray(kin(1).featLeg,{x},'all') , par.feats,'UniformOutput',false);
par.feats = cat(1, temp{:});

% trials
par.cond2use = [2 3]; % R and L hits, ~autowater

par.regularize = 1; % if 0, linear regression. if 1, ridge regression

par.lambdas = logspace(-3,3,20); % regularization parameters to search over. don't need to change unless you want to make it more fine-grained.

%% DECODING

close all

for sessix = 1:numel(meta)
    disp(['Decoding for session ' ' ' num2str(sessix) ' / ' num2str(numel(meta))])
    
    [X,Y,par] = preparePredictorsRegressors_v2(par, sessix, kin, regr,params);

    for ilambda = 1:numel(par.lambdas)
        lambda = par.lambdas(ilambda);
        cvmdl{ilambda} = fitrlinear(X.train,Y.train,'Learner','svm','KFold',par.nFolds,'Regularization','ridge','Lambda',lambda);
        loss(ilambda) = cvmdl{ilambda}.kfoldLoss;
        % cvpred{ilambda} = kfoldPredict(cvmdl{ilambda});
    end
    
    [~,bestmodel] = min(loss);
    mdl = cvmdl{bestmodel}; % we now have the best lambda, and they trained cvmodel with that lambda,
    % can predict test data using best cv mdl
    for i = 1:par.nFolds
        mdl_loss(i) = mdl.Trained{i}.loss(X.train,Y.train);
    end
    [~,bestmodel] = min(mdl_loss);
    testmdl = mdl.Trained{bestmodel}; % we now have the best lambda, and the trained model with that lambda,
    
    pred = testmdl.predict(X.test);   % use left-out test data for the prediction
    
    % Save the predictor coefficients from the testmdl
    loadings(:,sessix) = testmdl.Beta;  % Average the coefficients for each predictor term across folds; save these for each time point


    y = reshape(Y.test,Y.size.test(1),Y.size.test(2)); % original input data (standardized)
    yhat = reshape(pred,Y.size.test(1),Y.size.test(2)); % prediction
   


%     figure()
%     subplot(1,2,1); imagesc(y'); colorbar; subplot(1,2,2); imagesc(yhat');colorbar()
%     figure()
%     plot(yhat)
%     pause

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
%% Check how many test trials each session has 
nTestTrix = NaN(length(meta),2);
for sessix = 1:length(meta)
    nTestTrix(sessix,:)=[size(trueVals.Rhit{sessix},2),size(trueVals.Lhit{sessix},2)];
end
%% Baseline subtract CDTrialType
% Times that you want to use to baseline normalize CDTrialType
trialstart = mode(obj(1).bp.ev.bitStart)-mode(obj(1).bp.ev.(params(1).alignEvent));
start = find(obj(1).time>trialstart,1,'first');
samp = mode(obj(1).bp.ev.sample)-mode(obj(1).bp.ev.(params(1).alignEvent));
stop = find(obj(1).time<samp,1,'last');

cond2use = {'Rhit','Lhit'};
for sessix = 1:length(meta)
    for c = 1:length(cond2use)
        cond = cond2use{c};
        curr = modelpred.(cond){sessix};
        presampcurr = curr(start:stop,:);
        presampcurr = mean(presampcurr,1,'omitnan');
        presampcurr = mean(presampcurr,'omitnan');

        modelpred.(cond){sessix} = modelpred.(cond){sessix}-presampcurr;

        curr = trueVals.(cond){sessix};
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
%% Make heatmaps for a single session showing CDTrialType across trials and predicted CDTrialType
plotheatmap = 'yes';
sm = 20;
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
    for sessix = 1:length(meta)                                                                  % For each session...
        figure();
        cnt = 0;
        tempTrue = [];
        tempPred = [];
        l1 = size(trueVals.(cond2plot{1}){sessix},2)+0.5;
        % Combine the true values for CDTrialType and the model predicted
        % values across conditions
        % tempTrue = (time x [num left trials + num right trials])
        for c = 1:length(cond2plot)                                                 % For left and right trials...
            cond = cond2plot{c};
            currTrue = trueVals.(cond){sessix};                                     % Get the true single trial CDTrialType projections for that condition and session
            [~,sortix] = sort(mean(currTrue(start:stop,:),1,'omitnan'),'descend');  % Sort the true projections by average magnitude during the delay period
            tempTrue = [tempTrue,currTrue(:,sortix)];
            currPred = modelpred.(cond){sessix};                                    % Get the model predicted single trial CDTrialType projections
            tempPred = [tempPred,currPred(:,sortix)];
        end
        nTrials = size(tempTrue,2);                                                 % Total number of trials that are being plotted
        ax1 = subplot(1,2,1);                                                       % Plot true CDTrialType data on left subplot
        if strcmp(invertCD,'invert')
            imagesc(obj(sessix).time(par.timerange),1:nTrials,-1*tempTrue'); hold on                      % Heatmap of true data (sorted left trials will be on top, then a white line, then sorted right trials)
        else
            imagesc(obj(sessix).time(par.timerange),1:nTrials,tempTrue'); hold on;
        end
        line([obj(sessix).time(1),obj(sessix).time(end)],[l1,l1],'Color','black','LineStyle','--')
        

        ax2 = subplot(1,2,2);
        if strcmp(invertCD,'invert')
            imagesc(obj(sessix).time(par.timerange),1:nTrials,mySmooth(-1*tempPred,sm)'); hold on
        else
            imagesc(obj(sessix).time(par.timerange),1:nTrials,mySmooth(tempPred,sm)'); hold on
        end
% 
%         if strcmp(invertCD,'invert')
%             imagesc(obj(sessix).time(par.timerange),1:nTrials,-1*tempPred'); hold on
%         else
%             imagesc(obj(sessix).time(par.timerange),1:nTrials,tempPred'); hold on
%         end
        line([obj(sessix).time(1),obj(sessix).time(end)],[l1,l1],'Color','black','LineStyle','--')
        title(ax1,'CDTrialType - data')
        colorbar(ax1)
        clim(ax1,[-4 4])
        colormap(LeftRightDiverging_Colormap)
        xlabel(ax1,'Time from go cue (s)')
        xline(ax1,0,'k--','LineWidth',1)
        xline(ax1,-0.9,'k--','LineWidth',1)
        xline(ax1,-2.2,'k--','LineWidth',1)
        xlim(ax1,[-2.5 0])
        set(gca,'TickDir','out');
        
        title(ax2,'Model prediction')
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
end
%% Example plots by session for relating predicted and true CDTrialType
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
    tempR2 = Scatter_ModelPred_TrueCDTrialType(trueVals, modelpred, sessix, start, stop,meta,'invert');
    % Save R2 value for that session
    delR2_ALL = [delR2_ALL, tempR2];

    % Calculate averages and standard deviation for true CD and predicted CD  this session
    [avgCD,stdCD] = getAvgStd(trueVals,modelpred,sessix);

    colors = getColors();
    alph  = 0.2;

    %%% Plot an example session of CDlate prediction vs true value
    subplot(1,2,1)
    plotExampleCDTrialType_Pred(colors, obj, par, meta, avgCD, stdCD, sessix, trueVals,alph, tempR2,'invert');

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