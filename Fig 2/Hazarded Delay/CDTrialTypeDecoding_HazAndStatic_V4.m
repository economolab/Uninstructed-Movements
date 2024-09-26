% DECODING CDlate hazarded delay FROM ALL KINEMATIC FEATURES (ridge
% regression; regularization; train/test split)
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
addpath(genpath(fullfile(figpth,'Hazarded Delay')));
figpth = [basepth  '\Uninstructed-Movements\functions'];
addpath(genpath(fullfile(figpth,'HazDel funcs')));
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

params.tmin = -2.7;
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

% Haz delay params (delay lengths that were used in the session)
params.delay(1) = 0.3000;
params.delay(2) = 0.6000;
params.delay(3) = 1.2000;
params.delay(4) = 1.8000;
params.delay(5) = 2.4000;
params.delay(6) = 3.6000;
%% SPECIFY DATA TO LOAD

if strcmp(whichcomp,'LabPC')
    datapth = 'C:\Users\Jackie Birnbaum\Documents\Data';
elseif strcmp(whichcomp,'Laptop')
    datapth = 'C:\Users\Jackie\Documents\Grad School\Economo Lab';
end

meta = [];

% --- ALM ---
meta = loadJEB11_ALMVideo(meta,datapth);
meta = loadJEB12_ALMVideo(meta,datapth);
meta = loadJEB23_ALMVideo(meta,datapth);
meta = loadJEB24_ALMVideo(meta,datapth);

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
%% Load kinematic data
nSessions = numel(meta);
for sessix = 1:numel(meta)
    message = strcat('----Getting kinematic data for session',{' '},num2str(sessix), {' '},'out of',{' '},num2str(nSessions),'----');
    disp(message)
    kin(sessix) = getKinematics(obj(sessix), me(sessix), params(sessix));
end
%% For each session--Trials separated by delay length
% Get the trialIDs corresponding to each delay length
% Find the PSTH for R and L trials of each delay length
% Find the jaw velocities for R and L trials of each delay length
conditions = [2,3];
condfns = {'Rhit','Lhit'};
for sessix = 1:length(meta)
    del(sessix).delaylen = obj(sessix).bp.ev.goCue - obj(sessix).bp.ev.delay;       % Find the delay length for all trials
    del(sessix).del_trialid = getDelayTrix(params(sessix),conditions,del(sessix));  % Group the trials in each condition based on their delay length
    del(sessix).delPSTH = getPSTHbyDel(params(sessix),del(sessix),obj(sessix), condfns, conditions);             % Get avg PSTH for each delay length
end
%% Find CDLate the way that Inagaki et al., 2019, Nature does it (slight variation)
% cd late mode found using delay lengths of 0.3, 0.6, 1.2, and 1.8
dels4mode = [1,2,4];        % Relative to params.delay
conds4mode = [1,2];         % Relative to condfns
dels4proj = [2,3,4];        % Relative to params.delay
conds4proj = [1,2];         % Relative to condfns
cd_labels = 'late';
cd_epochs = 'goCue';
cd_times = [-0.6 -0.02]; % in seconds, relative to respective epochs

sm = 81;
for sessix = 1:length(meta)
    psth2use = [];
    % Get cond avg psth for the delay lengths that you want to use to calculate mode
    for c = 1:length(condfns)
        condpsth = del(sessix).delPSTH.trialPSTH.(condfns{c});
        temp = [];
        for dd = 1:length(dels4mode)
            currdel = dels4mode(dd);
            temp = cat(3,temp,condpsth{currdel});
        end
        avgcondps = mean(temp,3,'omitnan');
        psth2use = cat(3,psth2use,avgcondps);
    end
    del(sessix).cdlate_mode = calcCD_Haz(psth2use,cd_times,cd_epochs,cd_labels,del(sessix),obj(sessix),params(sessix));

    %%%% Project single trials onto the coding dimension that you
    %%%% specified
    nTrials = size(obj(sessix).trialdat,3);
    TrixProj = NaN(length(obj(sessix).time),nTrials);   % time x nTrials
    mode = del(sessix).cdlate_mode;
    for trix = 1:nTrials                                % For each trial...
        temp = obj(sessix).trialdat(:,:,trix);          % Get the PSTH for all cells on that trial
        TrixProj(:,trix) = mySmooth((temp*mode),sm);    % Project the trial PSTH onto the mode that you specified
    end
    del(sessix).singleProj = TrixProj;

    %%%% Condition-averaged projections onto coding dimension
    condproj = cell(1,length(dels4proj));
    temp = NaN(length(obj(sessix).time),length(conds4proj));
    for dd = 1:length(dels4proj)
        currdel = dels4proj(dd);
        for c = 1:length(conds4proj)
            currcond = conds4proj(c);
            condpsth = del(sessix).delPSTH.(condfns{currcond});
            currpsth = condpsth{currdel};
            temp(:,c) = mySmooth((currpsth*mode),sm);
        end
        condproj{dd} = temp;
    end
    del(sessix).condProj = condproj;
end
clearvars -except obj meta params me sav kin del condfns
% %% Sanity check for single trial projections onto CDs (can see on single trials that 
% % there is R and L selectivity)
% for sessix = 1:19
%     alltemp = [];
%     for ii = [2 3]
%         condtrix = params(sessix).trialid{ii};
%         temp1 = del(sessix).singleProj(:,condtrix);
%         alltemp = [alltemp, temp1];
%         plot(obj(1).time(1:271),mean(temp1(1:271,:),2,'omitnan')); hold on;
%     end
% %     imagesc(alltemp')
% %     colorbar
%     xlim([-2.5 0])
%     hold off
%     pause
% end
%% Predict CDTrialType from DLC features
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
stop = find(obj(1).time<go,1,'last');

par.timerange = start:stop;             % Trial start until go cue 
% par.timerange = 1:stop;               % Beginning of time-axis until go cue

par.timeaxis = obj(1).time(par.timerange);  % The times that correspond to the timerange used for decoding

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
par.dels2use = 1:4;

par.regularize = 0; % if 0, linear regression. if 1, ridge regression
%% DECODING

close all

for sessix = 1:numel(meta)
    disp(['Decoding for session ' ' ' num2str(sessix) ' / ' num2str(numel(meta))])
    
    [X,Y,delLength,par] = preparePredictorsRegressors_Haz(par, sessix, kin, del,params); 

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

    y = reshape(Y.train,Y.size.train(1),Y.size.train(2)); % original input data (standardized)
    yhat = reshape(pred,Y.size.train(1),Y.size.train(2)); % prediction
%     %%%%%%%%%%%% SANITY CHECK %%%%%%%%%%%%%%%%%%%%%
%     subplot(1,2,1)
%     plot(y)
%     subplot(1,2,2)
%     plot(yhat)
%     pause
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    for c = 1:length(par.cond2use)
        if c==1
            cond = 'Rhit';
        else
            cond = 'Lhit';
        end
        condix = par.condassign==c;
        for delix = 1:4
            switch delix
                case 1
                    ixs = delLength<0.4;
                case 2
                    ixs = delLength<0.7&delLength>0.5;
                case 3
                    ixs = delLength<1.3&delLength>1.1;
                case 4
                    ixs = delLength<1.9&delLength>1.7;
            end
            ixrange = find(condix&ixs);
            trueVals.(cond){sessix,delix} = y(:,ixrange);
            modelpred.(cond){sessix,delix} = yhat(:,ixrange);
        end
    end
end

disp('---FINISHED DECODING FOR ALL SESSIONS---')
%% Baseline subtract CDTrialType
basesub = 0;

% Times that you want to use to baseline normalize CDTrialType
% trialstart = mode(obj(1).bp.ev.bitStart)-mode(obj(1).bp.ev.(params(1).alignEvent));
% start = find(obj(1).time>trialstart,1,'first');
% samp = mode(obj(1).bp.ev.sample)-mode(obj(1).bp.ev.(params(1).alignEvent));
% stop = find(obj(1).time<samp,1,'last');

trialstart = mode(obj(1).bp.ev.bitStart)-mode(obj(1).bp.ev.(params(1).alignEvent));
start = find(par.timerange>trialstart,1,'first');
samp = mode(obj(1).bp.ev.sample)-mode(obj(1).bp.ev.(params(1).alignEvent));
stop = find(par.timerange<samp,1,'last');

if basesub==1
    cond2use = {'Rhit','Lhit'};
    for sessix = 1:length(meta)
        for c = cond2use
        for delix = 1:length(par.dels2use)
            curr = modelpred.(c){sessix,delix};
            presampcurr = curr(start:stop,:);
            presampcurr = mean(presampcurr,1,'omitnan');
            presampcurr = mean(presampcurr,'omitnan');

            modelpred.(c){sessix,delix} = modelpred.(c){sessix,delix}-presampcurr;

            curr = trueVals.(c){sessix,delix};
            presampcurr = curr(start:stop,:);
            presampcurr = mean(presampcurr,1,'omitnan');
            presampcurr = mean(presampcurr,'omitnan');

            trueVals.(c){sessix,delix} = trueVals.(c){sessix,delix}-presampcurr;
        end
        end
    end
end
clearvars -except avgloadings del kin meta modelpred obj par params trueVals
%% Sanity check to make sure that R/L trial assignments are being done properly in decoding
% sessix = 8;
% delix = 2;
% for sessix = 1:19
%     subplot(1,2,1)
%     imagesc(trueVals.Rhit{sessix,delix}'); colorbar
%     hold off
%     subplot(1,2,2)
%     imagesc(trueVals.Lhit{sessix,delix}'); colorbar
%     hold off;
%     pause
% end
%% Set any values before the trial start to be 0
start = 1;
gotime = 0;
samplength = 1.3;
presamp = 0.2;

condfns = {'Rhit','Lhit'};
dels = [0.3 0.6 1.2 1.8];
for sessix = 1:length(meta)
    for delix = 1:length(dels)
        trialstart = gotime-dels(delix)-samplength-presamp;
        stop = find(par.timerange<trialstart,1,'last');
        for cond = condfns
            modelpred.(cond{1}){sessix,delix}(start:stop,:) = 0;
            trueVals.(cond{1}){sessix,delix}(start:stop,:) = 0;
        end
    end
end
%%
%%%%%%%%%%%%%%%% SANITY CHECK %%%%%%%%%%%%%%%%%%%%%%%%%%
sm = 50;
delix = 3;
for sessix = 1:19
%     subplot(1,2,1)
    plot(par.timeaxis,mySmooth(trueVals.Lhit{sessix,delix}(:,1:10),sm),'Color',[1 0 0]); hold on;
    plot(par.timeaxis,mySmooth(modelpred.Lhit{sessix,delix}(:,1:10),sm),'Color',[1 0.5 0.5]); 

%     subplot(1,2,2)
    plot(par.timeaxis,mySmooth(trueVals.Rhit{sessix,delix}(:,1:10),sm),'Color',[0 0 1]); 
    plot(par.timeaxis,mySmooth(modelpred.Rhit{sessix,delix}(:,1:10),sm),'Color',[0.5 0.5 1]); hold off;
    pause
end
%% Make heatmaps for a single session showing CDTrialType across trials and predicted CDTrialType
sm = 30;
invertCD = 'invert';                    % 'Invert' or 'no' for whether or not you want to flip the sign of the CD projection

load('C:\Code\Uninstructed-Movements\LeftRightDiverging_Colormap.mat')

% Times that you want to use to sort CDTrialType
% delay = mode(obj(1).bp.ev.delay)-mode(obj(1).bp.ev.(params(1).alignEvent));
% start = find(obj(1).time>delay,1,'first');
% resp = mode(obj(1).bp.ev.goCue)-mode(obj(1).bp.ev.(params(1).alignEvent))-0.05;
% stop = find(obj(1).time<resp,1,'last');

delay = mode(obj(1).bp.ev.delay)-mode(obj(1).bp.ev.(params(1).alignEvent));
start = find(par.timeaxis>delay,1,'first');
resp = mode(obj(1).bp.ev.goCue)-mode(obj(1).bp.ev.(params(1).alignEvent))-0.05;
stop = find(par.timeaxis<resp,1,'last');

cond2plot = {'Lhit','Rhit'};
for sessix = 1:length(meta)                                                                  % For each session...
    figure();
    cnt = 0;
    tempTrue = [];
    tempPred = [];
    ll = [];
    % Combine the true values for CDTrialType and the model predicted values across conditions
    % tempTrue = (time x [num left trials + num right trials])
    for c = 1:length(cond2plot)                                                 % For left and right trials...
        cond = cond2plot{c};
        for delix = 1:length(par.dels2use)
            currTrue = trueVals.(cond){sessix,delix};                                     % Get the true single trial CDTrialType projections for that condition and session
            [~,sortix] = sort(mean(currTrue(start:stop,:),1,'omitnan'),'descend');  % Sort the true projections by average magnitude during the delay period
            tempTrue = [tempTrue,currTrue(:,sortix)];
            currPred = modelpred.(cond){sessix,delix};                                    % Get the model predicted single trial CDTrialType projections
            tempPred = [tempPred,currPred(:,sortix)];
            ll = [ll,size(tempTrue,2)];
        end
    end
    nTrials = size(tempTrue,2);                                                 % Total number of trials that are being plotted
    ax1 = subplot(1,2,1);                                                       % Plot true CDTrialType data on left subplot
    if strcmp(invertCD,'invert')
        imagesc(obj(sessix).time(par.timerange),1:nTrials,-1*tempTrue'); hold on                      % Heatmap of true data (sorted left trials will be on top, then a white line, then sorted right trials)
    else
        imagesc(obj(sessix).time(par.timerange),1:nTrials,tempTrue'); hold on;
    end

    cnt = 1;
    for lix = 1:length(ll)
        if lix==1||lix==5
            xx = -0.3;
        elseif lix==2||lix==6
            xx = -0.6;
        elseif lix==3||lix==7
            xx = -1.2;
        elseif lix==4||lix==8
            xx = -1.8;
        end
        line([xx,xx],[cnt,ll(lix)],'Color','black','LineStyle','--')
        line([xx-1.3,xx-1.3],[cnt,ll(lix)],'Color','black','LineStyle','--')
        cnt = ll(lix)+0.5;
    end
    line([obj(sessix).time(1),obj(sessix).time(end)],[ll(4),ll(4)],'Color','black','LineStyle','-')


    ax2 = subplot(1,2,2);
    if strcmp(invertCD,'invert')
        imagesc(obj(sessix).time(par.timerange),1:nTrials,mySmooth(-1*tempPred,sm+30)'); hold on
    else
        imagesc(obj(sessix).time(par.timerange),1:nTrials,mySmooth(tempPred,sm+30)'); hold on
    end
    cnt = 1;
    for lix = 1:length(ll)
        if lix==1||lix==5
            xx = -0.3;
        elseif lix==2||lix==6
            xx = -0.6;
        elseif lix==3||lix==7
            xx = -1.2;
        elseif lix==4||lix==8
            xx = -1.8;
        end
        line([xx,xx],[cnt,ll(lix)],'Color','black','LineStyle','--')
        line([xx-1.3,xx-1.3],[cnt,ll(lix)],'Color','black','LineStyle','--')
        cnt = ll(lix)+0.5;
    end
    line([obj(sessix).time(1),obj(sessix).time(end)],[ll(4),ll(4)],'Color','black','LineStyle','-')
    title(ax1,'CDTrialType - data')
    colorbar(ax1)
    clim(ax1,[-3.5 3.5])
    colormap(LeftRightDiverging_Colormap)
    xlabel(ax1,'Time from go cue (s)')
    %     xline(ax1,0,'k--','LineWidth',1)
    %     xline(ax1,-0.9,'k--','LineWidth',1)
    %     xline(ax1,-2.2,'k--','LineWidth',1)
    xlim(ax1,[-2.5 0])
    set(ax1,'color',0.25*[1 1 1]);
    set(gca,'TickDir','out');

    title(ax2,'Model prediction')
    xlabel(ax2,'Time from go cue (s)')
    %     xline(ax2,0,'k--','LineWidth',1)
    %     xline(ax2,-0.9,'k--','LineWidth',1)
    %     xline(ax2,-2.2,'k--','LineWidth',1)
    colorbar(ax2)
    colormap(LeftRightDiverging_Colormap)
    xlim(ax2,[-2.5 0])
    clim(ax2,[-1 1])
    set(ax2,'color',0.25*[1 1 1]);
    set(gca,'TickDir','out');

    sgtitle(['Example session:  ' meta(sessix).anm ' ' meta(sessix).date])
end
%% Example plots by session for relating predicted and true CDTrialType
delR2_ALL = [];

plotexample = 'no';

if strcmp(plotexample,'yes')
    plotrange = exsess;
else
    plotrange = 1:length(meta);
end

delix = 2;

delay = params(1).delay(delix) -mode(obj(1).bp.ev.(params(1).alignEvent));
start = find(par.timeaxis>delay,1,'first');
resp = mode(obj(1).bp.ev.goCue)-mode(obj(1).bp.ev.(params(1).alignEvent))-0.05;
stop = find(par.timeaxis<resp,1,'last');

colors = getColors();
alph  = 0.2;

for sessix = plotrange
    % Calculate averages and standard deviation for true CD and predicted CD  this session
    [avgCD,stdCD] = getAvgStd(trueVals,modelpred,sessix,delix);
    figure();

    %%% Plot a scatter plot for a single session of true CDlate and predicted CDlate for each trial
    %%% Each dot = an average value of CDlate during the delay period
    subplot(1,2,2)
    tempR2 = Scatter_ModelPred_TrueCDTrialType(trueVals, modelpred, sessix, start, stop,meta,'invert');

    %%% Plot an example session of CDlate prediction vs true value
    subplot(1,2,1)
    plotExampleCDTrialType_Pred_Haz(colors, obj, par, meta, avgCD, stdCD, sessix, trueVals,alph,'invert',delix);

    % Save R2 value for that session
    delR2_ALL = [delR2_ALL, tempR2];
end
%% Plot bar plot to show average R2 values across sessions
exsess = 1;

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
title(['CDChoice_RandDel: Ex session = ' meta(exsess).anm meta(exsess).date])
%% Print summary statistics
if par.regularize
    regtype = 'Ridge';
else
    regtype = 'none';
end

disp('---Summary statistics for CDTrialType prediction---')
disp(['Average R2 across all sessions (n = ' num2str(length(meta)) ' ) = ' num2str(mean(delR2_ALL))])
disp(['Standard deviation across all sessions = ' num2str(std(delR2_ALL))])
disp(['Regularization type: ' regtype])
t = datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss Z');
disp(t)

%% FUNCTIONS
function [avgCD,stdCD] = getAvgStd(trueVals,modelpred,sessix,delix)

avgCD.Rhit.true = mean(mySmooth(trueVals.Rhit{sessix,delix},81),2,'omitnan');      % Get average true CDlate for R and L hits for this session
avgCD.Lhit.true = mean(mySmooth(trueVals.Lhit{sessix,delix},81),2,'omitnan');
stdCD.Rhit.true = std(mySmooth(trueVals.Rhit{sessix,delix},81),0,2,'omitnan');     % Get standard deviation of true CDlate for R and L hits
stdCD.Lhit.true = std(mySmooth(trueVals.Lhit{sessix,delix},81),0,2,'omitnan');

modelpred.Rhit{sessix,delix} = fillmissing(modelpred.Rhit{sessix,delix},"nearest");
modelpred.Lhit{sessix,delix} = fillmissing(modelpred.Lhit{sessix,delix},"nearest");
infix = find(isinf(modelpred.Rhit{sessix,delix})); modelpred.Rhit{sessix,delix}(infix) = 0;
infix = find(isinf(modelpred.Lhit{sessix,delix})); modelpred.Lhit{sessix,delix}(infix) = 0;
avgCD.Rhit.pred = mean(mySmooth(modelpred.Rhit{sessix,delix},81),2,'omitnan');     % Get average predicted CDlate for R and L hits for this session
avgCD.Lhit.pred = mean(mySmooth(modelpred.Lhit{sessix,delix},81),2,'omitnan');
stdCD.Rhit.pred = std(mySmooth(modelpred.Rhit{sessix,delix},81),0,2,'omitnan');    % Get stdev of predicted CDlate for R and L hits for this session
stdCD.Lhit.pred = std(mySmooth(modelpred.Lhit{sessix,delix},81),0,2,'omitnan');
end