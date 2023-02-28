% Showing differences in movements across behavioral contexts
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
addpath(genpath(fullfile(otherpth,'Context Switching')));

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
%% Plot ME heatmap for example trials 
for sessix = 4%1:length(meta)
    date = meta(sessix).date;
    anm = meta(sessix).anm;

    feat = 'motion_energy';
    featix = find(strcmp(kin(sessix).featLeg,feat));
    ME = kin(sessix).dat(:,:,featix);
    cond2plot = 2:3;
    for c = 1:length(cond2plot)
        kin_by_cond.(feat){c} = ME(:,params(sessix).trialid{cond2plot(c)});         % Get the ME or move feature from the conditions that you want to plot
    end

    figure();
    numTrixPlot = 30;
    sm = 20;
    rangetoPlot = 1:numTrixPlot;
    for i=1:length(cond2plot)
        if i==1
            subplot(1,2,1);
        elseif i==2
            subplot(1,2,2);
        end
        imagesc(obj(1).time,1:numTrixPlot,mySmooth(kin_by_cond.(feat){i}(:,1:numTrixPlot),sm)'); colormap("linspecer");
        xlim([-2.7, 0])
        xlabel('Time before goCue (s)','FontSize',13)
        ylabel('Trials','FontSize',13)
        if i==1
            title('2AFC trials')
        elseif i==2
            title('AW trials')
        end
        c=colorbar;
        clim([0 0.75])
        ylabel(c,feat,'FontSize',12,'Rotation',90);
    end
    figtitle =  ['Example Trials from',anm,date];  % Name/title for session
    sgtitle(figtitle,'FontSize',16)
end
%% Plot kinematic features from example trials in 2AFC and AW
%'jaw_yvel_view1' ; 'nose_yvel_view1'
feat = 'nose_yvel_view1';                           % Specify the feature you want to look at
%feat = 'top_tongue_xdisp_view2';
featix = find(strcmp(kin(1).featLeg,feat));         % Find which index corresponds to that feature

cond = [2,3];                                       % Which conditions you want to plot the features for 
offset = 10;                                        % How offset you want the trace from each trial to be 
nTrials = 7;                                        % Number of trials from each condition that you want to plot
%plotExampleFeatTraces(meta,cond,params,kin,obj,offset,nTrials)

%% Find trials in which the animal switches between contexts
for sessix = 1:numel(meta)
    bp = obj(sessix).bp;
    [toAW_ix, toAFC_ix] = findSwitchTrials(bp);
    % toAW_ix = the first trial in an AW block; toAFC_ix = the first trial in a 2AFC block
    obj(sessix).toAW_ix = toAW_ix;  obj(sessix).toAFC_ix = toAFC_ix;
end

%% Find avg velocity for specified feature across trials in a given session on switch trials
nBufferTrials = 10;                              % Number of trials that you want to look at before and after switches
% Time period over which you want to average Feature Movement
start = find(obj(1).time>-3,1,'first');
samp = median(obj(1).bp.ev.sample)-median(obj(1).bp.ev.(params(1).alignEvent));
stop = find(obj(1).time<samp,1,'last');

%'jaw_yvel_view1' ; 'nose_yvel_view1'
feat = 'jaw_yvel_view1';                           % Specify the feature you want to look at
featix = find(strcmp(kin(1).featLeg,feat));         % Find which index corresponds to that feature

for sessix = 1:numel(meta)
    featkin = kin(sessix).dat(:,:,featix);
    
    % Find avg feature vel on 2AFC --> AWswitch trials (during the presample period)
    switchtype = 'toAW_ix';
    moveswitch(sessix).toAW_FeatVel = findFeatVel_SwitchAligned(nBufferTrials, obj, sessix, featkin,switchtype,start,stop);
   

    % Find avg feature vel on AW --> 2AFC switch trials (during the presample period)
    switchtype = 'toAFC_ix';
    moveswitch(sessix).to2AFC_FeatVel = findFeatVel_SwitchAligned(nBufferTrials, obj, sessix, featkin,switchtype,start,stop);
end

%% Plot switch-aligned Feature Velocity
% Concatenate all switch-aligned Feature Velocity from across sessions
toAW = []; to2AFC = [];
for sessix = 1:length(meta)
    toAW = [toAW; moveswitch(sessix).toAW_FeatVel];
    to2AFC = [to2AFC; moveswitch(sessix).to2AFC_FeatVel];
end
toAW = abs(toAW); to2AFC = abs(to2AFC);

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
ylabel([feat '(a.u.)'])

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
ylabel([feat '(a.u.)'])
end

function plotExampleFeatTraces(meta,cond,params,kin,obj,offset,nTrials)
for sessix = 1:length(meta)                         % For every session...
    figure();
    for c = 1:length(cond)                          % Go through each of the specified conditions...
        if c==1
            plottitle = '2AFC';
            col = 'black';
        else
            plottitle = 'AW';
            col = 'magenta';
        end
        temp = cond(c);
        condtrix = params(sessix).trialid{temp};    % Get the trials that correspond to this condition
        trix2plot = randsample(condtrix,nTrials);   % Randomly sample 'nTrials' number of trials from that condition to plot
        condfeat = kin(sessix).dat(:,trix2plot,featix); % Get the kinematic data from those trials for the desired feature
        subplot(1,2,c)
        for t = 1:size(condfeat,2)
            plot(obj(sessix).time,offset*t+condfeat(:,t),'Color',col,'LineWidth',2); hold on;   % Plot the movement trace from each of the trials that were randomly selected
        end
        title(plottitle)
        xlim([-3 0])
        xlabel(['Time from ' params(sessix).alignEvent ' (s)'])
        ylabel(feat)
    end
end
end