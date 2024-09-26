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
addpath(genpath(fullfile(utilspth,'figNP')));
figpth = [basepth  '\Uninstructed-Movements\Fig 2'];
addpath(genpath(fullfile(figpth,'funcs')));
addpath(genpath(fullfile(figpth,'Hazarded Delay')));
figpth = [basepth  '\Uninstructed-Movements\functions'];
addpath(genpath(fullfile(figpth,'HazDel funcs')));
figpth = [basepth  '\Uninstructed-Movements\Fig 5'];
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

% Haz delay params
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
%% Look at condition-averaged CDChoice across sessions
del2use = 2;        % Relative to dels4proj above
delLine = -1.2;
for cond = 1:length(condfns)
    temp = [];
    for sessix = 1:length(meta)
        curr = del(sessix).condProj{del2use}(:,cond);
        temp = [temp,curr];
    end
    avgCDAll.(condfns{cond}) = temp;
end

% Plot
figure();
colors = getColors();
cols = {colors.rhit, colors.lhit};
alph = 0.2;
for cond = 1:length(condfns)
    toplot = mean(avgCDAll.(condfns{cond}),2,'omitnan');
    err = 1.96*(std(avgCDAll.(condfns{cond}),0,2,'omitnan')./sqrt(length(meta)));
    ax = gca;
    shadedErrorBar(obj(sessix).time+1.2,toplot,err,{'Color',cols{cond},'LineWidth',2},alph,ax); hold on;
end