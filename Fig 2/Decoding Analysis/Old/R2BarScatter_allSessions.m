% DECODING CDlate FROM ALL KINEMATIC FEATURES
clear,clc,close all

whichcomp = 'Laptop';                                                % LabPC or Laptop

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

if strcmp(whichcomp,'LabPC')
    datapth = 'C:\Users\Jackie Birnbaum\Documents\Data';
elseif strcmp(whichcomp,'Laptop')
    datapth = 'C:\Users\Jackie\Documents\Grad School\Economo Lab';
end

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
regr = getCodingDimensions_2afc(obj,params,cond2use,cond2proj);

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
clc;

% params for decoding
rez.nFolds = 4;                                     % number of iterations (bootstrap)
rez.binSize = 30;                                   % Bin size to decode over (in milliseconds)
difTime = params(1).dt*1000;                        % Convert time-difference from s to ms
rez.dt = floor(rez.binSize / difTime);              % How many samples confer your bin size
rez.tm = obj(1).time(1:rez.dt:numel(obj(1).time));  % Create a new time axis over which to decode (based on the bin size that you want)
rez.numT = numel(rez.tm);                           % Number of time-points
rez.train = 1;                                      % fraction of trials to use for training (1-train for testing)

% match number of right and left hits, and right and left misses
cond2use = 2:5;                                     % Which conditions you want to include (R hit, L hit, R miss, L miss)
hitcond = [1 2];                                    % Which conditions out of cond2use are hit
misscond = [3 4];                                   % Which conditions out of cond2use are miss

% Do the decoding of CDTrialType from all DLC features
%%% trueVals = (1x1) struct with fields 'Rhit' and 'Lhit'.  Each of these fields is an (nSessions x 1) cell array.  
%%% Each cell contains (time x trials) array of the true CDTrialType values 
%%% modelpred is structured in the same way but this reflects the prediction of the multiple linear regression model
[trueVals, modelpred] = doCDTrialTypeDecoding_fromDLC(nSessions, kin, obj, cond2use, hitcond, misscond, regr,rez,params);

 disp('---FINISHED DECODING FOR ALL SESSIONS---')
%% Get R2 value for all sessions during delay
JEB13_delR2 = [];

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
colors = getColors();
load('C:\Users\Jackie\Documents\Grad School\Economo Lab\Code\Uninstructed-Movements\Fig 2\Decoding Analysis\delR2Vals_NoJEB13.mat')
delR2_ALL = [delR2, JEB13_delR2];                              % Group the R2 values for JEB13 and all other animals
anmNames = {'JEB6', 'JEB7', 'JEB7', 'EKH1','JGR2','JGR2','JGR3',};
JEB13_anmNames = {'JEB13','JEB13','JEB13','JEB13','JEB13','JEB13','JEB14','JEB14','JEB14','JEB14','JEB15','JEB15','JEB15','JEB15'};
anmNames_all = [anmNames,JEB13_anmNames];                               % Animal names for each session
nSessions = numel(anmNames_all);

uniqueAnm = unique(anmNames_all);

exsess = 3;                                                            % The index of the session that you want to be highlighted
markerSize = 60;
figure();
bar(mean(delR2_ALL),'FaceColor',colors.afc); hold on;         % Plot the average R2 value across all sessions
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
    scatter(1,delR2_ALL(sessix),markerSize,'filled',shape,'MarkerFaceColor',[0.65 0.65 0.65]); hold on;
end
scatter(1,delR2_ALL(exsess),markerSize,'filled','pentagram','black','MarkerEdgeColor','black')
legend([' ',anmNames_all])
ylim([0.4 1])
ax = gca;
ax.FontSize = 16;
title(['Ex session = Sesh ' num2str(exsess) '; Animal ' anmNames_all{exsess} ])