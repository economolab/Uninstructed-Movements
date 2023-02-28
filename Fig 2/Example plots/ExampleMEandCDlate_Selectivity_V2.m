% Example of session with selectivity in CDLate and Motion Energy
clear,clc,close all

% add paths
utilspth = 'C:\Users\Jackie\Documents\Grad School\Economo Lab\Code\Munib Uninstruct Move\uninstructedMovements_v2';
addpath(genpath(fullfile(utilspth,'DataLoadingScripts')));
addpath(genpath(fullfile(utilspth,'funcs')));
addpath(genpath(fullfile(utilspth,'utils')));
addpath(genpath(fullfile(utilspth,'fig3')));
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
meta = loadJEB6_ALMVideo(meta,datapth);
meta = loadJEB7_ALMVideo(meta,datapth);
meta = loadEKH1_ALMVideo(meta,datapth);
meta = loadEKH3_ALMVideo(meta,datapth);
meta = loadJGR2_ALMVideo(meta,datapth);
meta = loadJGR3_ALMVideo(meta,datapth);
%meta = loadJEB13_ALMVideo(meta,datapth);
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
orthogonalize = 'non-orthog';                                   % 'Orthogonalize' if you want to orthogonalize CDs to each other
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
%% Plot average ME and average CD-late for specified conditions (example session)
cond2use = [2,3];                                                   % Conditions you want to use
sessix = 6;                                                         % Session to use

feature = 'motion_energy';
featix = find(strcmp(kin(sessix).featLeg,feature));
ME = squeeze(kin(sessix).dat(:,:,featix));                          % Load single-trial ME values for this session
CDlate = regr(sessix).singleProj;                                   % Single-trial CDlate projections

sample = mode(obj(1).bp.ev.sample) - mode(obj(1).bp.ev.goCue);      % Timing of sample, delay and trialstart for plotting
delay = mode(obj(1).bp.ev.delay) - mode(obj(1).bp.ev.goCue);
trialStart = mode(obj(1).bp.ev.bitStart) - mode(obj(1).bp.ev.goCue);

colors = getColorsUpdated();
alph = 0.2;                                                                 % For opacity of confidence intervals
smooth = 21;                                                                % Smoothing average CDlate projections and ME
figure()    
for c = 1:length(cond2use)                                                  % For each condition
    if c==1                                                                 % Specify the color to plot in
        col = colors.rhit;
    else
        col = colors.lhit;
    end
    cond = cond2use(c);                                                     
    condtrix = params(sessix).trialid{cond};                                % Trials from this condition
    ntrix = length(condtrix);                                               % nTrials in condition
    condME = ME(:,condtrix);                                                % ME values for the trials in this condition
    toplotME = mySmooth(mean(condME,2,'omitnan'),smooth);                   % Take the average ME value for this condition
    errME = 1.96*(mySmooth(std(condME,0,2),21)/sqrt(ntrix));                % 95% confidence intervals 
    
    condCD = CDlate(:,condtrix);                                            % Do the same thing for CDlate
    toplotCD = mean(condCD,2,'omitnan');
    errCD = 1.96*(std(condCD,0,2)/sqrt(ntrix));

    subplot(2,1,1)                                                          % Plot the average ME for this condition, with error bars and trial lines 
    ax = gca;
    shadedErrorBar(regr(sessix).time,toplotME,errME,{'Color',col,'LineWidth',2},alph,ax); 
    hold on;
    title('Motion energy')
    xline(sample,'k--','LineWidth',1)
    xline(delay,'k--','LineWidth',1)
    xlim([-2.3 0])

    subplot(2,1,2)                                                          % Plot the average CDlate for this condition, with error bars and trial lines
    ax = gca;
    shadedErrorBar(regr(sessix).time,toplotCD,errCD,{'Color',col,'LineWidth',2},alph,ax); 
    hold on;
    set(gca, 'YDir','reverse')
    title('CD Late')
    xline(sample,'k--','LineWidth',1)
    xline(delay,'k--','LineWidth',1)
    xlim([-2.3 0])
    xlabel('Time from go cue (s)')
end