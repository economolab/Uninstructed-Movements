% same as st_elsayed.m but first remove CDs from full population single trial activity

clear,clc,close all

% add paths for data loading scripts, all fig funcs, and utils
utilspth = 'C:\Users\munib\Documents\Economo-Lab\code\uninstructedMovements_v2';
addpath(genpath(fullfile(utilspth,'DataLoadingScripts')));
addpath(genpath(fullfile(utilspth,'funcs')));
addpath(genpath(fullfile(utilspth,'utils')));
rmpath(genpath(fullfile(utilspth,'mc_stim/')))

% add paths for figure specific functions
addpath(genpath(pwd))

%% PARAMETERS
params.alignEvent          = 'goCue'; % 'jawOnset' 'goCue'  'moveOnset'  'firstLick'  'lastLick'

% time warping only operates on neural data for now.
% TODO: time warp for video and bpod data
params.timeWarp            = 0;  % piecewise linear time warping - each lick duration on each trial gets warped to median lick duration for that lick across trials
params.nLicks              = 20; % number of post go cue licks to calculate median lick duration for and warp individual trials to 

params.lowFR               = 1; % remove clusters with firing rates across all trials less than this val

% set conditions to calculate PSTHs for
params.condition(1)     = {'(hit|miss|no)'};                             % all trials
params.condition(end+1) = {'R&hit&~stim.enable&~autowater'};             % right hits, no stim, aw off
params.condition(end+1) = {'L&hit&~stim.enable&~autowater'};             % left hits, no stim, aw off
params.condition(end+1) = {'R&miss&~stim.enable&~autowater'};            % error right, no stim, aw off
params.condition(end+1) = {'L&miss&~stim.enable&~autowater'};            % error left, no stim, aw off

params.tmin = -1.6;
params.tmax = 2;
params.dt = 1/100;

% smooth with causal gaussian kernel
params.smooth = 15;

% cluster qualities to use
params.quality = {'all'}; % accepts any cell array of strings - special character 'all' returns clusters of any quality


params.traj_features = {{'tongue','left_tongue','right_tongue','jaw','trident','nose'},...
                        {'top_tongue','topleft_tongue','bottom_tongue','bottomleft_tongue','jaw','top_paw','bottom_paw','top_nostril','bottom_nostril'}};

params.feat_varToExplain = 80; % num factors for dim reduction of video features should explain this much variance

params.N_varToExplain = 80; % keep num dims that explains this much variance in neural data (when doing n/p)

params.advance_movement = 0.0;


%% SPECIFY DATA TO LOAD

datapth = '/Users/Munib/Documents/Economo-Lab/data/';

meta = [];

meta = loadJEB6_ALMVideo(meta,datapth);
meta = loadJEB7_ALMVideo(meta,datapth);
meta = loadEKH1_ALMVideo(meta,datapth);
meta = loadEKH3_ALMVideo(meta,datapth);
meta = loadJGR2_ALMVideo(meta,datapth);
meta = loadJGR3_ALMVideo(meta,datapth);
meta = loadJEB15_ALMVideo(meta,datapth);
% meta = loadJEB14_ALMVideo(meta,datapth);

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

%% Null and Potent Space - remove coding dimensions first

clearvars -except obj meta params me sav

cond2use = [2 3]; % left hit, right hit
for sessix = 1:numel(meta)
    trialdat = removeCodingDimensions(obj(sessix),params(sessix),cond2use);
    cond2use = [2 3 4 5]; % right hit, left hit, right miss, left miss
    rez(sessix) = singleTrial_elsayed_np(trialdat, obj(sessix), me(sessix), params(sessix), cond2use);
    cond2use = [1 2]; % right hits, left hits (corresponding to null/potent psths in rez)
    cd_null(sessix) = getCodingDimensions(rez(sessix).N_null_psth,trialdat,obj(sessix),params(sessix),cond2use);
    cd_potent(sessix) = getCodingDimensions(rez(sessix).N_potent_psth,trialdat,obj(sessix),params(sessix),cond2use);
end

cd_null_all = concatRezAcrossSessions(cd_null);
cd_potent_all = concatRezAcrossSessions(cd_potent);




%% plots

close all

sav = 0;


% -----------------------------------------------------------------------
% -- Null and Potent Space Single Trial Projections --
% -----------------------------------------------------------------------
 
% - projections showing move / quiet somehow
cond2use = [2 3]; % right hits, left hits
ndims = 3;
% plotSingleTrialNPHeatmaps(rez,params,me,ndims,cond2use)


% - how much variance in move and non-move time points
cond2use = [2 3]; % right hits, left hits
% plotVarianceInEpochs(rez,me,params,cond2use);                       

% - ve
% plotVarianceExplained_NP(rez);

% - ve over time
% plotVarianceExplained_NP_overTime(rez);


% -----------------------------------------------------------------------
% -- Null and Potent Space Trial-averaged Projections --
% -----------------------------------------------------------------------

ndims = 10; % how many n/p dimensions to plot, in order of variance explained
cond2plot = [3 4]; % right hit, left hit
% plot_NP_PSTH(rez,obj,params,ndims,cond2plot)


% -----------------------------------------------------------------------
% -- Coding Dimensions --
% -----------------------------------------------------------------------

% % - null
plotCDProj(cd_null_all,cd_null,sav)
% plotCDVarExp(cd_null_all,sav)
% plotSelectivity(cd_null_all,cd_null,sav)
% plotSelectivityExplained(cd_null_all,cd_null,sav)
% title('Null')

% % - potent
plotCDProj(cd_potent_all,cd_potent,sav)
% plotCDVarExp(cd_potent_all,sav)
% plotSelectivity(cd_potent_all,cd_potent,sav)
% plotSelectivityExplained(cd_potent_all,cd_potent,sav)
% title('Potent')


%% t=0 is the go cue, but only on trials where the animals were not moving PRIOR to the go cue
% same plots as plotSelectivityExplained 

% plotSelectivityExplained_GoCue_Transitions(rez,cd_potent_all,cd_potent,me,obj(1).time,sav)
findMovementBouts(rez,cd_potent_all,cd_potent,me,obj)

%% t=0 is transitions between non-movement and movement that do not coincide with the go cue
% same plots as plotSelectivityExplained 
































