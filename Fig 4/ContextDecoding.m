clear,clc,close all

% add paths for data loading scripts, all fig funcs, and utils
utilspth = 'C:\Users\Jackie\Documents\Grad School\Economo Lab\Code\Uninstructed-Movements\NullPotent_Context';
addpath(genpath(fullfile(utilspth,'DataLoadingScripts')));
addpath(genpath(fullfile(utilspth,'funcs')));
addpath(genpath(fullfile(utilspth,'utils')));
% rmpath(genpath(fullfile(utilspth,'fig1/')));
% rmpath(genpath(fullfile(utilspth,'mc_stim/')));

% add paths for figure specific functions
addpath(genpath(pwd))

%% PARAMETERS
params.alignEvent          = 'firstLick'; % 'jawOnset' 'goCue'  'moveOnset'  'firstLick'  'lastLick'

% time warping only operates on neural data for now.
% TODO: time warp for video and bpod data
params.timeWarp            = 0;  % piecewise linear time warping - each lick duration on each trial gets warped to median lick duration for that lick across trials
params.nLicks              = 20; % number of post go cue licks to calculate median lick duration for and warp individual trials to

params.lowFR               = 1; % remove clusters with firing rates across all trials less than this val

% set conditions to calculate PSTHs for
params.condition(1)     = {'(hit|miss|no)'};                             % all trials
params.condition(end+1) = {'hit&~stim.enable&~autowater'};               % all 2AFC hits, no stim
params.condition(end+1) = {'hit&~stim.enable&autowater'};                % all AW hits, no stim
params.condition(end+1) = {'miss&~stim.enable&~autowater'};              % error 2AFC, no stim, aw off
params.condition(end+1) = {'miss&~stim.enable&autowater'};              % error AW, no stim

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
meta = loadJEB14_ALMVideo(meta,datapth);
meta = loadJEB15_ALMVideo(meta,datapth);

% --- M1TJ ---
% meta = loadJEB14_M1TJVideo(meta,datapth);

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
%%
for sessix = 1:numel(meta)
    me(sessix) = loadMotionEnergy(obj(sessix), meta(sessix), params(sessix), datapth);
end
%% Decode ContextMode from DLC features and ME

clearvars -except obj meta params me sav

% -----------------------------------------------------------------------
% -- Curate Input Data --
% zscore single trial neural data (time*trials,neurons), for all trials
% -- Calculate null and potent spaces --
% null space from quiet time points
% potent space from moving time points
% -- Calculate coding directions from null and potent spaces --
% -----------------------------------------------------------------------

for sessix = 1:numel(meta)

    % -- input data
    trialdat_zscored = zscore_singleTrialNeuralData(obj(sessix).trialdat);

    % -- null and potent spaces
    cond2use = [2 3 4 5]; % 2AFC hit, AW hit, 2AFC miss, AW miss
    nullalltime = 0;      % use all time points to estimate null space if 1
    rez(sessix) = singleTrial_elsayed_np(trialdat_zscored, obj(sessix), me(sessix), params(sessix), cond2use, nullalltime);

    % -- coding dimensions
    cond2use = [1 2];          % 2AFC hits, AW hits (corresponding to null/potent psths in rez)
    cond2proj = [1:4];         % 2AFC hits, AW hits, 2AFC miss, AW miss (corresponding to null/potent psths in rez)
    cond2use_trialdat = [2 3]; % for calculating selectivity explained in full neural pop
    cd_null(sessix) = getCodingDimensions_Context(rez(sessix).N_null_psth,trialdat_zscored,obj(sessix),params(sessix),cond2use,cond2use_trialdat, cond2proj);
    cd_potent(sessix) = getCodingDimensions_Context(rez(sessix).N_potent_psth,trialdat_zscored,obj(sessix),params(sessix),cond2use,cond2use_trialdat, cond2proj);

end
