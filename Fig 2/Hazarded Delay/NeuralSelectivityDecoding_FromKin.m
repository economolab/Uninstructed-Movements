% From Inagaki 2019 Nature (Discrete Attractor Dynamics paper)
% "Neurons with significant selectivity (two-sided Wilcoxon rank sum test comparing spike counts in two correct trial types, P < 0.05) during the delay epochs were classified as selective cells. Selective cells were classified into lick-right preferring versus lick-left preferring, on the basis of their total spike counts during the delay epoch.
% For peri-stimulus time histograms (Figs. 5, 6, Extended Data Fig. 8k), only correct trials were included. For the peri-stimulus time histograms and selectivity of the 
% random delay task (Fig. 6, Extended Data Fig. 8k), only spikes before the go cue were pooled. Spikes were averaged over 100 ms with a 1-ms sliding window."

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
%% PARAMETERS
params.alignEvent          = 'delay'; % 'jawOnset' 'goCue'  'moveOnset'  'firstLick'  'lastLick'

% time warping only operates on neural data for now.
% TODO: time warp for video and bpod data
params.timeWarp            = 0;  % piecewise linear time warping - each lick duration on each trial gets warped to median lick duration for that lick across trials
params.nLicks              = 20; % number of post go cue licks to calculate median lick duration for and warp individual trials to

params.lowFR               = 1; % remove clusters with firing rates across all trials less than this val

% set conditions to calculate PSTHs for                            % all trials
params.condition(1) = {'R&hit&~stim.enable&~autowater&~early'};             % R 2AFC hits, no stim
params.condition(end+1) = {'L&hit&~stim.enable&~autowater&~early'};             % L 2AFC hits, no stim

params.tmin = -1.6;
params.tmax = 2.5;
params.dt = 1/200;

% smooth with causal gaussian kernel
params.smooth = 31;

% cluster qualities to use
params.quality = {'all'}; % accepts any cell array of strings - special character 'all' returns clusters of any quality

% vid features to use
params.traj_features = {{'tongue','left_tongue','right_tongue','jaw','trident','nose'},...
    {'top_tongue','topleft_tongue','bottom_tongue','bottomleft_tongue','jaw','top_nostril','bottom_nostril'}};

params.feat_varToExplain = 80; % num factors for dim reduction of video features should explain this much variance
params.N_varToExplain = 80; % keep num dims that explains this much variance in neural data (when doing n/p)
params.advance_movement = 0;

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
%% Find which cells are selective for trial-type during the delay period 
%clearvars -except obj meta params
smooth = 50;
modparams.quals2excl = {'Poor','Noisy'};
modparams.sm = 30;                                                    % Amount that you want to smooth PSTHs by
modparams.subTrials = 40;
for sessix = 1:length(meta)
    obj = findDelSelectiveCells(sessix, obj, meta, modparams, params); 
end
%% Calculate selectivity (spike rate difference) for all selective cells (on a given delay duration)
del2use = 1.2000;
cond2use = [1,2];
smooth = 200;
selectivity_All = [];
for sessix = 1:length(meta)
    selectivity_All = calcNeuralSelectivity_Haz(sessix, obj, del2use, params, smooth, cond2use);
end
%%
tempobj = obj(1); tempparams = params(1);
clearvars -except selectivity_All tempobj tempparams
obj = tempobj; params = tempparams;
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------------------------Static Delay-----------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% PARAMETERS
ctrlparams.alignEvent          = 'delay'; % 'jawOnset' 'goCue'  'moveOnset'  'firstLick'  'lastLick'

% time warping only operates on neural data for now.
% TODO: time warp for video and bpod data
ctrlparams.timeWarp            = 0;  % piecewise linear time warping - each lick duration on each trial gets warped to median lick duration for that lick across trials
ctrlparams.nLicks              = 20; % number of post go cue licks to calculate median lick duration for and warp individual trials to

ctrlparams.lowFR               = 1; % remove clusters with firing rates across all trials less than this val

% set conditions to calculate PSTHs for                            % all trials
ctrlparams.condition(1) = {'R&hit&~stim.enable&~autowater&~early'};             % R 2AFC hits, no stim
ctrlparams.condition(end+1) = {'L&hit&~stim.enable&~autowater&~early'};             % L 2AFC hits, no stim

ctrlparams.tmin = -1.6;
ctrlparams.tmax = 2.5;
ctrlparams.dt = 1/200;

% smooth with causal gaussian kernel
ctrlparams.smooth = 31;

% cluster qualities to use
ctrlparams.quality = {'all'}; % accepts any cell array of strings - special character 'all' returns clusters of any quality

% vid features to use
ctrlparams.traj_features = {{'tongue','left_tongue','right_tongue','jaw','trident','nose'},...
    {'top_tongue','topleft_tongue','bottom_tongue','bottomleft_tongue','jaw','top_paw','bottom_paw','top_nostril','bottom_nostril'}};

ctrlparams.feat_varToExplain = 80; % num factors for dim reduction of video features should explain this much variance
ctrlparams.N_varToExplain = 80; % keep num dims that explains this much variance in neural data (when doing n/p)
ctrlparams.advance_movement = 0;

datapth = 'C:\Users\Jackie\Documents\Grad School\Economo Lab';

ctrlmeta = [];

% --- ALM --- 
ctrlmeta = loadJEB6_ALMVideo(ctrlmeta,datapth);
ctrlmeta = loadJEB7_ALMVideo(ctrlmeta,datapth);
ctrlmeta = loadEKH1_ALMVideo(ctrlmeta,datapth);
ctrlmeta = loadEKH3_ALMVideo(ctrlmeta,datapth);
ctrlmeta = loadJGR2_ALMVideo(ctrlmeta,datapth);
ctrlmeta = loadJGR3_ALMVideo(ctrlmeta,datapth);
ctrlmeta = loadJEB13_ALMVideo(ctrlmeta,datapth);
ctrlmeta = loadJEB14_ALMVideo(ctrlmeta,datapth);
ctrlmeta = loadJEB15_ALMVideo(ctrlmeta,datapth);
ctrlmeta = loadJEB19_ALMVideo(ctrlmeta,datapth);

ctrlparams.probe = {ctrlmeta.probe}; % put probe numbers into params, one entry for element in ctrlmeta, just so i don't have to change code i've already written
%% LOAD DATA
% ----------------------------------------------
% -- Neural Data --
% obj (struct array) - one entry per session
% params (struct array) - one entry per session
% ----------------------------------------------
[ctrlobj,ctrlparams] = loadSessionData(ctrlmeta,ctrlparams);
%% Find which cells are selective for trial-type during the delay period 
smooth = 50;
modparams.quals2excl = {'Poor','Noisy'};
modparams.sm = 30;                                                    % Amount that you want to smooth PSTHs by
modparams.subTrials = 40;
for sessix = 1:length(ctrlmeta)
  ctrlobj = findSelectiveCells_Ctrl(sessix, ctrlobj, ctrlmeta, ctrlparams, modparams);
end
%% Calculate selectivity (spike rate difference) for all selective cells (on a given delay duration)
del2use = 0.9;
cond2use = [1,2];
smooth = 200;
selectivity_AllCtrl = [];
for sessix = 1:length(ctrlmeta)
    selectivity_AllCtrl = calcNeuralSelectivity_Haz(sessix, obj, del2use, params, smooth, cond2use);
end
