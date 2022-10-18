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
%% Null and Potent Space

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
%%
cd_null_all = concatRezAcrossSessions_Context(cd_null);
cd_potent_all = concatRezAcrossSessions_Context(cd_potent);

%% plots

close all

sav = 0;


% -----------------------------------------------------------------------
% -- Null and Potent Space Single Trial Projections --
% -----------------------------------------------------------------------

% % - projections showing move / quiet somehow
cond2use = [2 3]; % 2AFC hits, AW hits
ndims = 4; % top ndims variance explaining dimensions
% plotSingleTrialNPHeatmaps(rez,params,me,ndims,cond2use,meta);


% % - how much variance in move and non-move time points
cond2use = [2 3]; % 2AFC hits, AW hits
% plotVarianceInEpochs(rez,me,params,cond2use);

% % - ve
% plotVarianceExplained_NP(rez);

% % - ve over time (TODO)
% % % plotVarianceExplained_NP_overTime(rez);


% -----------------------------------------------------------------------
% -- Null and Potent Space Trial-averaged Projections --
% -----------------------------------------------------------------------

ndims = 5; % how many n/p dimensions to plot, in order of variance explained
cond2plot = [1 2]; % 2AFC hit, AW hit
plot_NP_PSTH(rez,obj,params,ndims,cond2plot,meta)


% -----------------------------------------------------------------------
% -- Coding Dimensions --
% -----------------------------------------------------------------------
plotmiss = 0;

titlestring = 'Null';
plotCDProj_Context(cd_null_all,cd_null,sav,titlestring,plotmiss)
% plotCDVarExp_Context(cd_null_all,sav,titlestring)
plotSelectivity_Context(cd_null_all,cd_null,sav,titlestring)
% plotSelectivityExplained_Context(cd_null_all,cd_null,sav,titlestring)

titlestring = 'Potent';
plotCDProj_Context(cd_potent_all,cd_potent,sav,titlestring,plotmiss)
% plotCDVarExp_Context(cd_potent_all,sav,titlestring)
plotSelectivity_Context(cd_potent_all,cd_potent,sav,titlestring)
% plotSelectivityExplained_Context(cd_potent_all,cd_potent,sav,titlestring)



%% t=0 is the go cue, but only on trials where the animals were not moving PRIOR to the go cue
% same plots as plotSelectivityExplained

%% t=0 is transitions between non-movement and movement that do not coincide with the go cue
% same plots as plotSelectivityExplained

stitch_dist = 0.025; % in seconds, stitch together movement bouts shorter than this
purge_dist = 0.1; % in seconds, remove move bouts shorter than this value, after stitching complete
tbout = 0.3; % move/non-move bout/transition required to be at least this long in seconds
[dat.mdat,dat.mdat_leg,dat.qdat,dat.qdat_leg,newme] = nonGC_moveTransitions(obj,me,params,stitch_dist,purge_dist,tbout);

%%
close all

% TODO: rewrite these functions to be more general

dim = 1; % dim to plot (most to least variance explained)

% plot_nonGC_moveTransitions_singleTrials(dat,obj,newme,rez,params,dim,meta)

% plot_nonGC_moveTransitions_singleTrials_v2(dat,obj,newme,rez,params,dim,meta) % random sampling of move to quiet bouts
% plot_nonGC_moveTransitions_singleTrials_v3(dat,obj,newme,rez,params,dim,meta) % random sampling of quiet to move bouts

ndims = 10;
% plot_nonGC_moveTransitions_singleTrials_v4(dat,obj,newme,rez,params,ndims,meta) % all bouts heat map, move to quiet 
% plot_nonGC_moveTransitions_singleTrials_v5(dat,obj,newme,rez,params,ndims,meta) % all bouts heat map, quiet to move

ndims = 10;
% plot_nonGC_moveTransitions_trialAvg(dat,obj,newme,rez,params,meta,ndims) % separate tiles for each dimension

% plot_nonGC_moveTransitions_trialAvg_v2(dat,obj,newme,rez,params,meta,ndims) % all dimensions plotted on same axis, 

% plot_nonGC_moveTransitions_trialAvg_v3(dat,obj,newme,rez,params,meta,ndims) % mean,stderr across dimensions

% plot_nonGC_moveTransitions_trialAvg_v4(dat,obj,newme,rez,params,meta,ndims) % var across dimensions, move to quiet
% plot_nonGC_moveTransitions_trialAvg_v5(dat,obj,newme,rez,params,meta,ndims) % var across dimensions, quiet to move

% plot_nonGC_moveTransitions_trialAvg_v6(dat,obj,newme,rez,params,meta,ndims) % sumsqmag across dimensions, move to quiet
% plot_nonGC_moveTransitions_trialAvg_v7(dat,obj,newme,rez,params,meta,ndims) % sumsqmag across dimensions, quiet to move



























