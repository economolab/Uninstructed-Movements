clear,clc,close all

% finds null/potent spaces for a single session based on elsayed method:


% * (most functions associated with this script are in the kinNullSpace directory) *

addpath(genpath(pwd))

%% PARAMETERS

% --SPECIFY WHICH ANIMAL AND SESSION TO LOAD
meta.anm = 'JEB7'; % 'JEB7'  'EKH3'  'JGR2'
meta.date = '2021-04-29'; % '2021-04-29'  '2021-08-11'  '2021-11-16'


% --SPECIFY PATH TO DATA HERE
if ispc
    meta.datapth = 'C:\Code\uninstructedMovements-Munib/data';
else
    meta.datapth = '/Users/Munib/Documents/Economo-Lab/data/';
end
% the data should be stored in the following structure:
% /params.datapth/
%  --- /DataObjects/
%  --- --- /animal_name_number/ (contains data_structure*.mat)
%  --- /lfads/
%  --- --- /input/ (contains anm_name_session_date*.mat/.h5)
%  --- --- /output/ (contains model_train/valid_anm_name_session_date*.h5)


% --SPECIFY PREPROCESSING PARAMETERS YOU WANT TO CHANGE BELOW THE FUNCTION CALL
% if loading lfads data, these default params will be overwritten by params
% used to generate lfads input data
params = getDefaultParams();
% params.probe = 2; % change default params.probe from '1' to '2'


% --SPECIFY METHOD OF DENOISING SINGLE TRIALS
% we need denoised, smooth single trial neural activity. We have two
% options:
% 1) load data that's ALREADY been passed through an lfads model
% 2) perform Factor Analysis on binned single trial data, followed by
%    smoothing
params.lfads_or_fa = 'lfads'; % 'lfads' or 'fa'
params.lfads_run = 'run3'; % 'run3' , leave empty to use most recent run
params.fcut_post_fa = 31; % if performing FA, cutoff freq to smooth rates and factors with a butterworth filter
params.feat_varToExplain = 80; % num factors for dim reduction of video features should explain this much variance
params.full_or_reduced = 'reduced'; % 'full'  or 'reduced' -- which data to use in regression
% using the full data will require another method, the system seems to be
% highly overfit as is
assert(strcmpi(params.full_or_reduced,'reduced'),'method to use full dimensional data doesnt exist, params.full_or_reduced should be set to `reduced`')

% --SPECIFY TIME POINTS AND LAG TO USE
params.prep = [-2.5 -0.05]; % initial and final time points (seconds) defining prep epoch, relative to alignevent
params.move = [-2.5 1.5];   % initial and final time points (seconds) defining move epoch, relative to alignevent
params.advance_movement = 0.025; % seconds, amount of time to advance movement data relative to neural data

%% NEURAL ACTIVITY

% getNeuralActivity() returns 4 main variables
% - dat: contains lfads/fa smoothed firing rates, factors, and trial numbers
%        dat.factors and dat.rates are size (time,factors/clusters,trials)
% - meta: session meta data
% - params: parameters used for preprocessing lfads input data
% - obj: preprocessed data obj
[meta,params,obj,dat] = getNeuralActivity(meta,params);
dat.factors = standardizeFactors(dat.factors);

% for i = 1:size(dat.factors,2)
%     figure(i)
% %     plot(obj.time,squeeze(dat.factors(:,i,1:numel(params.trialid{1}))))
%     imagesc(obj.time,1:numel(dat.trials),squeeze(dat.factors(:,i,:))');
% end
% 
% 
% % factor 8 shows differences within a condition, first half and second half
% figure;
% subplot(2,1,1)
% imagesc(obj.time,1:numel(params.trialid{1}),squeeze(dat.factors(:,8,1:numel(params.trialid{1})))');
% title('JEB7 4-29, Factor 8, Right Hit Trials')
% subplot(2,1,2)
% plot(obj.time,mean(temp(:,1:37),2),'LineWidth',2);
% hold on;
% plot(obj.time,mean(temp(:,38:75),2),'LineWidth',2)
% legend('First Half, avg right hit trials','Second Half, avg right hit trials')


%% MOTION ENERGY
% me is a struct with fields
% - data:       cell array (ntrials,1), with motion energy for each time point
% - moveThresh: a scalar - below this value of motion energy, animal is
%               said to be stationary
% - use:        0 or 1, indicating whether to use me.data. 1 if
%               a motionEnergy*.mat file is found, 0 if file not found
me = loadMotionEnergy(obj,meta,params,dat.trials); 

% To generate a motionEnergy*.mat file for a session, see https://github.com/economolab/videoAnalysisScripts/blob/main/motionEnergy.m


%% 

% We now have a struct, dat, that contains:
% - trials:        trial numbers in use
% - factors:       neural data that's been reduced using lfads or factor analysis (time,factors,trials)
% - rates:         denoised single trial neural data from lfads or smoothing (time,cells,trials)



%% elsayed method single trials

% TODO: variance explained of trial-averaged data

rates_or_factors = 'factors';

findcd = 0; % find CDs using elsayed method

rez = estimateQ(obj,dat,me,rates_or_factors,findcd);

rez = calVarExp_elsayed(rez);


%% plot projections

optimization_plots(rez,obj,dat,params);

f = figure(20);
pth = '/Users/Munib/Documents/Economo-Lab/code/uninstructedMovements/fig4/figs/elsayedNullSpace/';
fn = 'null_JEB7_2021-04-29_lfadsrun3';
mysavefig(f,pth,fn);

f = figure(21);
pth = '/Users/Munib/Documents/Economo-Lab/code/uninstructedMovements/fig4/figs/elsayedNullSpace/';
fn = 'potent_JEB7_2021-04-29_lfadsrun3';
mysavefig(f,pth,fn);


%% activity modes

% [latents,trials_by_type,modes] = activityModes_nullPotent(rez,dat,obj,params);
% TODO: var explained of trial averaged data, not single trials b/c the
% variability messes things up...

rez.N_null = rez.optim.N_null;
rez.N_potent = rez.optim.N_potent;

[latents,trials_by_type,modes] = cdNullSpace_elsayed(rez,dat,obj,params);

modes.varexp = activityModes_varexp(modes,rez);

plt.trial_types = [1 2];
plt.plot_mean = 0;
plt.colors = {[0 0.4470 0.7410], [0.6350 0.0780 0.1840]};
plotAllModes_nullPotent(obj,params,latents,trials_by_type,plt)
% 
% f = figure(211);
% pth = '/Users/Munib/Documents/Economo-Lab/code/uninstructedMovements/fig4/figs/elsayedNullSpace/';
% fn = 'nullCDs_JEB7_2021-04-29_lfadsrun3';
% mysavefig(f,pth,fn);
% 
% f = figure(214);
% pth = '/Users/Munib/Documents/Economo-Lab/code/uninstructedMovements/fig4/figs/elsayedNullSpace/';
% fn = 'potentCDs_JEB7_2021-04-29_lfadsrun3';
% mysavefig(f,pth,fn);



