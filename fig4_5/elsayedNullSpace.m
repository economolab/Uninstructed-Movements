clear,clc,close all


addpath(genpath(pwd))

%% PARAMETERS

% --SPECIFY WHICH ANIMAL AND SESSION TO LOAD
meta(1).anm = 'JEB7'; % 'JEB7'  'EKH3'  'JGR2'
meta(1).date = '2021-04-29'; % '2021-04-29'  '2021-08-11'  '2021-11-16'

meta(end+1).anm = 'JEB6'; % 'JEB7'  'EKH3'  'JGR2'
meta(end).date = '2021-04-18'; % '2021-04-29'  '2021-08-11'  '2021-11-16'

meta(end+1).anm = 'JGR2'; % 'JEB7'  'EKH3'  'JGR2'
meta(end).date = '2021-11-17'; % '2021-04-29'  '2021-08-11'  '2021-11-16'

meta = assignDataPath(meta);

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
dfparams = getDefaultParams();
% params.probe = 2; % change default params.probe from '1' to '2'


% --SPECIFY METHOD OF DENOISING SINGLE TRIALS
% we need denoised, smooth single trial neural activity. We have two
% options:
% 1) load data that's ALREADY been passed through an lfads model
% 2) perform Factor Analysis on binned single trial data, followed by
%    smoothing
dfparams.lfads_or_fa = 'lfads'; % 'lfads' or 'fa'
dfparams.lfads_run = 'run1'; % 'run3' , leave empty to use most recent run
dfparams.fcut_post_fa = 31; % if performing FA, cutoff freq to smooth rates and factors with a butterworth filter
dfparams.feat_varToExplain = 80; % num factors for dim reduction of video features should explain this much variance
dfparams.full_or_reduced = 'reduced'; % 'full'  or 'reduced' -- which data to use in regression
% using the full data will require another method, the system seems to be
% highly overfit as is
assert(strcmpi(dfparams.full_or_reduced,'reduced'),'method to use full dimensional data doesnt exist, params.full_or_reduced should be set to `reduced`')

% --SPECIFY TIME POINTS AND LAG TO USE
dfparams.prep = [-2.5 -0.05]; % initial and final time points (seconds) defining prep epoch, relative to alignevent
dfparams.move = [-2.5 1.5];   % initial and final time points (seconds) defining move epoch, relative to alignevent
dfparams.advance_movement = 0.025; % seconds, amount of time to advance movement data relative to neural data

%% NEURAL ACTIVITY

% getNeuralActivity() returns 4 main variables
% - dat: contains lfads/fa smoothed firing rates, factors, and trial numbers
%        dat.factors and dat.rates are size (time,factors/clusters,trials)
% - meta: session meta data
% - params: parameters used for preprocessing lfads input data
% - obj: preprocessed data obj

for i = 1:numel(meta)
    disp(['Loading data for ' meta(i).anm ' ' meta(i).date]);
    [meta(i),params(i),obj(i),dat(i)] = getNeuralActivity(meta(i),dfparams);
    gpfa(i) = getGPFAData(meta(i),'run1');
    me(i) = loadMotionEnergy(obj(i),meta(i),params(i),sort(dat(i).trials)); 
    
    disp('DONE');
    disp(' ');
end

%%

clearvars -except meta params obj dat gpfa dfparams me

%% MOTION ENERGY
% me is a struct with fields
% - data:       cell array (ntrials,1), with motion energy for each time point
% - moveThresh: a scalar - below this value of motion energy, animal is
%               said to be stationary
% - use:        0 or 1, indicating whether to use me.data. 1 if
%               a motionEnergy*.mat file is found, 0 if file not found

% To generate a motionEnergy*.mat file for a session, see https://github.com/economolab/videoAnalysisScripts/blob/main/motionEnergy.m


%% 

% We now have a struct, dat, that contains:
% - trials:        trial numbers in use
% - factors:       neural data that's been reduced using lfads or factor analysis (time,factors,trials)
% - rates:         denoised single trial neural data from lfads or smoothing (time,cells,trials)



%% elsayed method single trials
clear rez

rates_or_factors = 'factors';

for i = 1:numel(meta)
    rez(i) = estimateQ(obj(i),gpfa(i),me(i),params(i),rates_or_factors);
end

%% variance explained plots
close all

varexp_full = [rez(:).varexp_null ; rez(:).varexp_potent]';
f = figure; ax = axes(f);
vs = violinplot(varexp_full,{'Null','Potent'},...
    'EdgeColor',[1 1 1], 'ViolinAlpha',{0.2,1});
ylabel('Variance Explained of All Data')
ylim([0,1])
ax = gca;
ax.FontSize = 25;

% pth = '/Users/Munib/Documents/Economo-Lab/code/uninstructedMovements/fig4_5/figs/elsayedNullSpace';
% fn = 'run1_allSessions_varexp_allData';
% mysavefig(f,pth,fn);

varexp_sep = [rez(:).varexp_null_nonmove; rez(:).varexp_null_move ; rez(:).varexp_potent_move; rez(:).varexp_potent_nonmove]';
f = figure; ax = axes(f);
vs = violinplot(varexp_sep,{'Null, Non-move','Null, Move', 'Potent,Move','Potent, Non-move'},...
    'EdgeColor',[1 1 1], 'ViolinAlpha',{0.2,1});
ylabel('Variance Explained')
ylim([0,1])
ax = gca;
ax.FontSize = 25;

% pth = '/Users/Munib/Documents/Economo-Lab/code/uninstructedMovements/fig4_5/figs/elsayedNullSpace';
% fn = 'run1_allSessions_varexp_allSep';
% mysavefig(f,pth,fn);

%% plot projections

close all

cols = getColors();
clrs{1} = cols.rhit;
clrs{2} = cols.lhit;
for i = 1:numel(rez)
    %     optimization_plots(rez,obj,dat,params); % old
    
    temp = rez(i).optim.N_potent;
    
    f = figure;
    for dimix = 1:size(temp{1},2)
        ax = subplot(size(temp{1},2),1,dimix); hold on
        for j = 1:2
            plot(obj(1).time,temp{j}(:,dimix),'Color',clrs{j},'LineWidth',3)
        end
        xlim([obj(1).time(15),obj(1).time(end)])
        hold off;
    end
    
    temp = rez(i).optim.N_null;
    
    f = figure;
    for dimix = 1:size(temp{1},2)
        ax = subplot(size(temp{1},2),1,dimix); hold on;
        for j = 1:2
            plot(obj(1).time,temp{j}(:,dimix),'Color',clrs{j},'LineWidth',3)
        end
        xlim([obj(1).time(15),obj(1).time(end)])
        hold off;
    end
    
    
end



% f = figure(20);
% pth = '/Users/Munib/Documents/Economo-Lab/code/uninstructedMovements/fig4/figs/elsayedNullSpace/';
% fn = 'null_JEB7_2021-04-29_lfadsrun12';
% mysavefig(f,pth,fn);
% 
% f = figure(21);
% pth = '/Users/Munib/Documents/Economo-Lab/code/uninstructedMovements/fig4/figs/elsayedNullSpace/';
% fn = 'potent_JEB7_2021-04-29_lfadsrun12';
% mysavefig(f,pth,fn);


%% activity modes

% [latents,trials_by_type,modes] = activityModes_nullPotent(rez,dat,obj,params);
% TODO: var explained of trial averaged data, not single trials b/c the
% variability messes things up...

rez.N_null = rez.optim.N_null;
rez.N_potent = rez.optim.N_potent;

[latents,trials_by_type,modes] = cdNullSpace_elsayed(rez,dat,obj,params);

modes.varexp = activityModes_varexp(modes,rez);


%%
plt.trial_types = [1 2];
plt.plot_mean = 1;
plt.colors = {[0 0.4470 0.7410], [0.6350 0.0780 0.1840]};
plotAllModes_nullPotent(obj,params,latents,trials_by_type,plt)
% 
% f = figure(221);
% pth = '/Users/Munib/Documents/Economo-Lab/code/uninstructedMovements/fig4/figs/elsayedNullSpace/';
% fn = 'nullCDs_JEB7_2021-04-29_lfadsrun12';
% mysavefig(f,pth,fn);
% 
% f = figure(222);
% pth = '/Users/Munib/Documents/Economo-Lab/code/uninstructedMovements/fig4/figs/elsayedNullSpace/';
% fn = 'potentCDs_JEB7_2021-04-29_lfadsrun12';
% mysavefig(f,pth,fn);



