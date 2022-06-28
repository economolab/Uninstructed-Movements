clear,clc,close all


addpath(genpath(pwd))

%% PARAMETERS

% --SPECIFY WHICH ANIMAL AND SESSION TO LOAD
meta(1).anm = 'JEB7'; 
meta(1).date = '2021-04-29'; 

meta(end+1).anm = 'JEB7'; 
meta(end).date = '2021-04-30'; 

meta(end+1).anm = 'JEB6'; 
meta(end).date = '2021-04-18'; 

meta(end+1).anm = 'JGR2'; 
meta(end).date = '2021-11-16';

meta(end+1).anm = 'JGR2'; 
meta(end).date = '2021-11-17';

meta(end+1).anm = 'JGR3'; 
meta(end).date = '2021-11-18';

meta(end+1).anm = 'EKH3'; 
meta(end).date = '2021-08-11';

% meta(end+1).anm = 'EKH3'; 
% meta(end).date = '2021-08-07';

meta(end+1).anm = 'EKH1'; 
meta(end).date = '2021-08-07';




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

use = true(size(meta));
for i = 1:numel(meta)
    disp(['Loading data for ' meta(i).anm ' ' meta(i).date]);
%     [meta(i),params(i),obj(i),dat(i)] = getNeuralActivity(meta(i),dfparams);
    gpfa(i) = getGPFAData(meta(i),'run2');
    params(i) = gpfa(i).params;
    if isfield(gpfa(i).obj,'meta')
        gpfa(i).obj = rmfield(gpfa(i).obj,'meta');
    end
    obj(i) = gpfa(i).obj;
    me(i) = loadMotionEnergy(obj(i),meta(i),params(i),1:obj(i).bp.Ntrials); 
    if ~me(i).use
        use(i) = false;
    end
    disp('DONE');
    disp(' ');
end

gpfa = gpfa(use);
params = params(use);
meta = meta(use);
obj = obj(use);
me = me(use);

%%

clearvars -except meta params obj dat gpfa dfparams me



%% elsayed method single trials

clear rez

for i = 1:numel(meta)
    input_data = gpfa(i).gpfalatents; % dat(i).factors   dat(i).rates
    rez(i) = elsayedNullandPotentSpace(obj(i),input_data,me(i),params(i));
end

%% variance explained plots
close all

varexpPlots(rez);

%% plot projections

plotProjections(params,obj,rez)

%% activity modes

clear cdrez times

for i = 1:numel(rez)
    [cdrez(i),times] = cdNullSpace_elsayed(rez(i),obj(i),params(i));
end

rez = cdrez;

% modes.varexp = activityModes_varexp(modes,rez);

%% null space cds

nullSpaceCD(rez,obj,params,times)


%% potent space cds

potentSpaceCD(rez,obj,params,times)

