clear,clc,close all


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
params.advance_movement = 0.0; % seconds, amount of time to advance movement data relative to neural data

%% NEURAL ACTIVITY

% getNeuralActivity() returns 4 main variables
% - dat: contains lfads/fa smoothed firing rates, factors, and trial numbers
%        dat.factors and dat.rates are size (time,factors/clusters,trials)
% - meta: session meta data
% - params: parameters used for preprocessing lfads input data
% - obj: preprocessed data obj
[meta,params,obj,dat] = getNeuralActivity(meta,params);

%% lfads factors and gaussian smoothed single trials projected onto first nFactors pcs

lfads = dat.factors;

gauss = obj.trialdat(:,:,dat.trials);

sm = 31; winsize = sm*params.dt;
for i = 1:size(gauss,2) % for each clu
    gauss(:,i,:) = mySmooth(squeeze(gauss(:,i,:)),sm); % causal gaussian filter
end

temp = permute(gauss,[1 3 2]);
temp = reshape(temp,size(temp,1)*size(temp,2),size(temp,3));
[~,gauss] = pca(temp,'NumComponents',size(lfads,2));


%% plot one cluster, all trials, as a sanity check
% 
% figure; 
% subplot(211)
% plot(squeeze(lfads(:,1,:)));
% subplot(212)
% plot(squeeze(gauss(:,1,:)));

%% trials

clear condition
condition(1)     = {'R&hit&~stim.enable&~autowater&~early'};         % right hits, no stim, aw off
condition(end+1) = {'L&hit&~stim.enable&~autowater&~early'};         % left hits, no stim, aw off

trialNums = findTrials(obj, condition);

mask = ismember(dat.trials,trialNums{1});
trialid{1} = find(mask);
mask = ismember(dat.trials,trialNums{2});
trialid{2} = find(mask);



%% tongue angle

[kin,dat.featLeg] = getKinematicsFromVideo(obj,params,dat.trials);

[~,mask] = patternMatchCellArray(dat.featLeg,{'jaw_xvel_view1'},'all');
ix = find(mask);
ix = ix(2);
jawvel_y = kin(:,:,ix);

jawvel_y = fillmissing(jawvel_y,'nearest');

jawvel_y = mySmooth(jawvel_y,51);

%%

[yup,ylow] = envelope(jawvel_y,51,'rms');


%%

% n = 15;
% r = randsample(trialid{1},n);
% l = randsample(trialid{2},n);
% 
% % close all
% figure(1)
% % plot(jawvel_y(:,trialid{1}),'b'); hold on; plot(jawvel_y(:,trialid{2}),'r')
% plot(yup(:,r),'b'); hold on; 
% plot(yup(:,l),'r')
% plot(ylow(:,r),'b');
% plot(ylow(:,l),'r')


%% save data

lfads = permute(lfads,[1 3 2]);
lfads = reshape(lfads,size(lfads,1)*size(lfads,2),size(lfads,3));

% gauss = permute(gauss,[1 3 2]);
% gauss = reshape(gauss,size(gauss,1)*size(gauss,2),size(gauss,3));

yup = reshape(yup,size(yup,1)*size(yup,2),1);
ylow = reshape(ylow,size(ylow,1)*size(ylow,2),1);


save('jeb7_decoding','lfads','gauss','ylow','yup')













