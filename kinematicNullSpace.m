clear,clc,close all

% finds null/potent spaces for a single session based on the linear equation

%               V = WN, 

% where V is a video feature matrix and N is a matrix of 
% single trial neural activity. W is the transformation
% matrix between video features and neural activity . null(W) is a basis for
% the null space of W and colspace(W) is a basis for the row space of W (which
% is the potent space).

% * (all functions associated with this script are in the kinNullSpace directory) *

addpath(genpath(pwd))

%% PARAMETERS

% --SPECIFY WHICH ANIMAL AND SESSION TO LOAD
meta.anm = 'JEB7'; % 'JEB7'  'EKH3'  'JGR2'
meta.date = '2021-04-29'; % '2021-04-29'  '2021-08-11'  '2021-11-16'


% --SPECIFY PATH TO DATA HERE
meta.datapth = '/Users/Munib/Documents/Economo-Lab/data/';
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
params.lfads_run = ''; % 'run3' , leave empty to use most recent run
params.fcut_post_fa = 31; % if performing FA, cutoff freq to smooth rates and factors with a butterworth filter
params.feat_varToExplain = 90; % num factors for dim reduction of video features should explain this much variance
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

%% VIDEO FEATURES

% getKinematicsFromVideo() returns 2 variables
% - kin: feature matrix of size (time,features,trials). Features defined in
%         params.traj_features
% - featLeg: legend corresponding to features in kin (for 2nd dimension)
[kin,dat.featLeg] = getKinematicsFromVideo(obj,params,dat.trials);


% JAW ANGLE
jawAngle = getJawAngle(obj, dat.trials, params.alignEvent, params.advance_movement); % (time, trials)
dat.featLeg{end+1} = 'jaw_angle';

% create feature matrix, feats, and assign to a field in dat
dat.feats = cat(2,kin,reshape(jawAngle,size(jawAngle,1),1,size(jawAngle,2)));


% MOTION ENERGY
% me is a struct with fields
% - data:       cell array (ntrials,1), with motion energy for each time point
% - moveThresh: a scalar - below this value of motion energy, animal is
%               said to be stationary
% - use:        0 or 1, indicating whether to use me.data. 1 if
%               a motionEnergy*.mat file is found, 0 if file not found
me = loadMotionEnergy(obj,meta,params,dat.trials); 
if me.use
    dat.feats = cat(2,dat.feats,reshape(me.data,size(me.data,1),1,size(me.data,2))); 
    dat.featLeg{end+1} = 'motion_energy';
end
% To generate a motionEnergy*.mat file for a session, see https://github.com/economolab/videoAnalysisScripts/blob/main/motionEnergy.m


% DIMENSIONALITY REDUCTION
% many of the video features will be highly correlated, so we will perform PCA/FA
% on the matrix of features to reduce the dimensionality to a set of factors that
% best explain the movement captured by the video recordings
dat.feats_reduced = reduceDimensionVideoFeatures(dat.feats,params.feat_varToExplain,size(dat.factors,2));

disp('DONE CREATING FEATURE MATRIX AND REDUCED DIM FEATURE MATRIX')

% for i = 1:size(dat.feats_reduced,2)
%     figure(i)
% %     plot(obj.time,squeeze(dat.feats_reduced(:,i,1:numel(params.trialid{1}))))
% %     imagesc(obj.time,dat.trials,squeeze(dat.feats(:,i,:))'); title(dat.featLeg{i},'Interpreter','none')
%     imagesc(obj.time,dat.trials,squeeze(dat.feats_reduced(:,i,:))')
% %     colorbar
% %     caxis([-500 200])
% end

%% 

% We now have a struct, dat, that contains:
% - trials:        trial numbers in use
% - factors:       neural data that's been reduced using lfads or factor analysis (time,factors,trials)
% - rates:         denoised single trial neural data from lfads or smoothing (time,cells,trials)
% - featLeg:       cell array of strings describing the features in dat.feats (2nd dim)
% - feats:         displacements, velocity, jaw angle, and motion energy (time,feature,trials)
% - feats_reduced: pca reduced dat.feats (time,factors,trials);


%% REGRESSION

% here we pose the problem
%           
%               V = NW
%
% where V is a matrix of size (numVideoFeatures , time*trials)'
% and   N is a matrix of size (numFactors/numClusters , time*trials)'
% 
% We will estimate W with ridge regression (regularized linear regression).
% It is the linear transformation that takes neural activity to video
% activity
% The null space of W forms subspace that we will call the         NULL SPACE   of N
% The column space of W forms the subspace that we will call the   POTENT SPACE of N

rez = estimateW(dat,params,obj.time); % N,V are zscored neural activity and feature matrix


% V_ = N*W;
% ft = 1;
% figure; plot(V(1:1000,ft)); hold on; plot(V_(1:1000,ft))

    
%% NULL AND POTENT SPACE OF W

[rez.W_null,rez.W_potent,rez.N_null,rez.N_potent] = getNullPotentSpaces_SVD(rez.W',rez,dat); % transposing W b/c we found W by solving V'=N'*W, rather than V = WN


%% PREP TUNING

% this step is only valid if the number of null and potent dimensions are
% the same. For this to happen, you need R neural dimensions and R/2
% kinematic dimensions
if (size(rez.N,2)/2)==size(rez.V,2) 
    rez.tuning_ratio = preparatoryTuning(rez,obj.time,params);
end

%% 

% We now have a struct, rez, that contains results from regression (and
% data used)

%% random plots

% % for fa data
% figure; 
% plot(obj.time,squeeze(N_potent(:,1,1:numel(params.trialid{1}))),'b'); hold on
% trix = (numel(params.trialid{1})+1):((numel(params.trialid{1})+1) + numel(params.trialid{2}));
% plot(obj.time,squeeze(N_potent(:,1,trix)),'r')

% for lfads data
rhit = find(obj.bp.R & obj.bp.hit & ~obj.bp.autowater & ~obj.bp.early);
lhit = find(obj.bp.L & obj.bp.hit & ~obj.bp.autowater & ~obj.bp.early);
mask = ismember(dat.trials,rhit);
rhit = find(mask);
mask = ismember(dat.trials,lhit);
lhit = find(mask);

figure(1); sgtitle('Null Space Projections')
for i = 1:size(rez.N_null,2) % num null dims
    subplot(3,1,i); hold on
    for j = 1:numel(rhit)
        patchline(obj.time,squeeze(rez.N_null(:,i,rhit(j))),'EdgeColor','b','EdgeAlpha',0.4,'LineWidth',1.5);
    end
    for j = 1:numel(lhit)
        patchline(obj.time,squeeze(rez.N_null(:,i,lhit(j))),'EdgeColor','r','EdgeAlpha',0.4,'LineWidth',1.5);
    end
    xlim([obj.time(100) obj.time(end)])
    hold off
end

figure(2); sgtitle('Potent Space Projections')
for i = 1:size(rez.N_potent,2) % num null dims
    subplot(size(rez.N_potent,2),1,i); hold on
    for j = 1:numel(rhit)
        patchline(obj.time,squeeze(rez.N_potent(:,i,rhit(j))),'EdgeColor','b','EdgeAlpha',0.4,'LineWidth',1.5);
    end
    for j = 1:numel(lhit)
        patchline(obj.time,squeeze(rez.N_potent(:,i,lhit(j))),'EdgeColor','r','EdgeAlpha',0.4,'LineWidth',1.5);
    end
    xlim([obj.time(100) obj.time(end)])
    hold off
end


%% dPCA (TODO)








